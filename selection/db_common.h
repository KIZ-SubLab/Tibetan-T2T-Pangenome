#pragma once
// db_common.h - shared constants and helpers for the site-level null database.
//
// Statistical definition (fixed):
//   LOG_R  = log10(R_point  + 1e-4)
//   LOG_MU = log10(MU_point + 1e-12)
//   match  = |LOG_R_sim - LOG_R_obs| <= T  AND  |LOG_MU_sim - LOG_MU_obs| <= T
//   N_MATCH < MIN_MATCH -> insufficient match
//   P_SIM  = (N_GE + 1) / (N_MATCH + 1)

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <zlib.h>

namespace db {

constexpr double ER = 1e-4;
constexpr double EM = 1e-12;
constexpr int64_t MIN_MATCH = 10000;

// Two-level 2D grid over (LOG_R, LOG_MU).
constexpr double GRID_R_MIN = -4.5;
constexpr double GRID_R_MAX = 2.2;
constexpr double GRID_MU_MIN = -9.5;
constexpr double GRID_MU_MAX = -5.5;
constexpr int N_COARSE_R = 568;
constexpr int N_COARSE_MU = 60;
constexpr int N_FINE = 8;  // fine cells per coarse dimension
constexpr int64_t N_COARSE = N_COARSE_R * N_COARSE_MU;
constexpr int64_t N_LEAF = N_COARSE * N_FINE * N_FINE;  // 2,181,120
constexpr double DR = (GRID_R_MAX - GRID_R_MIN) / N_COARSE_R;
constexpr double DMU = (GRID_MU_MAX - GRID_MU_MIN) / N_COARSE_MU;
constexpr double FDR = DR / N_FINE;
constexpr double FDMU = DMU / N_FINE;

inline int coarse_r_idx(double logr) {
    double x = (logr - GRID_R_MIN) / DR;
    if (x < 0) return 0;
    int i = (int)x;
    return i >= N_COARSE_R ? N_COARSE_R - 1 : i;
}
inline int coarse_mu_idx(double logmu) {
    double x = (logmu - GRID_MU_MIN) / DMU;
    if (x < 0) return 0;
    int i = (int)x;
    return i >= N_COARSE_MU ? N_COARSE_MU - 1 : i;
}
inline int fine_r_idx(double logr, int cr) {
    double x = (logr - (GRID_R_MIN + cr * DR)) / FDR;
    if (x < 0) return 0;
    int i = (int)x;
    return i >= N_FINE ? N_FINE - 1 : i;
}
inline int fine_mu_idx(double logmu, int cmu) {
    double x = (logmu - (GRID_MU_MIN + cmu * DMU)) / FDMU;
    if (x < 0) return 0;
    int i = (int)x;
    return i >= N_FINE ? N_FINE - 1 : i;
}
inline int64_t leaf_id(double logr, double logmu) {
    int cr = coarse_r_idx(logr), cmu = coarse_mu_idx(logmu);
    int fr = fine_r_idx(logr, cr), fm = fine_mu_idx(logmu, cmu);
    return ((int64_t)cr * N_COARSE_MU + cmu) * (N_FINE * N_FINE) + fr * N_FINE + fm;
}

// Minimal .npy writer (row-major, little-endian, 1-D).
inline bool write_npy_u64(const std::string& path, const std::vector<uint64_t>& v) {
    std::string hdr = "{'descr': '<u8', 'fortran_order': False, 'shape': (" +
                      std::to_string(v.size()) + ",), }";
    size_t pad = (64 - ((10 + hdr.size() + 1) % 64)) % 64;
    hdr += std::string(pad, ' ');
    hdr += "\n";
    unsigned short hlen = (unsigned short)hdr.size();
    FILE* f = fopen(path.c_str(), "wb");
    if (!f) { perror(path.c_str()); return false; }
    unsigned char pre[10] = {0x93, 'N', 'U', 'M', 'P', 'Y', 1, 0, 0, 0};
    pre[8] = (unsigned char)(hlen & 0xFF);
    pre[9] = (unsigned char)(hlen >> 8);
    fwrite(pre, 1, 10, f);
    fwrite(hdr.data(), 1, hdr.size(), f);
    fwrite(v.data(), 8, v.size(), f);
    fclose(f);
    return true;
}

// Buffered gzip line reader; the returned pointer is valid until the next call.
struct GzLineReader {
    gzFile gz;
    char* buf;
    size_t cap;
    size_t lo = 0, len = 0;
    explicit GzLineReader(const std::string& path, size_t buf_cap = 1 << 26)
        : gz(gzopen(path.c_str(), "rb")), buf(new char[buf_cap]), cap(buf_cap) {}
    ~GzLineReader() { if (gz) gzclose(gz); delete[] buf; }
    const char* next() {
        for (;;) {
            for (size_t i = lo; i < len; i++) {
                if (buf[i] == '\n') {
                    buf[i] = '\0';
                    const char* r = buf + lo;
                    lo = i + 1;
                    return r;
                }
            }
            if (!gz) {
                if (lo >= len) return nullptr;
                buf[len] = '\0';
                const char* r = buf + lo;
                lo = len + 1;
                return r;
            }
            if (lo > 0) { memmove(buf, buf + lo, len - lo); len -= lo; lo = 0; }
            if (len >= cap - 1) { fprintf(stderr, "line too long\n"); exit(2); }
            int got = gzread(gz, buf + len, cap - 1 - len);
            if (got <= 0) { gzclose(gz); gz = nullptr; continue; }
            len += (size_t)got;
        }
    }
};

// Tab-separated field access (fields are tab-terminated, not NUL-terminated).
inline const char* fld(const char* line, int i) {
    for (int k = 0; k < i; k++) {
        while (*line && *line != '\t') line++;
        if (*line == '\t') line++;
    }
    return line;
}

}  // namespace db
