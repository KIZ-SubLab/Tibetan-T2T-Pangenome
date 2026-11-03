// 5_query.cpp - exact empirical P values against the null database.
// usage: 5_query <shard_dirs_csv> <queries> <out> <threads> [tol] [--full-cell-mode merged|context]

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <omp.h>

#include "db_common.h"
using namespace db;

static const double EPS = 1e-12;

struct DB {
    const double* pbs = nullptr;          // merged per-cell PBS
    const uint64_t* off = nullptr;
    const uint64_t* cnt = nullptr;
    int64_t n_leaf = 0;
    const double* ctx_r = nullptr;
    const double* ctx_m = nullptr;
    const uint64_t* ctx_off = nullptr;
    const uint64_t* ctx_cnt = nullptr;
    const double* pbs_flat = nullptr;
    const uint64_t* cell_off = nullptr;
    const uint32_t* cell_ctx = nullptr;
    int64_t n_ctx = 0;
};

static void* mmap_file(const std::string& path, size_t* size) {
    int fd = open(path.c_str(), O_RDONLY);
    if (fd < 0) { perror(path.c_str()); exit(2); }
    struct stat st;
    if (fstat(fd, &st) != 0) { perror("fstat"); exit(2); }
    *size = (size_t)st.st_size;
    void* p = mmap(nullptr, *size, PROT_READ, MAP_SHARED, fd, 0);
    if (p == MAP_FAILED) { perror("mmap"); exit(2); }
    close(fd);
    return p;
}

static const uint64_t* read_npy_u64(const std::string& path, int64_t* n) {
    size_t sz;
    unsigned char* p = (unsigned char*)mmap_file(path, &sz);
    if (sz < 16 || memcmp(p, "\x93NUMPY", 6) != 0) { fprintf(stderr, "bad npy\n"); exit(2); }
    unsigned short hlen = (unsigned short)(p[8] | (p[9] << 8));
    long long nn = 0;
    const char* sp = strstr((const char*)(p + 10), "'shape': (");
    if (!sp || sscanf(sp, "'shape': (%lld", &nn) != 1) {
        fprintf(stderr, "bad npy header\n"); exit(2);
    }
    if (n) *n = (int64_t)nn;
    return (const uint64_t*)(p + 10 + hlen);
}

static void load_db(const std::string& dir, DB& db) {
    size_t sz;
    db.pbs = (const double*)mmap_file(dir + "/pbs.bin", &sz);
    db.off = read_npy_u64(dir + "/leaf_offset.npy", &db.n_leaf);
    db.cnt = read_npy_u64(dir + "/leaf_count.npy", nullptr);
    int64_t n;
    db.ctx_r = (const double*)mmap_file(dir + "/ctx_logr.bin", &sz);
    db.n_ctx = (int64_t)(sz / 8);
    db.ctx_m = (const double*)mmap_file(dir + "/ctx_logmu.bin", &sz);
    db.ctx_off = (const uint64_t*)mmap_file(dir + "/ctx_offset.bin", &sz);
    db.ctx_cnt = (const uint64_t*)mmap_file(dir + "/ctx_count.bin", &sz);
    db.pbs_flat = (const double*)mmap_file(dir + "/pbs_flat.bin", &sz);
    db.cell_off = read_npy_u64(dir + "/cell_off.npy", &n);
    db.cell_ctx = (const uint32_t*)mmap_file(dir + "/cell_ctx.bin", &sz);
    fprintf(stderr, "shard %s: leaves=%lld contexts=%lld\n", dir.c_str(),
            (long long)db.n_leaf, (long long)db.n_ctx);
}

static inline void context_query(const DB& db, int64_t leaf,
                                 double lo_r, double hi_r, double lo_m, double hi_m,
                                 double pbs_obs, int64_t* nm, int64_t* ng) {
    uint64_t a = db.cell_off[leaf], b = db.cell_off[leaf + 1];
    for (uint64_t i = a; i < b; i++) {
        uint32_t c = db.cell_ctx[i];
        double r = db.ctx_r[c], m = db.ctx_m[c];
        if (r >= lo_r && r <= hi_r && m >= lo_m && m <= hi_m) {
            uint64_t cn = db.ctx_cnt[c], of = db.ctx_off[c];
            *nm += (int64_t)cn;
            const double* ps = db.pbs_flat + of;
            int64_t k = std::lower_bound(ps, ps + cn, pbs_obs) - ps;
            *ng += (int64_t)cn - k;
        }
    }
}

static void query_rect(const DB& db, bool merged_full, double tol,
                       double logr, double logmu, double pbs_obs,
                       int64_t* n_match, int64_t* n_ge) {
    double lo_r = logr - tol, hi_r = logr + tol, lo_m = logmu - tol, hi_m = logmu + tol;
    int64_t nm = 0, ng = 0;
    auto cr0 = coarse_r_idx(lo_r), cr1 = coarse_r_idx(hi_r);
    auto cm0 = coarse_mu_idx(lo_m), cm1 = coarse_mu_idx(hi_m);
    for (int cr = cr0; cr <= cr1; cr++) {
        double clo_r = GRID_R_MIN + cr * DR, chi_r = clo_r + DR;
        if (!(chi_r > lo_r && clo_r < hi_r)) continue;
        bool r_full = clo_r >= lo_r - EPS && chi_r <= hi_r + EPS;
        for (int cm = cm0; cm <= cm1; cm++) {
            double clo_m = GRID_MU_MIN + cm * DMU, chi_m = clo_m + DMU;
            if (!(chi_m > lo_m && clo_m < hi_m)) continue;
            bool m_full = clo_m >= lo_m - EPS && chi_m <= hi_m + EPS;
            bool coarse_full = r_full && m_full;
            int64_t base = ((int64_t)cr * N_COARSE_MU + cm) * (N_FINE * N_FINE);
            for (int fr = 0; fr < N_FINE; fr++) {
                double flo_r = clo_r + fr * FDR, fhi_r = flo_r + FDR;
                if (!coarse_full && !(fhi_r > lo_r && flo_r < hi_r)) continue;
                bool rf = flo_r >= lo_r - EPS && fhi_r <= hi_r + EPS;
                for (int fm = 0; fm < N_FINE; fm++) {
                    double flo_m = clo_m + fm * FDMU, fhi_m = flo_m + FDMU;
                    if (!coarse_full && !(fhi_m > lo_m && flo_m < hi_m)) continue;
                    int64_t leaf = base + fr * N_FINE + fm;
                    int64_t c = (int64_t)db.cnt[leaf];
                    if (!c) continue;
                    bool full = coarse_full || (rf &&
                                flo_m >= lo_m - EPS && fhi_m <= hi_m + EPS);
                    if (full && merged_full) {
                        int64_t o = (int64_t)db.off[leaf];
                        nm += c;
                        const double* ps = db.pbs + o;
                        int64_t k = std::lower_bound(ps, ps + c, pbs_obs) - ps;
                        ng += c - k;
                    } else {
                        context_query(db, leaf, lo_r, hi_r, lo_m, hi_m, pbs_obs, &nm, &ng);
                    }
                }
            }
        }
    }
    *n_match = nm;
    *n_ge = ng;
}

int main(int argc, char** argv) {
    if (argc < 5) {
        fprintf(stderr, "usage: %s <shard_dirs_csv> <queries> <out> <threads> "
                        "[tol] [--full-cell-mode merged|context]\n", argv[0]);
        return 2;
    }
    std::string dirs = argv[1], qpath = argv[2], opath = argv[3];
    int threads = atoi(argv[4]);
    double tol = 0.02;
    bool merged_full = true;
    for (int i = 5; i < argc; i++) {
        if (!strncmp(argv[i], "--full-cell-mode=", 17))
            merged_full = strncmp(argv[i] + 17, "context", 7) != 0;
        else if (argv[i][0] != '-')
            tol = strtod(argv[i], nullptr);
    }
    omp_set_num_threads(threads);

    std::vector<DB> dbs;
    size_t pos = 0;
    for (;;) {
        size_t comma = dirs.find(',', pos);
        std::string d = dirs.substr(pos, comma == std::string::npos
                                         ? std::string::npos : comma - pos);
        if (!d.empty()) { DB db; load_db(d, db); dbs.push_back(db); }
        if (comma == std::string::npos) break;
        pos = comma + 1;
    }

    std::vector<double> qr, qm, qp;
    {
        FILE* f = fopen(qpath.c_str(), "r");
        if (!f) { perror(qpath.c_str()); return 2; }
        char* line = nullptr;
        size_t cap = 0;
        while (getline(&line, &cap, f) > 0) {
            double a, b, c;
            if (sscanf(line, "%lf\t%lf\t%lf", &a, &b, &c) != 3) return 2;
            qr.push_back(a); qm.push_back(b); qp.push_back(c);
        }
        free(line);
        fclose(f);
    }
    int64_t Q = (int64_t)qr.size();
    std::vector<int64_t> order(Q), nmv(Q), ngv(Q);
    for (int64_t i = 0; i < Q; i++) order[i] = i;
    std::sort(order.begin(), order.end(), [&](int64_t a, int64_t b) {
        return qr[a] < qr[b] || (qr[a] == qr[b] && qm[a] < qm[b]);
    });

    double t0 = omp_get_wtime();
    #pragma omp parallel for schedule(dynamic, 256)
    for (int64_t oi = 0; oi < Q; oi++) {
        int64_t i = order[oi];
        int64_t nm = 0, ng = 0;
        for (const auto& db : dbs) {
            int64_t a, b;
            query_rect(db, merged_full, tol, qr[i], qm[i], qp[i], &a, &b);
            nm += a; ng += b;
        }
        nmv[i] = nm; ngv[i] = ng;
    }
    double wall = omp_get_wtime() - t0;

    FILE* fo = fopen(opath.c_str(), "w");
    if (!fo) { perror(opath.c_str()); return 2; }
    for (int64_t i = 0; i < Q; i++)
        fprintf(fo, "%d\t%.6g\t%lld\t%lld\n",
                nmv[i] >= MIN_MATCH ? 0 : -1, tol,
                (long long)nmv[i], (long long)ngv[i]);
    fclose(fo);
    fprintf(stderr, "%lld queries, %d threads, %.1f s -> %.1f queries/s\n",
            (long long)Q, threads, wall, Q / wall);
    return 0;
}
