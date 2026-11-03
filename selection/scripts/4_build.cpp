// 4_build.cpp - build the site-level null database for one shard.
// usage: 4_build <file_list> <out_dir> <work_dir> <threads> <map_dir> <map_suffix> <mutation_map>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <unistd.h>
#include <omp.h>

#include "db_common.h"
using namespace db;

static const int N_GROUP = 128;
static const size_t FLUSH_N = 65536;

struct __attribute__((packed)) Rec {
    double logr, logmu, pbs;
    uint32_t leaf;
};

struct Map {
    std::vector<int64_t> pos;
    std::vector<double> rate;
    std::vector<int64_t> bstart;
    std::vector<double> bmu;
};

static void load_maps(const std::string& dir, const std::string& suffix,
                      const std::string& mutation_map, Map maps[23]) {
    for (int c = 1; c <= 22; c++) {
        std::string path = dir + "/chr" + std::to_string(c) + suffix;
        FILE* f = fopen(path.c_str(), "r");
        if (!f) { fprintf(stderr, "cannot open %s\n", path.c_str()); exit(2); }
        char line[512];
        fgets(line, sizeof(line), f);
        while (fgets(line, sizeof(line), f)) {
            char* p = line;
            long long pos = strtoll(p, &p, 10);
            double rate = strtod(p, nullptr);
            maps[c].pos.push_back(pos);
            maps[c].rate.push_back(rate);
        }
        fclose(f);
    }
    GzLineReader rd(mutation_map);
    const char* line;
    rd.next();  // header
    while ((line = rd.next())) {
        const char* ch = fld(line, 0);
        if (strncmp(ch, "chr", 3) != 0) continue;
        int c = atoi(ch + 3);
        if (c < 1 || c > 22) continue;
        maps[c].bstart.push_back(strtoll(fld(line, 1), nullptr, 10));
        maps[c].bmu.push_back(strtod(fld(line, 9), nullptr));
    }
}

int main(int argc, char** argv) {
    if (argc < 8) {
        fprintf(stderr,
            "usage: %s <file_list> <out_dir> <work_dir> <threads> "
            "<recomb_map_dir> <recomb_map_suffix> <mutation_map>\n", argv[0]);
        return 2;
    }
    std::string listp = argv[1], out = argv[2], work = argv[3];
    int threads = atoi(argv[4]);
    std::string map_dir = argv[5], map_suffix = argv[6], mut_map = argv[7];
    omp_set_num_threads(threads);

    if (access((out + "/DONE").c_str(), F_OK) == 0) {
        printf("skip %s (DONE exists)\n", out.c_str());
        return 0;
    }
    std::string cmd = "mkdir -p " + out + " " + work;
    if (system(cmd.c_str()) != 0) return 2;

    Map maps[23];
    load_maps(map_dir, map_suffix, mut_map, maps);

    // File list: "<chrom>\t<path>" per line.
    std::vector<std::pair<int, std::string>> files;
    {
        FILE* f = fopen(listp.c_str(), "r");
        if (!f) { perror(listp.c_str()); return 2; }
        char lb[4096];
        while (fgets(lb, sizeof(lb), f)) {
            size_t L = strlen(lb);
            while (L && (lb[L - 1] == '\n' || lb[L - 1] == '\r')) lb[--L] = 0;
            if (!L) continue;
            char* tab = strchr(lb, '\t');
            if (!tab) { fprintf(stderr, "bad list line: %s\n", lb); return 2; }
            *tab = 0;
            int c = atoi(lb + 3);
            files.push_back({c, std::string(tab + 1)});
        }
        fclose(f);
    }
    printf("input files: %zu\n", files.size());

    // ---- pass 1: annotate + bucket by coarse cell ----
    std::vector<FILE*> gfile(threads * N_GROUP, nullptr);
    std::vector<std::vector<Rec>> gbuf(threads * N_GROUP);
    for (auto& v : gbuf) v.reserve(FLUSH_N);
    int64_t n_sites = 0, n_out = 0;
    double t0 = omp_get_wtime();

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int64_t local_sites = 0, local_out = 0;
        #pragma omp for schedule(dynamic, 1)
        for (size_t fi = 0; fi < files.size(); fi++) {
            int chrom = files[fi].first;
            Map& M = maps[chrom];
            GzLineReader rd(files[fi].second);
            const char* line;
            while ((line = rd.next())) {
                int64_t pos = strtoll(fld(line, 0), nullptr, 10);
                double pbs = strtod(fld(line, 1), nullptr);
                if (!std::isfinite(pbs)) continue;

                auto it = std::upper_bound(M.pos.begin(), M.pos.end(), pos);
                size_t gi = (it == M.pos.begin()) ? 1 : (size_t)(it - M.pos.begin());
                double r = M.rate[gi - 1];
                auto it2 = std::upper_bound(M.bstart.begin(), M.bstart.end(), pos);
                size_t bi = (it2 == M.bstart.begin()) ? 1 : (size_t)(it2 - M.bstart.begin());
                double mu = M.bmu[bi - 1];

                double logr = log10(r + ER), logmu = log10(mu + EM);
                if (logr < GRID_R_MIN || logr > GRID_R_MAX ||
                    logmu < GRID_MU_MIN || logmu > GRID_MU_MAX) {
                    local_out++;
                    continue;
                }
                int cr = coarse_r_idx(logr), cmu = coarse_mu_idx(logmu);
                int g = (int)(((int64_t)cr * N_COARSE_MU + cmu) * N_GROUP / N_COARSE);
                auto& buf = gbuf[tid * N_GROUP + g];
                buf.push_back({logr, logmu, pbs, (uint32_t)leaf_id(logr, logmu)});
                if (buf.size() >= FLUSH_N) {
                    FILE*& fh = gfile[tid * N_GROUP + g];
                    if (!fh) {
                        char p[2048];
                        snprintf(p, sizeof(p), "%s/g_t%02d_g%03d.bin", work.c_str(), tid, g);
                        fh = fopen(p, "wb");
                        if (!fh) { perror(p); exit(2); }
                    }
                    fwrite(buf.data(), sizeof(Rec), buf.size(), fh);
                    buf.clear();
                }
                local_sites++;
            }
        }
        for (int g = 0; g < N_GROUP; g++) {
            auto& buf = gbuf[tid * N_GROUP + g];
            if (buf.empty()) continue;
            FILE*& fh = gfile[tid * N_GROUP + g];
            if (!fh) {
                char p[2048];
                snprintf(p, sizeof(p), "%s/g_t%02d_g%03d.bin", work.c_str(), tid, g);
                fh = fopen(p, "wb");
            }
            fwrite(buf.data(), sizeof(Rec), buf.size(), fh);
            std::vector<Rec>().swap(buf);
        }
        #pragma omp atomic
        n_sites += local_sites;
        #pragma omp atomic
        n_out += local_out;
    }
    for (auto& f : gfile) if (f) fclose(f);
    printf("pass1: %lld sites, %lld out of grid\n", (long long)n_sites, (long long)n_out);

    // ---- pass 2: per group -> per leaf -> contexts + merged arrays ----
    FILE* f_lr = fopen((out + "/logr.bin").c_str(), "wb");
    FILE* f_lm = fopen((out + "/logmu.bin").c_str(), "wb");
    FILE* f_pb = fopen((out + "/pbs.bin").c_str(), "wb");
    FILE* f_cr = fopen((out + "/ctx_logr.bin").c_str(), "wb");
    FILE* f_cm = fopen((out + "/ctx_logmu.bin").c_str(), "wb");
    FILE* f_co = fopen((out + "/ctx_offset.bin").c_str(), "wb");
    FILE* f_cc = fopen((out + "/ctx_count.bin").c_str(), "wb");
    FILE* f_fl = fopen((out + "/pbs_flat.bin").c_str(), "wb");
    FILE* f_ce = fopen((out + "/cell_ctx.bin").c_str(), "wb");
    if (!f_lr || !f_lm || !f_pb || !f_cr || !f_cm || !f_co || !f_cc || !f_fl || !f_ce) {
        perror("output");
        return 2;
    }

    std::vector<uint64_t> leaf_cnt(N_LEAF, 0), cell_off(N_LEAF + 1, 0);
    uint64_t n_ctx = 0, flat_pos = 0, site_pos = 0;

    for (int g = 0; g < N_GROUP; g++) {
        std::vector<Rec> recs;
        for (int t = 0; t < threads; t++) {
            char p[2048];
            snprintf(p, sizeof(p), "%s/g_t%02d_g%03d.bin", work.c_str(), t, g);
            FILE* f = fopen(p, "rb");
            if (!f) continue;
            fseek(f, 0, SEEK_END);
            size_t n = (size_t)ftell(f) / sizeof(Rec);
            fseek(f, 0, SEEK_SET);
            size_t base = recs.size();
            recs.resize(base + n);
            if (n && fread(recs.data() + base, sizeof(Rec), n, f) != n) return 2;
            fclose(f);
            remove(p);
        }
        if (recs.empty()) continue;

        uint32_t lmin = N_LEAF, lmax = 0;
        for (const auto& r : recs) { lmin = std::min(lmin, r.leaf); lmax = std::max(lmax, r.leaf); }
        std::vector<uint64_t> cnt(lmax - lmin + 2, 0), cur(lmax - lmin + 2, 0);
        for (const auto& r : recs) cnt[r.leaf - lmin + 1]++;
        for (size_t i = 1; i < cnt.size(); i++) cnt[i] += cnt[i - 1];
        cur = cnt;
        std::vector<Rec> sorted(recs.size());
        for (const auto& r : recs) sorted[cur[r.leaf - lmin]++] = r;
        std::vector<Rec>().swap(recs);

        for (uint32_t lf = lmin; lf <= lmax; lf++) {
            uint64_t c = cnt[lf - lmin + 1] - cnt[lf - lmin];
            if (!c) continue;
            Rec* base = sorted.data() + cnt[lf - lmin];
            std::sort(base, base + c, [](const Rec& x, const Rec& y) {
                if (x.logr != y.logr) return x.logr < y.logr;
                if (x.logmu != y.logmu) return x.logmu < y.logmu;
                return x.pbs < y.pbs;
            });

            std::vector<double> cr_, cm_;
            std::vector<uint64_t> co_, cc_;
            uint64_t i = 0;
            while (i < c) {
                uint64_t j = i + 1;
                while (j < c && base[j].logr == base[i].logr && base[j].logmu == base[i].logmu) j++;
                cr_.push_back(base[i].logr);
                cm_.push_back(base[i].logmu);
                co_.push_back(flat_pos);
                cc_.push_back(j - i);
                for (uint64_t k = i; k < j; k++) fwrite(&base[k].pbs, 8, 1, f_fl);
                flat_pos += (j - i);
                i = j;
            }
            uint32_t nc = (uint32_t)cr_.size();
            fwrite(cr_.data(), 8, nc, f_cr);
            fwrite(cm_.data(), 8, nc, f_cm);
            fwrite(co_.data(), 8, nc, f_co);
            fwrite(cc_.data(), 8, nc, f_cc);
            for (uint32_t k = 0; k < nc; k++) {
                uint32_t id = (uint32_t)(n_ctx + k);
                fwrite(&id, 4, 1, f_ce);
            }
            n_ctx += nc;

            std::vector<Rec> merged(base, base + c);
            std::sort(merged.begin(), merged.end(),
                      [](const Rec& x, const Rec& y) { return x.pbs < y.pbs; });
            for (const auto& r : merged) {
                fwrite(&r.logr, 8, 1, f_lr);
                fwrite(&r.logmu, 8, 1, f_lm);
                fwrite(&r.pbs, 8, 1, f_pb);
            }
            site_pos += c;
            leaf_cnt[lf] = c;
            cell_off[lf + 1] = n_ctx;
        }
        std::vector<Rec>().swap(sorted);
    }
    fclose(f_lr); fclose(f_lm); fclose(f_pb); fclose(f_cr); fclose(f_cm);
    fclose(f_co); fclose(f_cc); fclose(f_fl); fclose(f_ce);

    std::vector<uint64_t> leaf_off(N_LEAF + 1, 0);
    for (int64_t l = 1; l < N_LEAF; l++) leaf_off[l] = leaf_off[l - 1] + leaf_cnt[l - 1];
    leaf_off[N_LEAF] = leaf_off[N_LEAF - 1] + leaf_cnt[N_LEAF - 1];

    std::vector<uint64_t> co(N_LEAF + 1, 0);
    uint64_t cur_ = 0;
    for (int64_t l = 0; l < N_LEAF; l++) {
        cur_ = cell_off[l + 1] ? cell_off[l + 1] : cur_;
        co[l + 1] = cur_;
    }
    write_npy_u64(out + "/leaf_offset.npy", leaf_off);
    write_npy_u64(out + "/leaf_count.npy", leaf_cnt);
    write_npy_u64(out + "/cell_off.npy", co);

    double wall = omp_get_wtime() - t0;
    FILE* js = fopen((out + "/stats.json").c_str(), "w");
    fprintf(js, "{\"n_site\": %llu, \"n_ctx\": %llu, \"wall_s\": %.1f}\n",
            (unsigned long long)site_pos, (unsigned long long)n_ctx, wall);
    fclose(js);
    FILE* d = fopen((out + "/DONE").c_str(), "w");
    fprintf(d, "ok\n");
    fclose(d);
    printf("DONE: sites=%llu contexts=%llu wall=%.1f s\n",
           (unsigned long long)site_pos, (unsigned long long)n_ctx, wall);
    return 0;
}
