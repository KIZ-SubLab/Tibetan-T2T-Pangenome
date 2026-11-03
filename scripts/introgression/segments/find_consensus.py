#!/usr/bin/env python3
"""
三方法共识筛选: 任意 overlap (>0) 即认为两方法支持, ≥2 方法检出则保留

输入: rm1kb-na-segments/{sprime,hmm,hmm_train}/
输出: rm1kb-na-segments/consensus/{sprime,hmm,hmm_train}/

SPrime: match.summary.txt (chr, seg, from, to, ...), 第一列无 chr 前缀
HMM:    decode.<个体>.diploid.txt (chrom, start, end, length, state, ...)
HMM_train: decode.<个体>.archaic.diploid.txt (同上)
"""

import sys, os, glob

# ============================================================
# 1. 加载片段 → {seg_id: (chr, start, end)}
# ============================================================

def load_sprime(path):
    """SPrime: 第一列无 chr 前缀, 加上前缀匹配 HMM"""
    bed = {}
    with open(path) as f:
        next(f)  # skip header
        for line in f:
            p = line.strip().split('\t')
            sid = f"sprime_{p[0]}_{p[1]}"          # chr + seg 保证唯一
            bed[sid] = ("chr" + p[0], int(p[2]), int(p[3]))
    return bed


def load_hmm(directory, prefix):
    """HMM/HMM_train 的 per-individual 文件"""
    bed = {}
    for fpath in sorted(glob.glob(f"{directory}/*.txt")):
        fname = os.path.basename(fpath)
        ind = (fname.replace("decode.", "")
                     .replace(".diploid.txt", "")
                     .replace(".archaic", ""))
        with open(fpath) as fh:
            next(fh)
            for i, line in enumerate(fh):
                p = line.strip().split('\t')
                sid = f"{prefix}_{ind}_{i}"
                bed[sid] = (p[0], int(p[1]), int(p[2]))
    return bed


# ============================================================
# 2. 按染色体建索引 → {chr: [(start, end, seg_id), ...]}
# ============================================================

def build_index(bed):
    idx = {}
    for sid, (c, s, e) in bed.items():
        idx.setdefault(c, []).append((s, e, sid))
    for c in idx:
        idx[c].sort(key=lambda x: x[0])   # 按 start 排序, 加速区间查询
    return idx


# ============================================================
# 3. 找被支持的片段: 与 idx_b 任意 overlap (>0) 即算支持
# ============================================================

def find_supported(bed_a, idx_b):
    supported = set()
    for sid, (c, s_a, e_a) in bed_a.items():
        if c not in idx_b:
            continue
        for s_b, e_b, sid_b in idx_b[c]:
            if e_b <= s_a:      # b 在 a 左边, 继续找下一个
                continue
            if s_b >= e_a:      # b 在 a 右边, 后续的更远, 停止
                break
            supported.add(sid)  # overlap > 0
            break
    return supported


# ============================================================
# 4. 筛选: 至少被另一个方法支持
# ============================================================

def consensus(sprime_bed, hmm_bed, htrain_bed,
              idx_sprime, idx_hmm, idx_htrain):

    sp_ok = (find_supported(sprime_bed, idx_hmm) |
             find_supported(sprime_bed, idx_htrain))
    hm_ok = (find_supported(hmm_bed, idx_sprime) |
             find_supported(hmm_bed, idx_htrain))
    ht_ok = (find_supported(htrain_bed, idx_sprime) |
             find_supported(htrain_bed, idx_hmm))
    return sp_ok, hm_ok, ht_ok


# ============================================================
# 5. 过滤原始文件并输出
# ============================================================

def filter_sprime(infile, outfile, trusted):
    with open(infile) as fin, open(outfile, "w") as fout:
        fout.write(next(fin))
        for line in fin:
            p = line.strip().split('\t')
            if f"sprime_{p[0]}_{p[1]}" in trusted:
                fout.write(line)


def filter_hmm(indir, outdir, prefix, trusted):
    os.makedirs(outdir, exist_ok=True)
    for fpath in sorted(glob.glob(f"{indir}/*.txt")):
        fname = os.path.basename(fpath)
        ind = (fname.replace("decode.", "")
                     .replace(".diploid.txt", "")
                     .replace(".archaic", ""))
        with open(fpath) as fin, open(f"{outdir}/{fname}", "w") as fout:
            fout.write(next(fin))
            for i, line in enumerate(fin):
                if f"{prefix}_{ind}_{i}" in trusted:
                    fout.write(line)


# ============================================================
# main
# ============================================================

def main():
    if len(sys.argv) < 2:
        sys.exit("Usage: python3 find_consensus.py <rm1kb-na-segments_dir>")

    base = sys.argv[1]
    out  = f"{base}/consensus"

    # 加载
    sprime  = load_sprime(f"{base}/sprime/match.summary.txt")
    hmm     = load_hmm(f"{base}/hmm", "hmm")
    htrain  = load_hmm(f"{base}/hmm_train", "htrain")

    print(f"SPrime:     {len(sprime):>6} segments")
    print(f"HMM:        {len(hmm):>6} segments")
    print(f"HMM_train:  {len(htrain):>6} segments")

    # 建索引
    idx_sprime = build_index(sprime)
    idx_hmm    = build_index(hmm)
    idx_htrain = build_index(htrain)

    # 共识筛选
    sp_ok, hm_ok, ht_ok = consensus(
        sprime, hmm, htrain, idx_sprime, idx_hmm, idx_htrain)

    print(f"\n=== >=2 方法支持 (any overlap > 0) ===")
    print(f"SPrime:     {len(sp_ok):>6} / {len(sprime):>6} "
          f"({100*len(sp_ok)/len(sprime):.1f}%)")
    print(f"HMM:        {len(hm_ok):>6} / {len(hmm):>6} "
          f"({100*len(hm_ok)/len(hmm):.1f}%)")
    print(f"HMM_train:  {len(ht_ok):>6} / {len(htrain):>6} "
          f"({100*len(ht_ok)/len(htrain):.1f}%)")

    # 输出
    filter_sprime(f"{base}/sprime/match.summary.txt",
                  f"{out}/sprime/match.summary.txt", sp_ok)
    filter_hmm(f"{base}/hmm", f"{out}/hmm", "hmm", hm_ok)
    filter_hmm(f"{base}/hmm_train", f"{out}/hmm_train", "htrain", ht_ok)

    # 统计
    n_s = sum(1 for _ in open(f"{out}/sprime/match.summary.txt")) - 1
    n_h = sum(sum(1 for _ in open(f"{out}/hmm/{os.path.basename(f)}")) - 1
              for f in glob.glob(f"{out}/hmm/*.txt"))
    n_t = sum(sum(1 for _ in open(f"{out}/hmm_train/{os.path.basename(f)}")) - 1
              for f in glob.glob(f"{out}/hmm_train/*.txt"))

    print(f"\n=== 最终 consensus 片段数 ===")
    print(f"SPrime:     {n_s}")
    print(f"HMM:        {n_h}")
    print(f"HMM_train:  {n_t}")
    print(f"输出: {out}/")


if __name__ == "__main__":
    main()
