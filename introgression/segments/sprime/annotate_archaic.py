#!/usr/bin/env python3
"""
替代 map_arch: 标注 SPrime score 文件中每个位点的古人类 allele 匹配情况

逻辑:
  SPrime ALLELE=0 → 渗入等位基因 = REF
  SPrime ALLELE=1 → 渗入等位基因 = ALT
  古人类 GT 包含渗入等位基因 → match
  古人类 GT 存在但不包含渗入等位基因 → mismatch
  位点不在 VCF 中(或被 mask 排除) → notcomp

用法 (测试):
  python3 annotate_archaic.py --chr 22 \
    --score tbn.yri.results/tbn.yri.22.score \
    --vcf archaic_filtered/Denisova25/22.vcf.gz \
    --mask mask_grch38/per_chr/nochr/22_mask.bed.gz \
    --tag Denisova25 \
    --out moss/test_py.chr22.mscore
"""

import argparse
import gzip
import sys
import os

def load_mask(mask_file):
    """加载 mask BED, 返回需要排除的 (chr, start, end) 集合"""
    masks = []
    with gzip.open(mask_file, 'rt') as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split('\t')
            masks.append((parts[0], int(parts[1]), int(parts[2])))
    return masks

def is_masked(chrom, pos, masks):
    """检查位点是否在 mask 区域内"""
    for m_chr, m_start, m_end in masks:
        if m_chr == chrom and m_start <= pos < m_end:
            return True
    return False

def load_archaic_variants(vcf_file, masks):
    """
    从 VCF 加载 variant 位点 (GT != 0/0 且 GT != ./.),
    跳过 mask 区域的位点。
    自动检测列数 (7列=无QUAL/FILTER/INFO, 10列=有)
    返回 {(chr, pos): gt_string}
    """
    variants = {}
    with gzip.open(vcf_file, 'rt') as f:
        col_layout = None  # 从 #CHROM 行自动检测
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                # 确定列结构: "FORMAT" 在第几列
                headers = line.strip().split('\t')
                if 'FORMAT' in headers:
                    fmt_col = headers.index('FORMAT')
                    sample_col = fmt_col + 1
                    col_layout = (fmt_col, sample_col)
                else:
                    # 默认: 7列格式 (FORMAT=5, SAMPLE=6 0-indexed)
                    col_layout = (8, 9)  # fallback
                continue
            if col_layout is None:
                continue

            parts = line.strip().split('\t')
            if len(parts) <= col_layout[1]:
                continue

            chrom = parts[0]
            pos = int(parts[1])
            alt = parts[4]

            # 跳过 mask 区域
            if is_masked(chrom, pos, masks):
                continue

            # 参考位点 (ALT=.) → 不需要存储, 默认就是 0/0
            if alt == '.':
                continue

            # 解析 GT
            fmt = parts[col_layout[0]].split(':')
            sample = parts[col_layout[1]].split(':')
            try:
                gt_idx = fmt.index('GT')
            except ValueError:
                continue
            gt = sample[gt_idx]

            # 只保留有 variant 的位点 (GT 含 1)
            if '1' in gt.replace('|', '/'):
                variants[(chrom, pos)] = gt

    return variants

def annotate_score(score_file, variants_list, tag_list, masks):
    """
    标注 score 文件的每个位点。
    variants_list: 每个古人类的 variant dict 列表
    tag_list: 对应的古人类标签列表
    """
    with open(score_file, 'r') as f:
        header = f.readline().strip().split('\t')
        # 输出 header
        out_header = header + tag_list
        yield out_header

        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue

            chrom = parts[0]
            pos = int(parts[1])
            ref = parts[3]
            alt = parts[4]
            allele = int(parts[6])  # 0=REF, 1=ALT

            # 渗入等位基因
            introgressed_allele = 0 if allele == 0 else 1

            # mask 区域 → 全部 notcomp
            if is_masked(chrom, pos, masks):
                annotations = ['notcomp'] * len(tag_list)
                yield parts + annotations
                continue

            annotations = []
            for variants in variants_list:
                if (chrom, pos) in variants:
                    gt = variants[(chrom, pos)]
                    # extract alleles from GT
                    gt_alleles = set()
                    for a in gt.replace('|', '/').split('/'):
                        if a.isdigit():
                            gt_alleles.add(int(a))

                    if introgressed_allele in gt_alleles:
                        annotations.append(tag_list[len(annotations)])
                    else:
                        annotations.append('mismatch')
                else:
                    # 不在 variant dict → archaic 是 0/0
                    if introgressed_allele == 0:
                        annotations.append(tag_list[len(annotations)])  # archaic REF = introgressed REF
                    else:
                        annotations.append('mismatch')  # archaic REF ≠ introgressed ALT

            yield parts + annotations


def main():
    parser = argparse.ArgumentParser(description='替代 map_arch 标注古人类 allele 匹配')
    parser.add_argument('--score', required=True, help='SPrime score 文件')
    parser.add_argument('--vcf', nargs='+', required=True, help='古人类 VCF (可多个)')
    parser.add_argument('--tag', nargs='+', required=True, help='古人类标签 (与 --vcf 一一对应)')
    parser.add_argument('--mask', required=True, help='mask BED 文件')
    parser.add_argument('--out', required=True, help='输出 mscore 文件')
    args = parser.parse_args()

    assert len(args.vcf) == len(args.tag), "--vcf 和 --tag 数量必须一致"

    print(f"加载 mask: {args.mask}", file=sys.stderr)
    masks = load_mask(args.mask)
    print(f"  mask 区域: {len(masks)}", file=sys.stderr)

    variants_list = []
    for vcf_file, tag in zip(args.vcf, args.tag):
        print(f"加载 {tag}: {vcf_file}", file=sys.stderr)
        variants = load_archaic_variants(vcf_file, masks)
        print(f"  variant 位点: {len(variants)}", file=sys.stderr)
        variants_list.append(variants)

    print(f"标注 score: {args.score}", file=sys.stderr)
    with open(args.out, 'w') as out:
        for row in annotate_score(args.score, variants_list, args.tag, masks):
            out.write('\t'.join(row) + '\n')

    print(f"输出: {args.out}", file=sys.stderr)


if __name__ == '__main__':
    main()
