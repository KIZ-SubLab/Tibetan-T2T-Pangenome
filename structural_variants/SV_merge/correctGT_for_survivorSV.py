#!/usr/bin/env python
import os
import sys
import gzip
from icecream import ic


def smart_open(file, opera):
    if opera == 'r':
        if os.path.exists(file) ==False:
            print("Can not open file {}".format(file))
            exit()
        else:
            if file.endswith(".gz"):
                out = gzip.open(file, 'rt')
            else:
                out = open(file, 'r')
    elif opera == 'w':
        if file.endswith(".gz"):
            out = gzip.open(file, 'wt')
        else:
            out = open(file, 'w')
    return out


def calculate_gt_score(gt):
    score = 0
    for char in gt:
        if char in '01':
            score += 2
        elif char == '/':
            score += 0.5
        elif char == '.':
            score -= 0.5
        elif char == '|':
            score += 1
    return score

def more_informative(gts):
    # 为每个 caller 计算 GT 分数
    gtss = {}
    for gt in gts:
        gtss[gt] = calculate_gt_score(gt)

    # 根据分数进行降序排序
    sorted_callers = sorted(gtss.items(), key=lambda item: item[1], reverse=True)
    return sorted_callers[0][0]


def parse_gt(info_format, info):
    items = info_format.split(":")
    gt_index = 0
    for index, content in enumerate(items):
        if content == "GT":
            gt_index = index
            break
    samples_gt = []
    for i in info:
        tmp = i.split(":")
        samples_gt.append(tmp[gt_index])
    return samples_gt

def correct_gt(input_file, output_file):
    out = smart_open(output_file, 'w')

    with smart_open(input_file, 'r') as file:
        for line in file:
            line = line.strip()
            tmp = line.split("\t")
            if line.startswith("##"):
                print(line, file=out)
            elif line.startswith("#CHROM"):
                new_out = "\t".join(tmp[:10])
                print(new_out, file=out)
            else:
                gt_info = parse_gt(tmp[8], tmp[9:12])
                #ic(gt_info)

                # 选择 GT 信息
                final_gt = more_informative(gt_info)
                # for haploid-resolved case, using "|"
                final_gt = final_gt.replace("/", "|")
                out_info = "\t".join(tmp[:8]) + "\tGT\t" + final_gt
                print(out_info, file=out)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("get consensus SV info from bcftools merged VCF")
        sys.exit(f"python3 {sys.argv[0]} *.vcf.gz output.vcf[.gz]")
    else:
        # 执行合并
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        correct_gt(input_file, output_file)

