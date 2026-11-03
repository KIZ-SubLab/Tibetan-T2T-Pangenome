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

def redundancy_removal(consensus_callers):
    r = set()
    for a,b in consensus_callers:
        r.add(a)
        r.add(b)
    return r

def similar_len(l1, l2):
    min_one = min(l1, l2)
    len_diff = abs(l1 - l2)
    if len_diff/min_one <= 0.1:
        return True
    else:
        return False


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
    gtss = {}
    for gt in gts:
        gtss[gt] = calculate_gt_score(gt)

    sorted_callers = sorted(gtss.keys(), key=lambda item: item[1], reverse=True)
    return sorted_callers[0]


def merge_sv_calls(input_file, output_file):
    out = smart_open(output_file, 'w')
    # VCF must be merge using this order
    input_order = {
        'pav': 0,
        'svimasm': 1,
        'lgvar': 2,
	    'pbsv':3,
	    'sniffles':4
    }

    with smart_open(input_file, 'r') as file:
        for line in file:
            line = line.strip()
            tmp = line.split("\t")
            if line.startswith("##"):
                print(line, file=out)
            elif line.startswith("#CHROM"):
                new_out = "\t".join(tmp[:10])
                print(new_out, file=out)
            elif ";" not in tmp[2]:
                continue
            else:
                results = {}

                # 处理 caller 信息
                chr_name, position = tmp[0], tmp[1]
                position = int(position) + 1
                caller_info = tmp[2].split(';')
                gt_info = tmp[9:14]
                caller_sv_info = {}
                svtypes = set()

                # 提取每个 caller 的信息
                for caller in caller_info:
                    caller_name, _, _, type_, length = caller.split("-")
                    caller_sv_info[caller_name] = {
                            'svtype': type_,
                            'svlen': int(length),
                            'pos': int(position),
                            'gt': gt_info[input_order[caller_name]]
                        }
                    svtypes.add(type_) # collect all svtype

                # 进行合并逻辑
                final_svtype = ''
                final_svlen = 0
                final_gt = ''
                consensus_callers = set()

                # 检查任何两个 caller 是否满足条件
                caller_names = list(caller_sv_info.keys())
                for i in range(len(caller_names)):
                    for j in range(i + 1, len(caller_names)):
                        caller1 = caller_names[i]
                        caller2 = caller_names[j]

                        # same svtype and length difference < 50 bp, so regard them as a cluster/same sv
                        if caller_sv_info[caller1]['svtype'] == caller_sv_info[caller2]['svtype']:
                            #if (abs(caller_sv_info[caller1]['svlen'] - caller_sv_info[caller2]['svlen']) < 50):
                            if similar_len(caller_sv_info[caller1]['svlen'], caller_sv_info[caller2]['svlen']):
                                # push this pair into a set
                                consensus_callers.add((caller1, caller2))
                # make pairs into single callers
                consensus_callers = redundancy_removal(consensus_callers)
                if len(consensus_callers) == 0:
                    continue
                final_supported_callers = "+".join(consensus_callers)

                # sort by svlen, then choose longest one
                # 根据 svlen 值进行逆向排序
                sorted_callers = sorted(caller_sv_info.keys(), key=lambda item: caller_sv_info[item]['svlen'], reverse=True)
                #ic(sorted_callers)
                final_svlen = caller_sv_info[sorted_callers[0]]['svlen']
                final_svtype = caller_sv_info[sorted_callers[0]]['svtype']


                consensus_callers_gt = []
                for caller in consensus_callers:
                    consensus_callers_gt.append(caller_sv_info[caller]['gt'])

                # 选择 GT 信息
                if len(consensus_callers_gt) == 0:
                    ic(chr_name, position)
                final_gt = more_informative(consensus_callers_gt)
                # for haploid-resolved case, using "|"
                final_gt = final_gt.replace("/", "|")
                # 如果有足够的条件则保存最终结果
                if len(consensus_callers) > 1:
                    results[f"{final_svtype}-{final_svlen}"] = final_gt
                    tmp[2] = f"{final_supported_callers}-{chr_name}-{position}-{final_svtype}-{final_svlen}"
                    tmp[9] = final_gt
                    out_info = "\t".join(tmp[:10])
                    print(out_info, file=out)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("get consensus SV info from bcftools merged VCF")
        sys.exit(f"python3 {sys.argv[0]} *.vcf.gz output.vcf[.gz]")
    else:
        # 执行合并
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        merge_sv_calls(input_file, output_file)

