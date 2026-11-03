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

def reduce_info(line):
    tmp = line.split("\t")
    info = tmp[7]
    atrri = {}
    a = info.split(";")
    for i in a:
        if "=" not in i:
            continue
        b = i.split("=")
        atrri[b[0]] = b[1]
    #SVTYPE=DEL;END=3696;SVLEN=-62
    tmp[7] = f"SVTYPE={atrri['SVTYPE']};END={atrri['END']};SVLEN={atrri['SVLEN']}"
    return "\t".join(tmp)



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

def vcf2dict(vcf):
    vcf_records = {}
    comp_INFO = []
    with smart_open(vcf, 'r') as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("##INFO="):
                comp_INFO.append(line)
            elif line.startswith("#"):
                continue
            else:
                tmp = line.split("\t")
                chrid, pos = tmp[0], tmp[1]
                reduced_line = reduce_info(line)
                if chrid not in vcf_records:
                    vcf_records.update({chrid: {pos: reduced_line}})
                else:
                    vcf_records[chrid].update({pos: reduced_line})
    return comp_INFO, vcf_records

def make_union(main_vcf, comp_vcf):
    comp_INFO, comp_records = vcf2dict(comp_vcf)
    header = []
    body = []
    with smart_open(main_vcf, 'r') as file:
        for line in file:
            line = line.strip()
            tmp = line.split("\t")
            if line.startswith("#"):
                header.append(line)
            else:
                body.append(line)
                chrid, pos = tmp[0],tmp[1]
                if chrid in comp_records:
                    if pos in comp_records[chrid]:
                        # del records that already in main vcf
                        del comp_records[chrid][pos]

    # print remaining records
    for c,p in comp_records.items():
        for pos in p:
            body.append(p[pos])

    print("\n".join(header))
    #print("\n".join(comp_INFO))
    print("\n".join(body))


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("add complementary vcf to main vcf")
        sys.exit(f"python3 {sys.argv[0]} main.vcf[.gz] complementary.vcf[.gz]")
    else:
        # 执行合并
        main_vcf = sys.argv[1]
        complementary_vcf = sys.argv[2]
        make_union(main_vcf, complementary_vcf)

