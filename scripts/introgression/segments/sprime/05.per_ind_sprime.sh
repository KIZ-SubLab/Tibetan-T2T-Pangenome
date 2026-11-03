#!/bin/bash
#SBATCH --job-name perIndSPrime
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --partition=wzhcnormal
#SBATCH -e perIndSPrime.err
#SBATCH -o perIndSPrime.log
#===============================================================================
# 个体级别 SPrime 渗入片段长度 (参考 a10086cao.sh 逻辑: addinfo.pl/addinfo2.pl + awk)
#
# 论文标准:
#   Neanderthal: archaic_allele=1, AltaiNean=match, Denisova=mismatch, Denisova25=mismatch,
#                NeaRate>0.6, DenRate<0.4, Den25Rate<0.4
#   Denisovan:   archaic_allele=1, AltaiNean=mismatch, (Denisova=match OR Denisova25=match),
#                NeaRate<0.3, (DenRate>0.4 OR Den25Rate>0.4)
#   Unknown:     archaic_allele=1, 其余
#===============================================================================

set -euo pipefail

SPRIME_DIR="/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/sprime"
MOSS="${SPRIME_DIR}/moss"
MATCH_SUMMARY="${SPRIME_DIR}/../output/rm1kb-na-segments/consensus/sprime/match.summary.txt"
VCF_DIR="/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/data/tbn/by_chr"
OUTDIR="${SPRIME_DIR}/per_ind"
ADINFO_PL="/work/home/songjiaxu/introgression/INDCAL/addinfo.pl"
ADINFO2_PL="/work/home/songjiaxu/introgression/INDCAL/addinfo2.pl"

mkdir -p "${OUTDIR}"/{Neanderthal,Denisovan,Unknown,WC}

INDIVIDUALS=(C101 C102 C103 C104 C105 C106 C107 C108
             TBN1 TBN2 TBN3 TBN4 TBN5 TBN6 TBN7 TBN8 TBN9
             TBN10 TBN11 TBN12 TBN13 TBN14 TBN15 TBN16 TBN17 TBN18 TBN19
             TBN20 TBN21 TBN23 TBN24 TBN25 TBN26 TBN27 TBN28)

# ============================================================
# addinfo.pl 输出列 (mscore行在前, clean GT在后, 1-based):
#   1-chr(无chr前缀)  2-pos  3-ID  4-REF_m  5-ALT_m
#   6-SEGMENT  7-ALLELE  8-SCORE  9-AltaiNean  10-Denisova  11-Denisova25
#   12-ref_c  13-alt_c  14-gt
#
# addinfo2.pl 追加的 segments 列 (从 join 输出):
#   15-seg_min_pos  16-seg_max_pos  17-n_variants
#   18-Nea_rate  19-Den_rate  20-Den25_rate
#
# 关键: 标签名本身 = match, "mismatch" = mismatch, "notcomp" = not comparable
# ============================================================

for sample in "${INDIVIDUALS[@]}"; do
    echo "========== ${sample} =========="

    for chr in {1..22}; do
        mscore="${MOSS}/out.chr${chr}.mscore"
        vcf="${VCF_DIR}/chr${chr}.vcf.gz"
        [ -f "$mscore" ] || continue
        [ -f "$vcf" ] || continue

        # 1. 提取样本 GT: chr pos ref alt gt
        bcftools view -s "${sample}" -c1 "${vcf}" 2>/dev/null | \
            grep -v '#' | awk -v c="$chr" '{print "chr"c, $2, $4, $5, $10}' \
            > "${sample}.chr${chr}.score.clean"

        [ -s "${sample}.chr${chr}.score.clean" ] || continue

        # 2. addinfo.pl: 按 POS 合并个体 GT 与 mscore 信息
        perl "$ADINFO_PL" "${sample}.chr${chr}.score.clean" "$mscore" \
            "out.chr${chr}.mscore.${sample}"

        # 3. 计算每个 SEGMENT 在个体中的范围 (min_pos, max_pos)
        awk '
        {
            seg = $6    # SEGMENT (mscore在前, addinfo.pl输出第6列)
            pos = $2    # POS
            if (seg in min_pos) {
                if (pos < min_pos[seg]) min_pos[seg] = pos
                if (pos > max_pos[seg]) max_pos[seg] = pos
            } else {
                min_pos[seg] = pos
                max_pos[seg] = pos
            }
        }
        END {
            for (s in min_pos)
                print s, min_pos[s], max_pos[s]
        }' "out.chr${chr}.mscore.${sample}" | \
            sort -k1,1n > "tmp.${sample}.chr${chr}.ranges"

        # 4. join match.summary + ranges → segments 文件
        awk -v chr="$chr" '$1 == chr' "$MATCH_SUMMARY" | \
            sort -k2,2n | \
            join -1 2 -2 1 --nocheck-order -o 1.1,1.2,2.2,2.3,1.5,1.6,1.7,1.8 \
                 - "tmp.${sample}.chr${chr}.ranges" | \
            sort -k2,2n > "${sample}.chr${chr}.segments"

        # 5. addinfo2.pl: 按 SEGMENT 追加片段信息
        perl "$ADINFO2_PL" "${sample}.chr${chr}.segments" \
            "out.chr${chr}.mscore.${sample}" \
            "out.chr${chr}.mscore.${sample}.addSeg"

        # 6. 按论文标准分类
        #    列索引: 7=ALLELE, 9=AltaiNean, 10=Denisova, 11=Denisova25
        #           15=min, 16=max, 18=Nea_rate, 19=Den_rate, 20=Den25_rate
        #    标签值: "AltaiNean"/"Denisova"/"Denisova25"=match, "mismatch"=mismatch
        #
        #    Neanderthal: archaic=1, matches Nea, mismatches both Denisovan
        awk -F'\t' '{
            if ($7==1 && $9=="AltaiNean" && $10=="mismatch" && $11=="mismatch" && $18>0.6 && $19<0.4 && $20<0.4)
                print
        }' "out.chr${chr}.mscore.${sample}.addSeg" > "${sample}.chr${chr}.Nea"

        #    Denisovan: archaic=1, mismatches Nea, matches at least one Denisovan
        awk -F'\t' '{
            if ($7==1 && $9=="mismatch" && ($10=="Denisova" || $11=="Denisova25") && $18<0.3 && ($19>0.4 || $20>0.4))
                print
        }' "out.chr${chr}.mscore.${sample}.addSeg" > "${sample}.chr${chr}.Den"

        #    Unknown: archaic=1, 其余
        awk -F'\t' '{
            if ($7==1 && $9=="mismatch" && $10=="mismatch" && $11=="mismatch")
                print
        }' "out.chr${chr}.mscore.${sample}.addSeg" > "${sample}.chr${chr}.Un"

        # 7. 计算总长度 (segment 坐标合并去重, count>29 的 segment 算入)
        for origin in Nea Den Un; do
            label=""
            case $origin in
                Nea) label="Neanderthal" ;;
                Den) label="Denisovan" ;;
                Un)  label="Unknown" ;;
            esac

            awk -F'\t' 'BEGIN { total=0 }
            {
                key = $15 "_" $16
                count[key]++
            }
            END {
                for (coord in count) {
                    if (count[coord] > 29) {
                        split(coord, pos, "_")
                        total += (pos[2] - pos[1])
                    }
                }
                printf "%.4f\n", (total / 1000000)
            }' "${sample}.chr${chr}.${origin}" >> "${OUTDIR}/WC/${sample}.allchr.${label}.rls"

            # 输出个体-染色体-片段详情
            awk -F'\t' -v origin="$label" '{
                key = $15 "_" $16
                seg_data[key]["min"] = $15
                seg_data[key]["max"] = $16
                seg_data[key]["count"]++
            }
            END {
                for (coord in seg_data) {
                    if (seg_data[coord]["count"] > 29) {
                        split(coord, pos, "_")
                        len = pos[2] - pos[1]
                        print pos[1], pos[2], len, seg_data[coord]["count"], origin
                    }
                }
            }' "${sample}.chr${chr}.${origin}" >> "${OUTDIR}/${label}/${sample}.txt"
        done

        # 清理中间文件
        rm -f "${sample}.chr${chr}.score.clean" \
              "out.chr${chr}.mscore.${sample}" \
              "out.chr${chr}.mscore.${sample}.addSeg" \
              "tmp.${sample}.chr${chr}.ranges" \
              "${sample}.chr${chr}.segments" \
              "${sample}.chr${chr}.Nea" "${sample}.chr${chr}.Den" "${sample}.chr${chr}.Un"
    done
done

# ============================================================
# 汇总
# ============================================================
echo ""
echo "=== 个体 SPrime 渗入片段长度 (Mb) ==="
echo -e "individual\tNeanderthal\tDenisovan\tUnknown" | tee "${OUTDIR}/per_ind_summary.txt"
for ind in "${INDIVIDUALS[@]}"; do
    n_nea=$(awk -F'[ \t]+' '{sum+=$3} END{printf "%.1f", sum/1e6+0}' "${OUTDIR}/Neanderthal/${ind}.txt" 2>/dev/null || echo "0.0")
    n_den=$(awk -F'[ \t]+' '{sum+=$3} END{printf "%.1f", sum/1e6+0}' "${OUTDIR}/Denisovan/${ind}.txt" 2>/dev/null || echo "0.0")
    n_unk=$(awk -F'[ \t]+' '{sum+=$3} END{printf "%.1f", sum/1e6+0}' "${OUTDIR}/Unknown/${ind}.txt" 2>/dev/null || echo "0.0")
    echo -e "$ind\t${n_nea}\t${n_den}\t${n_unk}" | tee -a "${OUTDIR}/per_ind_summary.txt"
done

echo "" | tee -a "${OUTDIR}/per_ind_summary.txt"
echo "=== 平均 ===" | tee -a "${OUTDIR}/per_ind_summary.txt"
awk '{nea+=$2; den+=$3; unk+=$4; n++} END{printf "平均\t%.1f\t%.1f\t%.1f\n", nea/n, den/n, unk/n}' "${OUTDIR}/per_ind_summary.txt" | tee -a "${OUTDIR}/per_ind_summary.txt"

echo ""
echo "$(date): 完成! 输出: ${OUTDIR}/"
