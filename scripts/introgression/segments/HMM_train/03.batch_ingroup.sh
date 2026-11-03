#!/bin/bash
#SBATCH --job-name ingroup2
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --partition=wzhcnormal
#SBATCH -e ingroup2.err
#SBATCH -o ingroup2.log
#===============================================================================
# 步骤3: hmmix create_ingroup — 用自己的 outgroup 过滤 TBN 个体
#
# 依赖: 01.create_outgroup.sh 产出 outgroup.txt
# 输出: obs.<个体>.txt × 35

set -euo pipefail

BASE_DIR="/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38"
OUTDIR="${BASE_DIR}/HMM_train/output"
mkdir -p "${OUTDIR}"

# 输入
TBN_VCF="${BASE_DIR}/data/tbn/by_chr/chr*.vcf.gz"
OUTGROUP="${OUTDIR}/outgroup.txt"
WEIGHTS="${BASE_DIR}/../hg38_strick_callability_mask.bed"
ANCESTRAL="${BASE_DIR}/../hg38_ancestral/homo_sapiens_ancestor_*.fa"

INDIVIDUALS=(
    C101 C102 C103 C104 C105 C106 C107 C108
    TBN1 TBN2 TBN3 TBN4 TBN5 TBN6 TBN7 TBN8 TBN9
    TBN10 TBN11 TBN12 TBN13 TBN14 TBN15 TBN16 TBN17 TBN18 TBN19
    TBN20 TBN21 TBN23 TBN24 TBN25 TBN26 TBN27 TBN28
)

echo "============================================"
echo "batch create_ingroup — ${#INDIVIDUALS[@]} 个个体 (用自己的 outgroup)"
echo "============================================"

SUCCESS=0
FAIL=0

for ind in "${INDIVIDUALS[@]}"; do
    OBS_FILE="${OUTDIR}/obs.${ind}.txt"

    if [ -f "${OBS_FILE}" ]; then
        echo "[SKIP] ${ind} — 已存在"
        ((SUCCESS++)) || true
        continue
    fi

    echo -n "[$(date +%H:%M:%S)] ${ind} ... "

    if hmmix create_ingroup \
        -ind="${ind}" \
        -vcf="${TBN_VCF}" \
        -outgroup="${OUTGROUP}" \
        -weights="${WEIGHTS}" \
        -ancestral="${ANCESTRAL}" \
        -out="${OUTDIR}/obs" 2>/dev/null; then
        n=$(wc -l < "${OBS_FILE}")
        echo "OK (${n} 行)"
        ((SUCCESS++)) || true
    else
        echo "FAILED"
        ((FAIL++)) || true
    fi
done

echo ""
echo "============================================"
echo "完成: 成功=${SUCCESS}, 失败=${FAIL}"
echo "============================================"
