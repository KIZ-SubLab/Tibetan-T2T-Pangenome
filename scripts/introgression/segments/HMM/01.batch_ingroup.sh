#!/bin/bash
#SBATCH --job-name ingroup
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --partition=wzhcnormal
#SBATCH -e ingroup.err
#SBATCH -o ingroup.log

set -euo pipefail

CONDA_BASE="/work/home/songjiaxu/anaconda3"
source "${CONDA_BASE}/bin/activate" "hmmix"

BASE_DIR="/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38"
TBN_VCF="${BASE_DIR}/data/tbn/by_chr/chr*.vcf.gz"

PREMADE_DIR="/work/home/songjiaxu/introgression/data_HMM_intro"
OUTGROUP="${PREMADE_DIR}/hg38_Outgroup_1000g_HGDP.txt"
WEIGHTS="${PREMADE_DIR}/hg38_strick_callability_mask.bed"

OUTDIR="${BASE_DIR}/HMM/output"
mkdir -p "${OUTDIR}"

INDIVIDUALS=(
    C101 C102 C103 C104 C105 C106 C107 C108
    TBN1 TBN2 TBN3 TBN4 TBN5 TBN6 TBN7 TBN8 TBN9
    TBN10 TBN11 TBN12 TBN13 TBN14 TBN15 TBN16 TBN17 TBN18 TBN19
    TBN20 TBN21 TBN23 TBN24 TBN25 TBN26 TBN27 TBN28
)

echo "============================================"
echo "批量 create_ingroup — ${#INDIVIDUALS[@]} 个个体"
echo "============================================"

SUCCESS=0
FAIL=0

for ind in "${INDIVIDUALS[@]}"; do
    OBS_FILE="${OUTDIR}/obs.${ind}.txt"

    if [ -f "${OBS_FILE}" ]; then
        echo "[SKIP] ${ind} — 文件已存在"
        ((SUCCESS++)) || true
        continue
    fi

    echo -n "[$(date +%H:%M:%S)] ${ind} ... "

    if hmmix create_ingroup \
        -ind="${ind}" \
        -vcf="${TBN_VCF}" \
        -outgroup="${OUTGROUP}" \
        -weights="${WEIGHTS}" \
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
