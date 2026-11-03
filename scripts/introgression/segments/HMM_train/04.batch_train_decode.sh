#!/bin/bash
#SBATCH --job-name train2
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --partition=wzhcnormal
#SBATCH -e train2.err
#SBATCH -o train2.log
#===============================================================================
# 步骤4+5: hmmix train + decode — 用自己的参数训练和解码
#
# 依赖: 03.batch_ingroup.sh 产出 obs 文件
# 输出: trained.<个体>.json + decode.<个体>.diploid.txt (含古人类标注)
#===============================================================================

set -euo pipefail

CONDA_BASE="/work/home/songjiaxu/anaconda3"
source "${CONDA_BASE}/bin/activate" "hmmix"

BASE_DIR="/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38"
OUTDIR="${BASE_DIR}/HMM_train/output"

WEIGHTS="${BASE_DIR}/../hg38_strick_callability_mask.bed"
MUTRATES="${OUTDIR}/mutationrate.bed"
INIT_PARAM="${BASE_DIR}/../Initialguesses.json"

# 古人类 VCF — per-chr 通配符 (hmmix 自动识别所有染色体)
ARCHAIC_VCF="/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/data/archaic/archaic_custom_bcf/individuals_highcov.*.bcf"

INDIVIDUALS=(
    C101 C102 C103 C104 C105 C106 C107 C108
    TBN1 TBN2 TBN3 TBN4 TBN5 TBN6 TBN7 TBN8 TBN9
    TBN10 TBN11 TBN12 TBN13 TBN14 TBN15 TBN16 TBN17 TBN18 TBN19
    TBN20 TBN21 TBN23 TBN24 TBN25 TBN26 TBN27 TBN28
)

echo "============================================"
echo "批量训练+解码 — ${#INDIVIDUALS[@]} 个个体 (diploid)"
echo "古人类 VCF (per-chr): ${ARCHAIC_VCF}"
echo "============================================"

for ind in "${INDIVIDUALS[@]}"; do
    OBS_FILE="${OUTDIR}/obs.${ind}.txt"
    TRAINED="${OUTDIR}/trained.${ind}.json"
    DECODE_OUT="${OUTDIR}/decode.${ind}"

    if [ ! -f "${OBS_FILE}" ]; then
        echo "[SKIP] ${ind} — 缺少 obs (先跑 03.batch_ingroup.sh)"
        continue
    fi

    echo ""
    echo "===== ${ind} ====="

    echo -n "[$(date +%H:%M:%S)] train ... "
    if [ -f "${TRAINED}" ]; then
        echo "SKIP (已存在)"
    else
        if hmmix train \
            -obs="${OBS_FILE}" \
            -weights="${WEIGHTS}" \
            -mutrates="${MUTRATES}" \
            -param="${INIT_PARAM}" \
            -out="${TRAINED}" 2>&1 | tail -1; then
            echo "OK"
        else
            echo "FAILED"
            continue
        fi
    fi

    # decode without archaic annotation (fallback)
    echo -n "[$(date +%H:%M:%S)] decode (no archaic) ... "
    if hmmix decode \
        -obs="${OBS_FILE}" \
        -weights="${WEIGHTS}" \
        -mutrates="${MUTRATES}" \
        -param="${TRAINED}" \
        -out="${DECODE_OUT}" \
        -extrainfo 2>/dev/null; then
        echo "OK"
        ls -lh "${DECODE_OUT}.diploid.txt"
    else
        echo "FAILED"
    fi

    # decode WITH archaic annotation
    echo -n "[$(date +%H:%M:%S)] decode (with archaic) ... "
    if hmmix decode \
        -obs="${OBS_FILE}" \
        -weights="${WEIGHTS}" \
        -mutrates="${MUTRATES}" \
        -param="${TRAINED}" \
        -out="${DECODE_OUT}.archaic" \
        -admixpop="${ARCHAIC_VCF}" \
        -extrainfo 2>/dev/null; then
        echo "OK"
        ls -lh "${DECODE_OUT}.archaic.diploid.txt"
    else
        echo "FAILED"
    fi
done

echo ""
echo "============================================"
echo "全部完成！"
echo "============================================"
