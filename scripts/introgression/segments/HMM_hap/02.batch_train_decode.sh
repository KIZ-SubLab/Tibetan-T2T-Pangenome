#!/bin/bash
#SBATCH --job-name train
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --partition=wzhcnormal
#SBATCH -e train.err
#SBATCH -o train.log

set -euo pipefail

CONDA_BASE="/work/home/songjiaxu/anaconda3"
source "${CONDA_BASE}/bin/activate" "hmmix"

PREMADE_DIR="/work/home/songjiaxu/introgression/data_HMM_intro"
WEIGHTS="${PREMADE_DIR}/hg38_strick_callability_mask.bed"
MUTRATES="${PREMADE_DIR}/hg38_mutationrate.bed"
INIT_PARAM="${PREMADE_DIR}/Initialguesses.json"
ARCHAIC_BCF="/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/data/archaic/merged_5ind_by_chr/individuals_highcov.*.bcf"
OUTDIR="/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/HMM/output_haploid"
CHROM="${1:-}"

INDIVIDUALS=(
    C101 C102 C103 C104 C105 C106 C107 C108
    TBN1 TBN2 TBN3 TBN4 TBN5 TBN6 TBN7 TBN8 TBN9
    TBN10 TBN11 TBN12 TBN13 TBN14 TBN15 TBN16 TBN17 TBN18 TBN19
    TBN20 TBN21 TBN23 TBN24 TBN25 TBN26 TBN27 TBN28
)

CHROM_ARG=""
[ -n "${CHROM}" ] && CHROM_ARG="-chrom=${CHROM}"

echo "============================================"
echo "批量训练+解码 — ${#INDIVIDUALS[@]} 个个体 (diploid)"
[ -n "${CHROM}" ] && echo "限定染色体: ${CHROM}"
echo "============================================"

for ind in "${INDIVIDUALS[@]}"; do
    OBS_FILE="${OUTDIR}/obs.${ind}.txt"
    TRAINED="${OUTDIR}/trained.${ind}.json"
    DECODE_OUT="${OUTDIR}/decode.${ind}"

    if [ ! -f "${OBS_FILE}" ]; then
        echo "[SKIP] ${ind} — 缺少 obs (先跑 01.batch_ingroup.sh)"
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
            -out="${TRAINED}" \
            -haploid \
            ${CHROM_ARG} 2>&1 | tail -1; then
            echo "OK"
        else
            echo "FAILED"
            continue
        fi
    fi

    echo -n "[$(date +%H:%M:%S)] decode ... "
    if hmmix decode \
        -obs="${OBS_FILE}" \
        -weights="${WEIGHTS}" \
        -mutrates="${MUTRATES}" \
        -param="${TRAINED}" \
        -out="${DECODE_OUT}" \
        -admixpop="${ARCHAIC_BCF}" \
        -extrainfo \
        -haploid \
        ${CHROM_ARG} 2>/dev/null; then
        echo "OK"
        ls -lh "${DECODE_OUT}.hap1.txt" "${DECODE_OUT}.hap2.txt"
    else
        echo "FAILED"
    fi
done

