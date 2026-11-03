#!/bin/bash
#SBATCH --job-name mutrate
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --partition=wzhcnormal
#SBATCH -e mutrate.err
#SBATCH -o mutrate.log

set -euo pipefail

CONDA_BASE="/work/home/songjiaxu/anaconda3"
source "${CONDA_BASE}/bin/activate" "hmmix"

BASE_DIR="/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38"
OUTDIR="${BASE_DIR}/HMM_train/output"

OUTGROUP="${OUTDIR}/outgroup.txt"
WEIGHTS="${BASE_DIR}/../hg38_strick_callability_mask.bed"
MUTRATES_OUT="${OUTDIR}/mutationrate.bed"

hmmix mutation_rate \
    -outgroup="${OUTGROUP}" \
    -weights="${WEIGHTS}" \
    -out="${MUTRATES_OUT}" \
    -window_size=1000000

echo "完成: ${MUTRATES_OUT}"
wc -l "${MUTRATES_OUT}"
head -3 "${MUTRATES_OUT}"
