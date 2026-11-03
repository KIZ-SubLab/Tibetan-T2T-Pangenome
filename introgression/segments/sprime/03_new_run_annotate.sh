#!/bin/bash
#SBATCH --job-name annotPy
#SBATCH --array=1-22%22
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --partition=wzhcnormal
#SBATCH -e tmp/annotPy_%a.err
#SBATCH -o tmp/annotPy_%a.log
#===============================================================================
# 用 Python 替代 map_arch, 标注 SPrime 渗入片段的古人类来源
# 3个古人类: Altai, Denisova, Denisova25
#===============================================================================

set -euo pipefail

chr="${SLURM_ARRAY_TASK_ID}"
DIR="/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38"

python3 ${DIR}/sprime/annotate_archaic.py \
    --score ${DIR}/sprime/tbn.yri.results/tbn.yri.${chr}.score \
    --vcf \
        ${DIR}/data/archaic_filtered/Altai/${chr}.vcf.gz \
        ${DIR}/data/archaic_filtered/Denisova/${chr}.vcf.gz \
        ${DIR}/data/archaic_filtered/Denisova25/${chr}.vcf.gz \
    --tag AltaiNean Denisova Denisova25 \
    --mask ${DIR}/../mask_grch38/per_chr/nochr/${chr}_mask.bed.gz \
    --out ${DIR}/sprime/moss/out.chr${chr}.mscore

echo "chr${chr}: $(awk -F'\t' 'NR>1{print $9,$10,$11}' ${DIR}/sprime/moss/out.chr${chr}.mscore | sort | uniq -c | sort -rn | head -5)"
