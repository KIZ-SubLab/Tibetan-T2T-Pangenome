#!/usr/bin/env bash
#SBATCH --job-name=introSV_EHH
#SBATCH --partition=wzhcnormal
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --time=04:00:00
#SBATCH --output=/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/output/rm1kb-na-segments/merged_fragments_v3/PBS_sv-snp/TMRCA_dating/introSV_fd_analysis/EHH/logs/04_EHH_%j.out
#SBATCH --error=/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/output/rm1kb-na-segments/merged_fragments_v3/PBS_sv-snp/TMRCA_dating/introSV_fd_analysis/EHH/logs/04_EHH_%j.err

set -euo pipefail
ROOT=/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/output/rm1kb-na-segments/merged_fragments_v3/PBS_sv-snp/TMRCA_dating/introSV_fd_analysis/EHH
source "$ROOT/config.sh"
python3 "$EHH_ROOT/scripts/run_ehh_all.py" \
  --manifest "$EHH_ROOT/metadata/group_manifest.tsv" \
  --proxy "$EHH_ROOT/results/highest_mappable_polarizable_proxy.tsv" \
  --clusters "$EHH_ROOT/metadata/phase_clusters.tsv" --phased-dir "$EHH_ROOT/vcf/phased" \
  --bcftools "$BCFTOOLS" --out-root "$EHH_ROOT" --window "$EHH_WINDOW_BP"
