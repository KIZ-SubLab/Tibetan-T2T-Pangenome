#!/usr/bin/env bash
#SBATCH --job-name=introSV_proxy
#SBATCH --partition=wzhcnormal
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --time=03:00:00
#SBATCH --output=/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/output/rm1kb-na-segments/merged_fragments_v3/PBS_sv-snp/TMRCA_dating/introSV_fd_analysis/EHH/logs/02_proxy_%j.out
#SBATCH --error=/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/output/rm1kb-na-segments/merged_fragments_v3/PBS_sv-snp/TMRCA_dating/introSV_fd_analysis/EHH/logs/02_proxy_%j.err

set -euo pipefail
ROOT=/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/output/rm1kb-na-segments/merged_fragments_v3/PBS_sv-snp/TMRCA_dating/introSV_fd_analysis/EHH
source "$ROOT/config.sh"

python3 "$EHH_ROOT/scripts/rank_ld_candidates.py" \
  --ld "$EHH_ROOT/results/all_SVs.CZA35.1Mb.ld.gz" \
  --ids "$EHH_ROOT/metadata/unique_original_sv_ids.txt" --snp-meta "$EHH_ROOT/metadata/snp_metadata.tsv" \
  --out "$EHH_ROOT/work/all_ranked_LD_candidates.tsv"
# EHH requires a fully observed core in all 35 diploid samples (AN=70).
python3 "$EHH_ROOT/scripts/filter_complete_candidates.py" \
  --input "$EHH_ROOT/work/all_ranked_LD_candidates.tsv" --required-an 70 \
  --output "$EHH_ROOT/work/all_ranked_LD_candidates.AN70.tsv"
awk 'BEGIN{FS=OFS="\t"} NR>1 && !seen[$3]++{print $4,$5-1,$5,$3,0,"+"}' \
  "$EHH_ROOT/work/all_ranked_LD_candidates.tsv" > "$EHH_ROOT/work/all_LD_candidates.CHM13.bed"
"$LIFTOVER" "$EHH_ROOT/work/all_LD_candidates.CHM13.bed" "$CHAIN" \
  "$EHH_ROOT/work/all_LD_candidates.GRCh38.bed" "$EHH_ROOT/work/all_LD_candidates.unmapped.bed"
python3 "$EHH_ROOT/scripts/select_valid_proxy.py" \
  --candidates "$EHH_ROOT/work/all_ranked_LD_candidates.AN70.tsv" \
  --ids "$EHH_ROOT/metadata/unique_original_sv_ids.txt" --sv-meta "$EHH_ROOT/metadata/sv_metadata.tsv" \
  --mapped "$EHH_ROOT/work/all_LD_candidates.GRCh38.bed" --unmapped "$EHH_ROOT/work/all_LD_candidates.unmapped.bed" \
  --hg38-ref "$HG38_REF" --ancestral-dir "$ANCESTRAL_DIR" \
  --out "$EHH_ROOT/results/highest_mappable_polarizable_proxy.tsv"

expected=$(grep -cve '^$' "$EHH_ROOT/metadata/unique_original_sv_ids.txt")
observed=$(awk -F '\t' 'NR>1 && $1=="ok"{n++}END{print n+0}' "$EHH_ROOT/results/highest_mappable_polarizable_proxy.tsv")
printf 'expected_unique_SVs\t%s\nproxy_status_ok\t%s\n' "$expected" "$observed" > "$EHH_ROOT/results/proxy_selection_counts.tsv"
test "$observed" -eq "$expected"
echo PROXY_SELECTION_COMPLETE
