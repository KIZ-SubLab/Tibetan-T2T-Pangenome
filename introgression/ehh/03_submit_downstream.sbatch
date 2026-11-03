#!/usr/bin/env bash
#SBATCH --job-name=introSV_submit_phase
#SBATCH --partition=wzhcnormal
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/output/rm1kb-na-segments/merged_fragments_v3/PBS_sv-snp/TMRCA_dating/introSV_fd_analysis/EHH/logs/03_submit_%j.out
#SBATCH --error=/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/output/rm1kb-na-segments/merged_fragments_v3/PBS_sv-snp/TMRCA_dating/introSV_fd_analysis/EHH/logs/03_submit_%j.err

set -euo pipefail
ROOT=/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/output/rm1kb-na-segments/merged_fragments_v3/PBS_sv-snp/TMRCA_dating/introSV_fd_analysis/EHH
source "$ROOT/config.sh"
python3 "$EHH_ROOT/scripts/make_phase_clusters.py" \
  --proxy "$EHH_ROOT/results/highest_mappable_polarizable_proxy.tsv" \
  --bcftools "$BCFTOOLS" --example-vcf "$SOURCE_ROOT/chr1_deepvariant_glnexus.vcf.gz" \
  --hg38-fai "${HG38_REF}.fai" --padding "$PHASE_PADDING_BP" \
  --out "$EHH_ROOT/metadata/phase_clusters.tsv"
n=$(awk 'END{print NR-1}' "$EHH_ROOT/metadata/phase_clusters.tsv")
test "$n" -gt 0
phase_job=$(sbatch --parsable --array="1-${n}%${PHASE_ARRAY_CONCURRENCY}" "$EHH_ROOT/scripts/03_phase_cluster.sbatch")
ehh_job=$(sbatch --parsable --dependency="afterok:${phase_job}" "$EHH_ROOT/scripts/04_plot_ehh.sbatch")
printf 'phase_array\t%s\t%s_tasks\nEHH_plot\t%s\n' "$phase_job" "$n" "$ehh_job" >> "$EHH_ROOT/job_ids.tsv"
echo "PHASE_JOB=$phase_job EHH_JOB=$ehh_job"
