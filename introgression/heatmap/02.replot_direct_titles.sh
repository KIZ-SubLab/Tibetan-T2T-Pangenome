#!/usr/bin/env bash
# Re-render direct-phased panels after shortening titles.  Existing SNP and
# haplotype tables are reused; this job does not rerun filtering or phasing.
# Required environment: HEAMAP_ROOT, MODE=den|nea|nea_ceu.
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3500M
#SBATCH --time=00:30:00

set -eo pipefail

: "${HEAMAP_ROOT:?HEAMAP_ROOT is required}"
: "${MODE:?MODE is required}"
: "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"

analysis_mode="${MODE%%_*}"
analysis_root="${HEAMAP_ROOT}/${analysis_mode}"
line="$(awk -F '\t' -v id="${SLURM_ARRAY_TASK_ID}" 'NR > 1 && $1 == id {print; exit}' "${analysis_root}/sv_manifest.tsv")"
if [[ -z "${line}" ]]; then
  printf 'No manifest row for task %s in %s\n' "${SLURM_ARRAY_TASK_ID}" "${analysis_root}" >&2
  exit 1
fi
IFS=$'\t' read -r task_id chr start end sv_type mid sv_id <<< "${line}"
case_dir="${analysis_root}/${sv_id}"

source /work/home/songjiaxu/anaconda3/etc/profile.d/conda.sh
conda activate plink
set -u

case "${MODE}" in
  den)
    Rscript "${analysis_root}/scripts/plot_direct_phased200.R" \
      "${case_dir}" "${mid}" "${sv_id}" Denisovan den 200
    ;;
  nea)
    Rscript "${analysis_root}/scripts/plot_direct_phased200.R" \
      "${case_dir}" "${mid}" "${sv_id}" Neanderthal nea 200
    ;;
  nea_ceu)
    Rscript "${analysis_root}/scripts_ceu35/plot_direct_ceu35.R" \
      "${case_dir}" "${case_dir}/direct_phased200/with_CEU35" \
      "${mid}" "${sv_id}" 200
    ;;
  *)
    printf 'Unsupported MODE: %s\n' "${MODE}" >&2
    exit 2
    ;;
esac

printf 'Replotted %s %s\n' "${MODE}" "${sv_id}"
