#!/usr/bin/env bash

# One entry point for Denisovan, Altai Neanderthal and Ghost SV heatmap batches.
# It validates arbitrary header-based case tables, optionally archives the
# previous enabled result trees, creates reproducible manifests, and submits
# base -> hetero200/direct_phased200 -> CEU35 -> summary Slurm dependencies.

set -eo pipefail

DEFAULT_ROOT="/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/output/rm1kb-na-segments/merged_fragments_v3/SV/denove_genotyping_SV_35tbn_5archic_dna/bam_archic/introgression_sv_heatmap/heamap"
DEFAULT_HEATBASE="/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/output/rm1kb-na-segments/merged_fragments_v3/SV/denove_genotyping_SV_35tbn_5archic_dna/bam_archic/introgression_sv_heatmap"
DEN_VCF_ROOT="/work/home/songjiaxu/introgression/fastq/den_new/snpad/all_chr_compressed"
NEA_VCF_ROOT="/work/home/songjiaxu/introgression/fastq/nea/Altai/snpad/all_chr_compressed"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HEATMAP_ROOT="${DEFAULT_ROOT}"
HEATBASE="${DEFAULT_HEATBASE}"
DEN_LIST=""
NEA_LIST=""
GHOST_LIST=""
ARCHIVE_EXISTING=0
ARCHIVE_LABEL=""
PREPARE_ONLY=0
INCLUDE_CEU=1
CEU_SEED=20260804
PANEL_N=200
WINDOW_BP=200000
MAX_WINDOW_BP=1600000
PARTITION=wzhcnormal
MAX_CONCURRENT=12
CPUS_PER_TASK=4
BASE_TIME=00:30:00
FOLLOWUP_TIME=00:30:00
CEU_TIME=01:00:00
SKIP_SUCCESS=0
CARRIER_COHERENT_IDS=""

usage() {
  cat <<'EOF'
Usage:
  bash run_all_sv_heatmaps.sh [options]

Required in normal use (defaults point to the current 326 tables):
  --den-list FILE              Denisovan case TSV; use 'none' to skip
  --nea-list FILE              Neanderthal case TSV; use 'none' to skip
  --ghost-list FILE            Ghost case TSV; use 'none' to skip

Output and safety:
  --heatmap-root DIR           Parent containing den/, nea/, ghost/, old/, pipeline/
  --archive-existing           Move enabled result directories into old/<label>/
  --archive-label LABEL        Stable archive name (default: pre_<run-id>)
  --resume                     Keep existing outputs and skip successful stages
  --prepare-only               Validate, archive/setup, but do not submit Slurm

Analysis controls:
  --skip-ceu                   Do not add CEU35 to Neanderthal/Ghost outputs
  --ceu-seed INT               Fixed CEU sampling seed (default 20260804)
  --panel-n INT                hetero/direct target SNP count (default 200)
  --window-bp INT              Initial bp on each side (default 200000)
  --max-window-bp INT          Ghost adaptive-window ceiling (default 1600000)
  --carrier-coherent-id ID     Explicitly enable result-informed special panel;
                               may be repeated (never enabled by default)

Slurm:
  --partition NAME             Default wzhcnormal
  --max-concurrent INT         Default 12 per population/stage
  --cpus INT                   Default 4

Examples:
  bash pipeline/run_all_sv_heatmaps.sh --archive-existing
  bash pipeline/run_all_sv_heatmaps.sh --den-list new.den.tsv \
    --nea-list new.nea.tsv --archive-existing --archive-label pre_new_cases
  bash pipeline/run_all_sv_heatmaps.sh --den-list none \
    --nea-list new.nea.tsv --resume
  bash pipeline/run_all_sv_heatmaps.sh --den-list none --nea-list none \
    --ghost-list new.ghost.tsv
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --den-list) DEN_LIST="$2"; shift 2 ;;
    --nea-list) NEA_LIST="$2"; shift 2 ;;
    --ghost-list) GHOST_LIST="$2"; shift 2 ;;
    --heatmap-root) HEATMAP_ROOT="$2"; shift 2 ;;
    --heatbase) HEATBASE="$2"; shift 2 ;;
    --archive-existing) ARCHIVE_EXISTING=1; shift ;;
    --archive-label) ARCHIVE_LABEL="$2"; shift 2 ;;
    --prepare-only) PREPARE_ONLY=1; shift ;;
    --skip-ceu) INCLUDE_CEU=0; shift ;;
    --ceu-seed) CEU_SEED="$2"; shift 2 ;;
    --panel-n) PANEL_N="$2"; shift 2 ;;
    --window-bp) WINDOW_BP="$2"; shift 2 ;;
    --max-window-bp) MAX_WINDOW_BP="$2"; shift 2 ;;
    --partition) PARTITION="$2"; shift 2 ;;
    --max-concurrent) MAX_CONCURRENT="$2"; shift 2 ;;
    --cpus) CPUS_PER_TASK="$2"; shift 2 ;;
    --resume) SKIP_SUCCESS=1; shift ;;
    --carrier-coherent-id)
      CARRIER_COHERENT_IDS="${CARRIER_COHERENT_IDS:+${CARRIER_COHERENT_IDS},}$2"
      shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) printf 'ERROR: unknown option: %s\n' "$1" >&2; usage >&2; exit 2 ;;
  esac
done

DEN_LIST="${DEN_LIST:-${HEATMAP_ROOT}/326.den_typing.tsv}"
NEA_LIST="${NEA_LIST:-${HEATMAP_ROOT}/326.nea_typing.tsv}"
GHOST_LIST="${GHOST_LIST:-${HEATMAP_ROOT}/326.ghost_typing.tsv}"
RUN_ID="$(date '+%Y%m%d_%H%M%S')"
ARCHIVE_LABEL="${ARCHIVE_LABEL:-pre_${RUN_ID}}"

for numeric in CEU_SEED PANEL_N WINDOW_BP MAX_WINDOW_BP MAX_CONCURRENT CPUS_PER_TASK; do
  value="${!numeric}"
  if [[ ! "${value}" =~ ^[1-9][0-9]*$ ]]; then
    printf 'ERROR: %s must be a positive integer, got %s\n' "${numeric}" "${value}" >&2
    exit 2
  fi
done
if [[ "${ARCHIVE_LABEL}" == */* || "${ARCHIVE_LABEL}" == "." || "${ARCHIVE_LABEL}" == ".." ]]; then
  printf 'ERROR: archive label must be one path component\n' >&2
  exit 2
fi

HEATMAP_ROOT="$(realpath -m -- "${HEATMAP_ROOT}")"
HEATBASE="$(realpath -m -- "${HEATBASE}")"
if [[ ! -d "${HEATMAP_ROOT}" || ! -d "${HEATBASE}" ]]; then
  printf 'ERROR: heatmap root or heatbase does not exist\n' >&2
  exit 1
fi
if [[ ! -r "${SCRIPT_DIR}/build_manifest.R" || ! -d "${SCRIPT_DIR}/modules" ]]; then
  printf 'ERROR: incomplete pipeline installation: %s\n' "${SCRIPT_DIR}" >&2
  exit 1
fi

active_lists=0
active_modes=()
for spec in "den:${DEN_LIST}" "nea:${NEA_LIST}" "ghost:${GHOST_LIST}"; do
  mode="${spec%%:*}"; list="${spec#*:}"
  [[ "${list}" == "none" ]] && continue
  active_lists=$((active_lists + 1))
  active_modes+=("${mode}")
  if [[ ! -r "${list}" ]]; then
    printf 'ERROR: cannot read %s case table: %s\n' "${mode}" "${list}" >&2
    exit 1
  fi
done
if (( active_lists == 0 )); then
  printf 'ERROR: all case tables are disabled\n' >&2
  exit 1
fi

source /work/home/songjiaxu/anaconda3/etc/profile.d/conda.sh
conda activate plink
set -u

for tool in Rscript bcftools plink plink2 SHAPEIT5_phase_common sbatch sha256sum; do
  command -v "${tool}" >/dev/null || {
    printf 'ERROR: required command unavailable: %s\n' "${tool}" >&2
    exit 1
  }
done
Rscript -e 'pkgs <- c("data.table", "pheatmap", "grid", "ragg"); stopifnot(all(vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)))'

if [[ "${ARCHIVE_EXISTING}" == "1" ]]; then
  # Moving these exact children is fast and recoverable because the archive is
  # on the same filesystem. Input case tables and pipeline/ remain in place.
  archive_dir="${HEATMAP_ROOT}/old/${ARCHIVE_LABEL}"
  if [[ -e "${archive_dir}" ]]; then
    printf 'ERROR: archive destination already exists: %s\n' "${archive_dir}" >&2
    exit 1
  fi
  mkdir -p "${archive_dir}"
  printf 'item\ttype\tbytes\tfile_count\toriginal_path\n' > "${archive_dir}/ARCHIVE_ITEMS.tsv"
  archive_count=0
  archive_names=("${active_modes[@]}")
  for legacy_doc in DEN_NEA_HAPLOTYPE_HEATMAP_WORKFLOW_GUIDE.md readme; do
    [[ -e "${HEATMAP_ROOT}/${legacy_doc}" ]] && archive_names+=("${legacy_doc}")
  done
  for name in "${archive_names[@]}"; do
    src="${HEATMAP_ROOT}/${name}"
    [[ -e "${src}" ]] || continue
    resolved="$(realpath -m -- "${src}")"
    case "${resolved}" in
      "${HEATMAP_ROOT}"/*) ;;
      *) printf 'ERROR: unsafe archive target: %s\n' "${resolved}" >&2; exit 1 ;;
    esac
    bytes="$(du -sb -- "${src}" | awk '{print $1}')"
    if [[ -d "${src}" ]]; then
      kind=directory
      files="$(find "${src}" -type f | wc -l)"
    else
      kind=file
      files=1
    fi
    printf '%s\t%s\t%s\t%s\t%s\n' "${name}" "${kind}" "${bytes}" "${files}" "${src}" \
      >> "${archive_dir}/ARCHIVE_ITEMS.tsv"
    mv -- "${src}" "${archive_dir}/${name}"
    archive_count=$((archive_count + 1))
  done
  {
    printf 'key\tvalue\n'
    printf 'archive_label\t%s\n' "${ARCHIVE_LABEL}"
    printf 'archived_at\t%s\n' "$(date -Is)"
    printf 'archive_path\t%s\n' "${archive_dir}"
    printf 'archived_items\t%s\n' "${archive_count}"
  } > "${archive_dir}/ARCHIVE_COMPLETE.tsv"
  printf 'Archived %d old top-level result items into %s\n' "${archive_count}" "${archive_dir}"
elif [[ "${SKIP_SUCCESS}" != "1" ]]; then
  for mode in "${active_modes[@]}"; do
    if [[ -e "${HEATMAP_ROOT}/${mode}" ]]; then
      printf 'ERROR: enabled result directory %s/ already exists. Use --archive-existing or --resume.\n' "${mode}" >&2
      exit 1
    fi
  done
fi

RUN_DIR="${HEATMAP_ROOT}/pipeline_runs/${RUN_ID}"
mkdir -p "${RUN_DIR}"
printf 'mode\tstage\tjob_id\tdependency\n' > "${RUN_DIR}/jobs.tsv"

install_mode() {
  local mode="$1" list="$2" label="$3" code="$4" vcf_root="$5" sample_name="$6" ceu="$7"
  local analysis_root="${HEATMAP_ROOT}/${mode}"
  mkdir -p "${analysis_root}/input" "${analysis_root}/scripts" \
    "${analysis_root}/logs/base" "${analysis_root}/logs/hetero200" \
    "${analysis_root}/logs/direct_phased200" "${analysis_root}/logs/summary"

  cases_copy="${analysis_root}/input/cases.tsv"
  if [[ "$(realpath -m -- "${list}")" != "$(realpath -m -- "${cases_copy}")" ]]; then
    cp -f -- "${list}" "${cases_copy}"
  fi
  cp -f -- "${SCRIPT_DIR}/modules/"* "${analysis_root}/scripts/"
  chmod +x "${analysis_root}/scripts/"*.sh
  Rscript "${SCRIPT_DIR}/build_manifest.R" \
    "${analysis_root}/input/cases.tsv" "${analysis_root}/sv_manifest.tsv" \
    "${analysis_root}/input/source_map.tsv" "${mode}"

  while IFS=$'\t' read -r task_id chr start end type mid sv_id; do
    [[ "${task_id}" == task_id ]] && continue
    mkdir -p "${analysis_root}/${sv_id}"
  done < "${analysis_root}/sv_manifest.tsv"

  if [[ "${ceu}" == "1" ]]; then
    mkdir -p "${analysis_root}/scripts_ceu35" "${analysis_root}/logs/ceu35"
    cp -f -- "${SCRIPT_DIR}/ceu35/"* "${analysis_root}/scripts_ceu35/"
    chmod +x "${analysis_root}/scripts_ceu35/"*.sh
    NEA_ROOT="${analysis_root}" CEU_ROOT="${HEATBASE}/1KGP/CEU" \
      CEU_SEED="${CEU_SEED}" bash "${analysis_root}/scripts_ceu35/prepare_ceu35.sh"
  fi

  local n_tasks input_hash
  n_tasks="$(awk 'END{print NR-1}' "${analysis_root}/sv_manifest.tsv")"
  input_hash="$(sha256sum "${analysis_root}/input/cases.tsv" | awk '{print $1}')"
  {
    printf 'key\tvalue\n'
    printf 'run_id\t%s\n' "${RUN_ID}"
    printf 'mode\t%s\n' "${mode}"
    printf 'archaic_label\t%s\n' "${label}"
    printf 'source_table\t%s\n' "$(realpath -m -- "${list}")"
    printf 'source_sha256\t%s\n' "${input_hash}"
    printf 'n_cases\t%s\n' "${n_tasks}"
    printf 'target_panel_n\t%s\n' "${PANEL_N}"
    printf 'window_bp_each_side\t%s\n' "${WINDOW_BP}"
    printf 'include_ceu35\t%s\n' "${ceu}"
    printf 'ceu_seed\t%s\n' "${CEU_SEED}"
    printf 'carrier_coherent_ids\t%s\n' "${CARRIER_COHERENT_IDS:-none}"
    printf 'created_at\t%s\n' "$(date -Is)"
  } > "${analysis_root}/RUN_METADATA.tsv"

  if [[ "${PREPARE_ONLY}" == "1" ]]; then
    printf 'Prepared %s (%s cases); Slurm submission skipped.\n' "${mode}" "${n_tasks}"
    return
  fi

  local common_export base_raw base_job het_raw het_job direct_raw direct_job
  local ceu_raw ceu_job summary_raw summary_job deps
  common_export="ALL,ANALYSIS_ROOT=${analysis_root},ARCHAIC_VCF_ROOT=${vcf_root},ARCHAIC_LABEL=${label},ARCHAIC_CODE=${code},ARCHAIC_SAMPLE_NAME=${sample_name},HEATBASE=${HEATBASE},WINDOW_BP=${WINDOW_BP},PANEL_N=${PANEL_N},SKIP_SUCCESS=${SKIP_SUCCESS},CARRIER_COHERENT_IDS=${CARRIER_COHERENT_IDS}"

  base_raw="$(sbatch --parsable --partition="${PARTITION}" \
    --job-name="${code}_base" --cpus-per-task="${CPUS_PER_TASK}" \
    --mem-per-cpu=3500M --time="${BASE_TIME}" \
    --array="1-${n_tasks}%${MAX_CONCURRENT}" \
    --output="${analysis_root}/logs/base/%A_%a.out" \
    --error="${analysis_root}/logs/base/%A_%a.err" \
    --export="${common_export}" "${analysis_root}/scripts/run_array.sh")"
  base_job="${base_raw%%;*}"
  printf '%s\tbase\t%s\tnone\n' "${mode}" "${base_job}" >> "${RUN_DIR}/jobs.tsv"

  het_raw="$(sbatch --parsable --partition="${PARTITION}" \
    --job-name="${code}_hetero200" --cpus-per-task=1 --mem-per-cpu=3500M \
    --time="${FOLLOWUP_TIME}" --dependency="afterok:${base_job}" \
    --array="1-${n_tasks}%${MAX_CONCURRENT}" \
    --output="${analysis_root}/logs/hetero200/%A_%a.out" \
    --error="${analysis_root}/logs/hetero200/%A_%a.err" \
    --export="${common_export}" "${analysis_root}/scripts/run_hetero200_array.sh")"
  het_job="${het_raw%%;*}"
  printf '%s\thetero200\t%s\tafterok:%s\n' "${mode}" "${het_job}" "${base_job}" >> "${RUN_DIR}/jobs.tsv"

  direct_raw="$(sbatch --parsable --partition="${PARTITION}" \
    --job-name="${code}_direct200" --cpus-per-task="${CPUS_PER_TASK}" \
    --mem-per-cpu=3500M --time="${FOLLOWUP_TIME}" \
    --dependency="afterok:${base_job}" --array="1-${n_tasks}%${MAX_CONCURRENT}" \
    --output="${analysis_root}/logs/direct_phased200/%A_%a.out" \
    --error="${analysis_root}/logs/direct_phased200/%A_%a.err" \
    --export="${common_export}" "${analysis_root}/scripts/run_direct_phased200_array.sh")"
  direct_job="${direct_raw%%;*}"
  printf '%s\tdirect_phased200\t%s\tafterok:%s\n' "${mode}" "${direct_job}" "${base_job}" >> "${RUN_DIR}/jobs.tsv"

  deps="${base_job}:${het_job}:${direct_job}"
  ceu_job=""
  if [[ "${ceu}" == "1" ]]; then
    ceu_raw="$(sbatch --parsable --partition="${PARTITION}" \
      --job-name="${code}_ceu35" --cpus-per-task="${CPUS_PER_TASK}" \
      --mem-per-cpu=3500M --time="${CEU_TIME}" \
      --dependency="afterok:${deps}" --array="1-${n_tasks}%${MAX_CONCURRENT}" \
      --output="${analysis_root}/logs/ceu35/%A_%a.out" \
      --error="${analysis_root}/logs/ceu35/%A_%a.err" \
      --export="ALL,NEA_ROOT=${analysis_root},HEATBASE=${HEATBASE},CEU_ROOT=${HEATBASE}/1KGP/CEU,CEU_SEED=${CEU_SEED},WINDOW_BP=${WINDOW_BP},SKIP_SUCCESS=${SKIP_SUCCESS},CARRIER_COHERENT_IDS=${CARRIER_COHERENT_IDS}" \
      "${analysis_root}/scripts_ceu35/run_ceu35_array.sh")"
    ceu_job="${ceu_raw%%;*}"
    printf '%s\tceu35\t%s\tafterok:%s\n' "${mode}" "${ceu_job}" "${deps}" >> "${RUN_DIR}/jobs.tsv"
    deps="${deps}:${ceu_job}"
  fi

  summary_raw="$(sbatch --parsable --partition="${PARTITION}" \
    --job-name="${code}_final_summary" --cpus-per-task=1 --mem-per-cpu=2000M \
    --time=00:15:00 --dependency="afterany:${deps}" \
    --output="${analysis_root}/logs/summary/%j.out" \
    --error="${analysis_root}/logs/summary/%j.err" \
    --export="ALL,ANALYSIS_ROOT=${analysis_root},ARCHAIC_LABEL=${label},ARCHAIC_CODE=${code},INCLUDE_CEU=${ceu}" \
    "${analysis_root}/scripts/run_all_summaries.sh")"
  summary_job="${summary_raw%%;*}"
  printf '%s\tfinal_summary\t%s\tafterany:%s\n' "${mode}" "${summary_job}" "${deps}" >> "${RUN_DIR}/jobs.tsv"
  printf 'Submitted %s: cases=%s base=%s hetero=%s direct=%s ceu=%s summary=%s\n' \
    "${mode}" "${n_tasks}" "${base_job}" "${het_job}" "${direct_job}" \
    "${ceu_job:-not_requested}" "${summary_job}"
}

install_ghost() {
  local list="$1"
  local analysis_root="${HEATMAP_ROOT}/ghost"
  mkdir -p "${analysis_root}/input" "${analysis_root}/scripts" \
    "${analysis_root}/logs/base" "${analysis_root}/logs/hetero200" \
    "${analysis_root}/logs/direct_phased200" "${analysis_root}/logs/ceu35" \
    "${analysis_root}/logs/summary"

  local cases_copy="${analysis_root}/input/cases.tsv"
  if [[ "$(realpath -m -- "${list}")" != "$(realpath -m -- "${cases_copy}")" ]]; then
    cp -f -- "${list}" "${cases_copy}"
  fi
  cp -f -- "${SCRIPT_DIR}/ghost/"* "${analysis_root}/scripts/"
  chmod +x "${analysis_root}/scripts/"*.sh
  Rscript "${SCRIPT_DIR}/build_manifest.R" \
    "${cases_copy}" "${analysis_root}/sv_manifest.tsv" \
    "${analysis_root}/input/source_map.tsv" ghost

  while IFS=$'\t' read -r task_id chr start end type mid sv_id; do
    [[ "${task_id}" == task_id ]] && continue
    mkdir -p "${analysis_root}/${sv_id}"
  done < "${analysis_root}/sv_manifest.tsv"

  if [[ "${INCLUDE_CEU}" == "1" ]]; then
    GHOST_ROOT="${analysis_root}" CEU_ROOT="${HEATBASE}/1KGP/CEU" \
      CEU_SEED="${CEU_SEED}" bash "${analysis_root}/scripts/prepare_ghost_ceu35.sh"
  fi

  local n_tasks input_hash
  n_tasks="$(awk 'END{print NR-1}' "${analysis_root}/sv_manifest.tsv")"
  input_hash="$(sha256sum "${cases_copy}" | awk '{print $1}')"
  {
    printf 'key\tvalue\n'
    printf 'run_id\t%s\n' "${RUN_ID}"
    printf 'mode\tghost\n'
    printf 'source_table\t%s\n' "$(realpath -m -- "${list}")"
    printf 'source_sha256\t%s\n' "${input_hash}"
    printf 'n_cases\t%s\n' "${n_tasks}"
    printf 'target_panel_n\t%s\n' "${PANEL_N}"
    printf 'initial_window_bp_each_side\t%s\n' "${WINDOW_BP}"
    printf 'max_window_bp_each_side\t%s\n' "${MAX_WINDOW_BP}"
    printf 'include_ceu35\t%s\n' "${INCLUDE_CEU}"
    printf 'ceu_seed\t%s\n' "${CEU_SEED}"
    printf 'selection_objective\tTBN_CHB_separation_only\n'
    printf 'ancient_samples\tDenisovan,Altai_Neanderthal\n'
    printf 'ancient_role\tdisplay_only\n'
    printf 'created_at\t%s\n' "$(date -Is)"
  } > "${analysis_root}/RUN_METADATA.tsv"

  if [[ "${PREPARE_ONLY}" == "1" ]]; then
    printf 'Prepared ghost (%s cases); Slurm submission skipped.\n' "${n_tasks}"
    return
  fi

  local common_export base_raw base_job het_raw het_job direct_raw direct_job
  local ceu_raw ceu_job summary_raw summary_job deps
  common_export="ALL,ANALYSIS_ROOT=${analysis_root},HEATBASE=${HEATBASE},DEN_VCF_ROOT=${DEN_VCF_ROOT},NEA_VCF_ROOT=${NEA_VCF_ROOT},WINDOW_BP=${WINDOW_BP},MAX_WINDOW_BP=${MAX_WINDOW_BP},PANEL_N=${PANEL_N},SKIP_SUCCESS=${SKIP_SUCCESS}"

  base_raw="$(sbatch --parsable --partition="${PARTITION}" \
    --job-name=ghost_base --cpus-per-task="${CPUS_PER_TASK}" \
    --mem-per-cpu=3500M --time=01:00:00 \
    --array="1-${n_tasks}%${MAX_CONCURRENT}" \
    --output="${analysis_root}/logs/base/%A_%a.out" \
    --error="${analysis_root}/logs/base/%A_%a.err" \
    --export="${common_export}" "${analysis_root}/scripts/run_ghost_array.sh")"
  base_job="${base_raw%%;*}"
  printf 'ghost\tbase\t%s\tnone\n' "${base_job}" >> "${RUN_DIR}/jobs.tsv"

  het_raw="$(sbatch --parsable --partition="${PARTITION}" \
    --job-name=ghost_hetero200 --cpus-per-task=1 --mem-per-cpu=3500M \
    --time=00:30:00 --dependency="afterok:${base_job}" \
    --array="1-${n_tasks}%${MAX_CONCURRENT}" \
    --output="${analysis_root}/logs/hetero200/%A_%a.out" \
    --error="${analysis_root}/logs/hetero200/%A_%a.err" \
    --export="${common_export}" "${analysis_root}/scripts/run_ghost_hetero_array.sh")"
  het_job="${het_raw%%;*}"
  printf 'ghost\thetero200\t%s\tafterok:%s\n' "${het_job}" "${base_job}" >> "${RUN_DIR}/jobs.tsv"

  direct_raw="$(sbatch --parsable --partition="${PARTITION}" \
    --job-name=ghost_direct200 --cpus-per-task="${CPUS_PER_TASK}" \
    --mem-per-cpu=3500M --time=00:45:00 --dependency="afterok:${base_job}" \
    --array="1-${n_tasks}%${MAX_CONCURRENT}" \
    --output="${analysis_root}/logs/direct_phased200/%A_%a.out" \
    --error="${analysis_root}/logs/direct_phased200/%A_%a.err" \
    --export="${common_export}" "${analysis_root}/scripts/run_ghost_direct_array.sh")"
  direct_job="${direct_raw%%;*}"
  printf 'ghost\tdirect_phased200\t%s\tafterok:%s\n' "${direct_job}" "${base_job}" >> "${RUN_DIR}/jobs.tsv"

  deps="${base_job}:${het_job}:${direct_job}"
  ceu_job=""
  if [[ "${INCLUDE_CEU}" == "1" ]]; then
    ceu_raw="$(sbatch --parsable --partition="${PARTITION}" \
      --job-name=ghost_ceu35 --cpus-per-task="${CPUS_PER_TASK}" \
      --mem-per-cpu=3500M --time=01:00:00 --dependency="afterok:${deps}" \
      --array="1-${n_tasks}%${MAX_CONCURRENT}" \
      --output="${analysis_root}/logs/ceu35/%A_%a.out" \
      --error="${analysis_root}/logs/ceu35/%A_%a.err" \
      --export="ALL,GHOST_ROOT=${analysis_root},HEATBASE=${HEATBASE},CEU_ROOT=${HEATBASE}/1KGP/CEU,CEU_SEED=${CEU_SEED},PANEL_N=${PANEL_N},SKIP_SUCCESS=${SKIP_SUCCESS}" \
      "${analysis_root}/scripts/run_ghost_ceu_array.sh")"
    ceu_job="${ceu_raw%%;*}"
    printf 'ghost\tceu35\t%s\tafterok:%s\n' "${ceu_job}" "${deps}" >> "${RUN_DIR}/jobs.tsv"
    deps="${deps}:${ceu_job}"
  fi

  summary_raw="$(sbatch --parsable --partition="${PARTITION}" \
    --job-name=ghost_final_summary --cpus-per-task=1 --mem-per-cpu=2000M \
    --time=00:20:00 --dependency="afterany:${deps}" \
    --output="${analysis_root}/logs/summary/%j.out" \
    --error="${analysis_root}/logs/summary/%j.err" \
    --export="ALL,GHOST_ROOT=${analysis_root}" \
    "${analysis_root}/scripts/run_ghost_summaries.sh")"
  summary_job="${summary_raw%%;*}"
  printf 'ghost\tfinal_summary\t%s\tafterany:%s\n' "${summary_job}" "${deps}" >> "${RUN_DIR}/jobs.tsv"
  printf 'Submitted ghost: cases=%s base=%s hetero=%s direct=%s ceu=%s summary=%s\n' \
    "${n_tasks}" "${base_job}" "${het_job}" "${direct_job}" \
    "${ceu_job:-not_requested}" "${summary_job}"
}

if [[ "${DEN_LIST}" != "none" ]]; then
  install_mode den "${DEN_LIST}" Denisovan den "${DEN_VCF_ROOT}" "" 0
fi
if [[ "${NEA_LIST}" != "none" ]]; then
  install_mode nea "${NEA_LIST}" Neanderthal nea "${NEA_VCF_ROOT}" Altai_Neanderthal "${INCLUDE_CEU}"
fi
if [[ "${GHOST_LIST}" != "none" ]]; then
  install_ghost "${GHOST_LIST}"
fi

{
  printf 'key\tvalue\n'
  printf 'run_id\t%s\n' "${RUN_ID}"
  printf 'prepared_at\t%s\n' "$(date -Is)"
  printf 'archive_existing\t%s\n' "${ARCHIVE_EXISTING}"
  printf 'archive_label\t%s\n' "${ARCHIVE_LABEL}"
  printf 'partition\t%s\n' "${PARTITION}"
  printf 'max_concurrent\t%s\n' "${MAX_CONCURRENT}"
  printf 'prepare_only\t%s\n' "${PREPARE_ONLY}"
} > "${RUN_DIR}/RUN.tsv"

printf 'Run metadata: %s\n' "${RUN_DIR}"
