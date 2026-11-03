#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: bash validate_config.sh CONFIG.env" >&2
  exit 2
fi
CONFIG=$(readlink -f "$1")
# shellcheck disable=SC1090
source "$CONFIG"

required=(BASE ANALYSIS_TAG OUTPUT_DIR SOURCE_VCF SAMPLE_LABEL CHROM REGION FOCAL_POS CACHE_TAG)
for name in "${required[@]}"; do
  [[ -n "${!name:-}" ]] || { echo "ERROR: missing config variable $name" >&2; exit 2; }
done
[[ "$ANALYSIS_TAG" =~ ^[A-Za-z0-9_.-]+$ ]] || { echo "ERROR: invalid ANALYSIS_TAG" >&2; exit 2; }
[[ "$SAMPLE_LABEL" =~ ^[A-Za-z0-9_.-]+$ ]] || { echo "ERROR: invalid SAMPLE_LABEL" >&2; exit 2; }
[[ "$OUTPUT_DIR" == "$BASE/runs/"* ]] || {
  echo "ERROR: OUTPUT_DIR must be inside BASE/runs" >&2; exit 2;
}
[[ "$CHROM" =~ ^chr([1-9]|1[0-9]|2[0-2])$ ]] || {
  echo "ERROR: CHROM must be an autosome written chr1..chr22" >&2; exit 2;
}
CHROM_NUM=${CHROM#chr}
[[ "$REGION" =~ ^([^:]+):([0-9]+)-([0-9]+)$ ]] || {
  echo "ERROR: REGION must look like chr1:231395000-232000000" >&2; exit 2;
}
REGION_CHROM=${BASH_REMATCH[1]}
REGION_START=${BASH_REMATCH[2]}
REGION_END=${BASH_REMATCH[3]}
[[ "$REGION_CHROM" == "$CHROM" ]] || { echo "ERROR: REGION chromosome differs from CHROM" >&2; exit 2; }
(( REGION_START >= 1 && REGION_END >= REGION_START )) || { echo "ERROR: invalid REGION coordinates" >&2; exit 2; }
(( FOCAL_POS >= REGION_START && FOCAL_POS <= REGION_END )) || {
  echo "ERROR: FOCAL_POS is outside REGION" >&2; exit 2;
}

BCFTOOLS=${BCFTOOLS:-/work/home/songjiaxu/bin/bcftools}
RELATE_ROOT=${RELATE_ROOT:-/work/home/songjiaxu/tools/relate}
CLUES_ROOT=${CLUES_ROOT:-/work/home/songjiaxu/tools/CLUES2-Companion/CLUES2}
CONDA_ENV=${CONDA_ENV:-omicverse-new}
GENETIC_MAP=${GENETIC_MAP:-}
ANCESTOR_FASTA=${ANCESTOR_FASTA:-}
if [[ -z "$GENETIC_MAP" ]]; then
  : "${GENETIC_MAP_TEMPLATE:?set GENETIC_MAP or GENETIC_MAP_TEMPLATE}"
  GENETIC_MAP=${GENETIC_MAP_TEMPLATE//@/$CHROM_NUM}
fi
if [[ -z "$ANCESTOR_FASTA" ]]; then
  : "${ANCESTOR_TEMPLATE:?set ANCESTOR_FASTA or ANCESTOR_TEMPLATE}"
  ANCESTOR_FASTA=${ANCESTOR_TEMPLATE//@/$CHROM_NUM}
fi

for path in "$CONFIG" "$SOURCE_VCF" "$GENETIC_MAP" "$ANCESTOR_FASTA"; do
  [[ -s "$path" ]] || { echo "ERROR: missing/empty file: $path" >&2; exit 2; }
done
[[ -s "$SOURCE_VCF.tbi" || -s "$SOURCE_VCF.csi" ]] || {
  echo "ERROR: SOURCE_VCF requires a .tbi or .csi index" >&2; exit 2;
}
[[ -x "$BCFTOOLS" ]] || { echo "ERROR: BCFTOOLS not executable: $BCFTOOLS" >&2; exit 2; }
"$BCFTOOLS" plugin -l 2>/dev/null | grep -Fx 'fill-tags' >/dev/null || {
  echo "ERROR: bcftools fill-tags plugin is unavailable" >&2; exit 2;
}
for exe in Relate RelateFileFormats RelateExtract RelateTreeView; do
  [[ -x "$RELATE_ROOT/bin/$exe" ]] || { echo "ERROR: missing RELATE executable $exe" >&2; exit 2; }
done
[[ -s "$CLUES_ROOT/inference.py" ]] || { echo "ERROR: missing CLUES2 inference.py" >&2; exit 2; }
command -v conda >/dev/null || { echo "ERROR: conda command is unavailable" >&2; exit 2; }
conda run -n "$CONDA_ENV" python -c \
  'import Bio, matplotlib, numpy, pandas, scipy; print("PYTHON_ENV_OK")' >/dev/null || {
  echo "ERROR: required Python packages are unavailable in conda env $CONDA_ENV" >&2; exit 2;
}
for py in deduplicate_haps_positions.py RelateToCLUES.fixed.py plot_clues2_trajectory.py plot_focal_tree.py validate_results.py; do
  [[ -s "$BASE/scripts/$py" ]] || { echo "ERROR: missing pipeline helper $py" >&2; exit 2; }
done
if [[ -n "${SAMPLES_FILE:-}" ]]; then
  [[ -s "$SAMPLES_FILE" ]] || { echo "ERROR: missing SAMPLES_FILE: $SAMPLES_FILE" >&2; exit 2; }
  missing=$(
    comm -23 <(awk 'NF{print $1}' "$SAMPLES_FILE" | sort -u) \
      <("$BCFTOOLS" query -l "$SOURCE_VCF" | sort -u) | wc -l
  )
  (( missing == 0 )) || { echo "ERROR: SAMPLES_FILE contains names absent from VCF" >&2; exit 2; }
fi

records=$(
  "$BCFTOOLS" query -r "$CHROM:$FOCAL_POS-$FOCAL_POS" \
    -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' "$SOURCE_VCF"
)
[[ -n "$records" ]] || { echo "ERROR: focal position absent from SOURCE_VCF" >&2; exit 2; }
filtered=$(printf '%s\n' "$records" | awk -F '\t' \
  -v id="${EXPECTED_FOCAL_ID:-}" -v ref="${EXPECTED_FOCAL_REF:-}" -v alt="${EXPECTED_FOCAL_ALT:-}" \
  '(!id || $3==id) && (!ref || $4==ref) && (!alt || $5==alt)')
[[ -n "$filtered" ]] || {
  echo "ERROR: focal records exist, but none matches EXPECTED_FOCAL_ID/REF/ALT" >&2
  printf '%s\n' "$records" >&2
  exit 2
}
if [[ $(printf '%s\n' "$filtered" | wc -l) -gt 1 ]]; then
  echo "ERROR: multiple focal records remain; set EXPECTED_FOCAL_ID/REF/ALT explicitly" >&2
  printf '%s\n' "$filtered" >&2
  exit 2
fi

sample_count=$("$BCFTOOLS" query -l "$SOURCE_VCF" | wc -l)
if [[ -n "${SAMPLES_FILE:-}" ]]; then
  sample_count=$(awk 'NF{print $1}' "$SAMPLES_FILE" | sort -u | wc -l)
fi
(( sample_count > 1 )) || { echo "ERROR: at least two diploid samples are required" >&2; exit 2; }
windows=${CLUES_WINDOWS:-"80kyr:2858:80000 160kyr:5715:160000"}
for spec in $windows; do
  [[ "$spec" =~ ^[A-Za-z0-9_.-]+:[1-9][0-9]*:[1-9][0-9]*$ ]] || {
    echo "ERROR: invalid CLUES_WINDOWS item: $spec" >&2; exit 2;
  }
done

printf 'CONFIG_OK\t%s\n' "$CONFIG"
printf 'analysis\t%s\nchromosome\t%s\nregion\t%s\nfocal\t%s:%s\n' \
  "$ANALYSIS_TAG" "$CHROM" "$REGION" "$CHROM" "$FOCAL_POS"
printf 'samples\t%s\nmap\t%s\nancestor\t%s\n' "$sample_count" "$GENETIC_MAP" "$ANCESTOR_FASTA"
printf 'focal_source_record\t%s\n' "$filtered"
