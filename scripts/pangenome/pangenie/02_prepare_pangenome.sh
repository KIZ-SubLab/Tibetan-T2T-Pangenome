#!/usr/bin/env bash
# Generate pangenome.vcf with the official PanGenie run-from-callset pipeline.
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${HERE}/config.sh"

IN="${WORKDIR}/input/panel.sorted.vcf"
OUTDIR="${WORKDIR}/results"
SNAKEFILE="${PANGENIE_REPO}/pipelines/run-from-callset/Snakefile"
CONFIG="${WORKDIR}/config.prepare.yaml"

[[ -s "${IN}" ]] || { echo "missing ${IN}; run 01_prepare_panel_vcf.sh first"; exit 1; }
[[ -s "${SNAKEFILE}" ]] || { echo "missing Snakefile: ${SNAKEFILE}"; exit 1; }
command -v snakemake >/dev/null
mkdir -p "${OUTDIR}" "${WORKDIR}/logs"

cat > "${CONFIG}" <<EOF
vcf: "${IN}"
reference: "${REFERENCE}"
reads: {}
pangenie_genotype: PanGenie
pangenie_index: PanGenie-index
outdir: "${OUTDIR}"
EOF

CONDA_OPT=""
if [[ "${USE_CONDA:-0}" == "1" ]]; then
    CONDA_OPT="--use-conda"
fi

# shellcheck disable=SC2086
snakemake \
    --snakefile "${SNAKEFILE}" \
    --configfile "${CONFIG}" \
    --cores "${THREADS}" \
    ${CONDA_OPT} \
    --rerun-incomplete \
    --printshellcmds \
    "${OUTDIR}/pangenome/pangenome.vcf"

echo "output: ${OUTDIR}/pangenome/pangenome.vcf"
