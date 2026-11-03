#!/bin/bash
# =============================================================================
# Run the enrichment pipeline in order.
#
# Usage:
#   cp config.template.sh config.sh   # fill in paths
#   bash run_pipeline.sh config.sh
#
# Steps 01 and 02 are skipped when preprocessed inputs are provided through
# GWAS_PAIRS_FROZEN and SV5KB_BACKGROUND_FROZEN, respectively.
# =============================================================================
set -euo pipefail

CONFIG="${1:-config.sh}"
# shellcheck source=/dev/null
source "${CONFIG}"

: "${OUTPUT_DIR:?OUTPUT_DIR must be set in the configuration file}"
N_PERM="${N_PERM:-100000}"
RANDOM_SEED="${RANDOM_SEED:-42}"
ALLOW_NETWORK="${ALLOW_NETWORK:-no}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "${OUTPUT_DIR}"

PAIRS="${GWAS_PAIRS_FROZEN:-${OUTPUT_DIR}/gwas/phenotype_gene_pairs.tsv}"
SV_BACKGROUND="${SV5KB_BACKGROUND_FROZEN:-${OUTPUT_DIR}/sv/SV5kb_background_GENCODEv38_proteinCoding.txt}"

# Validate that all required modules are present
for module in 01_prepare_gwas.py 02_prepare_sv_background.py \
              03_prepare_gene_sets.py 04_enrichment.py 05_union_modelB.py; do
    if [ ! -f "${HERE}/${module}" ]; then
        echo "Error: missing module ${HERE}/${module}" >&2
        exit 1
    fi
done

echo "==== [01] Preparing phenotype-gene associations ===="
if [ -z "${GWAS_PAIRS_FROZEN}" ]; then
    mkdir -p "${OUTPUT_DIR}/gwas"
    python3 "${HERE}/01_prepare_gwas.py" \
        --gwas "${GWAS_INPUT}" \
        --uri-labels "${GWAS_ONTOLOGY_MAP}" \
        --background "${GENCODE_PC_GENES}" \
        --out "${PAIRS}"
else
    echo "Using preprocessed phenotype-gene pairs: ${PAIRS}"
fi

echo "==== [02] Building SV-associated gene background ===="
if [ -z "${SV5KB_BACKGROUND_FROZEN}" ]; then
    mkdir -p "${OUTPUT_DIR}/sv"
    if [ -n "${SV5KB_INTERVALS}" ]; then
        python3 "${HERE}/02_prepare_sv_background.py" \
            --sv-regions "${SV_INPUT}" \
            --chm13-genes "${CHM13_GENE_ANNOTATION}" \
            --gencode-pc "${GENCODE_PC_GENES}" \
            --gencode-id-map "${GENCODE_GENE_ID_MAP}" \
            --outdir "${OUTPUT_DIR}/sv" \
            --sv-bed "${SV5KB_INTERVALS}" \
            --allow-network "${ALLOW_NETWORK}"
    else
        python3 "${HERE}/02_prepare_sv_background.py" \
            --sv-regions "${SV_INPUT}" \
            --chm13-genes "${CHM13_GENE_ANNOTATION}" \
            --gencode-pc "${GENCODE_PC_GENES}" \
            --gencode-id-map "${GENCODE_GENE_ID_MAP}" \
            --outdir "${OUTPUT_DIR}/sv" \
            --allow-network "${ALLOW_NETWORK}"
    fi
else
    echo "Using preprocessed SV background: ${SV_BACKGROUND}"
fi

echo "==== [03] Preparing candidate gene sets ===="
python3 "${HERE}/03_prepare_gene_sets.py" \
    --selection-input "${SELECTION_GENE_INPUT}" \
    --tesv-input "${TESV_GENE_INPUT}" \
    --hypoxia-input "${HYPOXIA_GENE_INPUT}" \
    --gencode-pc "${GENCODE_PC_GENES}" \
    --gencode-id-map "${GENCODE_GENE_ID_MAP}" \
    --sv-background "${SV_BACKGROUND}" \
    --outdir "${OUTPUT_DIR}/gene_sets" \
    --allow-network "${ALLOW_NETWORK}"
GENE_SETS="${OUTPUT_DIR}/gene_sets"

echo "==== [04a] GWAS enrichment: selection genes vs genome-wide background ===="
python3 "${HERE}/04_enrichment.py" --mode gwas \
    --target "${GENE_SETS}/selection_148.txt" \
    --background "${GENCODE_PC_GENES}" \
    --pairs "${PAIRS}" \
    --n-perm "${N_PERM}" --seed "${RANDOM_SEED}" \
    --restrict-to-background no \
    --outdir "${OUTPUT_DIR}/selection_gwas"

echo "==== [04b] Hypoxia enrichment: selection genes ===="
python3 "${HERE}/04_enrichment.py" --mode hypoxia \
    --target "${GENE_SETS}/selection_148.txt" \
    --background "${GENCODE_PC_GENES}" \
    --hypoxia "${GENE_SETS}/hypoxia_harmonized.txt" \
    --n-perm "${N_PERM}" --seed "${RANDOM_SEED}" \
    --outdir "${OUTPUT_DIR}/selection_hypoxia"

echo "==== [04c] GWAS enrichment: TE-SV genes vs SV-matched background ===="
python3 "${HERE}/04_enrichment.py" --mode gwas \
    --target "${GENE_SETS}/tesv_in_SV5kb_background.txt" \
    --background "${SV_BACKGROUND}" \
    --pairs "${PAIRS}" \
    --n-perm "${N_PERM}" --seed "${RANDOM_SEED}" \
    --restrict-to-background yes \
    --outdir "${OUTPUT_DIR}/tesv_gwas"

echo "==== [04d] Hypoxia enrichment: TE-SV genes ===="
python3 "${HERE}/04_enrichment.py" --mode hypoxia \
    --target "${GENE_SETS}/tesv_in_SV5kb_background.txt" \
    --background "${SV_BACKGROUND}" \
    --hypoxia "${GENE_SETS}/hypoxia_harmonized.txt" \
    --n-perm "${N_PERM}" --seed "${RANDOM_SEED}" \
    --outdir "${OUTPUT_DIR}/tesv_hypoxia"

echo "==== [05] Union enrichment with stratified Model B permutations ===="
python3 "${HERE}/05_union_modelB.py" \
    --union "${GENE_SETS}/union297.txt" \
    --background-b "${SV_BACKGROUND}" \
    --genome "${GENCODE_PC_GENES}" \
    --pairs "${PAIRS}" \
    --hypoxia "${GENE_SETS}/hypoxia_harmonized.txt" \
    --n-perm "${N_PERM}" --seed "${RANDOM_SEED}" \
    --outdir "${OUTPUT_DIR}/union_modelB"

echo "==== Analysis completed: ${OUTPUT_DIR} ===="
