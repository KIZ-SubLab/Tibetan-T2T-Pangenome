# =============================================================================
# Configuration template for the enrichment pipeline.
# Copy to config.sh and fill in paths for your system.
# All modules receive inputs through these variables (see README.md).
# =============================================================================

# ---- GWAS Catalog ----
# Raw GWAS Catalog "All associations" full-list TSV.
GWAS_INPUT=""
# Trait URI to label dictionary (TSV with header: TRAIT_URI<TAB>TRAIT_LABEL).
GWAS_ONTOLOGY_MAP=""

# ---- Annotation ----
# GENCODE v38 protein-coding gene symbols (one per line).
GENCODE_PC_GENES=""
# GENCODE v38 identifier lookup (TSV without header:
#   gene_id<TAB>gene_name<TAB>gene_type).
GENCODE_GENE_ID_MAP=""
# CHM13 v2.0 RefSeq Liftoff v5.2 gene-level annotation (TSV with columns
#   CHROM START_1BASED END_1BASED STRAND GENE_ID GENE_NAME GENE_BIOTYPE ...).
CHM13_GENE_ANNOTATION=""

# ---- Structural variants ----
# Merged SV +/- 5 kb intervals (1-based inclusive TSV with columns
#   CHROM START_1BASED END_1BASED SV_ID SVTYPE SV_POS SV_END).
SV_INPUT=""
# Optional precomputed SV +/- 5 kb intervals in BED format.
SV5KB_INTERVALS=""

# ---- Candidate gene lists ----
# One gene symbol per line.
SELECTION_GENE_INPUT=""
TESV_GENE_INPUT=""
HYPOXIA_GENE_INPUT=""

# ---- Optional preprocessed inputs ----
# Phenotype-gene pairs, e.g. the 1,459-phenotype set (skips step 01).
GWAS_PAIRS_FROZEN=""
# SV +/- 5 kb matched background gene list (skips step 02).
SV5KB_BACKGROUND_FROZEN=""

# ---- Run parameters ----
OUTPUT_DIR=""
N_PERM=100000
RANDOM_SEED=42
# "yes" enables Ensembl REST lookups for unambiguous gene symbol mapping.
ALLOW_NETWORK="no"
