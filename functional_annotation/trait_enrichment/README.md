# GWAS trait and hypoxia gene enrichment

## Overview

This workflow supports the functional enrichment analyses in Figure 4. It tests
candidate genes for enrichment among genes associated with GWAS traits and in a
curated hypoxia-related gene set. The analyses use SNP-based selection genes,
TESV candidate genes, and their union. Permutations use a gene background matched
to each analysis.

## Requirements

- Python 3 (with `pandas`, `numpy`, `scipy`)
- `bedtools` (used for interval intersection)
- Optional network access to the Ensembl REST API for gene identifier mapping
  (`ALLOW_NETWORK=yes`)

## Input files

### Primary inputs

| Variable | Description | Format |
|---|---|---|
| `GWAS_INPUT` | GWAS Catalog "All associations" full-list TSV | TSV with the standard 38-column GWAS Catalog header |
| `GWAS_ONTOLOGY_MAP` | Trait URI to trait label dictionary | TSV: `TRAIT_URI` `TRAIT_LABEL` (header required) |
| `GENCODE_PC_GENES` | GENCODE v38 protein-coding gene names | One gene symbol per line |
| `GENCODE_GENE_ID_MAP` | GENCODE v38 gene identifier lookup | TSV, no header: `gene_id` `gene_name` `gene_type` |
| `CHM13_GENE_ANNOTATION` | CHM13 v2.0 RefSeq Liftoff v5.2 gene-level annotation | TSV: `CHROM` `START_1BASED` `END_1BASED` `STRAND` `GENE_ID` `GENE_NAME` `GENE_BIOTYPE` |
| `SV_INPUT` | Merged SV ±5 kb intervals | TSV: `CHROM` `START_1BASED` `END_1BASED` `SV_ID` `SVTYPE` `SV_POS` `SV_END` |
| `SELECTION_GENE_INPUT` | Selection candidate gene symbols | One symbol per line |
| `TESV_GENE_INPUT` | TE–SV candidate gene symbols | One symbol per line |
| `HYPOXIA_GENE_INPUT` | Curated hypoxia-related gene symbols | One symbol per line |

### Preprocessed inputs (optional)

These may be supplied directly if the upstream preprocessing has already been
performed. When provided, the corresponding preparation step is skipped.

| Variable | Description |
|---|---|
| `GWAS_PAIRS_FROZEN` | Phenotype–gene pairs (`TRAIT_URI`, `TRAIT_LABEL`, `GENE`), e.g. the 1,459-phenotype set used in the analysis |
| `SV5KB_BACKGROUND_FROZEN` | Precomputed SV ±5 kb matched background gene list (e.g. 8,933 genes) |
| `SV5KB_INTERVALS` | Precomputed SV ±5 kb intervals in BED format |

## Analysis

**GWAS phenotype gene sets.** Associations with `INTERGENIC = 0` and a non-empty
`MAPPED_GENE` are retained. Associations are filtered at `PVALUE_MLOG >= 7.30103`
(P ≤ 5 × 10⁻⁸), falling back to `P-VALUE <= 5e-8` when `PVALUE_MLOG` is missing.
`MAPPED_TRAIT_URI` is expanded into one row per ontology term, and `MAPPED_GENE`
is split on commas. Gene symbols are matched exactly to the GENCODE v38
protein-coding universe; pairs are deduplicated on
(`TRAIT_URI`, `TRAIT_LABEL`, `GENE`) and phenotypes with fewer than 20 distinct
protein-coding genes are removed.

**SV-associated gene background.** Gene-level intervals from the CHM13 v2.0
RefSeq Liftoff v5.2 annotation are intersected with the merged SV ±5 kb intervals
(at least 1 bp overlap; `bedtools intersect -wa -wb`). Gene symbols are
harmonized to GENCODE v38 protein-coding symbols by exact match, by three
Ensembl-ID-verified symbol updates (ZFTA→C11orf95, SMIM45→LINC00634, PALS2→MPP6),
and, if enabled, by unambiguous mapping through a single stable Ensembl gene ID.
No fuzzy matching is performed.

**Candidate gene sets.** Selection and TE–SV symbols are harmonized under the
same rules. TE–SV genes outside the SV ±5 kb matched background are excluded from
the matched-background analyses while remaining eligible for union-based
analyses that use the genome-wide background. The union gene set is
Selection ∪ TE–SV (within the SV background).

**Enrichment.** Selection genes are tested against the genome-wide protein-coding
background; TE–SV genes are tested against the SV ±5 kb matched background, with
phenotype gene sets first restricted to that background and re-filtered at ≥20
genes. The union gene set is tested with the stratified Model B permutation
described below.

**Union Model B permutation.** Each permutation samples 167 genes without
replacement from the SV ±5 kb background, excludes them from the genome-wide
background, samples 130 additional genes, and combines both sets into 297 unique
genes. The same random gene set is used to compute overlaps with all GWAS
phenotypes and with the hypoxia gene set.

## Usage

```bash
cp config.template.sh config.sh   # fill in input and output paths
bash run_pipeline.sh config.sh
```

Individual steps can be run independently, for example:

```bash
python3 04_enrichment.py --mode gwas \
    --target gene_sets/selection_148.txt \
    --background /path/to/gencode_v38_protein_coding_genes.txt \
    --pairs gwas/phenotype_gene_pairs.tsv \
    --n-perm 100000 --seed 42 \
    --outdir results/selection_gwas
```

## Output files

| Output | Description |
|---|---|
| `GWAS_enrichment_full.tsv` | One row per tested phenotype: observed and expected overlap, fold enrichment, empirical P, BH-FDR, significance flag |
| `GWAS_significant_phenotypes.tsv` | Phenotypes with FDR < 0.05 and fold enrichment > 1 |
| `Hypoxia_enrichment_statistics.tsv` | Observed/expected overlap, fold enrichment, empirical P and overlapping genes for the hypoxia test |
| `Hypoxia_permutation_null_distribution.tsv` | Null overlap counts from the permutation procedure |
| `Hypoxia_overlap_genes.txt` | Genes shared between the target set and the hypoxia gene set |

## Statistical methods

- Background sampling is performed without replacement.
- Random numbers are generated with `numpy.random.default_rng(seed)`.
- `N_PERM = 100,000`; `seed = 42` by default.
- Expected overlap is the permutation mean.
- Fold enrichment = observed overlap / expected overlap.
- Empirical P = (1 + number of null overlaps ≥ observed) / (N_PERM + 1).
- GWAS phenotypes are corrected with the Benjamini–Hochberg procedure over all
  phenotypes tested in that analysis; significance requires FDR < 0.05 and
  fold enrichment > 1.
- The hypoxia test of the union gene set is a single prespecified test and is not
  corrected for multiple comparisons.
- Union Model B sampling: 167 genes from the SV ±5 kb background, then 130 genes
  from the remaining genome-wide background, yielding 297 unique genes per
  permutation; GWAS and hypoxia overlaps are computed from the same random sets.

Default analysis sizes are reported in the pipeline output and are not fixed in
code; the scripts accept any gene lists and backgrounds of the appropriate format.

## References

- GWAS Catalog (EMBL-EBI): https://www.ebi.ac.uk/gwas/
- GENCODE: https://www.gencodegenes.org/
- CHM13 gene annotation (RefSeq Liftoff v5.2): UCSC/RefSeq annotation resources
- Ensembl REST API (gene identifier lookup): https://rest.ensembl.org/
- bedtools: https://bedtools.readthedocs.io/
- pandas, NumPy, SciPy: standard scientific Python packages
