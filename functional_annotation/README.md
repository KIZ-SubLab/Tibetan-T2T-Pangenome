# Functional annotation and enrichment

This workflow annotates Tibetan-enriched structural variants (TESVs) with
ENCODE4 candidate cis-regulatory elements (cCREs). It models population and time
effects in HUVEC ATAC-seq and RNA-seq data and tests enrichment against matched
structural variants. Run the commands from this directory.

The [trait and hypoxia gene enrichment workflow](trait_enrichment/) analyzes selection and TESV candidate genes for Figure 4 using GWAS trait annotations and a curated hypoxia-related gene set.

## Contents

| Script | Purpose |
| --- | --- |
| `lift_ccre.py` | UCSC liftOver of registry intervals with accession/class preservation |
| `annotate_sv326_with_encode4_ccre.py` | Direct and nearby SV-cCRE annotation |
| `annotate_sv326_ccre_tissue_zscores.py` | Accession-based joins to two-header Z-score matrices |
| `annotate_sv326_ccre_target_genes.py` | Accession-based 3D-chromatin, CRISPR and eQTL target links |
| `fit_additive_population_model.R` | TMM, limma-voom and individual-blocked additive model |
| `prepare_regions.py` | Genomewide regions, independent of selected SVs |
| `matched_enrichment.py` | Stratified background sampling, empirical P values and QC |
| `test_*.py` | Synthetic regression and end-to-end tests |

## Environment

```bash
conda env create -f environment.yml
conda activate tesv-annotation
python -m unittest discover -v
Rscript test_population_model.R
```

The environment specification uses the recorded analysis package versions; it is
not a lockfile and has not been solved locally. UCSC liftOver is external. Cooler
is needed only for the optional historical Hi-C region preparation path.

## Input contracts

- SV CSV: `SV_ID` or `svID_5caller`, e.g. `chr2-46472155-46475559-DEL-3404`.
  IDs use 1-based inclusive reference coordinates. `PBS` is required for the
  enrichment selected set. Additional selection and frequency columns are kept.
- Universe CSV: `KEY`, `SV_ID` or `svID_5caller`; `AF_TIB`, `AF_HAN` in [0,1].
  Every selected ID must be present exactly once; fuzzy matching is prohibited.
- Lifted cCRE BED: six tab-separated columns, chromosome, start, end, rDHS
  accession, cCRE accession, class. Coordinates are BED0 half-open in CHM13v2.0.
- Count matrices: tab-separated, `peak_id`/`feature_id` for ATAC or
  `gene_id`/`feature_id` for RNA, plus 50 sample columns `A-1-1` through `W-5-5`.
  A = Tibetan, W = Han; final suffix 1--5 = 0 h, 6 h, 1 d, 3 d, 5 d.
  Use integer raw ATAC fragments and unrounded RSEM expected RNA counts, not TPM/RPKM.
- Raw V4 gene-link files: accession, gene ID, gene symbol, biotype followed by
  evidence fields. Supply all three collections for the manuscript analysis.
- Z-score matrices: experiment header, biosample header, then accession-indexed
  numeric records. Do not remove the second header.

## 1. cCRE annotation

```bash
python lift_ccre.py --bed data/GRCh38-cCREs.bed \
  --chain data/grch38-chm13v2.chain --outdir results/liftover

python annotate_sv326_with_encode4_ccre.py --sv-csv data/SV326.csv \
  --ccre-bed results/liftover/ccre_t2t.bed --nearby-bp 1000 \
  --ins-overlap-mode span --outdir results/ccre

python annotate_sv326_ccre_tissue_zscores.py \
  --ccre-hits results/ccre/sv326_encode4_ccre_hits.tsv \
  --zscore-dir data/zscores --zscore-files ENCFF833CSA.txt.gz ENCFF560UZM.txt.gz \
  --z-threshold 1.64 --keep-all-scores --outdir results/tissue

python annotate_sv326_ccre_target_genes.py \
  --ccre-hits results/ccre/sv326_encode4_ccre_hits.tsv \
  --gene-link-dir data/Human_Gene_Links --outdir results/targets
```

liftOver uses `-minMatch=0.95` without `-multiple`. Source accessions/classes are
restored after mapping BED4 record IDs; ENCODE's six columns are not standard
BED6 score/strand fields. Returned spans, including internal chain gaps, are used.
Direct overlap requires at least one shared base. An INS with equal endpoints
occupies one reference base; `span` otherwise uses reported coordinates, not
inserted sequence length. `--ins-overlap-mode point` is a sensitivity analysis.
Nearby hits are disjoint from direct hits and include the exact distance cutoff.

Z >= 1.64 is descriptive reference activity, not tissue specificity or an SV
effect test. Raw biosample names and scores are retained; keyword tissue groups
are only summaries. Missing entries are not measured zero activity. Public gene
links nominate targets; they do not demonstrate deletion-genotype eQTL effects.

## 2. Population and hypoxia tests

```bash
Rscript fit_additive_population_model.R --matrix data/atac_raw_counts.tsv \
  --assay atac --outdir results/population_atac
Rscript fit_additive_population_model.R --matrix data/rna_expected_counts.tsv \
  --assay rna --outdir results/population_rna
```

Filter CPM >= 1 in >= 5 samples; TMM normalization; limma-voom with categorical
`population + time`, `duplicateCorrelation(individual)` and robust empirical-Bayes
moderation. Population uses a moderated t-test; time uses a moderated F-test of
the four post-baseline coefficients. BH correction is over all retained features
within each assay/test. The biological sample size is ten individuals, not fifty
independent replicates. No population-by-time interaction or genotype model is
included. Primary outputs: `population_Tibetan_vs_Han.tsv`, `time_omnibus.tsv`,
`normalized_logCPM.tsv.gz`, `model_qc.tsv`, and `sessionInfo.txt`.

## 3. Matched background enrichment

```bash
python prepare_regions.py --ccre-bed results/liftover/ccre_t2t.bed \
  --rna-time results/population_rna/time_omnibus.tsv \
  --gene-link data/Human_Gene_Links/V4-hg38.Gene-Links.3D-Chromatin.txt \
  --gene-link data/Human_Gene_Links/V4-hg38.Gene-Links.CRISPR.txt \
  --gene-link data/Human_Gene_Links/V4-hg38.Gene-Links.eQTLs.txt \
  --outdir results/regions

python matched_enrichment.py --selected data/SV326.csv \
  --universe data/universe_sv.csv --manifest results/regions/region_manifest.tsv \
  --expected-selected 326 --n-perm 10000 --seed 20260928 \
  --ins-mode span --outdir results/enrichment
```

The selected set excludes rows without numeric PBS (the appended VNTRs in the
original input); the expected-count check must pass. VNTR annotation is supported,
but these VNTRs are not part of the 326-SV matched enrichment test.
Controls exclude selected IDs. Match SV type, right-closed length bins and frequency
bins; the proxy is `min((AF_TIB+AF_HAN)/2, 1-(AF_TIB+AF_HAN)/2)` by default.
Draw 326 distinct controls per iteration without replacement; controls can recur
between iterations. Inadequate strata fail rather than relaxing matching.
The same draws are used across endpoints. Count each SV once per endpoint.
Enrichment P = `(1 + number of null counts >= observed)/(B+1)`; fold enrichment =
observed/null mean. Inspect matching balance; bin matching is not exact continuous
matching and does not control every genomic confounder.

The two manuscript endpoints are `cCRE_any_direct` and `cCRE_hypoxia_gene_direct`.
The latter uses any linked gene with time-test FDR < 0.05 and maximum absolute
post-baseline log2FC >= 0.5, without tissue restriction. Use full genomewide RNA
results and gene links, not links already filtered to TESVs. Gene IDs/symbols must
match the RNA identifiers; unresolved/conflicting matches are reported in QC.

**Historical multiplicity:** region preparation retains the original nine-endpoint
manifest and optional ATAC/Hi-C preparation code to preserve the original test-family
definitions. The cCRE-only example reports other endpoints as MISSING and the run
as PARTIAL. This is expected, not evidence of zero overlap. Its two nominal
empirical P values are usable, but BH values are not the original full-family
results when other tests are missing (missing tests count as P=1). To reproduce
historical BH values, provide all original genomewide inputs or use the original
complete region manifest/BEDs with `--require-complete`. Do not reinterpret two
displayed endpoints as a new two-test correction family.

No Hi-C population-effect result is needed for the cCRE-only commands. The retained
optional legacy O/E path is descriptive, not a differential-loop significance test.
This package starts from quantified matrices and candidate regions; FASTQ alignment,
peak calling and FitHiC discovery are not implemented here.

## Validation and release scope

See `REVIEW.md` for changes, checks and limitations. Supply research data separately;
this folder intentionally contains neither a data-source table nor private data.
Use a fresh output directory for each run to avoid confusing stale and new results.
