# Code review and validation

## Scope

Reviewed publication snapshot of the three analysis modules described in the
Methods. No plotting code or source-data tables are included. No original source
scripts or biological results were modified. The package is analysis-only, not
a raw-read processing workflow.

## Changes in this snapshot

- Removed embedded enrichment plots and R MDS/volcano outputs.
- Replaced cluster-specific default paths with relative input paths.
- Preserved selection matching, interval rules, random seed, sampling algorithm
  and original endpoint families in enrichment. Historical optional preparation
  paths remain to allow complete-family FDR reproduction; they are not needed by
  the cCRE-only example.
- Corrected single-base non-INS coordinates to one BED base rather than two.
- Included nearby hits at the exact specified gap; `nearby-bp=0` is direct-only.
- Rejected duplicate SV IDs, invalid coordinates and fractional BED coordinates.
- R: rejected missing counts instead of converting them to zero, duplicate feature
  IDs instead of renaming them, invalid sample designs and zero-size libraries.
  ATAC must contain raw integer counts; RSEM expected RNA counts may be fractional.
- Added R `sessionInfo.txt` output and a server-side synthetic model test.
- Added known Z-score matrix assay labels and header/row-width checks.
- Left unobserved tissue summary values missing rather than filling with zero.
- Made missing requested gene-link files fail explicitly. Missing optional
  selection metadata remain NA; empty target results no longer break report reads.
- Added UCSC liftOver wrapper: BED4 record IDs preserve BED3+3 attributes,
  minMatch=0.95, no multiple mappings, unmapped-record/QC outputs.

These validation and interval-edge fixes can change behavior on malformed data,
single-base non-INS records or exact-window-edge nearby hits. They do not constitute
a rerun or replacement of the manuscript results; use the original result tables
until a full-data rerun is compared.

## Checks

Local result: all 18 Python tests passed; eight Python sources parsed and all six
analysis command-line entry points passed `--help` checks.

- Python source syntax and command-line help checks.
- Synthetic cCRE overlap -> two-header Z-score -> target-gene workflow,
  including no matched target genes.
- Direct/nearby boundaries, nested intervals, one-base INS/non-INS, accession-alias
  deduplication, invalid coordinates and frequency rejection.
- Selection/background separation, unique within-draw controls, reproducibility,
  insufficient-stratum rejection, BH calculations and all-endpoint synthetic run.
- A 10,000-draw synthetic test of the empirical-P calculation.
- Static scan for Chinese text, user/cluster absolute paths and plotting imports.

## Limitations

Rscript and UCSC liftOver were not available on the local validation machine.
The R model and real-genome liftOver have not been executed for this snapshot.
Run `Rscript test_population_model.R` in the server environment. The conda
environment specification has not been freshly installed/solved. Python tests use
the local runtime rather than the exact server package build. Full biological
data, real Cooler extraction, and original numerical results were not rerun.

The retained keyword-based tissue grouping is descriptive, not a validated tissue
ontology. Public gene links are not evidence of a causal SV effect. Gene-symbol
and Ensembl ID agreement must be checked in `gene_link_scan_qc.tsv`; ambiguous
matches are excluded. Multiple linked genes may qualify one cCRE, but overlap
counts remain one per SV per endpoint.
