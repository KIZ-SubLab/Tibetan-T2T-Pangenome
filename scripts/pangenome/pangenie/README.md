# PanGenie genotyping pipeline

Minimal command-line scripts for genotyping short-read samples with PanGenie from a
phased, multi-sample reference-panel VCF.

The pipeline is:

```
panel VCF -> pangenome VCF -> PanGenie index -> per-sample genotyping -> merged VCF
```

## Requirements

- PanGenie v4.x (tested with v4.1.1 SIF and v4.2.1)
- Singularity/Apptainer (or install PanGenie natively)
- `bcftools`, `samtools`, `bgzip`/`tabix`
- `snakemake` and `python3` for the official `run-from-callset` preparation pipeline
- The official PanGenie repository (for `pipelines/run-from-callset`)

## Inputs

- Reference genome FASTA
- Fully phased, multi-sample panel VCF
- Paired-end FASTQ files for each sample
- Sample sheet: `sample<TAB>R1<TAB>R2` (one sample per line)

## Configuration

Edit `config.sh` and set the paths and thread count:

```bash
PANGENIE_SIF=/path/to/pangenie.sif
PANGENIE_REPO=/path/to/pangenie
REFERENCE=/path/to/reference.fa
PANEL_VCF=/path/to/panel.phased.vcf.gz
WORKDIR=/path/to/pangenie_run
SAMPLE_SHEET=/path/to/samples.tsv
THREADS=16
```

`config.sh` also defines `pangenie_exec`, which runs PanGenie through the SIF. If
PanGenie is installed natively, replace it with `pangenie_exec() { "$@"; }`.

## Steps

```bash
# 1. Clean the panel VCF (symbolic DUP/INS, symbolic DEL, missing phase) and sort
bash 01_prepare_panel_vcf.sh

# 2. Generate pangenome.vcf with the official PanGenie preparation pipeline
bash 02_prepare_pangenome.sh

# 3. Build the PanGenie index
bash 03_build_index.sh

# 4. Genotype one sample
bash 04_genotype_sample.sh SAMPLE /path/to/R1.fastq.gz /path/to/R2.fastq.gz

#    or all samples in the sample sheet
bash 05_batch_genotype.sh

# 5. Compress and merge all per-sample VCFs
bash 06_merge_vcfs.sh

# 6. Optional: convert bubble genotypes to biallelic variant genotypes
bash 07_convert_to_biallelic.sh
```

## Outputs

- `results/pangenome/pangenome.vcf`
- `results/index/pangenie_*`
- `typing/<sample>/<sample>_genotyping.vcf`
- `results/genotypes/all_samples.merged.vcf.gz`
- `results/genotypes/all_samples.biallelic.vcf.gz` (optional step)

## Notes

- `PanGenie-index` in v4.1.1 supports `-t` only (not `-j`).
- `02_prepare_pangenome.sh` uses the official pipeline's `Snakefile`; add
  `--use-conda` if you want Snakemake to build the rule environments.
- `01_prepare_panel_vcf.sh` removes symbolic records other than `<DEL>` and resolves
  symbolic `<DEL>` records against the reference. If your callset is already
  sequence-resolved and fully phased, this step can be skipped.
- `04_genotype_sample.sh` converts FASTQ to a combined FASTA because PanGenie is
  run in FASTA mode here; paired-end order does not need to be preserved.
- For large cohorts, `06_merge_vcfs.sh` merges per-sample VCFs in groups of 250
  and then merges the group files.
- `07_convert_to_biallelic.sh` uses `convert-to-biallelic.py` from the official
  repository and the `callset-biallelic.vcf.gz` produced by the preparation
  pipeline.
