# Tandem repeat analysis

These scripts support the analysis of tandem repeats as complex structural variants in Figure 5. The workflow genotypes repeat loci from PacBio HiFi reads, compares repeat copy-number distributions between populations, and annotates repeat motifs at selected loci.

| Script | Analysis |
| --- | --- |
| [01-trgt.sh](01-trgt.sh) | Genotype tandem repeats with TRGT using T2T-CHM13 v2.0 and a VAMOS motif catalog. |
| [02-TRDS.cal.r](02-TRDS.cal.r) | Calculate the tandem repeat disparity score (TRDS) between population copy-number distributions using the Wasserstein metric. |
| [03-Rst_cal.r](03-Rst_cal.r) | Estimate per-locus R_ST across repeated equal-size subsamples with missingness and heterozygosity filters. |
| [04-vampire.sh](04-vampire.sh) | Annotate repeat motifs in locus-specific FASTA files with VAMPIRE. |

TRGT uses the reference FASTA, the repeat catalog, and a HiFi alignment BAM for each sample. The population analyses use locus-by-haplotype copy-number matrices. Set the input paths and population column ranges in each R script to match your cohort. The motif annotation script takes per-locus directories with the sample sequences and motif FASTA. Install TRGT, VAMPIRE, R, and the R packages data.table and waddR.
