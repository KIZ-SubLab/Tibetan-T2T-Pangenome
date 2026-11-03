# SV calling pipeline

SV discovery from PacBio HiFi reads and haplotype-resolved assemblies against
T2T-CHM13 v2.0.

| Script | Input | Caller |
|---|---|---|
| `01_sniffles2.sh` | HiFi reads | minimap2 + Sniffles2 v2.3.3 |
| `02_pbsv.sh` | HiFi reads | pbmm2 + PBSV v2.9.0 |
| `03_svim_asm.sh` | hap1/hap2 assemblies | minimap2 (asm5) + SVIM-asm v1.0.3 |
| `04_pav.sh` | hap1/hap2 assemblies | PAV v1.1.2 |

## Requirements

- minimap2 v2.26-r1175
- samtools v1.19.2, bcftools
- Sniffles2 v2.3.3
- pbsv v2.9.0 (includes pbmm2)
- SVIM-asm v1.0.3
- PAV v1.1.2

## Inputs

- Reference genome: T2T-CHM13 v2.0 FASTA
- HiFi reads: adapter-filtered FASTQ
- Assemblies: haplotype-resolved FASTA (hap1 = paternal, hap2 = maternal)
- Mask BED: centromeres, telomeres and their ±5-Mb flanking regions in
  reference coordinates

## Usage

Edit the variable block at the top of each script (paths, sample name,
threads), then run:

```bash
bash 01_sniffles2.sh
bash 02_pbsv.sh
bash 03_svim_asm.sh
bash 04_pav.sh
```

Each script performs alignment (read-based callers only), SV calling, and
filtering in one run.

## Filtering

Calls are filtered to keep:

- SVs of at least 50 bp (and at most 1 Mbp downstream)
- precise breakpoints (read-based: `IMPRECISE` removed; SVIM-asm: `FILTER=PASS`)
- variants outside centromeres, telomeres and their ±5-Mb flanks

## Output

Per-sample filtered VCF: `<sample>.sniffles.clean.vcf`, `<sample>.pbsv.clean.vcf`,
`<sample>.svimasm.clean.vcf`, `<sample>.pav.clean.vcf`.
