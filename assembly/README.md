# Genome assembly pipeline

Trio-binned, ultra-long-read-assisted de novo genome assembly with
[hifiasm](https://github.com/chhylp123/hifiasm), followed by long-read
polishing with [nextPolish2](https://github.com/Nextomics/nextPolish2).
The scripts reproduce the assembly workflow used for the
telomere-to-telomere (T2T) genomes in this study.

## Workflow

```
Illumina PE reads ──fastp──> cleaned reads ──yak k19/k21/k31──> k-mer databases
PacBio HiFi BAM ──HiFiAdapterFilt──> filtered HiFi reads
ONT pass reads ──NanoFilt (>= 100 kb)──> ultra-long ONT reads

paternal k31 k-mers ─┐
maternal k31 k-mers ─┤
                     ▼
      hifiasm (trio phasing + --ul ONT) -> phased hap1 / hap2 contigs
                     ▼
   minimap2 (HiFi) -> nextPolish2 (k21 + k31) -> polished hap1 / hap2
```

The assembly is phased with parental Illumina k-mer sets, while ONT reads
>= 100 kb are used as ultra-long reads to resolve repeats and extend contigs
towards telomeres. The two haplotypes are then polished independently with
HiFi reads and the child's short-read k-mer sets.

## Input data

```
${PROJECT_ROOT}/data_tmp/
├── NGS/${sample}/${sample}_1.fq.gz        # Illumina paired-end reads
├── NGS/${sample}/${sample}_2.fq.gz
├── HIFI/${sample}/${sample}.hifi_reads.bam
└── ONT/${sample}/${sample}.pass.fastq.gz
```

`${sample}` follows the sample IDs in the sample sheet. Paths are configured
in `config.sh` and can be overridden with environment variables.

### Sample sheet

`config.sh` points to a whitespace- or tab-delimited sample sheet with four
columns. Lines starting with `#` are ignored.

```
# family  paternal  maternal  child
TBN1      T2T007    T2T008    T2T009
TBN2      T2T010    T2T011    T2T012
```

The `family` ID is used as the assembly prefix; `paternal`, `maternal` and
`child` are sample IDs used to locate the raw data.

## Steps

| Script | Purpose | Key parameters |
|---|---|---|
| `01_ngs_qc_kmer.sh` | fastp trimming; yak k-mer counting | fastp `-u 30 -q 20`; yak `-k19/-k21/-k31 -b37` |
| `02_hifi_filter_qc.sh` | HiFi adapter/contaminant filtering | HiFiAdapterFilt `-l 44 -m 97`; seqkit stats, samtools stats |
| `03_ont_filter_qc.sh` | Ultra-long ONT read selection | NanoFilt `-l 100000`; NanoStat |
| `04_hifiasm_trio.sh` | Phased assembly | `hifiasm -1 pat.k31.yak -2 mat.k31.yak --ul ONT_100k` |
| `05_assembly_stats.sh` | Contiguity and inheritance statistics | seqkit stats; yak trioeval; optional unimap dot plot |
| `06_nextpolish2.sh` | Haplotype-aware polishing | minimap2 `-ax map-hifi`; nextPolish2 with k21 + k31 |

k-mer sets are generated once per sample:

| Sample | k | Purpose |
|---|---|---|
| parents | 19 | `yak trioeval` (SE/HE) in QC |
| parents | 31 | hifiasm trio phasing (`-1` / `-2`) |
| child | 21 | `yak qv` and nextPolish2 |
| child | 31 | nextPolish2 |

## Usage

```bash
cd assembly
# edit paths, sample sheet location and threads
vim config.sh

bash 01_ngs_qc_kmer.sh
bash 02_hifi_filter_qc.sh
bash 03_ont_filter_qc.sh
bash 04_hifiasm_trio.sh
bash 05_assembly_stats.sh
bash 06_nextpolish2.sh
```

The scripts are ordinary Bash scripts and make no assumptions about a
scheduler. On a cluster, submit each script inside the allocation prepared for
your environment, for example with your own `sbatch` wrapper. All steps are
single-sample loops over the sample sheet and process each trio independently.

## Output

```
01_NGS_QC/${sample}/${sample}_1.clean.fq.gz
01_NGS_QC/${sample}/${sample}_2.clean.fq.gz
01_NGS_QC/${sample}/${sample}.k19.yak | .k21.yak | .k31.yak

02_hifi_filter_QC/${child}/${child}.hifi_reads.filt.fastq.gz
02_hifi_filter_QC/${child}/${child}.hifi.seqkit.stats
02_hifi_filter_QC/${child}/${child}.hifi.samtools.stats

03_ont_filter_QC/${child}/${child}.pass_100k.fastq.gz
03_ont_filter_QC/${child}/${child}.pass.stat
03_ont_filter_QC/${child}/${child}.pass_100k.stat

04_hifiasm/${child}/${family}_hifiasm_trio_ont.dip.hap1_ctg.fa
04_hifiasm/${child}/${family}_hifiasm_trio_ont.dip.hap2_ctg.fa
04_hifiasm/${child}/${family}_hifiasm_trio_ont.dip.hap{1,2}_ctg.trioeval
04_hifiasm/${child}/${family}_hifiasm_trio_ont.dip.hap{1,2}_ctg.seqkit.stats

05_nextpolish2/${child}/${family}.hap1.np2.fa
05_nextpolish2/${child}/${family}.hap2.np2.fa
```

The polished `*.hap{1,2}.np2.fa` files are the final assemblies used for
quality control and downstream analyses.

## Notes

- Genome size is assumed to be ~3.1 Gb (human). Adjust the ONT length cutoff
  (100 kb) and any resource-related settings for other species.
- `THREADS` in `config.sh` controls the parallelism of each step;
  `NP2_THREADS` controls nextPolish2 separately (default 5).
- `05_assembly_stats.sh` produces dot plots only when `REFERENCE_FASTA` and
  `UNIMAP_DOTPLOT` are both set in `config.sh`.
- Record the exact tool versions used in your run (`tool --version`) alongside
  the generated results. The commands are otherwise self-contained.
