# Genome assembly quality control

Quality control of the final, nextPolish2-polished haplotype assemblies.
The scripts reproduce the QV, error-detection and reference-based evaluation
used in this study.

## Evaluated assemblies

```
${ASSEMBLY_DIR}/${child}/${family}.hap1.np2.fa
${ASSEMBLY_DIR}/${child}/${family}.hap2.np2.fa
```

By default `ASSEMBLY_DIR` points to `05_nextpolish2` produced by the
`genome_assembly` pipeline. The sample sheet is the same four-column file used
there (`family, paternal, maternal, child`).

Required intermediate files from the assembly pipeline:

```
${NGS_QC_DIR}/${child}/${child}.k21.yak     # QV
${NGS_QC_DIR}/${paternal}/${paternal}.k19.yak
${NGS_QC_DIR}/${maternal}/${maternal}.k19.yak
${HIFI_QC_DIR}/${child}/${child}.hifi_reads.filt.fastq.gz
```

## Steps

| Script | Metric | Tool |
|---|---|---|
| `01_yak_qv_trioeval.sh` | Consensus QV (k21); SE and HE (k19) | yak qv, yak trioeval |
| `02_inspector.sh` | Base- and structure-level errors | inspector |
| `03_quast.sh` | Contiguity, completeness and accuracy vs CHM13v2 / GRCh38 | QUAST |
| `04_summary.sh` | Combines all metrics into one table | awk |

### Metric definitions

- **QV**: `yak qv -p <child k21> <assembly>`; the reported value is the third
  field of the `QV` line (consistent with the original analysis).
- **SE / HE**: `yak trioeval -p <paternal k19> <maternal k19> <assembly>`;
  the W rate is reported as SE and the H rate as HE.
- **inspector**: `inspector.py --datatype hifi` with the filtered HiFi reads;
  the summary reports contig counts, N50, mapping rate, structural errors and
  its own QV estimate.
- **QUAST**: both haplotypes are evaluated together against CHM13v2 and
  GRCh38 with `--large --est-ref-size 3100000000 --no-icarus`; N50, NG50 and
  genome fraction are extracted from `report.txt`.

## Usage

```bash
cd QC
# edit paths and references
vim config.sh

bash 01_yak_qv_trioeval.sh
bash 02_inspector.sh
bash 03_quast.sh
bash 04_summary.sh
```

`inspector.py` should be run in its own environment (it is distributed as a
Python package with third-party dependencies). Set
`RUN_INSPECTOR_CORRECT=true` in `config.sh` to additionally run
`inspector-correct.py` and re-evaluate the corrected assembly.

## Output

```
QC_results/${child}/${family}.${hap}.np2.qv.k21.txt
QC_results/${child}/${family}.${hap}.np2.trioeval.k19.txt
QC_results/${child}/${family}.${hap}.np2.inspector/summary_statistics
QC_results/${child}/${family}.np2.chm13.quast/report.txt
QC_results/${child}/${family}.np2.hg38.quast/report.txt

QC_results/T2T_QC_summary.tsv
```

`T2T_QC_summary.tsv` contains one row per haplotype with QV, SE, HE,
inspector statistics and QUAST N50 / NG50 / genome fraction against both
references.
