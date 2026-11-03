# Minigraph-Cactus pangenome

Build a chromosome-level pangenome graph from ~940 haplotype assemblies with
[Minigraph-Cactus](https://github.com/ComparativeGenomicsToolkit/cactus), then
filter the VCF with [vcfbub](https://github.com/pangenome/vcfbub) and normalise
it with [vcfwave](https://github.com/vcflib/vcflib).

## Requirements

cactus >= 3.x (`cactus-pangenome`, `vcfbub`), vcflib `vcfwave`, `bcftools`,
`bgzip`, `singularity`.

## Steps

1. **seqfile** — one per chromosome, `name<TAB>fasta`, reference sequence first,
   then every haplotype.
2. `01_mc.sh` — `cactus-pangenome` per chromosome (`--filter 97`, `--vcf clip`,
   `--mapCores 16`).
3. `02_vcfbub.sh` — `vcfbub -l 0 -a 100000`.
4. `03_vcfwave.sh` — split each VCF into 5 Mb windows, `vcfwave -I 1000` per
   window, concatenate, sort and index per chromosome.

Edit the paths and resource settings at the top of each script, then run them
in order.

## Output

```
vcfbub/<chr>.vcfbub.vcf.gz
wave/final/<chr>.wave.vcf.gz (+ .tbi)
```

## Notes

- vcfwave can run for days without output on very large alleles (tens of kb) or
  on telomere/centromere repeats. Remove those records before running it.
- Do not start several vcfwave processes in parallel inside one container; they
  can deadlock. Parallelise by splitting the window list across jobs instead.
- The final `bcftools sort` is required: vcfwave can emit a position rollback at
  a window boundary, which leaves the concatenated VCF un-sortable/un-indexable.
