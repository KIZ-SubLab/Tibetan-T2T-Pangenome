# Site-level PBS simulation

Empirical P values for single-variant PBS from a neutral *msprime* simulation.
Each simulated site is annotated with its own point recombination (`R`) and
mutation (`MU`) rates; observed and simulated sites are matched on
`(LOG_R, LOG_MU)`.

## Requirements

- Python: `msprime`, `scikit-allel==1.3.13`, `numpy`
- `bcftools`
- C++17 compiler with OpenMP and zlib
- Inputs: gmap recombination map per chromosome, 25-kb mutation-rate table,
  sample-ID lists for TIB / HAN / EUR, observed VCF

## Statistics

```
FST  = a / (a + b + c)                    (Weir & Cockerham 1984)
T    = -ln(1 - FST)
PBS  = (T_TIB-HAN + T_TIB-EUR - T_HAN-EUR) / 2
LOG_R  = log10(R_point  + 1e-4)
LOG_MU = log10(MU_point + 1e-12)
match  = |LOG_R_sim - LOG_R_obs| <= tol  and  |LOG_MU_sim - LOG_MU_obs| <= tol
P_SIM  = (N_GE + 1) / (N_MATCH + 1)        (N_MATCH >= 10,000 required)
```

BH-FDR is applied to the `PBS > 0` and valid sites only.

## Steps

```bash
# 1. simulate one chromosome / replicate
python scripts/1_simulate.py --chrom chr22 --rep 1 \
    --recomb-map maps/chr22.gmap --mutation-map maps/mutation_25kb.tsv.gz \
    --out vcf/rep0001.chr22.vcf.gz

# 2. per-site PBS (add --no-clean for simulated VCFs)
python scripts/2_site_pbs.py --vcf vcf/rep0001.chr22.vcf.gz --chrom chr22 --no-clean \
    --tib TIB.samples --han HAN.samples --eur EUR.samples \
    --out sites/rep0001.chr22.pbs.tsv.gz

# 3. annotate observed sites with point rates -> queries + index
python scripts/3_annotate.py --sites observed/chr22.pbs.tsv.gz --chrom chr22 \
    --recomb-map maps/chr22.gmap --mutation-map maps/mutation_25kb.tsv.gz \
    --out queries/chr22.tsv

# 4. build the null database (one shard = a group of replicates x 22 chromosomes)
g++ -O2 -std=c++17 -fopenmp -o bin/4_build scripts/4_build.cpp -I. -lz
echo -e "chr22\tsites/rep0001.chr22.pbs.tsv.gz" > files.txt
./bin/4_build files.txt shards/SHARD_0001 work 16 maps _noMask.scaled.gmap \
    maps/mutation_25kb.tsv.gz

# 5. query
g++ -O2 -std=c++17 -fopenmp -o bin/5_query scripts/5_query.cpp -I.
./bin/5_query shards/SHARD_0001 queries/chr22.tsv results/q.out 32 0.02 \
    --full-cell-mode merged

# 6. aggregate and BH
python scripts/6_aggregate.py --queries queries/chr22.tsv \
    --index queries/chr22.tsv.index --shards 'results/shard_*.out' --out final.tsv
```

## Output

```
CHROM POS PBS LOG_R LOG_MU N_MATCH_TOTAL N_GE_TOTAL P_SIM_PLUS1 P_BH_PBSpos MATCH_VALID
```

`MATCH_VALID=0` and `NA` P values mark sites with fewer than 10,000 simulated
matches. `NA` in `P_BH_PBSpos` means the site was outside the BH family
(`PBS <= 0` or invalid), not that its adjusted P equals 1.

## Key parameters

| | |
|---|---|
| Grid | 568 x 60 coarse, 8 x 8 fine (2,181,120 cells) |
| LOG_R range / LOG_MU range | [-4.5, 2.2] / [-9.5, -5.5] |
| Matching tolerance | 0.02 (both dimensions, fixed) |
| Minimum matches | 10,000 |
| Replicates | 10,000 per chromosome |

## Note

The simulation produces SNPs only. Indels are matched against the SNP
site-level null using their own point rates.
