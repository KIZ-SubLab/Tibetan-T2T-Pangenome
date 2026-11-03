#!/usr/bin/env python
"""Simulate one chromosome for one replicate under a neutral demographic model.

Output: a VCF containing biallelic, segregating, single-mutation sites.

Example
-------
python 01_simulate_chromosome.py --chrom chr22 --rep 1 \
    --recomb-map maps/chr22_noMask.scaled.gmap \
    --mutation-map maps/mutation_rate_25kb.tsv.gz \
    --out vcf/rep0001.chr22.vcf.gz
"""
import argparse
import gzip
import math

import numpy as np

# Population sizes and split times (generations, generation time 29 years).
GEN_TIME = 29
NE_TIB, NE_CHB, NE_EAS, NE_OOA, NE_ANC = 36260, 50870, 10206, 2831, 11293
T_TIB_CHB = 16871 / GEN_TIME
T_EAS_CEU = 46000 / GEN_TIME
T_OOA = 119000 / GEN_TIME
NE_CEU_PRESENT = 2512 * math.exp(0.0016 * T_EAS_CEU)
CEU_GROWTH = 0.0016

# Sample sizes (diploid individuals).
N_TIB, N_CHB, N_CEU = 35, 102, 38

DEFAULT_SEED = 20260825


def build_demography(msprime):
    d = msprime.Demography()
    d.add_population(name="TIB", initial_size=NE_TIB)
    d.add_population(name="CHB", initial_size=NE_CHB)
    d.add_population(name="CEU", initial_size=NE_CEU_PRESENT, growth_rate=CEU_GROWTH)
    d.add_population(name="EAS_ANC", initial_size=NE_EAS)
    d.add_population(name="OOA", initial_size=NE_OOA)
    d.add_population_split(time=T_TIB_CHB, derived=["TIB", "CHB"], ancestral="EAS_ANC")
    d.add_population_split(time=T_EAS_CEU, derived=["EAS_ANC", "CEU"], ancestral="OOA")
    d.add_population_parameters_change(time=T_OOA, initial_size=NE_ANC, population="OOA")
    d.sort_events()
    return d


def load_gmap(path, msprime):
    """Convert a gmap (pos, cM/Mb, cM) to an msprime.RateMap (per bp/gen)."""
    pos, rate = [], []
    with open(path) as f:
        f.readline()
        for line in f:
            a = line.split()
            pos.append(float(a[0]))
            rate.append(float(a[1]) * 1e-8)  # cM/Mb -> per bp per generation
    return msprime.RateMap(position=pos, rate=rate[:-1])


def load_mutation_map(path, chrom, seq_len, msprime):
    """Convert a 25-kb mutation-rate table to an msprime.RateMap.

    The RateMap is clipped to `seq_len` (the sequence length of the
    recombination RateMap), as msprime requires both maps to span the same
    interval exactly.
    """
    pos, rate = [0.0], []
    with gzip.open(path, "rt") as f:
        f.readline()
        for line in f:
            a = line.rstrip("\n").split("\t")
            if a[0] != chrom:
                continue
            start, end, mu = int(a[1]), int(a[2]), float(a[9])
            if start >= seq_len:
                break
            pos.append(float(min(end, seq_len)))
            rate.append(mu)
    if not rate:
        raise ValueError(f"no mutation-map bins for {chrom}")
    if pos[-1] < seq_len:          # extend the last rate to the sequence end
        pos.append(float(seq_len))
        rate.append(rate[-1])
    return msprime.RateMap(position=pos, rate=rate)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--rep", type=int, required=True)
    ap.add_argument("--recomb-map", required=True)
    ap.add_argument("--mutation-map", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=DEFAULT_SEED)
    args = ap.parse_args()

    import msprime
    import allel

    chr_num = int(args.chrom.replace("chr", ""))
    ss = np.random.SeedSequence([args.seed, args.rep, chr_num])
    ancestry_seed, mutation_seed = (int(x.generate_state(1, dtype=np.uint32)[0])
                                    for x in ss.spawn(2))

    demography = build_demography(msprime)
    recomb_map = load_gmap(args.recomb_map, msprime)
    mut_map = load_mutation_map(args.mutation_map, args.chrom,
                                recomb_map.sequence_length, msprime)

    ts = msprime.sim_ancestry(
        samples=[
            msprime.SampleSet(N_TIB, population="TIB", ploidy=2),
            msprime.SampleSet(N_CHB, population="CHB", ploidy=2),
            msprime.SampleSet(N_CEU, population="CEU", ploidy=2),
        ],
        demography=demography,
        recombination_rate=recomb_map,
        random_seed=ancestry_seed,
    )
    mts = msprime.sim_mutations(
        ts,
        rate=mut_map,
        model=msprime.BinaryMutationModel(),
        discrete_genome=True,
        random_seed=mutation_seed,
    )

    # Keep single-mutation sites that are segregating in the sampled individuals.
    sites = list(mts.sites())
    keep_site = np.array([len(s.mutations) == 1 for s in sites])
    positions = np.array([s.position for s in sites], dtype=np.int64)[keep_site]
    gm = mts.genotype_matrix().T           # (n_haplotypes, n_sites)
    gm = gm[:, keep_site]
    ac = gm.sum(axis=0)
    seg = (ac > 0) & (ac < gm.shape[0])
    positions, gm = positions[seg], gm[:, seg]

    n_hap, n_sites = gm.shape
    n_ind = n_hap // 2
    hap0, hap1 = gm[0::2, :], gm[1::2, :]
    codes = hap0 + hap1
    lookup = np.array(["0/0", "0/1", "1/1"])

    sample_names = ([f"TIB{i:03d}" for i in range(1, N_TIB + 1)]
                    + [f"HAN{i:03d}" for i in range(1, N_CHB + 1)]
                    + [f"EUR{i:03d}" for i in range(1, N_CEU + 1)])

    with gzip.open(args.out, "wt") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                + "\t".join(sample_names) + "\n")
        chunk = 50000
        for i0 in range(0, n_sites, chunk):
            i1 = min(n_sites, i0 + chunk)
            gs = lookup[codes[:, i0:i1]].T
            for i in range(i0, i1):
                f.write(f"{args.chrom}\t{positions[i]}\t.\tA\tG\t.\tPASS\t.\tGT\t"
                        + "\t".join(gs[i - i0]) + "\n")

    print(f"{args.chrom} rep {args.rep}: {n_sites} segregating sites -> {args.out}")


if __name__ == "__main__":
    main()
