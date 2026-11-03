#!/usr/bin/env bash
# Minigraph-Cactus: build one pangenome graph per chromosome.
set -euo pipefail

REF=REF
SIF=/path/to/cactus.sif
OUT=mc
CORES=64
MEM=960Gi

for CHR in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 \
           chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY; do
    SEQFILE=$OUT/$CHR/$CHR.seqfile    # name<TAB>fasta, reference first

    singularity exec "$SIF" cactus-pangenome "$OUT/$CHR/js" "$SEQFILE" \
        --batchSystem single_machine \
        --maxCores "$CORES" --maxMemory "$MEM" \
        --retryCount 3 --restart \
        --outName "$CHR-mc" --outDir "$OUT/$CHR/graph" \
        --reference "$REF" --vcfReference "$REF" --refContigs "$CHR" \
        --filter 97 --giraffe filter --gbz clip filter --gfa clip \
        --vcf clip --chrom-vg clip filter \
        --mapCores 16
done
