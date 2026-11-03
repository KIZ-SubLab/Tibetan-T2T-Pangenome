trgt genotype \
    --genome chm13v2.0.fa \
    --repeats CHM13.VAMOS.motif.bed \
    --reads ${sample}.hifi.bam \
    --output-prefix ${sample}_chm13_hifi_VAMOS \
    --sample-name ${sample}_chm13_hifi \
    -k XY \
    -t 36

