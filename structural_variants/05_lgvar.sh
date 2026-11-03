SAMPLE_LIST="./sample.txt"
REF=/path/to/reference.fa
PAT=/path/to/sample.paternal.fasta
MAT=/path/to/sample.maternal.fasta
PAF=/path/to/sample.aln.paf

run_lgvar(){
    SAMPLE="$1"
    mkdir -p ./${SAMPLE}
    cd ./${SAMPLE}
    singularity exec ../../lgvar_v1.4.1.sif /LGVAR/LGVAR run \
        -r ${REF} \
        -q1 ${MAT}/${SAMPLE}_hap1.fasta \
        -q2 ${PAT}/${SAMPLE}_hap2.fasta \
        -p1 ${PAF}/${SAMPLE}_hap1_chm13.paf \
        -p2 ${PAF}/${SAMPLE}_hap2_chm13.paf \
        -m sensitive \
        -dv 0.05 \
        -d 500 \
        -s ${SAMPLE} \
        -cp1 hap1.tsv \
        -cp2 hap2.tsv
    cd ..
}

export -f run_lgvar
cat ${SAMPLE_LIST} | xargs -I {} -P 8 bash -c 'run_lgvar "$@"' bash {}
