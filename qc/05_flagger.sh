#1. create a whole-genome BED file
cd "${sample}" || exit 1
ln -s "${sample}.fasta" ./
samtools faidx "${sample}.fasta"
cat "${sample}.fasta.fai" | awk '{print $1"\t0\t"$2}' > whole_genome.bed
    
#2.convert bam to cov
# Put the path to the whole-genome bed file in a json file
echo "{" > annotations_path.json
echo \"whole_genome\" : \"whole_genome.bed\" >> annotations_path.json
echo "}" >> annotations_path.json


# Convert bam to cov.gz with bam2cov program
singularity run flagger.v1.2.0.sif bam2cov \
       	--bam "${sample}.aln.sorted.bam" \
	--output ./coverage_file.cov.gz \
	--annotationJson ./annotations_path.json \
	--threads 16 \
	--baselineAnnotation whole_genome

#3.run hmm_flagger
mkdir -p ./hmm_flagger_outputs

singularity run flagger.v1.2.0.sif hmm_flagger \
	--input ./coverage_file.cov.gz \
	--outputDir ./hmm_flagger_outputs \
	--alphaTsv alpha_optimum_trunc_exp_gaussian_w_16000_n_50.HiFi_DC_1.2_DEC_2024.v1.1.0.tsv \
	--labelNames Err,Dup,Hap,Col \
	--threads 16

