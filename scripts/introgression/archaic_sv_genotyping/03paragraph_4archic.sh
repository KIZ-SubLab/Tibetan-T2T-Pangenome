#!/bin/bash
mkdir -p bam_archic;cd bam_archic
for i in `cat /work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/output/rm1kb-na-segments/merged_fragments_v3/SV/denove_genotyping_SV_35tbn_5archic_dna/archic.sample.txt`
do
mkdir -p ${i}
echo "#!/bin/bash" > ${i}/command.del.sh
echo "#SBATCH --job-name=${i}.del.paragraph" >> ${i}/command.del.sh
echo "#SBATCH --nodes=1" >> ${i}/command.del.sh
echo "#SBATCH --ntasks=1" >> ${i}/command.del.sh
echo "#SBATCH --cpus-per-task=32" >> ${i}/command.del.sh
echo "#SBATCH --partition=wzhcnormal" >> ${i}/command.del.sh
echo "#SBATCH -o ${i}.del.out" >> ${i}/command.del.sh
echo "echo \"${i}\" > target" >> ${i}/command.del.sh
echo "" >> ${i}/command.del.sh
echo "perl /work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/output/rm1kb-na-segments/merged_fragments_v3/SV/denove_genotyping_SV_35tbn_5archic_dna/make.sampletxt.pl /work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/output/rm1kb-na-segments/merged_fragments_v3/SV/denove_genotyping_SV_35tbn_5archic_dna/archic.depth.txt target ${i}.sample.txt" >> ${i}/command.del.sh
echo "" >> ${i}/command.del.sh
echo "## DEL typing" >> ${i}/command.del.sh

echo "#!/bin/bash" > ${i}/command.ins.sh
echo "#SBATCH --job-name=${i}.ins.paragraph" >> ${i}/command.ins.sh
echo "#SBATCH --nodes=1" >> ${i}/command.ins.sh
echo "#SBATCH --ntasks=1" >> ${i}/command.ins.sh
echo "#SBATCH --cpus-per-task=32" >> ${i}/command.ins.sh
echo "#SBATCH --partition=wzhcnormal" >> ${i}/command.ins.sh
echo "#SBATCH -o ${i}.ins.out" >> ${i}/command.ins.sh
echo "time python /work/home/songjiaxu/singlecells/tbn_han_fastq/1001TBN_hg19vcf/38tbn_NGS/mapping_results/scripts/paragraph/bin/multigrmpy.py -i /work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/output/rm1kb-na-segments/merged_fragments_v3/SV/denove_genotyping_SV_35tbn_5archic_dna/5caller_paragraph_ready.vcf.gz -m ${i}.sample.txt -r /work/home/songjiaxu/introgression/fastq/ref/chm13v2.0.fa -o ${i}.ins -t 32" >> ${i}/command.ins.sh
done


