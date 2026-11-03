#!/bin/bash
#SBATCH --job-name run.sprime
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --partition=wzhcnormal
#SBATCH -e run.sprime.err
#SBATCH -o run.sprime.log

bcftools view \
-r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 \
-Oz \
-o TBN35_YRI108.grch38.merge.00.exMiss.noChr.autosome.vcf.gz \
/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/sprime/TBN35_YRI108.grch38.merge.00.exMiss.noChr.vcf.gz
tabix -p vcf TBN35_YRI108.grch38.merge.00.exMiss.noChr.autosome.vcf.gz

outgroup=/work/home/songjiaxu/introgression/data_t2t/outgroup.txt
gt=TBN35_YRI108.grch38.merge.00.exMiss.noChr.autosome.vcf.gz

for chr in {1..22};
do
map=/work/home/songjiaxu/introgression/data_t2t/plink.GRCh38.map/plink.all.GRCh38.map
out=tbn.yri.${chr}
java -jar /work/home/songjiaxu/introgression/data_t2t/sprime.jar gt=${gt} outgroup=${outgroup} map=${map} out=${out} chrom=${chr}
done
