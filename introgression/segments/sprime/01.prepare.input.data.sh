#!/bin/bash
#SBATCH --job-name 01prepare.input.data
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --partition=wzhcnormal
#SBATCH -e 01.prepare.input.data.err
#SBATCH -o 01.prepare.input.data.log

vcf1=/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/data/tbn/TBN35.snps.unphased.vcf.gz
ovcf2=/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/data/yri/YRI108.grch38.biSNV.vcf.gz
bcftools merge ${vcf1} ${ovcf2} --missing-to-ref -Oz -o TBN35_YRI108.grch38.merge.00.vcf.gz
bcftools view -h TBN35_YRI108.grch38.merge.00.vcf.gz > header.txt
bcftools view -H TBN35_YRI108.grch38.merge.00.vcf.gz \
| grep -v '\./' \
| grep -v '/\.' \
| grep -v '\.|' \
| grep -v '|\.' \
> body.txt

cat header.txt body.txt \
| bgzip -c > TBN35_YRI108.grch38.merge.00.exMiss.vcf.gz
tabix -p vcf TBN35_YRI108.grch38.merge.00.exMiss.vcf.gz

zcat TBN35_YRI108.grch38.merge.00.exMiss.vcf.gz | grep -v -E 'random|chrUn|alt|HLA|EBV' | sed 's/chr//g' | bgzip > TBN35_YRI108.grch38.merge.00.exMiss.noChr.vcf.gz
tabix -p vcf TBN35_YRI108.grch38.merge.00.exMiss.noChr.vcf.gz
