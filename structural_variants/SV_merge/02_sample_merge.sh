ls /path/to/5caller2ref/01-merge/C*/mix_merge/bcftools_plus_survivor.sort.vcf |sort > vcf.list
wc -l vcf.list
[ -s vcf.list ] || exit -1

#date
caller_merge='mix_merge'
popu_merge='survivor_1oAll'
indi_support=1
target_region='chm13v2_telo_cent.complement.bed'
SURVIVOR merge vcf.list 200 $indi_support 1 0 0 50 ${caller_merge}.sv.${popu_merge}.vcf
echo "Merge done..."


# sort vcf and bgzip, index
bcftools sort -T ./tmp ${caller_merge}.sv.${popu_merge}.vcf -o ${caller_merge}.sv.${popu_merge}.sort.vcf
bgzip -c ${caller_merge}.sv.${popu_merge}.sort.vcf > ${caller_merge}.sv.${popu_merge}.sort.vcf.gz
tabix ${caller_merge}.sv.${popu_merge}.sort.vcf.gz
echo "sort done..."

# stat
SURVIVOR stats CN1v1.0_${cohort}_${caller_merge}.sv.all.sort.vcf 50 -1 -1 CN1v1.0_${cohort}_${caller_merge}.sv.all.sort.vcf.gz.sum

# exclude records in centromere and telomere
bcftools view -T $target_region --threads 8 -o ${caller_merge}.sv.${popu_merge}.sort.ctRM.vcf ${caller_merge}.sv.${popu_merge}.sort.vcf.gz
bgzip -c ${caller_merge}.sv.${popu_merge}.sort.ctRM.vcf > ${caller_merge}.sv.${popu_merge}.sort.ctRM.vcf.gz
tabix ${caller_merge}.sv.${popu_merge}.sort.ctRM.vcf.gz
truvari collapse -r 500 -p 0 -P 0.5 -s 50 -S 100000 -k common -i ${caller_merge}.sv.${popu_merge}.sort.ctRM.vcf.gz -o ${caller_merge}.sv.${popu_merge}.sort.ctRM.truvari_collapse.vcf -c ${caller_merge}.sv.${popu_merge}.sort.ctRM.vcf.collapsed

# exclude INV 
grep -v 'INV' ${caller_merge}.sv.${popu_merge}.sort.ctRM.truvari_collapse.vcf > ${caller_merge}.sv.${popu_merge}.sort.ctRM.truvari_collapse.invRM.vcf
grep -v '^chr[XYM]' ${caller_merge}.sv.${popu_merge}.sort.ctRM.truvari_collapse.invRM.vcf > ${caller_merge}.sv.${popu_merge}.sort.ctRM.truvari_collapse.invRM.autosome.vcf
date


