#!/bin/bash

if [[ $# != 1   ]];then
	echo "Usage: $0 list"
	exit 1
fi

list=$1
wk=`pwd`
# get data, current_done
asm_caller=/path/to/5caller2ref/asm2ref/
sv=/path/to/5caller2ref/hifi2ref/
bin=/path/to/5caller2ref/bin

cat $list|grep -v '^#' | while read id
do
	[ -d $id ] || mkdir $id
	cd $id

	[ -d bcftools ] || mkdir bcftools
	cd bcftools
	bcftools merge -m none -o merge.bcftools.vcf ${asm_caller}/pav/pav_${id}.vcf.gz ${asm_caller}/svimasm/svimasm_${id}.vcf.gz ${asm_caller}/lgvar/lgvar_${id}.vcf.gz $sv/pbsv/pbsv_${id}.vcf.gz $sv/sniffles/sniffles_${id}.vcf.gz --force-samples
	# 3 samples -> one sample, choose best GT
	python3 $bin/consensusSV_from_bcftoolsMerge.py merge.bcftools.vcf merge.bcftools.uniform.vcf && echo "$id bcftools merge done"
	cd ..

	[ -d survivor ] || mkdir survivor
	cd survivor
	# merge 2 out 5 callers
	for tool in pav svimasm lgvar;do
		echo $asm_caller/$tool/${tool}_${id}.vcf >> vcf.list
	done
	for tool in pbsv sniffles;do
		echo $sv/$tool/${tool}_${id}.vcf >> vcf.list
	done
	SURVIVOR merge vcf.list 200 2 1 0 0 50 merge.survivor.vcf && echo "$id done"
	python3 $bin/correctGT_for_survivorSV.py merge.survivor.vcf merge.survivor.uniform.vcf
	cd ..

	[ -d mix_merge ] || mkdir mix_merge
	cd mix_merge
	python3 $bin/unionSVsets.py ../bcftools/merge.bcftools.uniform.vcf ../survivor/merge.survivor.uniform.vcf > bcftools_plus_survivor.vcf
	bcftools sort bcftools_plus_survivor.vcf -o bcftools_plus_survivor.sort.vcf

	cd ..

	cd $wk
done

