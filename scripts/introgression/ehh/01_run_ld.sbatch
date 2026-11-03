#!/usr/bin/env bash
#SBATCH --job-name=introSV_LD
#SBATCH --partition=wzhcnormal
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/output/rm1kb-na-segments/merged_fragments_v3/PBS_sv-snp/TMRCA_dating/introSV_fd_analysis/EHH/logs/01_LD_%j.out
#SBATCH --error=/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/output/rm1kb-na-segments/merged_fragments_v3/PBS_sv-snp/TMRCA_dating/introSV_fd_analysis/EHH/logs/01_LD_%j.err

set -euo pipefail
ROOT=/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/output/rm1kb-na-segments/merged_fragments_v3/PBS_sv-snp/TMRCA_dating/introSV_fd_analysis/EHH
source "$ROOT/config.sh"
mkdir -p "$EHH_ROOT"/{logs,metadata,work,vcf,plink,results}

cut -f1,3 "$SAMPLE_MAP" | tail -n +2 > "$EHH_ROOT/metadata/deepvariant_to_sv326.reheader.tsv"
cut -f3 "$SAMPLE_MAP" | tail -n +2 > "$EHH_ROOT/metadata/cza35.sv326.samples.txt"

"$BCFTOOLS" view --threads "$SLURM_CPUS_PER_TASK" -S "$EHH_ROOT/metadata/cza35.sv326.samples.txt" \
  -i "ID=@$EHH_ROOT/metadata/unique_original_sv_ids.txt" -Ou "$SV_VCF" \
  | "$BCFTOOLS" +fill-tags -Oz -o "$EHH_ROOT/vcf/target_SVs.CZA35.CHM13v2.vcf.gz" -- -t AC,AN,AF
"$BCFTOOLS" index -f -t "$EHH_ROOT/vcf/target_SVs.CZA35.CHM13v2.vcf.gz"
expected=$(grep -cve '^$' "$EHH_ROOT/metadata/unique_original_sv_ids.txt")
observed=$("$BCFTOOLS" index -n "$EHH_ROOT/vcf/target_SVs.CZA35.CHM13v2.vcf.gz")
test "$observed" -eq "$expected"

"$BCFTOOLS" query -f '%ID\t%CHROM\t%POS\n' "$EHH_ROOT/vcf/target_SVs.CZA35.CHM13v2.vcf.gz" \
  | awk -v w="$LD_WINDOW_BP" 'BEGIN{OFS="\t"}{s=$3-w;if(s<1)s=1;print $2,s,$3+w}' \
  | sort -k1,1V -k2,2n \
  | awk 'BEGIN{OFS="\t"} NR==1{c=$1;s=$2;e=$3;next} $1==c && $2<=e+1{if($3>e)e=$3;next} {print c,s,e;c=$1;s=$2;e=$3} END{if(NR)print c,s,e}' \
  > "$EHH_ROOT/metadata/ld_windows_1Mb.tsv"

: > "$EHH_ROOT/work/snp_parts.list"
idx=0
while read -r chrom start end; do
  idx=$((idx+1)); n=${chrom#chr}; src=${SNP_VCF_PATTERN//\{CHR\}/$n}
  part="$EHH_ROOT/vcf/part_$(printf '%03d' "$idx").${chrom}.${start}-${end}.CZA35.CHM13v2.snps.vcf.gz"
  tmp="$EHH_ROOT/work/part_$(printf '%03d' "$idx").subset.vcf.gz"
  renamed="$EHH_ROOT/work/part_$(printf '%03d' "$idx").renamed.vcf.gz"
  "$BCFTOOLS" view --threads "$SLURM_CPUS_PER_TASK" -r "${chrom}:${start}-${end}" -m2 -M2 -v snps -Oz -o "$tmp" "$src"
  "$BCFTOOLS" reheader -s "$EHH_ROOT/metadata/deepvariant_to_sv326.reheader.tsv" -o "$renamed" "$tmp"
  "$BCFTOOLS" annotate --threads "$SLURM_CPUS_PER_TASK" --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -Ou "$renamed" \
    | "$BCFTOOLS" +fill-tags -Oz -o "$part" -- -t AC,AN,AF
  "$BCFTOOLS" index -f -t "$part"
  printf '%s\n' "$part" >> "$EHH_ROOT/work/snp_parts.list"
done < "$EHH_ROOT/metadata/ld_windows_1Mb.tsv"

"$BCFTOOLS" concat --threads "$SLURM_CPUS_PER_TASK" -a -f "$EHH_ROOT/work/snp_parts.list" \
  -Oz -o "$EHH_ROOT/vcf/all_case_windows.CZA35.CHM13v2.biallelic_snps.vcf.gz"
"$BCFTOOLS" index -f -t "$EHH_ROOT/vcf/all_case_windows.CZA35.CHM13v2.biallelic_snps.vcf.gz"

"$BCFTOOLS" query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/AF\n' \
  "$EHH_ROOT/vcf/all_case_windows.CZA35.CHM13v2.biallelic_snps.vcf.gz" \
  | awk 'BEGIN{OFS="\t";print "id","chrom","pos","ref","alt","ac","an","af"}{print}' > "$EHH_ROOT/metadata/snp_metadata.tsv"
"$BCFTOOLS" query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/AF\n' \
  "$EHH_ROOT/vcf/target_SVs.CZA35.CHM13v2.vcf.gz" \
  | awk 'BEGIN{OFS="\t";print "id","chrom","pos","ref","alt","ac","an","af"}{print}' > "$EHH_ROOT/metadata/sv_metadata.tsv"

"$BCFTOOLS" view -Ov "$EHH_ROOT/vcf/target_SVs.CZA35.CHM13v2.vcf.gz" \
  | python3 "$EHH_ROOT/scripts/make_sv_ld_vcf.py" \
  | "$BCFTOOLS" view -Oz -o "$EHH_ROOT/vcf/target_SVs.CZA35.CHM13v2.LDcoded.vcf.gz"
"$BCFTOOLS" index -f -t "$EHH_ROOT/vcf/target_SVs.CZA35.CHM13v2.LDcoded.vcf.gz"

"$PLINK" --vcf "$EHH_ROOT/vcf/all_case_windows.CZA35.CHM13v2.biallelic_snps.vcf.gz" \
  --vcf-half-call missing --memory 22000 --threads "$SLURM_CPUS_PER_TASK" --double-id --allow-extra-chr \
  --make-bed --out "$EHH_ROOT/plink/cza35.snps"
"$PLINK" --vcf "$EHH_ROOT/vcf/target_SVs.CZA35.CHM13v2.LDcoded.vcf.gz" \
  --vcf-half-call missing --memory 22000 --threads "$SLURM_CPUS_PER_TASK" --double-id --allow-extra-chr \
  --make-bed --out "$EHH_ROOT/plink/cza35.svs"
"$PLINK" --bfile "$EHH_ROOT/plink/cza35.snps" \
  --bmerge "$EHH_ROOT/plink/cza35.svs.bed" "$EHH_ROOT/plink/cza35.svs.bim" "$EHH_ROOT/plink/cza35.svs.fam" \
  --memory 22000 --threads "$SLURM_CPUS_PER_TASK" --allow-extra-chr --make-bed --out "$EHH_ROOT/plink/cza35.snps_svs"
"$PLINK" --bfile "$EHH_ROOT/plink/cza35.snps_svs" --allow-extra-chr --memory 22000 --threads "$SLURM_CPUS_PER_TASK" \
  --ld-snp-list "$EHH_ROOT/metadata/unique_original_sv_ids.txt" --r2 gz --ld-window-kb 1000 \
  --ld-window 999999 --ld-window-r2 0 --out "$EHH_ROOT/results/all_SVs.CZA35.1Mb"

echo LD_COMPLETE
