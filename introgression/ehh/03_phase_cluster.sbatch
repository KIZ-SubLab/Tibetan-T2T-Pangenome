#!/usr/bin/env bash
#SBATCH --job-name=introSV_phase
#SBATCH --partition=wzhcnormal
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=00:10:00
#SBATCH --output=/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/output/rm1kb-na-segments/merged_fragments_v3/PBS_sv-snp/TMRCA_dating/introSV_fd_analysis/EHH/logs/03_phase_%A_%a.out
#SBATCH --error=/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/output/rm1kb-na-segments/merged_fragments_v3/PBS_sv-snp/TMRCA_dating/introSV_fd_analysis/EHH/logs/03_phase_%A_%a.err

set -euo pipefail
ROOT=/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/output/rm1kb-na-segments/merged_fragments_v3/PBS_sv-snp/TMRCA_dating/introSV_fd_analysis/EHH
source "$ROOT/config.sh"
manifest="$EHH_ROOT/metadata/phase_clusters.tsv"
row=$(awk -F '\t' -v i="$SLURM_ARRAY_TASK_ID" 'NR>1 && $1==i{print;exit}' "$manifest")
test -n "$row"
IFS=$'\t' read -r idx cluster chrom cs ce hs he proxy_ids sv_ids <<< "$row"
n=${chrom#chr}
src=${SNP_VCF_PATTERN//\{CHR\}/$n}
map=${GENETIC_MAP_PATTERN//\{CHR\}/$n}
ref=${REFERENCE_PATTERN//\{CHR\}/$n}
work="$EHH_ROOT/work/phase/$cluster"
out="$EHH_ROOT/vcf/phased/${cluster}.CZA35.GRCh38.phased.vcf.gz"
mkdir -p "$work" "$EHH_ROOT/vcf/phased"

"$BCFTOOLS" view -r "${chrom}:${cs}-${ce}" -m2 -M2 -v snps -Ov "$src" \
  | python3 "$EHH_ROOT/scripts/strip_to_gt.py" --sample-map "$EHH_ROOT/metadata/deepvariant_to_sv326.reheader.tsv" \
  | "$BCFTOOLS" annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -Ou \
  | "$BCFTOOLS" +fill-tags -Oz -o "$work/input.CHM13.gt.vcf.gz" -- -t AC,AN,AF
"$BCFTOOLS" index -f -t "$work/input.CHM13.gt.vcf.gz"

export JAVA_HOME=/work/home/songjiaxu/anaconda3/envs/gatk/lib/jvm
export PATH=/work/home/songjiaxu/anaconda3/envs/gatk/bin:${PATH}
"$GATK" --java-options "-Xmx8g" LiftoverVcf --INPUT "$work/input.CHM13.gt.vcf.gz" --OUTPUT "$work/lifted.raw.vcf.gz" \
  --CHAIN "$CHAIN" --REJECT "$work/lifted.rejected.vcf.gz" --REFERENCE_SEQUENCE "$HG38_REF" \
  --RECOVER_SWAPPED_REF_ALT true --WRITE_ORIGINAL_POSITION true --WARN_ON_MISSING_CONTIG true --CREATE_INDEX true

"$BCFTOOLS" view -r "${chrom}:${hs}-${he}" -m2 -M2 -v snps -Ou "$work/lifted.raw.vcf.gz" \
  | "$BCFTOOLS" norm -d all -Ou \
  | "$BCFTOOLS" +fill-tags -Ou -- -t AC,AN,AF,F_MISSING \
  | "$BCFTOOLS" view -i 'INFO/AC>0 && INFO/AC<INFO/AN && INFO/F_MISSING=0' -Ob -o "$work/target.unphased.bcf"
"$BCFTOOLS" index -f "$work/target.unphased.bcf"
"$BCFTOOLS" view -r "${chrom}:${hs}-${he}" -Ou "$ref" \
  | "$BCFTOOLS" +fill-tags -Ob -o "$work/reference.with_tags.bcf" -- -t AC,AN,AF
"$BCFTOOLS" index -f "$work/reference.with_tags.bcf"

for id in ${proxy_ids//,/ }; do
  test "$("$BCFTOOLS" query -i "ID==\"$id\"" -f '%ID\n' "$work/target.unphased.bcf" | wc -l)" -eq 1
done

"$SHAPEIT_COMMON" --input "$work/target.unphased.bcf" --reference "$work/reference.with_tags.bcf" --map "$map" \
  --region "${chrom}:${hs}-${he}" --output "$work/common.phased.bcf" --thread "$SLURM_CPUS_PER_TASK" \
  --seed "$SHAPEIT_SEED" --log "$EHH_ROOT/logs/${cluster}.phase_common.log"
"$BCFTOOLS" index -f "$work/common.phased.bcf"
"$SHAPEIT_RARE" --input "$work/target.unphased.bcf" --input-region "${chrom}:${hs}-${he}" \
  --scaffold "$work/common.phased.bcf" --scaffold-region "${chrom}:${hs}-${he}" --map "$map" \
  --output "$work/full.phased.bcf" --thread "$SLURM_CPUS_PER_TASK" --seed "$SHAPEIT_SEED" \
  --log "$EHH_ROOT/logs/${cluster}.phase_rare.log"
"$BCFTOOLS" index -f "$work/full.phased.bcf"
"$BCFTOOLS" view -Oz -o "$out" "$work/full.phased.bcf"
"$BCFTOOLS" index -f -t "$out"
for id in ${proxy_ids//,/ }; do
  "$BCFTOOLS" query -i "ID==\"$id\"" -f '%ID[\t%GT]\n' "$out" \
    | awk 'BEGIN{ok=0}{ok=1;for(i=2;i<=NF;i++)if($i!~/^[01]\|[01]$/)exit 2}END{if(!ok)exit 3}'
done
printf '%s\t%s\t%s\t%s\t%s\n' "$cluster" "$chrom" \
  "$("$BCFTOOLS" index -n "$work/target.unphased.bcf")" \
  "$("$BCFTOOLS" index -n "$work/common.phased.bcf")" \
  "$("$BCFTOOLS" index -n "$work/full.phased.bcf")" > "$EHH_ROOT/results/${cluster}.phase_QC.tsv"
