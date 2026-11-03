#!/usr/bin/env bash
# Read-only source audit for one ghost SV. Temporary regional VCFs are written
# only under the compute-node scratch directory and removed on exit.
# Required: GHOST_TABLE, REPORT_DIR, HEATBASE, DEN_VCF_ROOT, NEA_VCF_ROOT.
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2500M
#SBATCH --time=00:20:00

set -eo pipefail
: "${GHOST_TABLE:?GHOST_TABLE is required}"
: "${REPORT_DIR:?REPORT_DIR is required}"
: "${HEATBASE:?HEATBASE is required}"
: "${DEN_VCF_ROOT:?DEN_VCF_ROOT is required}"
: "${NEA_VCF_ROOT:?NEA_VCF_ROOT is required}"
: "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"

line="$(awk -F '\t' -v n="${SLURM_ARRAY_TASK_ID}" 'NR == n + 1 {print; exit}' "${GHOST_TABLE}")"
[[ -n "${line}" ]] || { printf 'Missing input row %s\n' "${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
IFS=$'\t' read -r original_id chr start end sv_type rest <<< "${line}"
mid=$(((start + end) / 2))
left=$((mid - 200000)); (( left < 1 )) && left=1
right=$((mid + 200000))
region="${chr}:${left}-${right}"

source /work/home/songjiaxu/anaconda3/etc/profile.d/conda.sh
conda activate plink
set -u

tmp_base="${SLURM_TMPDIR:-/tmp}"
tmp_dir="$(mktemp -d "${tmp_base%/}/ghost_check.${SLURM_ARRAY_TASK_ID}.XXXXXX")"
cleanup_tmp() {
  case "${tmp_dir}" in
    "${tmp_base%/}"/ghost_check."${SLURM_ARRAY_TASK_ID}".*)
      rm -f -- \
        "${tmp_dir}/chb.vcf.gz" "${tmp_dir}/chb.vcf.gz.tbi" \
        "${tmp_dir}/tbn.vcf.gz" "${tmp_dir}/tbn.vcf.gz.tbi" \
        "${tmp_dir}/den.vcf.gz" "${tmp_dir}/den.vcf.gz.tbi" \
        "${tmp_dir}/nea.vcf.gz" "${tmp_dir}/nea.vcf.gz.tbi"
      rmdir -- "${tmp_dir}" 2>/dev/null || true
      ;;
  esac
}
trap cleanup_tmp EXIT

chb_src="${HEATBASE}/1KGP/CHB/1KGP.${chr}.CHB.biallelic_snp.vcf.gz"
tbn_src="${HEATBASE}/TBN/TBN35.chm13.biSNV.vcf.gz"
den_src="${DEN_VCF_ROOT}/${chr}.vcf.gz"
nea_src="${NEA_VCF_ROOT}/${chr}.vcf.gz"

bcftools view -r "${region}" -m2 -M2 -v snps -Oz -o "${tmp_dir}/chb.vcf.gz" "${chb_src}"
bcftools view -r "${region}" -m2 -M2 -v snps -Oz -o "${tmp_dir}/tbn.vcf.gz" "${tbn_src}"
bcftools view -r "${region}" -m2 -M2 -v snps \
  -i 'GT!="mis" && FMT/DP>=2 && FMT/GQ>=20' \
  -Oz -o "${tmp_dir}/den.vcf.gz" "${den_src}"
bcftools view -r "${region}" -m2 -M2 -v snps \
  -i 'GT!="mis" && FMT/DP>=2 && FMT/GQ>=20' \
  -Oz -o "${tmp_dir}/nea.vcf.gz" "${nea_src}"
for f in "${tmp_dir}"/*.vcf.gz; do bcftools index -f -t "${f}"; done

modern_exact="$({ bcftools isec -n=2 -w1 -c none \
  "${tmp_dir}/chb.vcf.gz" "${tmp_dir}/tbn.vcf.gz" -Ou \
  | bcftools view -m2 -M2 -v snps -H | wc -l; } | tr -d ' ')"
fourway_callable="$({ bcftools isec -n=4 -w1 -c none \
  "${tmp_dir}/chb.vcf.gz" "${tmp_dir}/tbn.vcf.gz" \
  "${tmp_dir}/den.vcf.gz" "${tmp_dir}/nea.vcf.gz" -Ou \
  | bcftools view -m2 -M2 -v snps -H | wc -l; } | tr -d ' ')"
den_callable="$(bcftools index -n "${tmp_dir}/den.vcf.gz")"
nea_callable="$(bcftools index -n "${tmp_dir}/nea.vcf.gz")"

mkdir -p "${REPORT_DIR}/parts"
{
  printf 'task_id\toriginal_id\tchr\tstart\tend\tmid\tregion\tmodern_exact\tden_callable\tnea_callable\tfourway_callable\n'
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${SLURM_ARRAY_TASK_ID}" "${original_id}" "${chr}" "${start}" "${end}" \
    "${mid}" "${region}" "${modern_exact}" "${den_callable}" \
    "${nea_callable}" "${fourway_callable}"
} > "${REPORT_DIR}/parts/${SLURM_ARRAY_TASK_ID}.tsv"
