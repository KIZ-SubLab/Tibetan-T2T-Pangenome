#!/bin/bash
#SBATCH --job-name=prepare_paragraph_chm13
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=wzhcnormal
#SBATCH --output=prepare_paragraph_chm13.%j.out

set -euo pipefail

BASE=/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/output/rm1kb-na-segments/merged_fragments_v3/SV/denove_genotyping_SV_35tbn_5archic_dna
INPUT="${BASE}/5caller_new_clean_renamed.sorted.vcf.gz"
REFERENCE=/work/home/songjiaxu/introgression/fastq/ref/chm13v2.0.fa
PYTHON=/work/home/songjiaxu/anaconda3/envs/gatk/bin/python
SCRIPT="${BASE}/prepare_paragraph_chm13.py"

COMBINED="${BASE}/5caller_paragraph_chm13_fixed.ins_del.vcf.gz"
INS="${BASE}/5caller_paragraph_chm13_fixed.ins.vcf.gz"
DEL="${BASE}/5caller_paragraph_chm13_fixed.del.vcf.gz"
STATS="${BASE}/5caller_paragraph_chm13_fixed.stats.tsv"

for required in "$INPUT" "$REFERENCE" "${REFERENCE}.fai" "$SCRIPT"; do
    if [[ ! -s "$required" ]]; then
        echo "ERROR: missing or empty input: $required" >&2
        exit 1
    fi
done

"$PYTHON" "$SCRIPT" \
    --input "$INPUT" \
    --reference "$REFERENCE" \
    --output "$COMBINED" \
    --stats "$STATS"

bcftools index --threads 8 -f -t "$COMBINED"

bcftools view --threads 8 -i 'INFO/SVTYPE="INS"' -Oz -o "$INS" "$COMBINED"
bcftools index --threads 8 -f -t "$INS"

bcftools view --threads 8 -i 'INFO/SVTYPE="DEL"' -Oz -o "$DEL" "$COMBINED"
bcftools index --threads 8 -f -t "$DEL"

# Strict checks: any failure stops the job.
bcftools index -n "$COMBINED"
bcftools index -n "$INS"
bcftools index -n "$DEL"

if bcftools norm --check-ref w -f "$REFERENCE" "$COMBINED" -Ou 2>&1 \
    >/dev/null | grep -q '^REF_MISMATCH'; then
    echo "ERROR: REF mismatch remains in $COMBINED" >&2
    exit 1
fi

if bcftools query -f '%REF\t%ALT\n' "$COMBINED" \
    | awk 'substr($1,1,1) != substr($2,1,1) {bad++} END {exit(bad>0)}'; then
    :
else
    echo "ERROR: REF/ALT padding-base mismatch remains in $COMBINED" >&2
    exit 1
fi

echo "Finished successfully."
echo "Combined: $COMBINED"
echo "INS:      $INS"
echo "DEL:      $DEL"
echo "Stats:    $STATS"
