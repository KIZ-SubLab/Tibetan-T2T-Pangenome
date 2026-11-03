#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: bash run_relate_clues.sh CONFIG.env" >&2
  exit 2
fi
CONFIG=$(readlink -f "$1")
# shellcheck disable=SC1090
source "$CONFIG"

bash "$BASE/scripts/validate_config.sh" "$CONFIG"

PIPELINE_VERSION=1.0
BCFTOOLS=${BCFTOOLS:-/work/home/songjiaxu/bin/bcftools}
RELATE_ROOT=${RELATE_ROOT:-/work/home/songjiaxu/tools/relate}
RELATE_BIN=$RELATE_ROOT/bin
CLUES_ROOT=${CLUES_ROOT:-/work/home/songjiaxu/tools/CLUES2-Companion/CLUES2}
CONDA_ENV=${CONDA_ENV:-omicverse-new}
MU=${MU:-1.25e-8}
INITIAL_NE=${INITIAL_NE:-10000}
GENERATION_YEARS=${GENERATION_YEARS:-28}
SEED=${SEED:-1}
BRANCH_SAMPLES=${BRANCH_SAMPLES:-100}
CLUES_DF=${CLUES_DF:-200}
CLUES_INTEGRATION_POINTS=${CLUES_INTEGRATION_POINTS:-40}
CLUES_H=${CLUES_H:-0.5}
CLUES_CI=${CLUES_CI:-0.95}
CLUES_WINDOWS=${CLUES_WINDOWS:-"80kyr:2858:80000 160kyr:5715:160000"}
RELATE_MEMORY_GB=${RELATE_MEMORY_GB:-110}
THREADS=${SLURM_CPUS_PER_TASK:-${THREADS:-20}}
RESUME=${RESUME:-NO}

CHROM_NUM=${CHROM#chr}
[[ "$REGION" =~ ^([^:]+):([0-9]+)-([0-9]+)$ ]]
REGION_START=${BASH_REMATCH[2]}
REGION_END=${BASH_REMATCH[3]}
GENETIC_MAP=${GENETIC_MAP:-${GENETIC_MAP_TEMPLATE//@/$CHROM_NUM}}
ANCESTOR_FASTA=${ANCESTOR_FASTA:-${ANCESTOR_TEMPLATE//@/$CHROM_NUM}}

mkdir -p "$OUTPUT_DIR/logs"
if [[ -s "$OUTPUT_DIR/RUN_COMPLETE" ]]; then
  echo "ERROR: completed run already exists: $OUTPUT_DIR" >&2
  exit 2
fi
if [[ -s "$OUTPUT_DIR/run_config_snapshot.env" ]]; then
  [[ "$RESUME" == YES ]] || {
    echo "ERROR: incomplete run directory exists; set RESUME=YES after inspection or use a new ANALYSIS_TAG" >&2
    exit 2
  }
  diff -q <(grep -v '^RESUME=' "$CONFIG") \
    <(grep -v '^RESUME=' "$OUTPUT_DIR/run_config_snapshot.env") >/dev/null || {
    echo "ERROR: resume config differs from the saved snapshot (apart from RESUME)" >&2
    exit 2
  }
else
  cp "$CONFIG" "$OUTPUT_DIR/run_config_snapshot.env"
fi
mkdir -p "$OUTPUT_DIR/inputs" "$OUTPUT_DIR/results" "$OUTPUT_DIR/results/qc"
exec > >(tee -a "$OUTPUT_DIR/logs/pipeline.log") 2> >(tee -a "$OUTPUT_DIR/logs/pipeline.log" >&2)

LOCK_DIR=
on_exit() {
  local rc=$?
  if [[ -n "$LOCK_DIR" && -d "$LOCK_DIR" ]]; then
    rmdir "$LOCK_DIR" 2>/dev/null || true
  fi
  if (( rc != 0 )); then
    printf 'FAILED\t%s\texit=%s\n' "$(date -Is)" "$rc" >> "$OUTPUT_DIR/RUN_FAILED.log"
  fi
}
trap on_exit EXIT

echo "[$(date -Is)] RELATE/CLUES pipeline v$PIPELINE_VERSION"
echo "analysis=$ANALYSIS_TAG target=$REGION focal=$CHROM:$FOCAL_POS"

input_signature=$(
  {
    printf '%s\n' "$PIPELINE_VERSION" "$SOURCE_VCF" "$CHROM" "$SAMPLE_LABEL" \
      "$GENETIC_MAP" "$ANCESTOR_FASTA" "$MU" "$INITIAL_NE" "$SEED" "$RELATE_MEMORY_GB"
    stat -Lc '%n\t%s\t%Y' "$SOURCE_VCF" "$GENETIC_MAP" "$ANCESTOR_FASTA"
    if [[ -n "${SAMPLES_FILE:-}" ]]; then
      stat -Lc '%n\t%s\t%Y' "$SAMPLES_FILE"
      awk 'NF{print $1}' "$SAMPLES_FILE" | sort -u
    else
      "$BCFTOOLS" query -l "$SOURCE_VCF"
    fi
  } | sha256sum | awk '{print $1}'
)
safe_cache_tag=$(printf '%s' "$CACHE_TAG" | tr -c 'A-Za-z0-9_.-' '_')
CACHE_DIR=$BASE/cache/${safe_cache_tag}_${input_signature:0:12}
mkdir -p "$BASE/cache"

INPUT_CACHE=$CACHE_DIR/input
RELATE_CACHE=$CACHE_DIR/relate
DEMOGRAPHY_CACHE=$CACHE_DIR/demography
QC_CACHE=$CACHE_DIR/qc
VCF_PREFIX=$INPUT_CACHE/cohort.$CHROM.qc
QC_VCF=$VCF_PREFIX.vcf.gz
RAW_HAPS=$INPUT_CACHE/raw.haps
RAW_SAMPLE=$INPUT_CACHE/raw.sample
POLARIZED_HAPS=$INPUT_CACHE/hap.haps
UNIQUE_HAPS=$INPUT_CACHE/hap.unique.haps
MAP_FILE=$INPUT_CACHE/$CHROM.map
POPLABELS=$INPUT_CACHE/cohort.haploid.poplabels
RELATE_PREFIX=$RELATE_CACHE/cohort.$CHROM.relate
COAL_PREFIX=$DEMOGRAPHY_CACHE/cohort_haploid

if [[ ! -s "$CACHE_DIR/CACHE_COMPLETE" ]]; then
  mkdir -p "$CACHE_DIR"
  LOCK_DIR=$CACHE_DIR/.build_lock
  if ! mkdir "$LOCK_DIR" 2>/dev/null; then
    echo "ERROR: another job is building this cache: $CACHE_DIR" >&2
    exit 3
  fi
  mkdir -p "$INPUT_CACHE" "$RELATE_CACHE" "$DEMOGRAPHY_CACHE" "$QC_CACHE"
  if [[ -s "$CACHE_DIR/cache_signature.sha256" ]]; then
    [[ $(cat "$CACHE_DIR/cache_signature.sha256") == "$input_signature" ]] || {
      echo "ERROR: existing incomplete cache has a different signature" >&2
      exit 3
    }
  else
    printf '%s\n' "$input_signature" > "$CACHE_DIR/cache_signature.sha256"
  fi

  echo "[$(date -Is)] Building chromosome cache: $CACHE_DIR"
  if [[ ! -s "$QC_VCF" ]]; then
    if [[ -n "${SAMPLES_FILE:-}" ]]; then
      N_EXPECTED_SAMPLES=$(awk 'NF{print $1}' "$SAMPLES_FILE" | sort -u | wc -l)
    else
      N_EXPECTED_SAMPLES=$("$BCFTOOLS" query -l "$SOURCE_VCF" | wc -l)
    fi
    N_EXPECTED_HAPLOTYPES=$((2 * N_EXPECTED_SAMPLES))
    view_args=(-m2 -M2 -v snps -r "$CHROM")
    if [[ -n "${SAMPLES_FILE:-}" ]]; then
      view_args+=(-S "$SAMPLES_FILE")
    fi
    "$BCFTOOLS" view "${view_args[@]}" -Ou "$SOURCE_VCF" \
      | "$BCFTOOLS" +fill-tags -Ou -- -t AC,AN \
      | "$BCFTOOLS" view -i "AN=$N_EXPECTED_HAPLOTYPES && AC>0 && AC<AN" -Oz -o "$QC_VCF"
    "$BCFTOOLS" index -t "$QC_VCF"
  fi
  "$BCFTOOLS" query -l "$QC_VCF" > "$INPUT_CACHE/source.samples.txt"
  N_SAMPLES=$(wc -l < "$INPUT_CACHE/source.samples.txt")
  N_HAPLOTYPES=$((2 * N_SAMPLES))
  (( N_SAMPLES > 1 )) || { echo "ERROR: fewer than two samples after filtering" >&2; exit 3; }

  if "$BCFTOOLS" query -f '[%GT\n]' "$QC_VCF" \
      | awk 'index($0,"/"){found=1} END{exit !found}'; then
    echo "ERROR: unphased genotype remains after filtering" >&2
    exit 3
  fi
  missing_count=$("$BCFTOOLS" query -f '[%GT\n]' "$QC_VCF" | awk 'index($0,"."){n++} END{print n+0}')
  (( missing_count == 0 )) || { echo "ERROR: missing genotype remains after filtering" >&2; exit 3; }

  if [[ ! -s "$RAW_HAPS" || ! -s "$RAW_SAMPLE" ]]; then
    "$RELATE_BIN/RelateFileFormats" --mode ConvertFromVcf -i "$VCF_PREFIX" \
      --haps "$RAW_HAPS" --sample "$RAW_SAMPLE"
  fi
  if [[ ! -s "$INPUT_CACHE/ancestor.$CHROM.fa" ]]; then
    awk -v chrom="$CHROM" 'NR==1{print ">"chrom;next}{print}' "$ANCESTOR_FASTA" \
      > "$INPUT_CACHE/ancestor.$CHROM.fa"
  fi
  if [[ ! -s "$POLARIZED_HAPS" ]]; then
    "$RELATE_BIN/RelateFileFormats" --mode FlipHapsUsingAncestor \
      --haps "$RAW_HAPS" --sample "$RAW_SAMPLE" \
      --ancestor "$INPUT_CACHE/ancestor.$CHROM.fa" -o "$INPUT_CACHE/hap"
  fi
  if [[ ! -s "$UNIQUE_HAPS" ]]; then
    python3 "$BASE/scripts/deduplicate_haps_positions.py" \
      --input "$POLARIZED_HAPS" --output "$UNIQUE_HAPS" \
      --audit "$QC_CACHE/duplicate_position_resolution.tsv"
  fi
  awk 'NR>1 && $3<=previous{exit 1}{previous=$3}' "$UNIQUE_HAPS" || {
    echo "ERROR: polarized haps positions are not strictly increasing" >&2; exit 3;
  }
  if [[ ! -s "$MAP_FILE" ]]; then
    awk 'NF>=4{print $1,$2,$4}' "$GENETIC_MAP" > "$MAP_FILE"
  fi
  [[ $(wc -l < "$MAP_FILE") -gt 1 ]] || { echo "ERROR: converted genetic map is empty" >&2; exit 3; }
  if [[ ! -s "$POPLABELS" ]]; then
    awk -v pop="$SAMPLE_LABEL" 'BEGIN{print "sample population group sex"} NR>2{print $1,pop,pop,1;print $1,pop,pop,1}' \
      "$RAW_SAMPLE" > "$POPLABELS"
  fi
  [[ $(( $(wc -l < "$POPLABELS") - 1 )) -eq "$N_HAPLOTYPES" ]] || {
    echo "ERROR: poplabels/haplotype count mismatch" >&2; exit 3;
  }

  {
    printf 'field\tvalue\n'
    printf 'source_vcf\t%s\nchromosome\t%s\nsamples\t%s\nhaplotypes\t%s\n' \
      "$SOURCE_VCF" "$CHROM" "$N_SAMPLES" "$N_HAPLOTYPES"
    printf 'complete_polymorphic_snvs\t%s\n' "$("$BCFTOOLS" view -H "$QC_VCF" | wc -l)"
    printf 'polarized_before_coordinate_dedup\t%s\n' "$(wc -l < "$POLARIZED_HAPS")"
    printf 'polarized_after_coordinate_dedup\t%s\n' "$(wc -l < "$UNIQUE_HAPS")"
    printf 'duplicate_records_removed\t%s\n' "$(( $(wc -l < "$QC_CACHE/duplicate_position_resolution.tsv") - 1 ))"
  } > "$QC_CACHE/input_qc.tsv"

  if [[ ! -s "$RELATE_PREFIX.anc" || ! -s "$RELATE_PREFIX.mut" ]]; then
    cd "$RELATE_CACHE"
    set +e
    "$RELATE_BIN/Relate" --mode All --haps "$UNIQUE_HAPS" --sample "$RAW_SAMPLE" \
      --map "$MAP_FILE" -m "$MU" -N "$INITIAL_NE" --memory "$RELATE_MEMORY_GB" \
      --seed "$SEED" -o "cohort.$CHROM.relate"
    relate_rc=$?
    set -e
    if [[ ! -s "$RELATE_PREFIX.anc" || ! -s "$RELATE_PREFIX.mut" ]]; then
      echo "ERROR: RELATE failed with status $relate_rc and final anc/mut are absent" >&2
      exit 3
    fi
    if (( relate_rc != 0 )); then
      echo "WARNING: RELATE returned $relate_rc but final anc/mut exist"
    fi
  fi
  if [[ ! -s "$COAL_PREFIX.coal" ]]; then
    "$RELATE_ROOT/scripts/EstimatePopulationSize/EstimatePopulationSize.sh" \
      -i "$RELATE_PREFIX" -m "$MU" --poplabels "$POPLABELS" \
      --threshold 0 --noplot --seed "$SEED" -o "$COAL_PREFIX"
  fi
  [[ -s "$COAL_PREFIX.coal" ]] || { echo "ERROR: population-history .coal is absent" >&2; exit 3; }
  printf 'COMPLETE\t%s\n' "$(date -Is)" > "$CACHE_DIR/CACHE_COMPLETE"
  rmdir "$LOCK_DIR"
  LOCK_DIR=
else
  echo "[$(date -Is)] Reusing chromosome cache: $CACHE_DIR"
  [[ $(cat "$CACHE_DIR/cache_signature.sha256") == "$input_signature" ]] || {
    echo "ERROR: cache signature mismatch" >&2; exit 3;
  }
fi

N_SAMPLES=$(wc -l < "$INPUT_CACHE/source.samples.txt")
N_HAPLOTYPES=$((2 * N_SAMPLES))
for file in "$QC_VCF" "$UNIQUE_HAPS" "$RAW_SAMPLE" "$RELATE_PREFIX.anc" "$RELATE_PREFIX.mut" "$COAL_PREFIX.coal"; do
  [[ -s "$file" ]] || { echo "ERROR: incomplete cache file: $file" >&2; exit 4; }
done

echo "[$(date -Is)] Extracting target region and focal site"
if [[ ! -s "$OUTPUT_DIR/inputs/target_region.qc.vcf.gz" ]]; then
  "$BCFTOOLS" view -r "$REGION" -Oz -o "$OUTPUT_DIR/inputs/target_region.qc.vcf.gz" "$QC_VCF"
  "$BCFTOOLS" index -t "$OUTPUT_DIR/inputs/target_region.qc.vcf.gz"
fi
[[ $("$BCFTOOLS" view -H "$OUTPUT_DIR/inputs/target_region.qc.vcf.gz" | wc -l) -gt 0 ]] || {
  echo "ERROR: target region contains no usable SNPs" >&2; exit 4;
}

HAP_FOCAL=$(awk -v p="$FOCAL_POS" '$3==p{print;n++} END{if(n!=1)exit 1}' "$UNIQUE_HAPS") || {
  echo "ERROR: focal position is absent or duplicated after ancestral polarization/deduplication" >&2
  exit 4
}
read -r _hchr HAP_ID _hpos HAP_ANC HAP_DER _states <<< "$HAP_FOCAL"
if [[ -n "${EXPECTED_FOCAL_ID:-}" && "$HAP_ID" != "$EXPECTED_FOCAL_ID" ]]; then
  echo "ERROR: retained focal ID $HAP_ID differs from EXPECTED_FOCAL_ID=$EXPECTED_FOCAL_ID" >&2
  exit 4
fi

SOURCE_RECORDS=$("$BCFTOOLS" query -r "$CHROM:$FOCAL_POS-$FOCAL_POS" \
  -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' "$QC_VCF")
SOURCE_RECORD=$(printf '%s\n' "$SOURCE_RECORDS" | awk -F '\t' -v id="$HAP_ID" '$3==id{print;n++} END{if(n!=1)exit 1}') || {
  echo "ERROR: cannot match retained haps focal record back to QC VCF" >&2; exit 4;
}
printf '%s\n' "$SOURCE_RECORD" > "$OUTPUT_DIR/inputs/focal_record.source.tsv"
read -r SRC_CHROM SRC_POS SRC_ID SRC_REF SRC_ALT _src_gt <<< "$SOURCE_RECORD"
if [[ -n "${EXPECTED_FOCAL_REF:-}" && "$SRC_REF" != "$EXPECTED_FOCAL_REF" ]]; then
  echo "ERROR: source REF $SRC_REF differs from EXPECTED_FOCAL_REF=$EXPECTED_FOCAL_REF" >&2; exit 4
fi
if [[ -n "${EXPECTED_FOCAL_ALT:-}" && "$SRC_ALT" != "$EXPECTED_FOCAL_ALT" ]]; then
  echo "ERROR: source ALT $SRC_ALT differs from EXPECTED_FOCAL_ALT=$EXPECTED_FOCAL_ALT" >&2; exit 4
fi

RESULTS=$OUTPUT_DIR/results
mkdir -p "$RESULTS/local_genealogy" "$RESULTS/branch_samples" "$RESULTS/figures" "$RESULTS/tree"
LOCAL_PREFIX=$RESULTS/local_genealogy/region
if [[ ! -s "$LOCAL_PREFIX.anc" || ! -s "$LOCAL_PREFIX.mut" ]]; then
  "$RELATE_BIN/RelateExtract" --mode AncMutForSubregion \
    --first_bp "$REGION_START" --last_bp "$REGION_END" \
    --anc "$RELATE_PREFIX.anc" --mut "$RELATE_PREFIX.mut" -o "$LOCAL_PREFIX"
fi

RELATE_FOCAL=$(awk -F';' -v p="$FOCAL_POS" '$2==p{print;n++} END{if(n!=1)exit 1}' "$RELATE_PREFIX.mut") || {
  echo "ERROR: focal mutation is absent or duplicated in RELATE .mut" >&2; exit 4;
}
head -n 1 "$RELATE_PREFIX.mut" > "$RESULTS/focal_mutation_record.tsv"
printf '%s\n' "$RELATE_FOCAL" >> "$RESULTS/focal_mutation_record.tsv"

awk -v p="$FOCAL_POS" '$3==p{for(i=6;i<=NF;i++)print $i}' "$UNIQUE_HAPS" \
  > "$RESULTS/derived_haplotypes.txt"
[[ $(wc -l < "$RESULTS/derived_haplotypes.txt") -eq "$N_HAPLOTYPES" ]] || {
  echo "ERROR: focal derived haplotype vector length mismatch" >&2; exit 4;
}
DERIVED_COUNT=$(awk '{s+=$1} END{print s+0}' "$RESULTS/derived_haplotypes.txt")
DERIVED_FREQ=$(awk -v n="$DERIVED_COUNT" -v d="$N_HAPLOTYPES" 'BEGIN{printf "%.8f",n/d}')
(( DERIVED_COUNT > 0 && DERIVED_COUNT < N_HAPLOTYPES )) || {
  echo "ERROR: focal derived allele is not polymorphic" >&2; exit 4;
}

IFS=';' read -r SNP_INDEX MUT_POS MUT_DIST MUT_ID TREE_INDEX BRANCH_INDEX NOT_MAPPING IS_FLIPPED AGE_BEGIN AGE_END ALLELES _ <<< "$RELATE_FOCAL"
IFS='/' read -r MUT_ANC MUT_DER <<< "$ALLELES"
{
  printf 'chromosome\tposition\tsource_vcf_id\tsource_ref\tsource_alt\tancestral\tderived\tderived_count\tderived_frequency\tn_samples\tn_haplotypes\tage_begin_generations\tage_end_generations\tage_begin_years\tage_end_years\n'
  awk -v c="$CHROM" -v p="$FOCAL_POS" -v id="$MUT_ID" -v sr="$SRC_REF" -v sa="$SRC_ALT" \
    -v a="$MUT_ANC" -v d="$MUT_DER" -v dc="$DERIVED_COUNT" -v f="$DERIVED_FREQ" \
    -v ns="$N_SAMPLES" -v nh="$N_HAPLOTYPES" -v lo="$AGE_BEGIN" -v hi="$AGE_END" -v g="$GENERATION_YEARS" \
    'BEGIN{OFS="\t";print c,p,id,sr,sa,a,d,dc,f,ns,nh,lo,hi,lo*g,hi*g}'
} > "$RESULTS/focal_summary.tsv"

BRANCH_PREFIX=$RESULTS/branch_samples/focal
if [[ ! -s "$BRANCH_PREFIX.newick" ]]; then
  "$RELATE_ROOT/scripts/SampleBranchLengths/SampleBranchLengths.sh" \
    -i "$RELATE_PREFIX" -m "$MU" --coal "$COAL_PREFIX.coal" \
    --num_samples "$BRANCH_SAMPLES" --first_bp "$FOCAL_POS" --last_bp "$FOCAL_POS" \
    --format n --seed "$SEED" -o "$BRANCH_PREFIX"
fi
if [[ ! -s "$RESULTS/focal_times.txt" ]]; then
  conda run -n "$CONDA_ENV" python "$BASE/scripts/RelateToCLUES.fixed.py" \
    --RelateSamples "$BRANCH_PREFIX.newick" \
    --DerivedFile "$RESULTS/derived_haplotypes.txt" --out "$RESULTS/focal"
fi

printf 'window\tt_cutoff_generations\tmax_plot_years\tlogLR\tneg_log10_p\tepoch_start\tepoch_end\tselection_mle\tci_lower\tci_upper\n' \
  > "$RESULTS/result_summary.tsv"
window_labels=()
for spec in $CLUES_WINDOWS; do
  IFS=':' read -r label cutoff max_years <<< "$spec"
  window_labels+=("$label")
  out_dir=$RESULTS/clues2_$label
  mkdir -p "$out_dir"
  prefix=$out_dir/focal
  if [[ ! -s "${prefix}_freqs.txt" || ! -s "${prefix}_post.txt" ]]; then
    conda run -n "$CONDA_ENV" python "$CLUES_ROOT/inference.py" \
      --times "$RESULTS/focal_times.txt" --popFreq "$DERIVED_FREQ" --coal "$COAL_PREFIX.coal" \
      --tCutoff "$cutoff" --CI "$CLUES_CI" --df "$CLUES_DF" \
      --integration_points "$CLUES_INTEGRATION_POINTS" --h "$CLUES_H" --out "$prefix"
  fi
  figure=$RESULTS/figures/focal_CLUES2_$label
  conda run -n "$CONDA_ENV" python "$BASE/scripts/plot_clues2_trajectory.py" \
    --freqs "${prefix}_freqs.txt" --post "${prefix}_post.txt" --out "$figure" \
    --max-years "$max_years" --generation-years "$GENERATION_YEARS" \
    --title "$SAMPLE_LABEL | $MUT_ID | $CHROM:$FOCAL_POS | derived $MUT_DER" \
    --derived-label "$MUT_DER" --sample-label "$SAMPLE_LABEL" \
    --n-samples "$N_SAMPLES" --n-haplotypes "$N_HAPLOTYPES"
  inference=$(awk 'NR==2{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' "${prefix}_inference.txt")
  ci=$(awk 'NR==2{print $2"\t"$3}' "${prefix}_CI.txt")
  printf '%s\t%s\t%s\t%s\t%s\n' "$label" "$cutoff" "$max_years" "$inference" "$ci" \
    >> "$RESULTS/result_summary.tsv"
done

PLOT_PREFIX=$RESULTS/tree/focal_tree
"$RELATE_BIN/RelateTreeView" --mode TreeView --anc "$RELATE_PREFIX.anc" \
  --mut "$RELATE_PREFIX.mut" --snp_of_interest "$SNP_INDEX" -o "$PLOT_PREFIX"
set +e
"$RELATE_BIN/RelateTreeView" --mode BranchesBelowMutation --anc "$RELATE_PREFIX.anc" \
  --mut "$RELATE_PREFIX.mut" --snp_of_interest "$SNP_INDEX" -o "$PLOT_PREFIX.focal"
branch_rc=$?
set -e
if (( branch_rc != 0 )) || [[ ! -s "$PLOT_PREFIX.focal.plotcoords.mut" ]]; then
  printf 'branchID\n' > "$PLOT_PREFIX.focal.plotcoords.mut"
fi
"$RELATE_BIN/RelateTreeView" --mode MutationsOnBranches --anc "$RELATE_PREFIX.anc" \
  --mut "$RELATE_PREFIX.mut" --haps "$UNIQUE_HAPS" --sample "$RAW_SAMPLE" \
  --snp_of_interest "$SNP_INDEX" -o "$PLOT_PREFIX" || true
conda run -n "$CONDA_ENV" python "$BASE/scripts/plot_focal_tree.py" \
  --plotcoords "$PLOT_PREFIX.plotcoords" \
  --mutation-branches "$PLOT_PREFIX.focal.plotcoords.mut" \
  --mut "$RELATE_PREFIX.mut" --chrom "$CHROM" --position "$FOCAL_POS" \
  --years-per-generation "$GENERATION_YEARS" --sample-label "$SAMPLE_LABEL" --out "$PLOT_PREFIX"

cp "$QC_CACHE/input_qc.tsv" "$RESULTS/input_qc.tsv"
{
  printf 'field\tvalue\n'
  printf 'pipeline_version\t%s\nanalysis_tag\t%s\nregion\t%s\nfocal\t%s:%s\n' \
    "$PIPELINE_VERSION" "$ANALYSIS_TAG" "$REGION" "$CHROM" "$FOCAL_POS"
  printf 'chromosome_cache\t%s\ncache_signature\t%s\nsource_vcf\t%s\n' \
    "$CACHE_DIR" "$input_signature" "$SOURCE_VCF"
  printf 'genetic_map\t%s\nancestor_fasta\t%s\nbranch_samples\t%s\n' \
    "$GENETIC_MAP" "$ANCESTOR_FASTA" "$BRANCH_SAMPLES"
  printf 'mu\t%s\ninitial_ne\t%s\ngeneration_years\t%s\nclues_windows\t%s\n' \
    "$MU" "$INITIAL_NE" "$GENERATION_YEARS" "$CLUES_WINDOWS"
} > "$RESULTS/run_manifest.tsv"

conda run -n "$CONDA_ENV" python "$BASE/scripts/validate_results.py" \
  --output-dir "$OUTPUT_DIR" --windows "${window_labels[@]}" \
  --expected-df "$CLUES_DF" --expected-haplotypes "$N_HAPLOTYPES"
if grep -q $'\tFAIL\t' "$RESULTS/validation_report.tsv"; then
  echo "ERROR: validation report contains FAIL" >&2
  exit 5
fi

(
  cd "$OUTPUT_DIR"
  find results -type f ! -name final_checksums.sha256 -print0 \
    | sort -z | xargs -0 sha256sum > results/final_checksums.sha256
)
printf 'COMPLETE\t%s\n' "$(date -Is)" > "$OUTPUT_DIR/RUN_COMPLETE"
echo "DONE: $OUTPUT_DIR"
