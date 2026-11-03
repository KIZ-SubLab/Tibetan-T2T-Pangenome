#!/usr/bin/env bash
#SBATCH --job-name=trace_ultimate
#SBATCH --partition=wzhcnormal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=02:00:00
#SBATCH --chdir=/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/output/rm1kb-na-segments/merged_fragments_v3/TRACE_find.super.archic
#SBATCH --output=logs/trace_ultimate_%j.out
#SBATCH --error=logs/trace_ultimate_%j.err

# TRACE/ghost/super-archaic one-command Slurm workflow for the 35-individual data set.
#
# User entry point:
#   sbatch 00_TRACE_ULTIMATE.sbatch
#
# This file is both the dispatcher and every worker.  The initial job submits a
# dependency DAG.  After the SINGER manifest exists, a second dispatcher reads
# its exact task count and submits the compute-heavy continuation.  Completed,
# validated products are reused; the script never deletes scientific outputs.

set -euo pipefail

ROOT=/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/output/rm1kb-na-segments/merged_fragments_v3/TRACE_find.super.archic
SELF="$ROOT/00_TRACE_ULTIMATE.sbatch"
PARTITION="${PARTITION:-wzhcnormal}"
SINGER_MAX="${SINGER_MAX:-400}"
MERGE_MAX="${MERGE_MAX:-100}"
EXTRACT_MAX="${EXTRACT_MAX:-400}"
TRACE_MAX="${TRACE_MAX:-70}"
GHOST_MAX="${GHOST_MAX:-22}"
DRY_RUN="${DRY_RUN:-0}"
STAGE="${PIPELINE_STAGE:-MASTER}"

TRACE_ENV=/work/home/songjiaxu/anaconda3/envs/trace
GATK_ENV=/work/home/songjiaxu/anaconda3/envs/gatk
export PATH="$TRACE_ENV/bin:$GATK_ENV/bin:/work/home/songjiaxu/bin:$PATH"
TRACE_PY="$TRACE_ENV/bin/python"
BCFTOOLS="$GATK_ENV/bin/bcftools"
REF_DIR=/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/data/ref_hg38
ANC_DIR=/work/home/songjiaxu/introgression/data_HMM_intro/hg38_ancestral
ARCHAIC_DIR=/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/data/archaic
YRI_DIR=/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/data/yri/by_chr
STRICT_DIR="$ROOT/resources/accessibility_mask/1000G_GRCh38_StrictMask_20160622/by_chrom"
GLOBAL_MANIFEST="$ROOT/singer/manifests/all_autosomes.global_chunks.tsv"
ULTIMATE="$ROOT/ultimate"

mkdir -p "$ROOT/logs" "$ULTIMATE/logs" "$ULTIMATE/submissions" "$ULTIMATE/checkpoints"

log() {
    printf '[%s] [%s] %s\n' "$(date -Is)" "$STAGE" "$*"
}

die() {
    log "ERROR: $*" >&2
    exit 1
}

need_file() {
    [[ -s "$1" ]] || die "missing or empty file: $1"
}

need_exec() {
    [[ -x "$1" ]] || die "missing executable: $1"
}

checkpoint() {
    local name="$1"
    {
        printf 'stage\t%s\n' "$name"
        printf 'status\tPASS\n'
        printf 'completed\t%s\n' "$(date -Is)"
        printf 'job_id\t%s\n' "${SLURM_JOB_ID:-manual}"
        printf 'array_task_id\t%s\n' "${SLURM_ARRAY_TASK_ID:-NA}"
    } > "$ULTIMATE/checkpoints/${name}.PASS.tsv"
}

strip_job_id() {
    local value="$1"
    printf '%s\n' "${value%%;*}"
}

submit_stage() {
    local stage="$1" name="$2" dependency="$3" array_spec="$4"
    local cpus="$5" mem="$6" walltime="$7"
    local -a opts
    opts=(
        --parsable
        --partition="$PARTITION"
        --job-name="$name"
        --nodes=1
        --ntasks=1
        --cpus-per-task="$cpus"
        --mem="$mem"
        --time="$walltime"
        --output="$ULTIMATE/logs/${name}_%A_%a.out"
        --error="$ULTIMATE/logs/${name}_%A_%a.err"
        --export="ALL,PIPELINE_STAGE=${stage},ULTIMATE_RUN_ID=${ULTIMATE_RUN_ID}"
    )
    [[ -n "$dependency" ]] && opts+=(--dependency="afterok:${dependency}")
    [[ -n "$array_spec" ]] && opts+=(--array="$array_spec")
    if [[ "$DRY_RUN" == 1 ]]; then
        printf 'DRYRUN_%s' "$stage"
        return 0
    fi
    strip_job_id "$(sbatch "${opts[@]}" "$SELF")"
}

record_job() {
    local stage="$1" job_id="$2" dependency="$3" array_spec="$4"
    printf '%s\t%s\t%s\t%s\n' "$stage" "$job_id" "${dependency:-NONE}" "${array_spec:-NONE}" \
        >> "$RUN_DIR/submitted_jobs.tsv"
    log "submitted $stage: job=$job_id dependency=${dependency:-NONE} array=${array_spec:-NONE}"
}

run_master() {
    ULTIMATE_RUN_ID="${ULTIMATE_RUN_ID:-$(date +%Y%m%dT%H%M%S)_${SLURM_JOB_ID:-manual}}"
    export ULTIMATE_RUN_ID
    RUN_DIR="$ULTIMATE/submissions/$ULTIMATE_RUN_ID"
    mkdir -p "$RUN_DIR"
    printf 'stage\tjob_id\tdependency\tarray\n' > "$RUN_DIR/submitted_jobs.tsv"
    printf '%s\n' "$ULTIMATE_RUN_ID" > "$ULTIMATE/LATEST_RUN_ID"
    ln -sfn "$RUN_DIR" "$ULTIMATE/latest"

    log "building initial dependency chain; run_id=$ULTIMATE_RUN_ID"
    local j_preflight j_prep j_lift j_input j_manifest j_dispatch
    j_preflight="$(submit_stage PREFLIGHT tr_u_preflight '' '' 1 2G 02:00:00)"
    record_job PREFLIGHT "$j_preflight" '' ''
    j_prep="$(submit_stage PREP tr_u_prep "$j_preflight" '' 4 8G 12:00:00)"
    record_job PREP "$j_prep" "$j_preflight" ''
    j_lift="$(submit_stage LIFTOVER tr_u_lift "$j_prep" "1-22%22" 4 8G 12:00:00)"
    record_job LIFTOVER "$j_lift" "$j_prep" '1-22%22'
    j_input="$(submit_stage SINGER_INPUT tr_u_sinput "$j_lift" "1-22%22" 1 2G 04:00:00)"
    record_job SINGER_INPUT "$j_input" "$j_lift" '1-22%22'
    j_manifest="$(submit_stage MANIFEST tr_u_manifest "$j_input" '' 2 8G 04:00:00)"
    record_job MANIFEST "$j_manifest" "$j_input" ''
    j_dispatch="$(submit_stage DISPATCH_AFTER_MANIFEST tr_u_dispatch "$j_manifest" '' 1 2G 02:00:00)"
    record_job DISPATCH_AFTER_MANIFEST "$j_dispatch" "$j_manifest" ''

    if [[ "$DRY_RUN" == 1 ]]; then
        log "dry-run complete; heavy stages are submitted by DISPATCH_AFTER_MANIFEST"
    else
        log "initial DAG submitted. Later job IDs will be appended by job $j_dispatch"
        log "job ledger: $RUN_DIR/submitted_jobs.tsv"
    fi
}

run_dispatch_after_manifest() {
    RUN_DIR="$ULTIMATE/submissions/${ULTIMATE_RUN_ID:?ULTIMATE_RUN_ID missing}"
    need_file "$GLOBAL_MANIFEST"
    local chunks
    chunks="$(( $(wc -l < "$GLOBAL_MANIFEST") - 1 ))"
    [[ "$chunks" -gt 0 ]] || die "global SINGER manifest has no tasks"
    log "submitting continuation for $chunks exact SINGER chunks"

    local j_singer j_vsinger j_merge j_vmerge j_extract j_vextract
    local j_infer j_summary j_trace j_ghost j_gagg j_gsnp j_rg2 j_post j_final
    j_singer="$(submit_stage SINGER tr_u_singer '' "0-$((chunks - 1))%${SINGER_MAX}" 2 4G 54:00:00)"
    record_job SINGER "$j_singer" '' "0-$((chunks - 1))%${SINGER_MAX}"
    j_vsinger="$(submit_stage VERIFY_SINGER tr_u_vsinger "$j_singer" '' 1 4G 04:00:00)"
    record_job VERIFY_SINGER "$j_vsinger" "$j_singer" ''
    j_merge="$(submit_stage MERGE tr_u_merge "$j_vsinger" "0-1099%${MERGE_MAX}" 5 10G 24:00:00)"
    record_job MERGE "$j_merge" "$j_vsinger" "0-1099%${MERGE_MAX}"
    j_vmerge="$(submit_stage VERIFY_MERGE tr_u_vmerge "$j_merge" '' 1 4G 04:00:00)"
    record_job VERIFY_MERGE "$j_vmerge" "$j_merge" ''
    j_extract="$(submit_stage EXTRACT tr_u_extract "$j_vmerge" "0-1099%${EXTRACT_MAX}" 8 16G 18:00:00)"
    record_job EXTRACT "$j_extract" "$j_vmerge" "0-1099%${EXTRACT_MAX}"
    j_vextract="$(submit_stage VERIFY_EXTRACT tr_u_vextract "$j_extract" '' 1 4G 04:00:00)"
    record_job VERIFY_EXTRACT "$j_vextract" "$j_extract" ''
    j_infer="$(submit_stage INFER tr_u_infer "$j_vextract" "0-69%${TRACE_MAX}" 4 8G 36:00:00)"
    record_job INFER "$j_infer" "$j_vextract" "0-69%${TRACE_MAX}"
    j_summary="$(submit_stage SUMMARIZE tr_u_summary "$j_infer" "0-69%${TRACE_MAX}" 2 4G 12:00:00)"
    record_job SUMMARIZE "$j_summary" "$j_infer" "0-69%${TRACE_MAX}"
    j_trace="$(submit_stage AGGREGATE_TRACE tr_u_trace "$j_summary" '' 4 8G 12:00:00)"
    record_job AGGREGATE_TRACE "$j_trace" "$j_summary" ''
    j_ghost="$(submit_stage GHOST_CHR tr_u_ghost "$j_trace" "1-22%${GHOST_MAX}" 4 12G 24:00:00)"
    record_job GHOST_CHR "$j_ghost" "$j_trace" "1-22%${GHOST_MAX}"
    j_gagg="$(submit_stage GHOST_AGGREGATE tr_u_gagg "$j_ghost" '' 2 8G 04:00:00)"
    record_job GHOST_AGGREGATE "$j_gagg" "$j_ghost" ''
    j_gsnp="$(submit_stage GHOST_SNPS tr_u_gsnp "$j_gagg" '' 2 12G 12:00:00)"
    record_job GHOST_SNPS "$j_gsnp" "$j_gagg" ''
    j_rg2="$(submit_stage RULEGROUP2_SUPERDEEP tr_u_rg2 "$j_gagg" '' 4 12G 24:00:00)"
    record_job RULEGROUP2_SUPERDEEP "$j_rg2" "$j_gagg" ''
    j_post="$(submit_stage SUPER_POST tr_u_spost "$j_rg2" '' 4 16G 24:00:00)"
    record_job SUPER_POST "$j_post" "$j_rg2" ''
    j_final="$(submit_stage FINAL_VALIDATE tr_u_final "${j_gsnp}:${j_post}" '' 2 8G 06:00:00)"
    record_job FINAL_VALIDATE "$j_final" "${j_gsnp}:${j_post}" ''
    log "full DAG submitted; terminal validation job=$j_final"
}

run_preflight() {
    need_exec "$TRACE_PY"
    need_exec "$BCFTOOLS"
    need_exec "$GATK_ENV/bin/gatk"
    need_exec "$TRACE_ENV/bin/trace-extract"
    need_exec "$TRACE_ENV/bin/trace-infer"
    need_exec "$TRACE_ENV/bin/trace-summarize"
    need_exec "$ROOT/software/SINGER-0.1.8/singer-0.1.8-beta-linux-x86_64/singer_master"
    command -v sbatch >/dev/null || die "sbatch is not on PATH"
    need_exec "$TRACE_ENV/bin/bedtools"
    need_file "$ROOT/vcf/TBN35-v2-mc-chm13.vcfbub.a100k.wave.vcf.gz"
    need_file "$ROOT/vcf/TBN35-v2-mc-chm13.vcfbub.a100k.wave.vcf.gz.tbi"
    need_file "$ROOT/vcf/chm13v2-grch38.chain"
    for chr in $(seq 1 22); do
        need_file "$REF_DIR/chr${chr}.fa"
        need_file "$ANC_DIR/homo_sapiens_ancestor_${chr}.fa"
        need_file "$STRICT_DIR/strictmask_chr${chr}.bed"
        need_file "$ROOT/resources/genetic_maps/GRCh38_1KGP/trace_hapmap4/genetic_map_hg38_chr${chr}.txt"
        need_file "$ARCHAIC_DIR/individuals_highcov.${chr}.bcf"
        [[ -s "$ARCHAIC_DIR/individuals_highcov.${chr}.bcf.csi" || \
           -s "$ARCHAIC_DIR/individuals_highcov.${chr}.bcf.tbi" ]] \
            || die "archaic BCF index missing for chr${chr}"
        need_file "$YRI_DIR/chr${chr}.vcf.gz"
        [[ -s "$YRI_DIR/chr${chr}.vcf.gz.tbi" || -s "$YRI_DIR/chr${chr}.vcf.gz.csi" ]] \
            || die "YRI VCF index missing for chr${chr}"
    done
    for script in \
        01_filter_missing5_fill0.sbatch 02_liftover_biallelic_snp.sbatch \
        03_prepare_singer_vcf.sbatch generate_singer_chunk_manifests.py \
        04_singer_chunks.sbatch 06_merge_singer_arg.sbatch 07_trace_extract.sbatch \
        08_trace_infer.sbatch 09_trace_summarize.sbatch 10_aggregate_trace.sbatch \
        prepare_trace_analysis_masks.sh make_sample_node_map.py \
        ghost_classification/run_ghost_chr.sbatch ghost_classification/run_aggregate.sbatch \
        ghost_classification/scripts/23_export_rg1_ghost_snps.py \
        super_archaic/19_run_paper_rulegroup2.sbatch \
        super_archaic/20_map_strict_DEN_to_host_background.py \
        super_archaic/21_plot_chr22_strict_DEN_case.py \
        super_archaic/22_export_strict_superarchaic_SNPs.py; do
        need_file "$ROOT/$script"
    done
    {
        printf 'tool\tversion_or_path\n'
        printf 'pipeline\tTRACE_ULTIMATE_v1\n'
        printf 'root\t%s\n' "$ROOT"
        printf 'trace_python\t%s\n' "$($TRACE_PY --version 2>&1)"
        printf 'bcftools\t%s\n' "$($BCFTOOLS --version | head -1)"
        printf 'gatk\t%s\n' "$($GATK_ENV/bin/gatk --version 2>&1 | tr '\n' ' ')"
        printf 'singer\t%s\n' "$ROOT/software/SINGER-0.1.8/singer-0.1.8-beta-linux-x86_64/singer_master"
        printf 'checked\t%s\n' "$(date -Is)"
    } > "$ULTIMATE/software_and_inputs.tsv"
    checkpoint PREFLIGHT
}

run_prep() {
    local out="$ROOT/preprocess/TBN35-v2-mc-chm13.vcfbub.a100k.wave.35ind.miss5pct.fill0.vcf.gz"
    local manifest="${out%.vcf.gz}.manifest.txt"
    if [[ -s "$out" && -s "$out.tbi" && -s "$manifest" ]]; then
        [[ "$($BCFTOOLS query -l "$out" | wc -l)" -eq 35 ]] || die "prepared VCF sample count is not 35"
        [[ "$($BCFTOOLS index -n "$out")" -gt 0 ]] || die "prepared VCF index reports no records"
        grep -qx 'output_samples=35' "$manifest" || die "prepared VCF manifest sample QA missing"
        grep -q '^fill_rule=replace_each_missing_GT_allele_dot_with_reference_allele_0$' "$manifest" \
            || die "prepared VCF fill-rule manifest mismatch"
        log "validated existing prepared VCF"
    else
        [[ ! -e "$out" && ! -e "$out.tbi" ]] \
            || die "partial prepared output exists; inspect it before rerunning: $out"
        bash "$ROOT/01_filter_missing5_fill0.sbatch"
    fi
    checkpoint PREP
}

run_liftover() {
    bash "$ROOT/02_liftover_biallelic_snp.sbatch"
    checkpoint "LIFTOVER_chr${SLURM_ARRAY_TASK_ID}"
}

run_singer_input() {
    bash "$ROOT/03_prepare_singer_vcf.sbatch"
    checkpoint "SINGER_INPUT_chr${SLURM_ARRAY_TASK_ID}"
}

run_manifest() {
    "$TRACE_PY" "$ROOT/generate_singer_chunk_manifests.py" "$ROOT" --chunk-size 2000000 --min-snps 10
    bash "$ROOT/prepare_trace_analysis_masks.sh"
    "$TRACE_PY" - "$ROOT" "$GLOBAL_MANIFEST" <<'PY'
import csv, os, sys
from pathlib import Path
root, out = Path(sys.argv[1]), Path(sys.argv[2])
tmp = out.with_name(out.name + f'.tmp.{os.getpid()}')
global_task = 0
with tmp.open('w', newline='') as fo:
    w = csv.writer(fo, delimiter='\t', lineterminator='\n')
    w.writerow(['global_task_id', 'chromosome', 'local_task_id', 'start', 'end',
                'snp_count', 'vcf_prefix', 'output_prefix'])
    for chrom_i in range(1, 23):
        p = root / 'singer/manifests' / f'chr{chrom_i}.chunks.tsv'
        with p.open(newline='') as fi:
            for row in csv.DictReader(fi, delimiter='\t'):
                w.writerow([global_task, row['chromosome'], row['task_id'], row['start'],
                            row['end'], row['snp_count'], row['vcf_prefix'], row['output_prefix']])
                global_task += 1
if global_task < 1000:
    raise RuntimeError(f'Unexpectedly few SINGER chunks: {global_task}')
os.replace(tmp, out)
print(f'global SINGER tasks={global_task}')
PY
    checkpoint MANIFEST
}

run_singer() {
    need_file "$GLOBAL_MANIFEST"
    local global_task row chromosome local_task
    global_task="${SLURM_ARRAY_TASK_ID:?array task missing}"
    row="$(awk -F '\t' -v n="$((global_task + 2))" 'NR==n {print; exit}' "$GLOBAL_MANIFEST")"
    [[ -n "$row" ]] || die "no global manifest row for task $global_task"
    row="${row%$'\r'}"
    IFS=$'\t' read -r _ chromosome local_task _ <<< "$row"
    export CHR="${chromosome#chr}"
    export SLURM_ARRAY_TASK_ID="$local_task"
    bash "$ROOT/04_singer_chunks.sbatch"
}

run_verify_singer() {
    local expected observed
    expected="$(( $(wc -l < "$GLOBAL_MANIFEST") - 1 ))"
    observed="$(find "$ROOT/singer/chunks" -type f -name completed.tsv | wc -l)"
    [[ "$observed" -eq "$expected" ]] || die "completed SINGER chunks $observed != $expected"
    while IFS=$'\t' read -r _ _ _ _ _ _ _ prefix; do
        [[ "$prefix" == output_prefix ]] && continue
        need_file "$(dirname "$prefix")/completed.tsv"
        for posterior in 250 299; do
            need_file "${prefix}_nodes_${posterior}.txt"
            need_file "${prefix}_branches_${posterior}.txt"
            need_file "${prefix}_muts_${posterior}.txt"
        done
    done < "$GLOBAL_MANIFEST"
    {
        printf 'metric\tvalue\nexpected_chunks\t%s\ncompleted_chunks\t%s\n' "$expected" "$observed"
        printf 'retained_posteriors\t250-299\nstatus\tPASS\ncompleted\t%s\n' "$(date -Is)"
    } > "$ROOT/singer/singer_completion_manifest.tsv"
    checkpoint VERIFY_SINGER
}

map_posterior_task() {
    local global_task="${SLURM_ARRAY_TASK_ID:?array task missing}"
    export CHR="$((global_task / 50 + 1))"
    export SLURM_ARRAY_TASK_ID="$((global_task % 50))"
}

run_merge() {
    map_posterior_task
    bash "$ROOT/06_merge_singer_arg.sbatch"
}

run_verify_merge() {
    local tsz manifests
    tsz="$(find "$ROOT/singer/merged" -type f -name '*.tsz' | wc -l)"
    manifests="$(find "$ROOT/singer/merged" -type f -name '*.tsz.manifest.json' | wc -l)"
    [[ "$tsz" -eq 1100 && "$manifests" -eq 1100 ]] \
        || die "merged ARG mismatch: tsz=$tsz manifests=$manifests expected=1100"
    "$TRACE_PY" "$ROOT/make_sample_node_map.py" \
        "$ROOT/singer/input/chr1/chr1.grch38.biallelic_snp.ancestral_ref.vcf" \
        "$ROOT/trace/sample_node_map.tsv"
    [[ "$(wc -l < "$ROOT/trace/sample_node_map.tsv")" -eq 71 ]] || die "node map QA failed"
    checkpoint VERIFY_MERGE
}

run_extract() {
    map_posterior_task
    bash "$ROOT/07_trace_extract.sbatch"
}

run_verify_extract() {
    local npz manifests list chrom posterior path
    npz="$(find "$ROOT/trace/extract" -type f -name '*.npz' | wc -l)"
    manifests="$(find "$ROOT/trace/extract" -type f -name '*.manifest.tsv' | wc -l)"
    [[ "$npz" -eq 1100 && "$manifests" -eq 1100 ]] \
        || die "TRACE extraction mismatch: npz=$npz manifests=$manifests expected=1100"
    mkdir -p "$ROOT/trace/data_files"
    for chrom in $(seq 1 22); do
        list="$ROOT/trace/data_files/chr${chrom}.posterior_npz.txt"
        : > "$list"
        for posterior in $(seq 250 299); do
            path="$ROOT/trace/extract/chr${chrom}/chr${chrom}.t15000.strictmask.posterior_${posterior}.all70.npz"
            need_file "$path"
            printf '%s\n' "$path" >> "$list"
        done
        [[ "$(wc -l < "$list")" -eq 50 ]] || die "posterior list QA failed for chr${chrom}"
    done
    checkpoint VERIFY_EXTRACT
}

run_infer() {
    bash "$ROOT/08_trace_infer.sbatch"
}

run_summarize() {
    local hap="${SLURM_ARRAY_TASK_ID:?array task missing}" pad dir
    pad="$(printf '%02d' "$hap")"
    dir="$ROOT/trace/summary/hap_${pad}"
    if [[ -s "$dir/completed.tsv" && "$(find "$dir" -maxdepth 1 -name '*.summary.txt' | wc -l)" -eq 22 ]]; then
        log "validated existing TRACE summary for haplotype $hap"
    else
        bash "$ROOT/09_trace_summarize.sbatch"
    fi
}

run_aggregate_trace() {
    bash "$ROOT/10_aggregate_trace.sbatch"
    need_file "$ROOT/trace/final/TRACE_introgressed_tracts.all_haplotypes.tsv.gz"
    need_file "$ROOT/trace/final/TRACE_final_QA.json"
    need_file "$ROOT/trace/final/TRACE_introgressed_tracts_overview.pdf"
    checkpoint AGGREGATE_TRACE
}

run_ghost_chr() {
    bash "$ROOT/ghost_classification/run_ghost_chr.sbatch"
}

run_ghost_aggregate() {
    bash "$ROOT/ghost_classification/run_aggregate.sbatch"
    "$TRACE_PY" "$ROOT/ghost_classification/scripts/validate_final.py" \
        --root "$ROOT" \
        --output-dir "$ROOT/ghost_classification"
    need_file "$ROOT/ghost_classification/final/TRACE_archaic_classification.all_autosomes.tsv.gz"
    need_file "$ROOT/ghost_classification/final/VALIDATION.json"
    checkpoint GHOST_AGGREGATE
}

run_ghost_snps() {
    "$TRACE_PY" "$ROOT/ghost_classification/scripts/23_export_rg1_ghost_snps.py" --root "$ROOT"
    need_file "$ROOT/ghost_classification/final/sites_RG1_Ghost/SUCCESS"
    need_file "$ROOT/ghost_classification/final/sites_RG1_Ghost/TRACE_Ghost.ND00_candidate_SNPs.unique_sites.vcf.gz"
    need_file "$ROOT/ghost_classification/final/sites_RG1_Ghost/TRACE_Ghost.ND00_candidate_SNPs.unique_sites.vcf.gz.tbi"
    checkpoint GHOST_SNPS
}

run_rulegroup2_superdeep() {
    bash "$ROOT/super_archaic/19_run_paper_rulegroup2.sbatch"
    need_file "$ROOT/super_archaic/paper_rulegroup2/backgrounds/TRACE_archaic_classification.RuleGroup2.all_autosomes.tsv.gz"
    need_file "$ROOT/super_archaic/paper_rulegroup2/final_paper_t2/validation.json"
    need_file "$ROOT/super_archaic/paper_rulegroup2/final_strict_DEN/summary.json"
    checkpoint RULEGROUP2_SUPERDEEP
}

run_super_post() {
    "$TRACE_PY" "$ROOT/super_archaic/20_map_strict_DEN_to_host_background.py" "$ROOT"
    "$TRACE_PY" "$ROOT/super_archaic/21_plot_chr22_strict_DEN_case.py" "$ROOT"
    "$TRACE_PY" "$ROOT/super_archaic/22_export_strict_superarchaic_SNPs.py" "$ROOT"
    need_file "$ROOT/super_archaic/paper_rulegroup2/final_strict_DEN/TRACE_strict_superarchaic_DEN.green_with_orange_host.per_haplotype.tsv.gz"
    need_file "$ROOT/super_archaic/paper_rulegroup2/final_strict_DEN/sites/TRACE_strict_superarchaic_DEN.SNPs.summary.json"
    need_file "$ROOT/super_archaic/paper_rulegroup2/final_strict_DEN/figures/TRACE_strict_superarchaic_DEN.chr22_47p8Mb.case_7samples.pdf"
    checkpoint SUPER_POST
}

run_final_validate() {
    local final="$ULTIMATE/final"
    mkdir -p "$final"
    "$TRACE_PY" - "$ROOT" "$final" <<'PY'
import csv, gzip, json, os, sys
from collections import Counter
from pathlib import Path

root, out = Path(sys.argv[1]), Path(sys.argv[2])
required = {
    'TRACE tracts': root/'trace/final/TRACE_introgressed_tracts.all_haplotypes.tsv.gz',
    'TRACE QA': root/'trace/final/TRACE_final_QA.json',
    'RG1 classification': root/'ghost_classification/final/TRACE_archaic_classification.all_autosomes.tsv.gz',
    'RG1 validation': root/'ghost_classification/final/VALIDATION.json',
    'Ghost SNP summary': root/'ghost_classification/final/sites_RG1_Ghost/TRACE_Ghost_SNP_export.summary.json',
    'Ghost SNP VCF': root/'ghost_classification/final/sites_RG1_Ghost/TRACE_Ghost.ND00_candidate_SNPs.unique_sites.vcf.gz',
    'RG2 classification': root/'super_archaic/paper_rulegroup2/backgrounds/TRACE_archaic_classification.RuleGroup2.all_autosomes.tsv.gz',
    'super-deep validation': root/'super_archaic/paper_rulegroup2/final_paper_t2/validation.json',
    'strict DEN summary': root/'super_archaic/paper_rulegroup2/final_strict_DEN/summary.json',
    'strict DEN SNP summary': root/'super_archaic/paper_rulegroup2/final_strict_DEN/sites/TRACE_strict_superarchaic_DEN.SNPs.summary.json',
    'chr22 case QA': root/'super_archaic/paper_rulegroup2/final_strict_DEN/figures/TRACE_strict_superarchaic_DEN.chr22_47p8Mb.case_7samples.qa.json',
}
for label, path in required.items():
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f'missing {label}: {path}')

def count_gz_rows(path):
    with gzip.open(path, 'rt') as f:
        return max(0, sum(1 for _ in f)-1)

rg1 = Counter()
with gzip.open(required['RG1 classification'], 'rt') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        rg1[row['classification']] += 1
rg2 = Counter()
with gzip.open(required['RG2 classification'], 'rt') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        rg2[row['classification']] += 1
ghost = json.loads(required['Ghost SNP summary'].read_text())
strict = json.loads(required['strict DEN summary'].read_text())
superdeep = json.loads(required['super-deep validation'].read_text())
summary = {
    'status': 'PASS',
    'pipeline': 'TRACE_ULTIMATE_v1',
    'coordinate_system': 'GRCh38 autosomes',
    'samples': 35,
    'haplotypes': 70,
    'singer_chunks': sum(1 for _ in open(root/'singer/manifests/all_autosomes.global_chunks.tsv'))-1,
    'singer_retained_posteriors': '250-299 (50)',
    'TRACE_introgressed_tract_records': count_gz_rows(required['TRACE tracts']),
    'RuleGroup1_classification_counts': dict(rg1),
    'RuleGroup2_classification_counts': dict(rg2),
    'Ghost_coordinate_union_bp': ghost['ghost_coordinate_union_bp'],
    'Ghost_ND00_candidate_unique_SNPs': ghost['nd00_candidate_unique_sites'],
    'strict_DEN_superarchaic_haplotype_records': strict['strict_DEN_segments_per_haplotype'],
    'strict_DEN_superarchaic_synthetic_segments': strict['strict_DEN_segments_synthetic'],
    'strict_DEN_superarchaic_union_bp': strict['strict_DEN_superarchaic_union_bp'],
    'superdeep_validation_status': superdeep.get('status'),
    'validation_errors': [],
    'excluded_deprecated_branches': [
        't31500 HMM rerun (deprecated; super-deep is detected from t15000 local t2)',
        'old Rule Group 1 strict-superarchaic outputs',
        'smoke/example and repair-only jobs',
        'CMS/locuszoom intersection outside TRACE core (allele-direction audit required)',
    ],
}
(out/'PIPELINE_VALIDATION.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False)+'\n')
print(json.dumps(summary, indent=2, ensure_ascii=False))
PY

    while IFS=$'\t' read -r link target; do
        [[ "$link" == link_name ]] && continue
        need_file "$ROOT/$target"
        ln -sfn "$ROOT/$target" "$final/$link"
    done <<'EOF'
link_name	target
01_TRACE_tracts.tsv.gz	trace/final/TRACE_introgressed_tracts.all_haplotypes.tsv.gz
02_TRACE_overview.pdf	trace/final/TRACE_introgressed_tracts_overview.pdf
03_RG1_all_classifications.tsv.gz	ghost_classification/final/TRACE_archaic_classification.all_autosomes.tsv.gz
04_RG1_burden_per_individual.tsv	ghost_classification/final/TRACE_archaic_burden.per_individual.tsv
05_Ghost_ND00_candidate_SNPs.vcf.gz	ghost_classification/final/sites_RG1_Ghost/TRACE_Ghost.ND00_candidate_SNPs.unique_sites.vcf.gz
05_Ghost_ND00_candidate_SNPs.vcf.gz.tbi	ghost_classification/final/sites_RG1_Ghost/TRACE_Ghost.ND00_candidate_SNPs.unique_sites.vcf.gz.tbi
06_Ghost_ND00_candidate_SNPs.tsv.gz	ghost_classification/final/sites_RG1_Ghost/TRACE_Ghost.ND00_candidate_SNPs.unique_sites.tsv.gz
07_RG2_all_classifications.tsv.gz	super_archaic/paper_rulegroup2/backgrounds/TRACE_archaic_classification.RuleGroup2.all_autosomes.tsv.gz
08_superdeep_NEA_DEN_per_haplotype.tsv.gz	super_archaic/paper_rulegroup2/final_paper_t2/TRACE_superdeep_t2_gt31500.per_haplotype.tsv.gz
09_superdeep_NEA_DEN_overall.tsv	super_archaic/paper_rulegroup2/final_paper_t2/TRACE_superdeep_t2.synthetic_overall.tsv
10_Figure4C_like.pdf	super_archaic/paper_rulegroup2/final_paper_t2/figures/Figure4C_like_t2_gt31500_NEA_DEN.pdf
11_strict_DEN_superarchaic.tsv.gz	super_archaic/paper_rulegroup2/final_strict_DEN/TRACE_strict_superarchaic_DEN.per_haplotype.tsv.gz
12_strict_DEN_orange_hosts.tsv.gz	super_archaic/paper_rulegroup2/final_strict_DEN/TRACE_strict_superarchaic_DEN.orange_host_intervals.per_haplotype.tsv.gz
13_strict_DEN_superarchaic_SNPs.vcf.gz	super_archaic/paper_rulegroup2/final_strict_DEN/sites/TRACE_strict_superarchaic_DEN.SNPs.branch_mutations.unique_sites.vcf.gz
13_strict_DEN_superarchaic_SNPs.vcf.gz.tbi	super_archaic/paper_rulegroup2/final_strict_DEN/sites/TRACE_strict_superarchaic_DEN.SNPs.branch_mutations.unique_sites.vcf.gz.tbi
14_chr22_case_7samples.pdf	super_archaic/paper_rulegroup2/final_strict_DEN/figures/TRACE_strict_superarchaic_DEN.chr22_47p8Mb.case_7samples.pdf
15_methods_guide.md	docs/GUIDE_TRACE_SUPERARCHAIC_BEGINNER.zh-CN.md
EOF

    printf 'category\tfile_or_directory\tmeaning\n' > "$final/RESULTS_INDEX.tsv"
    printf 'TRACE\ttrace/final/\tTRACE t=15000片段、个体/单倍型负担与总览图\n' >> "$final/RESULTS_INDEX.tsv"
    printf 'RG1\tghost_classification/final/\t普通NEA/DEN/Ghost Rule Group 1来源分类\n' >> "$final/RESULTS_INDEX.tsv"
    printf 'Ghost_SNP\tghost_classification/final/sites_RG1_Ghost/\tGhost区域SNP及ND00候选分支突变；二者不可混称\n' >> "$final/RESULTS_INDEX.tsv"
    printf 'RG2\tsuper_archaic/paper_rulegroup2/backgrounds/\t加入YRI后、用于super-archaic背景的论文Rule Group 2分类\n' >> "$final/RESULTS_INDEX.tsv"
    printf 'Superdeep\tsuper_archaic/paper_rulegroup2/final_paper_t2/\tNEA/DEN背景内t2>31500的super-deep候选和Figure 4C类图\n' >> "$final/RESULTS_INDEX.tsv"
    printf 'Strict_DEN\tsuper_archaic/paper_rulegroup2/final_strict_DEN/\t严格DEN super-archaic、橙色宿主、位点与chr22案例图\n' >> "$final/RESULTS_INDEX.tsv"
    printf 'Validation\tultimate/final/PIPELINE_VALIDATION.json\t终极动态计数和验收状态\n' >> "$final/RESULTS_INDEX.tsv"
    printf 'PASS\n' > "$final/SUCCESS"
    checkpoint FINAL_VALIDATE
    log "all requested analyses completed and validated: $final"
}

case "$STAGE" in
    MASTER) run_master ;;
    DISPATCH_AFTER_MANIFEST) run_dispatch_after_manifest ;;
    PREFLIGHT) run_preflight ;;
    PREP) run_prep ;;
    LIFTOVER) run_liftover ;;
    SINGER_INPUT) run_singer_input ;;
    MANIFEST) run_manifest ;;
    SINGER) run_singer ;;
    VERIFY_SINGER) run_verify_singer ;;
    MERGE) run_merge ;;
    VERIFY_MERGE) run_verify_merge ;;
    EXTRACT) run_extract ;;
    VERIFY_EXTRACT) run_verify_extract ;;
    INFER) run_infer ;;
    SUMMARIZE) run_summarize ;;
    AGGREGATE_TRACE) run_aggregate_trace ;;
    GHOST_CHR) run_ghost_chr ;;
    GHOST_AGGREGATE) run_ghost_aggregate ;;
    GHOST_SNPS) run_ghost_snps ;;
    RULEGROUP2_SUPERDEEP) run_rulegroup2_superdeep ;;
    SUPER_POST) run_super_post ;;
    FINAL_VALIDATE) run_final_validate ;;
    *) die "unknown PIPELINE_STAGE=$STAGE" ;;
esac
