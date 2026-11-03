#!/usr/bin/env bash
set -euo pipefail

WORK_ROOT=/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/output/rm1kb-na-segments/merged_fragments_v3/PBS_sv-snp/TMRCA_dating/locuszoom_cms_fd_daf
RUN_NAME=chr2_45900000_46500000_DEN_PBS_fD
RUN_ROOT="$WORK_ROOT/$RUN_NAME"
PIPELINE="$WORK_ROOT/reusable_pipeline"
PBS_FILE="$WORK_ROOT/snp.all.sites.pbs.PBS_gt0.tsv.gz"
FD_WINDOWS=/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/output/rm1kb-na-segments/merged_fragments_v3/PBS_sv-snp/TMRCA_dating/EPAS1/fd_analysis/comparisons/DEN/results/windows_100SNP_step50.tsv
REFERENCE_CACHE=/work/home/songjiaxu/jianxin_cms/caculate_results/cms_corr_4_parament/daf_SGPL1_single/relate_time_395/96/north_locuszoom_cms/_reference_cache
CHROM=chr2
START=45900000
END=46500000

if [[ -e "$RUN_ROOT" ]]; then
  echo "Refusing to overwrite existing run directory: $RUN_ROOT" >&2
  exit 2
fi
for file in "$PBS_FILE" "$FD_WINDOWS" "$REFERENCE_CACHE/gencode.v50.basic.annotation.gtf.gz" "$REFERENCE_CACHE/encodeCcreCombined_ENCODE3.bb"; do
  [[ -s "$file" ]] || { echo "Missing required input: $file" >&2; exit 2; }
done

mkdir -p "$RUN_ROOT"/{data,figures,metadata,qa,refs,scripts}
ln -s "$REFERENCE_CACHE/gencode.v50.basic.annotation.gtf.gz" "$RUN_ROOT/refs/gencode.v50.basic.annotation.gtf.gz"
ln -s "$REFERENCE_CACHE/encodeCcreCombined_ENCODE3.bb" "$RUN_ROOT/refs/encodeCcreCombined_ENCODE3.bb"
cp "$PIPELINE"/*.py "$RUN_ROOT/scripts/"
cp "$PIPELINE/run_locuszoom_cms_fd_daf.sh" "$RUN_ROOT/scripts/"

python3 "$RUN_ROOT/scripts/01_prepare_pbs_fd_tables.py" "$RUN_ROOT" \
  --pbs-file "$PBS_FILE" --fd-windows "$FD_WINDOWS" \
  --gtf "$RUN_ROOT/refs/gencode.v50.basic.annotation.gtf.gz" \
  --chrom "$CHROM" --start "$START" --end "$END" \
  > "$RUN_ROOT/metadata/prepare.stdout.log" 2> "$RUN_ROOT/metadata/prepare.stderr.log"

python3 "$RUN_ROOT/scripts/fetch_ucsc_layered_h3k27ac.py" "$RUN_ROOT" \
  --chrom "$CHROM" --start "$START" --end "$END" \
  > "$RUN_ROOT/metadata/h3k27ac.stdout.log" 2> "$RUN_ROOT/metadata/h3k27ac.stderr.log"

python3 "$RUN_ROOT/scripts/02_plot_pbs_fd_locus.py" "$RUN_ROOT" \
  --chrom "$CHROM" --start "$START" --end "$END" \
  > "$RUN_ROOT/metadata/plot.stdout.log" 2> "$RUN_ROOT/metadata/plot.stderr.log"

python3 "$RUN_ROOT/scripts/03_final_qa.py" "$RUN_ROOT" \
  --chrom "$CHROM" --start "$START" --end "$END" \
  > "$RUN_ROOT/qa/qa.stdout.log" 2> "$RUN_ROOT/qa/qa.stderr.log"

sha256sum "$PBS_FILE" "$FD_WINDOWS" "$REFERENCE_CACHE/gencode.v50.basic.annotation.gtf.gz" \
  "$REFERENCE_CACHE/encodeCcreCombined_ENCODE3.bb" > "$RUN_ROOT/metadata/input_sha256.txt"

cat > "$RUN_ROOT/README.zh-CN.md" <<'EOF'
# chr2:45,900,000–46,500,000 PBS–DEN fD 区域整合图

主图以 PBS 为纵坐标；同一个 PBS 值同时控制点颜色（区域低值为灰色，区域高值为橙色）。背景蓝色柱为已完成 DEN 分析中的 `fD(EAS,TBN;Denisovan,YRI)`，采用 100-SNP 窗口和 50-SNP 步长。原始 fD 分析从 46,030,415 bp 开始，图中更靠左的区段不填补、不外推。

下方依次保留 ENCODE3 combined cCRE、UCSC hg38 Layered H3K27ac 七细胞系和 GENCODE v50 basic 基因轨道。输出包括可编辑 SVG、PDF、300 dpi PNG 和 600 dpi LZW TIFF。

复现：`bash reusable_pipeline/run_locuszoom_cms_fd_daf.sh`。入口脚本拒绝覆盖已有结果目录；如需重跑，请使用新的 RUN_NAME。
EOF

echo "$RUN_ROOT"
