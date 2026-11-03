#!/bin/bash
#SBATCH --job-name plotSPrime
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --partition=wzhcnormal
#SBATCH -e plot.err
#SBATCH -o plot.log
#===============================================================================
# 06.plot.sh — 汇总 SPrime 渗入片段的古人类匹配率，画二维密度图
#
# 输入: moss/out.chr{1..22}.mscore (Python annotate_archaic.py 输出)
# 输出: match.summary.txt + 轮廓图 png
#===============================================================================

set -euo pipefail

SPRIME_DIR="/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/sprime"
MOSSDIR="${SPRIME_DIR}/moss"
OUTDIR="${SPRIME_DIR}/plot"
mkdir -p "${OUTDIR}"

# Step 1: 汇总每个渗入片段的匹配率
echo "=== Step 1: score_summary.r ==="
Rscript "${SPRIME_DIR}/score_summary.r" \
    "${MOSSDIR}" \
    "${OUTDIR}/match.summary.txt"

echo ""
echo "=== 汇总结果前几行 ==="
head -5 "${OUTDIR}/match.summary.txt" 2>/dev/null || echo "(文件为空, 可能还没有 mscore 文件)"

# Step 2: 画二维密度图
echo ""
echo "=== Step 2: plot_contour.r ==="
if [ -f "${OUTDIR}/match.summary.txt" ]; then
    Rscript "${SPRIME_DIR}/plot_contour.r" \
        "${OUTDIR}/match.summary.txt" \
        "${OUTDIR}/tbn.contour"
    echo ""
    ls -lh "${OUTDIR}/"*.png 2>/dev/null || echo "(无图片)"
else
    echo "跳过: 无 match.summary.txt"
fi

echo ""
echo "=== 完成 ==="
