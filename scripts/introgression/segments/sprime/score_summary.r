#!/usr/bin/env Rscript
#===============================================================================
# score_summary.r — 汇总每个渗入片段的古人类匹配率
#
# 输入: moss/ 目录下所有 .mscore 文件
# 输出: match.summary.txt
#
# 逻辑:
#   对每个 SEGMENT(渗入片段), 统计各古人类的可比对数(notcomp排除)
#   匹配率 = 1 - mismatch/(match+mismatch)
#   需 ≥10 可比对位点, 否则 NA
#===============================================================================

args = commandArgs(trailingOnly=TRUE)
if(length(args) < 2) {
    cat("Usage: Rscript score_summary.r <mscore_dir> <out_file>\n")
    quit(status=1)
}

mscore_dir = args[1]
out_file    = args[2]

minn = 10  # 最少可比对位点数

# tag: 该列的标签名，match 可以是 "match" 或 tag 名本身
calmatchrate = function(x, tag, minn) {
    is_match    = (x == "match") | (x == tag)
    is_mismatch = (x == "mismatch")
    total = sum(is_match | is_mismatch)
    if(total >= minn) {
        mis = sum(is_mismatch)
        return(round(1 - mis/total, digits=4))
    } else {
        return(NA)
    }
}

files = list.files(path=mscore_dir, pattern="*.mscore$", full.names=TRUE)
cat(sprintf("找到 %d 个 mscore 文件\n", length(files)))

out = c()
for(f in files) {
    da = read.csv(file=f, header=TRUE, sep="\t")
    if(nrow(da) == 0) next

    chr = as.character(da[1, 1])
    seglist = unique(da$SEGMENT)
    n_cols = ncol(da)
    archaic_tags = colnames(da)[9:n_cols]

    for(seg in seglist) {
        idx = da$SEGMENT == seg
        n_variants = sum(idx)
        if(n_variants < minn) next  # 整个片段位点太少, 跳过

        seg_pos = da$POS[idx]
        x = c(chr, seg, min(seg_pos), max(seg_pos), n_variants)

        for(j in 9:n_cols) {
            tag = colnames(da)[j]
            x = c(x, calmatchrate(da[idx, j], tag, minn))
        }
        out = rbind(out, x)
    }
}

if(is.null(out) || nrow(out) == 0) {
    cat("没有有效片段!\n")
    quit(status=0)
}

colnames(out) = c("chr", "seg", "from", "to", "n_variants", archaic_tags)
out = out[order(out[,1], as.numeric(out[,3])), ]
write.table(out, out_file, quote=FALSE, row.names=FALSE, sep="\t")

cat(sprintf("输出: %s (%d 片段)\n", out_file, nrow(out)))
cat(sprintf("列名: %s\n", paste(colnames(out), collapse=", ")))
