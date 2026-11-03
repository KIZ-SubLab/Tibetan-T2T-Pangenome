#!/usr/bin/env Rscript
#===============================================================================
# plot_contour.r — 古人类匹配率二维密度图
#
# 输入:  score_summary.r 输出的 match.summary.txt
# 输出: 两张图:
#   tbn.contour_Altai_vs_Denisova.png   (尼人 vs 丹人)
#   tbn.contour_Denisova_vs_Denisova25.png  (两个丹人)
#===============================================================================

library(MASS)

args = commandArgs(trailingOnly=TRUE)
if(length(args) < 2) {
    cat("Usage: Rscript plot_contour.r <summary_file> <output_prefix>\n")
    quit(status=1)
}

infile  = args[1]
outpref = args[2]

da = read.csv(file=infile, header=TRUE, sep="\t")
cat(sprintf("读取 %d 个片段, 列: %s\n", nrow(da), paste(colnames(da), collapse=", ")))

# 提取各古人类匹配率列
# 新 pipeline 列名: AltaiNean, Denisova, Denisova25
# 旧 pipeline 兼容: AltaiNean, AltaiDeni
col_names = colnames(da)

# 找尼人列和丹人列
nean_col = grep("AltaiNean|Nean|nean", col_names)[1]
den_col  = grep("Denisova$|AltaiDeni|Deni", col_names)[1]
den25_col = grep("Denisova25", col_names)[1]

cat(sprintf("尼人列: %s (%d), 丹人列: %s (%d), 丹人25列: %s (%d)\n",
    col_names[nean_col], nean_col,
    col_names[den_col], den_col,
    col_names[den25_col], den25_col))

plot_contour_pair = function(x, y, xlab, ylab, filename) {
    ok = !is.na(x) & !is.na(y)
    x = x[ok]; y = y[ok]
    if(length(x) < 10) {
        cat(sprintf("跳过 %s: 有效点太少 (%d)\n", filename, length(x)))
        return()
    }

    png(filename, width=1100, height=1100, res=250)
    contour(
        kde2d(x, y, n=100, lims=c(0,1,0,1)),
        levels=c(seq(0.3, 0.9, 0.1), seq(1, 10, 1)),
        xaxs="i", yaxs="i", xlab=xlab, ylab=ylab,
        main="", las=1, cex.lab=1.3, cex.axis=1.2,
        cex.main=1.5, lty=5, labcex=0.8
    )
    grid(lty=3, col="gray")
    dev.off()
    cat(sprintf("输出: %s (%d 点)\n", filename, length(x)))
}

# 图1: 尼人 vs 丹人 (Denisova)
if(!is.na(nean_col) && !is.na(den_col)) {
    plot_contour_pair(
        da[, nean_col], da[, den_col],
        paste("Match to Neanderthal (", col_names[nean_col], ")", sep=""),
        paste("Match to Denisovan (", col_names[den_col], ")", sep=""),
        paste(outpref, "_Nean_vs_Denisova.png", sep="")
    )
}

# 图2: 尼人 vs 丹人25
if(!is.na(nean_col) && !is.na(den25_col)) {
    plot_contour_pair(
        da[, nean_col], da[, den25_col],
        paste("Match to Neanderthal (", col_names[nean_col], ")", sep=""),
        paste("Match to Denisova25", sep=""),
        paste(outpref, "_Nean_vs_Denisova25.png", sep="")
    )
}

# 图3: 两个丹人互相对比
if(!is.na(den_col) && !is.na(den25_col)) {
    plot_contour_pair(
        da[, den_col], da[, den25_col],
        paste("Match to Denisova", sep=""),
        "Match to Denisova25",
        paste(outpref, "_Denisova_vs_Denisova25.png", sep="")
    )
}

cat("全部完成!\n")
