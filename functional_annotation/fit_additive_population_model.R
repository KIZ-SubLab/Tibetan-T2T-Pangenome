#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

parse_args <- function(x) {
  out <- list(min_cpm = 1, min_samples = 5L, fdr = 0.05)
  i <- 1L
  while (i <= length(x)) {
    if (i == length(x)) stop("Missing value for ", x[[i]])
    key <- x[[i]]; value <- x[[i + 1L]]
    if (key == "--matrix") out$matrix <- value
    else if (key == "--assay") out$assay <- value
    else if (key == "--outdir") out$outdir <- value
    else if (key == "--min-cpm") out$min_cpm <- as.numeric(value)
    else if (key == "--min-samples") out$min_samples <- as.integer(value)
    else if (key == "--fdr") out$fdr <- as.numeric(value)
    else stop("Unknown option: ", key)
    i <- i + 2L
  }
  required <- c("matrix", "assay", "outdir")
  missing <- required[!vapply(required, function(k) !is.null(out[[k]]), logical(1))]
  if (length(missing)) stop("Missing arguments: ", paste(missing, collapse = ", "))
  if (!out$assay %in% c("atac", "rna")) stop("--assay must be atac or rna")
  if (!is.finite(out$min_cpm) || out$min_cpm < 0 ||
      is.na(out$min_samples) || out$min_samples < 1 || out$min_samples > 50 ||
      !is.finite(out$fdr) || out$fdr <= 0 || out$fdr > 1) stop("Invalid threshold")
  out
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(opt$matrix)) stop("Count matrix missing: ", opt$matrix)

tab <- read.delim(opt$matrix, check.names = FALSE, stringsAsFactors = FALSE,
                  quote = "", comment.char = "")
sample_cols <- grep("^[AW]-[0-9]+-[1-5]$", names(tab), value = TRUE)
if (length(sample_cols) != 50L) stop("Expected 50 A/W samples; found ", length(sample_cols))
if (anyDuplicated(names(tab))) stop("Duplicate input column names")

time_order <- c("0h", "6h", "1d", "3d", "5d")
parts <- strsplit(sample_cols, "-", fixed = TRUE)
meta <- data.frame(
  sample = sample_cols,
  population_code = vapply(parts, `[[`, character(1), 1L),
  individual = vapply(parts, function(z) paste(z[[1]], z[[2]], sep = "-"), character(1)),
  time = time_order[as.integer(vapply(parts, `[[`, character(1), 3L))],
  stringsAsFactors = FALSE
)
meta$population <- factor(ifelse(meta$population_code == "A", "Tibetan", "Han"),
                          levels = c("Han", "Tibetan"))
meta$time <- factor(meta$time, levels = time_order)
if (!all(table(meta$population, meta$time) == 5L) ||
    !all(table(meta$individual) == 5L) || length(unique(meta$individual)) != 10L) {
  stop("Expected five individuals per population and five timepoints per individual")
}
if (!all(table(meta$individual, meta$time) == 1L)) {
  stop("Each individual must have exactly one observation per timepoint")
}
write.table(meta, file.path(opt$outdir, "sample_metadata.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)

feature_candidates <- if (opt$assay == "atac") c("peak_id", "feature_id") else c("gene_id", "feature_id")
feature_id_col <- feature_candidates[feature_candidates %in% names(tab)][1]
if (is.na(feature_id_col)) stop("No feature ID column found")
feature_cols <- setdiff(names(tab), sample_cols)
feature_meta <- tab[, feature_cols, drop = FALSE]
feature_ids <- as.character(tab[[feature_id_col]])
if (anyNA(feature_ids) || any(!nzchar(trimws(feature_ids))) || anyDuplicated(feature_ids)) {
  stop("Feature IDs must be nonmissing and unique; do not silently rename genes or peaks")
}
feature_meta$feature_id <- feature_ids

counts <- as.matrix(tab[, sample_cols, drop = FALSE])
storage.mode(counts) <- "numeric"
if (any(!is.finite(counts))) stop("Counts contain missing or nonfinite values")
if (any(counts < 0)) stop("Counts contain negative values")
if (any(colSums(counts) <= 0)) stop("Every sample must have positive library size")
if (opt$assay == "atac" && any(counts != floor(counts))) {
  stop("ATAC input must contain raw integer fragment counts, not CPM/RPKM")
}
rownames(counts) <- feature_ids

# The primary question is the average Tibetan-Han difference across all five
# timepoints. Time is a categorical adjustment; no interaction is fitted.
design <- model.matrix(~ population + time, data = meta)
rownames(design) <- meta$sample
population_coef <- "populationTibetan"
time_coefs <- grep("^time", colnames(design))
if (!population_coef %in% colnames(design) || length(time_coefs) != 4L) {
  stop("Unexpected additive design columns: ", paste(colnames(design), collapse = ","))
}

y <- DGEList(counts = counts)
keep <- rowSums(cpm(y) >= opt$min_cpm) >= opt$min_samples
if (!any(keep)) stop("No features passed CPM filtering")
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")
kept_meta <- feature_meta[match(rownames(y), feature_meta$feature_id), , drop = FALSE]

v0 <- voom(y, design, plot = FALSE)
corfit <- duplicateCorrelation(v0, design, block = meta$individual)
rho <- corfit$consensus.correlation
if (!is.finite(rho)) stop("duplicateCorrelation returned a nonfinite value")
v <- voom(y, design, block = meta$individual, correlation = rho, plot = FALSE)
fit <- eBayes(lmFit(v, design, block = meta$individual, correlation = rho), robust = TRUE)

population <- topTable(fit, coef = population_coef, number = Inf,
                       sort.by = "none", adjust.method = "BH")
population$feature_id <- rownames(population)
population$effect_definition <- "Tibetan_minus_Han_average_across_five_timepoints"
population$significant_FDR_005 <- population$adj.P.Val <= opt$fdr
population <- merge(kept_meta, population, by = "feature_id", all.y = TRUE, sort = FALSE)
write.table(population, file.path(opt$outdir, "population_Tibetan_vs_Han.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

time_test <- topTable(fit, coef = time_coefs, number = Inf,
                      sort.by = "none", adjust.method = "BH")
time_test$feature_id <- rownames(time_test)
time_test$effect_definition <- "common_time_effect_adjusted_for_population"
time_test$significant_FDR_005 <- time_test$adj.P.Val <= opt$fdr
time_test <- merge(kept_meta, time_test, by = "feature_id", all.y = TRUE, sort = FALSE)
write.table(time_test, file.path(opt$outdir, "time_omnibus.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

normalized <- data.frame(kept_meta, v$E[rownames(y), sample_cols, drop = FALSE],
                         check.names = FALSE)
con <- gzfile(file.path(opt$outdir, "normalized_logCPM.tsv.gz"), "wt")
write.table(normalized, con, sep = "\t", row.names = FALSE, quote = FALSE)
close(con)

qc <- data.frame(
  key = c("assay", "input_features", "retained_features", "samples", "individuals",
          "normalization", "model", "population_estimand", "time_encoding",
          "individual_block", "duplicate_correlation", "min_cpm", "min_samples",
          "population_nominal_P005", "population_FDR005", "time_omnibus_FDR005"),
  value = c(opt$assay, nrow(counts), nrow(y), ncol(y), length(unique(meta$individual)),
            "edgeR_TMM", "limma_voom_duplicateCorrelation",
            "average_Tibetan_minus_Han_across_time", "categorical_5_levels", "yes",
            rho, opt$min_cpm, opt$min_samples,
            sum(population$P.Value < .05, na.rm = TRUE),
            sum(population$adj.P.Val <= opt$fdr, na.rm = TRUE),
            sum(time_test$adj.P.Val <= opt$fdr, na.rm = TRUE)),
  stringsAsFactors = FALSE
)
write.table(qc, file.path(opt$outdir, "model_qc.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
writeLines(capture.output(sessionInfo()), file.path(opt$outdir, "sessionInfo.txt"))

message("[INFO] ", opt$assay, " retained=", nrow(y), "; rho=", signif(rho, 4),
        "; population_FDR005=", sum(population$adj.P.Val <= opt$fdr, na.rm = TRUE))
