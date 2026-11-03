#!/usr/bin/env Rscript
# Run from the package root after installing the R dependencies.
set.seed(17)
samples <- unlist(lapply(c("A", "W"), function(pop) unlist(lapply(1:5, function(i)
  paste(pop, i, 1:5, sep = "-")))))
individual <- sub("-[1-5]$", "", samples)
time <- as.integer(sub(".*-", "", samples))
donors <- unique(individual)
baseline <- exp(rnorm(300, 4, 0.7))
donor_effect <- matrix(rnorm(300 * 10, 0, 0.25), nrow = 300)
mu <- sapply(seq_along(samples), function(j)
  baseline * exp(donor_effect[, match(individual[j], donors)] + 0.1 * (time[j]-1)))
counts <- matrix(rnbinom(length(mu), mu = mu, size = 20), nrow = 300)
tab <- data.frame(peak_id = paste0("p", 1:300), counts, check.names = FALSE)
names(tab)[-1] <- samples
tmp <- tempfile("tesv_r_test_")
dir.create(tmp)
write.table(tab, file.path(tmp, "counts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
status <- system2(file.path(R.home("bin"), "Rscript"), c("fit_additive_population_model.R",
  "--matrix", shQuote(file.path(tmp, "counts.tsv")), "--assay", "atac",
  "--outdir", shQuote(file.path(tmp, "out"))))
stopifnot(status == 0)
for (file in c("population_Tibetan_vs_Han.tsv", "time_omnibus.tsv")) {
  x <- read.delim(file.path(tmp, "out", file))
  stopifnot(nrow(x) > 0, all(is.finite(x$P.Value)),
            all(x$P.Value >= 0 & x$P.Value <= 1),
            isTRUE(all.equal(x$adj.P.Val, p.adjust(x$P.Value, "BH"))))
}
message("Population-model smoke test passed. Temporary outputs: ", tmp)
