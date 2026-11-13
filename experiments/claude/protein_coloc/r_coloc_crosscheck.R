# Runs R coloc::coloc.abf on per-pair windows exported by r_coloc_crosscheck.py.
# Usage: Rscript r_coloc_crosscheck.R <manifest.csv> <out.csv>
# manifest columns: pair, path, p1, p2, p12, s
# Each window CSV has columns snp, BETA_TRAIT, SE_TRAIT, BETA_PROTEIN, SE_PROTEIN.
#
# Checked against coloc 5.2.3 source (process.dataset / approx.bf.estimates): when beta and
# varbeta are given, type "cc" uses sd.prior = 0.2 and ignores s and MAF; type "quant" with sdY
# given uses sd.prior = 0.15 * sdY. These match ColocAbf's CaseControlTraitModel (W = 0.2) and
# QuantitativeTraitModel(sd_y = 1) (W = 0.15).
suppressPackageStartupMessages(library(coloc))
args <- commandArgs(trailingOnly = TRUE)
manifest <- read.csv(args[1], stringsAsFactors = FALSE)
rows <- vector("list", nrow(manifest))
for (i in seq_len(nrow(manifest))) {
  m <- manifest[i, ]
  d <- read.csv(m$path, stringsAsFactors = FALSE)
  stopifnot(!anyDuplicated(d$snp))
  invisible(capture.output(res <- suppressWarnings(suppressMessages(coloc.abf(
    dataset1 = list(beta = d$BETA_TRAIT, varbeta = d$SE_TRAIT^2, type = "cc", s = m$s, snp = d$snp),
    dataset2 = list(beta = d$BETA_PROTEIN, varbeta = d$SE_PROTEIN^2, type = "quant", sdY = 1, snp = d$snp),
    p1 = m$p1, p2 = m$p2, p12 = m$p12
  )))))
  s <- res$summary
  rows[[i]] <- data.frame(
    pair = m$pair, nsnps = s[["nsnps"]],
    R_PP_H0 = s[["PP.H0.abf"]], R_PP_H1 = s[["PP.H1.abf"]], R_PP_H2 = s[["PP.H2.abf"]],
    R_PP_H3 = s[["PP.H3.abf"]], R_PP_H4 = s[["PP.H4.abf"]]
  )
}
write.csv(do.call(rbind, rows), args[2], row.names = FALSE)
