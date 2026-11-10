# Taxa driving the ordination axes ---------------------------------------------------
#   PCA  : squared loadings (unit-norm eigenvectors) x 100
#   CA   : column weight x squared loading x 100
#   PCoA : covariance between transformed abundances and axis scores,
#          contribution = share of the total absolute covariance
#   NMDS / UMAP : Spearman correlation with the axes (scale-free axes),
#          contribution = share of the total absolute correlation

compute_drivers <- function(res, counts, axes = c(1, 2)) {
  engine <- res$method$engine
  x <- res$x
  ax <- res$scores[, axes, drop = FALSE]

  if (engine == "pca") {
    score <- res$loadings[, axes, drop = FALSE]
    contrib <- score^2 * 100
  } else if (engine == "ca") {
    score <- res$loadings[, axes, drop = FALSE]
    contrib <- res$weights * score^2 * 100
  } else {
    score <- if (engine == "pcoa") stats::cov(x, ax) else suppressWarnings(stats::cor(x, ax, method = "spearman"))
    score[is.na(score)] <- 0
    contrib <- sweep(abs(score), 2, colSums(abs(score)), "/") * 100
  }

  rel <- counts / rowSums(counts)
  lab <- colnames(res$scores)[axes]
  out <- data.frame(
    taxon = colnames(x),
    score_x = score[, 1],
    score_y = score[, 2],
    contrib_x = contrib[, 1],
    contrib_y = contrib[, 2],
    mean_rel_pct = colMeans(rel)[colnames(x)] * 100,
    prevalence_pct = colMeans(counts > 0)[colnames(x)] * 100,
    row.names = NULL,
    check.names = FALSE
  )
  out <- out[order(out$contrib_x + out$contrib_y, decreasing = TRUE), ]
  names(out)[2:5] <- c(
    paste(lab[1], "score"), paste(lab[2], "score"),
    paste(lab[1], "contrib. %"), paste(lab[2], "contrib. %")
  )
  out
}
