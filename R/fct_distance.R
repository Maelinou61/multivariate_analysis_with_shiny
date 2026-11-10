# Distances (input: transformed samples x taxa matrix) ----------------------------

compute_distance <- function(x, method, pseudocount = 1) {
  switch(method,
    bray = vegan::vegdist(x, "bray"),
    jaccard = vegan::vegdist(x, "jaccard", binary = TRUE),
    euclidean = stats::dist(x),
    chisq = vegan::vegdist(x, "chisq"),
    # Aitchison = Euclidean distance between CLR-transformed samples (on raw counts)
    aitchison = vegan::vegdist(x, "aitchison", pseudocount = pseudocount),
    # Robust Aitchison = Euclidean distance on vegan's robust CLR (zeros imputed);
    # falls back to zeros kept at 0 if the imputation fails on degenerate data
    robust.aitchison = tryCatch(
      vegan::vegdist(x, "robust.aitchison"),
      error = function(e) stats::dist(rclr_no_impute(x))
    ),
    stop("Unknown distance: ", method)
  )
}

# Keeps the samples flagged TRUE in `keep`
subset_dist <- function(d, keep) {
  if (all(keep)) return(d)
  stats::as.dist(as.matrix(d)[keep, keep, drop = FALSE])
}
