# Distances (input: transformed samples x taxa matrix) ----------------------------

compute_distance <- function(x, method) {
  switch(method,
    bray = vegan::vegdist(x, "bray"),
    jaccard = vegan::vegdist(x, "jaccard", binary = TRUE),
    euclidean = stats::dist(x),
    chisq = vegan::vegdist(x, "chisq"),
    stop("Unknown distance: ", method)
  )
}

# Keeps the samples flagged TRUE in `keep`
subset_dist <- function(d, keep) {
  if (all(keep)) return(d)
  stats::as.dist(as.matrix(d)[keep, keep, drop = FALSE])
}
