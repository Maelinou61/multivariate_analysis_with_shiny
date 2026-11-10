# Ordination engines --------------------------------------------------------------
# Every engine returns the same structure:
#   scores      : matrix, samples x axes
#   var_expl    : % variance per axis (NULL for NMDS / UMAP)
#   loadings    : taxa x axes (PCA / CA only, NULL otherwise)
#   weights     : taxon weights (CA only)
#   stat        : named list of extra values (e.g. NMDS stress)

# Large problems use a truncated eigen/SVD decomposition (RSpectra)
FULL_DECOMP_MAX <- 2000

ord_pcoa <- function(d, k = 10) {
  n <- attr(d, "Size")
  k <- min(k, n - 2)
  # Gower double-centring of -D^2/2
  A <- -0.5 * as.matrix(d)^2
  G <- A - rowMeans(A) - rep(colMeans(A), each = n) + mean(A)
  rm(A)
  if (n <= FULL_DECOMP_MAX) {
    e <- eigen(G, symmetric = TRUE)
    total <- sum(e$values[e$values > 0])
    vals <- e$values[seq_len(k)]
    vecs <- e$vectors[, seq_len(k), drop = FALSE]
  } else {
    e <- RSpectra::eigs_sym(G, k, which = "LA")
    total <- sum(diag(G))
    vals <- e$values
    vecs <- e$vectors
  }
  pos <- vals > 0
  vals <- vals[pos]
  scores <- vecs[, pos, drop = FALSE] %*% diag(sqrt(vals), nrow = length(vals))
  rownames(scores) <- labels(d)
  list(scores = scores, var_expl = vals / total * 100)
}

ord_pca <- function(x, k = 10) {
  xc <- scale(x, center = TRUE, scale = FALSE)
  total <- sum(xc^2)
  k <- min(k, nrow(x) - 1, ncol(x))
  if (min(dim(xc)) <= FULL_DECOMP_MAX) {
    s <- svd(xc, nu = k, nv = k)
    s$d <- s$d[seq_len(k)]
  } else {
    s <- RSpectra::svds(xc, k)
  }
  scores <- s$u %*% diag(s$d, nrow = k)
  rownames(scores) <- rownames(x)
  rownames(s$v) <- colnames(x)
  list(scores = scores, var_expl = s$d^2 / total * 100, loadings = s$v)
}

ord_ca <- function(x, k = 10) {
  fit <- vegan::cca(x)
  k <- min(k, length(fit$CA$eig))
  eig <- fit$CA$eig
  sc <- vegan::scores(fit, display = "sites", choices = seq_len(k), scaling = 1)
  list(
    scores = as.matrix(sc),
    var_expl = (eig / sum(eig) * 100)[seq_len(k)],
    loadings = fit$CA$v[, seq_len(k), drop = FALSE],
    weights = fit$colsum
  )
}

# NMDS iterates on all pairwise ranks: hours beyond a few thousand samples
NMDS_MAX <- 1500

ord_nmds <- function(d, k = 2, seed = 1) {
  n <- attr(d, "Size")
  if (n > NMDS_MAX) {
    stop(sprintf(
      "NMDS is limited to %s samples (you have %s). Subset the samples or use PCoA/UMAP.",
      fmt_n(NMDS_MAX), fmt_n(n)
    ))
  }
  set.seed(seed)
  fit <- vegan::metaMDS(d, k = k, trymax = 20, autotransform = FALSE, trace = FALSE)
  list(scores = as.matrix(fit$points), stat = list(stress = fit$stress))
}

ord_umap <- function(input, n_neighbors = 15, min_dist = 0.1, seed = 1) {
  n <- if (inherits(input, "dist")) attr(input, "Size") else nrow(input)
  n_neighbors <- min(n_neighbors, n - 1)
  set.seed(seed)
  sc <- uwot::umap(input, n_neighbors = n_neighbors, min_dist = min_dist, n_components = 2)
  rownames(sc) <- if (inherits(input, "dist")) labels(input) else rownames(input)
  list(scores = sc, stat = list(n_neighbors = n_neighbors, min_dist = min_dist))
}

# Runs one registry entry from filtered counts.
# Returns the standard structure plus: method, x (transformed), d (or NULL), log.
run_analysis <- function(counts, method_id, params = list()) {
  m <- METHODS[[method_id]]
  if (is.null(m)) stop("Unknown analysis: ", method_id)
  seed <- params$seed %||% 1

  x <- transform_counts(counts, m$transform, params$pseudocount %||% 1)

  need_dist <- m$engine %in% ENGINES_NEED_DIST ||
    (m$engine == "umap" && m$distance != "euclidean")
  d <- if (need_dist) compute_distance(x, m$distance) else NULL

  res <- switch(m$engine,
    pcoa = ord_pcoa(d),
    pca = ord_pca(x),
    ca = ord_ca(x),
    nmds = ord_nmds(d, seed = seed),
    umap = ord_umap(
      if (is.null(d)) x else d,
      n_neighbors = params$n_neighbors %||% 15,
      min_dist = params$min_dist %||% 0.1,
      seed = seed
    )
  )
  colnames(res$scores) <- paste0(m$prefix, seq_len(ncol(res$scores)))

  res$method_id <- method_id
  res$method <- m
  res$x <- x
  res$d <- d
  res$log <- list(
    log_entry("Transform", TRANSFORM_LABELS[[m$transform]], "transformation"),
    log_entry("Distance", m$dist_label, "distance"),
    log_entry(
      "Ordination", m$prefix,
      caption = if (!is.null(res$stat$stress)) sprintf("stress %.3f", res$stat$stress)
                else if (!is.null(res$var_expl)) sprintf("%.1f%% on 2 axes", sum(res$var_expl[1:2]))
                else "embedding",
      tooltip = m$label
    )
  )
  res
}

# Distance used for PERMANOVA (computed lazily for engines that did not need it)
analysis_distance <- function(res) {
  if (!is.null(res$d)) return(res$d)
  compute_distance(res$x, res$method$distance)
}
