# Small simulated data set: 3 groups with different composition
make_toy <- function(n_per = 10, n_taxa = 25, seed = 42) {
  set.seed(seed)
  groups <- rep(c("A", "B", "C"), each = n_per)
  base <- matrix(rexp(3 * n_taxa), 3)
  counts <- t(sapply(seq_along(groups), function(i) {
    rpois(n_taxa, 20 * base[match(groups[i], c("A", "B", "C")), ])
  }))
  rownames(counts) <- sprintf("S%02d", seq_along(groups))
  colnames(counts) <- sprintf("T%02d", seq_len(n_taxa))
  meta <- data.frame(
    run = rownames(counts),
    group = groups,
    depth = rowSums(counts),
    block = rep(c("x", "y"), length.out = length(groups)),
    multi = rep(c("a|b", "a"), length.out = length(groups)),
    row.names = rownames(counts)
  )
  list(counts = counts, meta = meta)
}
