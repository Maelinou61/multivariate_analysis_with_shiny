# Sample and taxon filters ------------------------------------------------------

# Keeps samples whose `column` value is in `values` (exact match).
filter_by_metadata <- function(counts, meta, column = NULL, values = NULL) {
  if (is.null(column) || column == "" || length(values) == 0) {
    return(list(counts = counts, meta = meta, log = list()))
  }
  keep <- !is.na(meta[[column]]) & as.character(meta[[column]]) %in% values
  if (sum(keep) < 3) stop("Fewer than 3 samples left after the metadata filter.")
  list(
    counts = counts[keep, , drop = FALSE],
    meta = meta[keep, , drop = FALSE],
    log = list(log_dims(
      "Subset", counts[keep, , drop = FALSE],
      sprintf("%s ∈ {%s}", column, paste(values, collapse = ", "))
    ))
  )
}

# Abundance filters:
#   min_depth      : minimum total count per sample (samples with 0 are always removed)
#   min_prevalence : minimum fraction of samples in which a taxon is present (0-1)
#   min_count      : minimum total count of a taxon across samples
filter_by_abundance <- function(counts, meta, min_depth = 0, min_prevalence = 0, min_count = 0) {
  n0 <- dim(counts)

  depth <- rowSums(counts)
  keep_s <- depth >= max(min_depth, 1)
  counts <- counts[keep_s, , drop = FALSE]

  prev <- colSums(counts > 0) / nrow(counts)
  keep_t <- prev >= min_prevalence & colSums(counts) >= max(min_count, 1)
  counts <- counts[, keep_t, drop = FALSE]

  # Samples may become empty once taxa are removed
  keep_s2 <- rowSums(counts) > 0
  counts <- counts[keep_s2, , drop = FALSE]
  meta <- meta[rownames(counts), , drop = FALSE]

  if (nrow(counts) < 3 || ncol(counts) < 2) {
    stop("The filters are too strict: fewer than 3 samples or 2 taxa remain.")
  }

  tip <- sprintf(
    "−%s samples (depth < %s or empty), −%s taxa (prevalence < %s%% or total < %s)",
    fmt_n(n0[1] - nrow(counts)), fmt_n(max(min_depth, 1)),
    fmt_n(n0[2] - ncol(counts)), format(min_prevalence * 100), fmt_n(max(min_count, 1))
  )
  list(counts = counts, meta = meta, log = list(log_dims("Filtered", counts, tip)))
}

apply_filters <- function(counts, meta, settings) {
  a <- filter_by_metadata(counts, meta, settings$meta_column, settings$meta_values)
  b <- filter_by_abundance(
    a$counts, a$meta,
    min_depth = settings$min_depth,
    min_prevalence = settings$min_prevalence,
    min_count = settings$min_count
  )
  list(counts = b$counts, meta = b$meta, log = c(a$log, b$log))
}
