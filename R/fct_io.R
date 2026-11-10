# Reading and matching input files -------------------------------------------
# Conventions used everywhere downstream:
#   counts : numeric matrix, samples in ROWS, taxa in COLUMNS (vegan style)
#   meta   : data.frame, one row per sample, rownames = sample IDs,
#            rows in the same order as rownames(counts)

read_table_auto <- function(path) {
  data.table::fread(
    path,
    data.table = FALSE, check.names = FALSE,
    na.strings = c("", "NA"), showProgress = FALSE
  )
}

# Removes stray quotes, e.g. 'Ackermannviridae' -> Ackermannviridae
clean_ids <- function(x) gsub("^['\"]+|['\"]+$", "", trimws(as.character(x)))

# Reads an abundance table with taxa in rows and samples in columns
# (first column = taxon IDs) and returns a samples x taxa matrix.
read_abundance <- function(path) {
  df <- read_table_auto(path)
  if (ncol(df) < 3) {
    stop("The abundance matrix needs a taxon column and at least 2 sample columns.")
  }
  taxa <- clean_ids(df[[1]])
  df <- df[-1]

  non_numeric <- names(df)[!vapply(df, is.numeric, logical(1))]
  if (length(non_numeric)) {
    stop(sprintf(
      "%d sample column(s) are not numeric, e.g. %s.",
      length(non_numeric), paste(head(non_numeric, 3), collapse = ", ")
    ))
  }

  m <- as.matrix(df)
  rm(df)
  m[is.na(m)] <- 0
  if (any(m < 0)) stop("The abundance matrix contains negative values.")

  colnames(m) <- clean_ids(colnames(m))
  if (anyDuplicated(colnames(m))) stop("Duplicated sample IDs in the abundance matrix.")

  # Duplicated taxon names are summed
  if (anyDuplicated(taxa)) {
    m <- rowsum(m, taxa, reorder = FALSE)
  } else {
    rownames(m) <- taxa
  }
  t(m)
}

read_metadata <- function(path) {
  meta <- read_table_auto(path)
  names(meta)[names(meta) == ""] <- "V1"
  meta[] <- lapply(meta, function(col) {
    if (is.character(col)) {
      col <- trimws(col)
      col[col == ""] <- NA
    }
    col
  })
  meta
}

# Column of `meta` whose values best match the sample IDs
guess_id_column <- function(meta, sample_ids) {
  overlap <- vapply(meta, function(col) sum(clean_ids(col) %in% sample_ids), numeric(1))
  if (max(overlap) == 0) return(NA_character_)
  names(meta)[which.max(overlap)]
}

RANKS <- c("Phylum", "Class", "Order", "Family", "Genus", "Species", "vOTU", "Other")

guess_rank <- function(filename) {
  f <- tolower(basename(filename))
  hits <- c(
    vOTU = "votu", Species = "species", Genus = "genus", Family = "family",
    Order = "order", Class = "class", Phylum = "phylum"
  )
  for (rank in names(hits)) {
    if (grepl(paste0("(^|[^a-z])", hits[[rank]], "([^a-z]|$)"), f)) return(rank)
  }
  "Other"
}

rank_plural <- function(rank) {
  switch(rank,
    Phylum = "phyla", Class = "classes", Order = "orders", Family = "families",
    Genus = "genera", Species = "species", vOTU = "vOTUs", "taxa"
  )
}

# Keeps samples present in both files, aligns metadata to the matrix rows.
match_samples <- function(counts, meta, id_col) {
  ids <- clean_ids(meta[[id_col]])
  dup <- duplicated(ids) | is.na(ids)
  meta <- meta[!dup, , drop = FALSE]
  ids <- ids[!dup]

  keep <- intersect(rownames(counts), ids)
  if (length(keep) < 3) {
    stop(sprintf(
      "Only %d sample(s) are shared by the matrix and the metadata column '%s'.",
      length(keep), id_col
    ))
  }

  only_matrix <- setdiff(rownames(counts), ids)
  only_meta <- setdiff(ids, rownames(counts))

  meta <- meta[match(keep, ids), , drop = FALSE]
  rownames(meta) <- keep
  counts <- counts[keep, , drop = FALSE]

  tip <- sprintf(
    "%s sample(s) only in the matrix, %s only in the metadata%s",
    fmt_n(length(only_matrix)), fmt_n(length(only_meta)),
    if (any(dup)) sprintf(", %s duplicated/empty IDs ignored", fmt_n(sum(dup))) else ""
  )

  list(
    counts = counts,
    meta = meta,
    log = list(log_dims("Matched", counts, tip))
  )
}
