# Which metadata columns can be used where --------------------------------------

# A column is a row index if it is integer-like and equal to 0..n-1 or 1..n,
# or if it has no real name (e.g. pandas' unnamed index column).
is_index_column <- function(name, col) {
  if (name %in% c("V1", "...1", "X", "Unnamed: 0")) return(TRUE)
  if (!is.numeric(col) || anyNA(col)) return(FALSE)
  s <- sort(col)
  n <- length(col)
  isTRUE(all.equal(s, 0:(n - 1))) || isTRUE(all.equal(s, seq_len(n)))
}

has_multi_values <- function(col, sep = "|") {
  is.character(col) && any(grepl(sep, col, fixed = TRUE), na.rm = TRUE)
}

# Columns usable for colour / shape / PERMANOVA.
# Returns list(categorical = <names>, numeric = <names>).
plot_columns <- function(meta, id_col = NULL, max_levels = 30) {
  categorical <- character(0)
  numeric <- character(0)
  n <- nrow(meta)
  for (nm in setdiff(names(meta), id_col)) {
    col <- meta[[nm]]
    if (is_index_column(nm, col) || has_multi_values(col)) next
    n_lev <- length(unique(stats::na.omit(col)))
    if (n_lev < 2) next
    if (is.numeric(col)) {
      numeric <- c(numeric, nm)
    } else if (n_lev <= max_levels && n_lev < n) {
      categorical <- c(categorical, nm)
    }
  }
  list(categorical = categorical, numeric = numeric)
}

# Columns usable to subset samples (exact match on values, '|' fields allowed).
filter_columns <- function(meta, id_col = NULL, max_levels = 100) {
  keep <- vapply(setdiff(names(meta), id_col), function(nm) {
    col <- meta[[nm]]
    n_lev <- length(unique(stats::na.omit(col)))
    !is_index_column(nm, col) && n_lev >= 2 && n_lev <= max_levels
  }, logical(1))
  names(keep)[keep]
}
