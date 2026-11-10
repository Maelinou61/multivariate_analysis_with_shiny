# Transformations (input: samples x taxa counts) ---------------------------------

TRANSFORM_LABELS <- c(
  none = "Raw counts",
  total = "Relative abundance",
  pa = "Presence/absence",
  hellinger = "Hellinger",
  clr = "CLR",
  rclr = "Robust CLR"
)

clr <- function(x, pseudocount = 1) {
  lx <- log(x + pseudocount)
  lx - rowMeans(lx)
}

# Robust CLR: log-ratio to the geometric mean of the non-zero values only.
# vegan (>= 2.6-6) then fills the zeros by low-rank matrix completion (OptSpace),
# as in the original notebook. If that fails (tiny or degenerate matrices),
# zeros are left at 0.
rclr <- function(x) {
  out <- tryCatch(
    unclass(vegan::decostand(x, "rclr")),
    error = function(e) NULL
  )
  if (!is.null(out)) {
    attr(out, "parameters") <- NULL
    dimnames(out) <- dimnames(x)
    return(out)
  }
  rclr_no_impute(x)
}

rclr_no_impute <- function(x) {
  lx <- log(x)
  lx[!is.finite(lx)] <- NA
  out <- lx - rowMeans(lx, na.rm = TRUE)
  out[is.na(out)] <- 0
  out
}

transform_counts <- function(counts, method, pseudocount = 1) {
  x <- switch(method,
    none = counts * 1.0,
    total = vegan::decostand(counts, "total"),
    pa = vegan::decostand(counts, "pa"),
    hellinger = vegan::decostand(counts, "hellinger"),
    clr = clr(counts, pseudocount),
    rclr = rclr(counts),
    stop("Unknown transformation: ", method)
  )
  attr(x, "decostand") <- NULL
  x
}
