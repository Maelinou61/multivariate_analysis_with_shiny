# Processing log -----------------------------------------------------------
# Every pipeline step returns one or more log entries. The summary strip
# renders them in order: a big value, a small caption, and a hover tooltip.

log_entry <- function(step, value, caption = NULL, tooltip = NULL) {
  list(step = step, value = value, caption = caption, tooltip = tooltip)
}

fmt_n <- function(x) format(x, big.mark = " ", scientific = FALSE, trim = TRUE)

# Entry describing the size of a samples x taxa matrix
log_dims <- function(step, counts, tooltip = NULL) {
  log_entry(
    step,
    value = sprintf("%s × %s", fmt_n(nrow(counts)), fmt_n(ncol(counts))),
    caption = "samples × taxa",
    tooltip = tooltip
  )
}
