# Plots ------------------------------------------------------------------------------

# Categorical colours, assigned in this fixed order. Validated for colour-blind
# separation; beyond 8 groups the rarest ones are folded into "Other".
PALETTE <- c(
  "#2a78d6", "#eb6834", "#1baf7a", "#eda100",
  "#e87ba4", "#008300", "#6250d6", "#e34948"
)
OTHER_COLOUR <- "#c9c8c3"
SEQ_LOW <- "#cfe0f6"
SEQ_HIGH <- "#123c6e"

axis_label <- function(res, i) {
  nm <- colnames(res$scores)[i]
  if (!is.null(res$var_expl) && i <= length(res$var_expl)) {
    sprintf("%s (%.1f%%)", nm, res$var_expl[i])
  } else {
    nm
  }
}

plot_title <- function(res, rank) {
  sprintf("%s · %s level", res$method$label, rank)
}

plot_subtitle <- function(res, rank, colour = NULL, highlight = NULL) {
  parts <- c(
    sprintf("%s samples · %s %s", fmt_n(nrow(res$scores)), fmt_n(ncol(res$x)), rank_plural(rank))
  )
  if (!is.null(res$stat$stress)) parts <- c(parts, sprintf("stress %.3f", res$stat$stress))
  if (length(highlight)) parts <- c(parts, sprintf("highlighted %s: %s", colour, paste(highlight, collapse = ", ")))
  paste(parts, collapse = " · ")
}

# Groups beyond the palette size are folded into "Other" (least frequent first)
fold_levels <- function(v, max_levels = length(PALETTE)) {
  v <- as.character(v)
  v[is.na(v)] <- "NA"
  tab <- sort(table(v), decreasing = TRUE)
  if (length(tab) <= max_levels) return(factor(v, levels = names(tab)))
  keep <- names(tab)[seq_len(max_levels - 1)]
  v[!v %in% keep] <- "Other"
  factor(v, levels = c(keep, "Other"))
}

# Data used by the ordination plot: scores joined to metadata
plot_data <- function(res, meta, axes = c(1, 2)) {
  df <- data.frame(
    sample = rownames(res$scores),
    x = res$scores[, axes[1]],
    y = res$scores[, axes[2]],
    stringsAsFactors = FALSE
  )
  cbind(df, meta[df$sample, , drop = FALSE])
}

make_ordination_plot <- function(res, meta, rank, axes = c(1, 2),
                                 colour = NULL, shape = NULL, highlight = NULL,
                                 ellipse = FALSE, size = 2, alpha = 0.8) {
  df <- plot_data(res, meta, axes)
  if (identical(colour, "")) colour <- NULL
  if (identical(shape, "")) shape <- NULL

  is_num <- !is.null(colour) && is.numeric(df[[colour]])
  if (!is.null(colour) && !is_num) {
    if (length(highlight)) {
      v <- as.character(df[[colour]])
      v[is.na(v) | !v %in% highlight] <- "Other"
      df$.colour <- factor(v, levels = c(intersect(highlight, v), "Other"))
    } else {
      df$.colour <- fold_levels(df[[colour]])
    }
  }
  if (!is.null(shape)) df$.shape <- fold_levels(df[[shape]], 6)

  df$.text <- paste0(
    "<b>", df$sample, "</b>",
    if (!is.null(colour)) paste0("<br>", colour, ": ", df[[colour]]) else "",
    if (!is.null(shape)) paste0("<br>", shape, ": ", df[[shape]]) else ""
  )
  if (is_num) df$.colour <- df[[colour]]
  aes_pts <- ggplot2::aes(x = x, y = y, text = .text)
  if (!is.null(colour)) aes_pts[["colour"]] <- ggplot2::aes(colour = .colour)[["colour"]]
  if (!is.null(shape)) aes_pts[["shape"]] <- ggplot2::aes(shape = .shape)[["shape"]]

  # "Other"/NA points go in a first layer so the groups of interest are drawn on
  # top (plotly draws one trace per layer and level, in layer order)
  background <- !is_num & !is.null(colour) & df$.colour %in% c("Other", "NA")
  if (is.null(colour) || is_num) background <- rep(FALSE, nrow(df))

  p <- ggplot2::ggplot(df) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
    ggplot2::geom_point(aes_pts, data = df[background, ], size = size, alpha = alpha) |> suppressWarnings() +
    ggplot2::geom_point(aes_pts, data = df[!background, ], size = size, alpha = alpha) |> suppressWarnings()

  if (is_num) {
    p <- p + ggplot2::scale_colour_gradient(low = SEQ_LOW, high = SEQ_HIGH, name = colour)
  } else if (!is.null(colour)) {
    lv <- levels(df$.colour)
    cols <- stats::setNames(PALETTE[seq_along(lv)], lv)
    cols[names(cols) %in% c("Other", "NA")] <- OTHER_COLOUR
    p <- p + ggplot2::scale_colour_manual(values = cols, name = colour)
    if (ellipse) {
      ell <- df[!df$.colour %in% c("Other", "NA"), ]
      ok <- names(which(table(ell$.colour) >= 4))
      p <- p + ggplot2::stat_ellipse(
        data = ell[ell$.colour %in% ok, ],
        ggplot2::aes(x = x, y = y, colour = .colour),
        linewidth = 0.6, show.legend = FALSE
      )
    }
  }
  if (!is.null(shape)) p <- p + ggplot2::scale_shape_discrete(name = shape)

  p +
    ggplot2::labs(
      x = axis_label(res, axes[1]), y = axis_label(res, axes[2]),
      title = plot_title(res, rank),
      subtitle = plot_subtitle(res, rank, colour, highlight)
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey94"),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey35")
    )
}

make_scree_plot <- function(res, n = 10) {
  v <- res$var_expl[seq_len(min(n, length(res$var_expl)))]
  df <- data.frame(axis = factor(colnames(res$scores)[seq_along(v)], levels = colnames(res$scores)[seq_along(v)]), pct = v)
  ggplot2::ggplot(df, ggplot2::aes(axis, pct)) +
    ggplot2::geom_col(fill = PALETTE[1], width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f", pct)), vjust = -0.4, size = 3.5, colour = "grey25") +
    ggplot2::labs(x = NULL, y = "Variance explained (%)", title = res$method$label) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.1))) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank())
}
