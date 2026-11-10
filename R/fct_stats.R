# PERMANOVA + dispersion test ------------------------------------------------------
# Always run on the analysis distance matrix, never on the ordination scores.
#
# vegan::adonis2 permutes the full n x n matrix at every permutation, which takes
# tens of minutes for several thousand samples. For a single term the same
# pseudo-F can be computed by permuting the group labels instead:
#   SS_total  = sum(d_ij^2) / n
#   SS_within = sum_g sum_{i,j in g} d_ij^2 / n_g        (categorical)
#   SS_model  = tr(H G H) with G the Gower-centred matrix (numeric, H = hat matrix)
# Results are identical to adonis2 (see tests/testthat/test-stats.R).

permanova_core <- function(D2, group, perms, progress = NULL) {
  n <- nrow(D2)
  ss_total <- sum(D2) / (2 * n)

  if (is.factor(group)) {
    df_model <- nlevels(group) - 1
    ss_within <- function(g) {
      Z <- stats::model.matrix(~ g - 1)
      ng <- colSums(Z)
      sum(colSums(Z * (D2 %*% Z)) / ng) / 2
    }
    ss_model_fun <- function(g) ss_total - ss_within(g)
  } else {
    df_model <- 1
    A <- -0.5 * D2
    G <- A - rowMeans(A) - rep(colMeans(A), each = n) + mean(A)
    ss_model_fun <- function(g) {
      v <- g - mean(g)
      as.numeric(crossprod(v, G %*% v)) / sum(v^2)
    }
  }
  df_res <- n - df_model - 1
  f_of <- function(ss) (ss / df_model) / ((ss_total - ss) / df_res)

  ss_obs <- ss_model_fun(group)
  f_obs <- f_of(ss_obs)
  n_perm <- nrow(perms)
  f_perm <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    f_perm[i] <- f_of(ss_model_fun(group[perms[i, ]]))
    if (!is.null(progress) && i %% 10 == 0) progress(i / n_perm)
  }
  p <- (sum(f_perm >= f_obs - sqrt(.Machine$double.eps)) + 1) / (nrow(perms) + 1)

  data.frame(
    Df = c(df_model, df_res, n - 1),
    SumOfSqs = c(ss_obs, ss_total - ss_obs, ss_total),
    R2 = c(ss_obs, ss_total - ss_obs, ss_total) / ss_total,
    F = c(f_obs, NA, NA),
    `Pr(>F)` = c(p, NA, NA),
    row.names = c("Model", "Residual", "Total"),
    check.names = FALSE
  )
}

# Distance of each sample to its group centroid (pseudo-Euclidean space,
# same as vegan::betadisper(type = "centroid")), then a one-way ANOVA.
dispersion_test <- function(D2, group) {
  Z <- stats::model.matrix(~ group - 1)
  ng <- colSums(Z)
  to_group <- (D2 %*% Z) / rep(ng, each = nrow(D2))       # mean d^2 to each group
  within <- colSums(Z * (D2 %*% Z)) / (2 * ng^2)           # group "spread" term
  gi <- as.integer(group)
  z <- to_group[cbind(seq_along(gi), gi)] - within[gi]
  disp <- sqrt(abs(z))
  stats::anova(stats::lm(disp ~ group))[["Pr(>F)"]][1]
}

run_permanova <- function(d, meta, variable, permutations = 999, strata = NULL,
                          seed = 1, progress = NULL) {
  has_strata <- !is.null(strata) && strata != ""
  group <- meta[[variable]]
  keep <- !is.na(group)
  if (has_strata) keep <- keep & !is.na(meta[[strata]])

  group <- group[keep]
  if (!is.numeric(group)) group <- droplevels(factor(group))
  if (is.factor(group) && nlevels(group) < 2) stop("The grouping variable needs at least 2 groups.")
  if (is.numeric(group) && stats::var(group) == 0) stop("The variable is constant.")

  D2 <- as.matrix(subset_dist(d, keep))^2
  ctrl <- permute::how(nperm = permutations)
  if (has_strata) permute::setBlocks(ctrl) <- factor(meta[[strata]][keep])
  set.seed(seed)
  perms <- permute::shuffleSet(length(group), control = ctrl)

  tab <- permanova_core(D2, group, perms, progress)
  rownames(tab)[1] <- variable
  disp_p <- if (is.factor(group)) dispersion_test(D2, group) else NA_real_

  list(
    table = tab,
    variable = variable,
    r2 = tab$R2[1],
    f = tab$F[1],
    p = tab[["Pr(>F)"]][1],
    disp_p = disp_p,
    n = sum(keep),
    n_groups = if (is.factor(group)) nlevels(group) else NA_integer_,
    permutations = nrow(perms),
    strata = strata
  )
}

fmt_p <- function(p) {
  if (is.na(p)) return("–")
  if (p < 0.001) "p < 0.001" else sprintf("p = %.3f", p)
}

permanova_log <- function(res) {
  log_entry(
    "PERMANOVA",
    sprintf("R² %.3f", res$r2),
    fmt_p(res$p),
    tooltip = sprintf(
      "%s, n = %s, %d permutations%s. Dispersion test: %s",
      res$variable, fmt_n(res$n), res$permutations,
      if (!is.null(res$strata) && res$strata != "") paste0(", strata = ", res$strata) else "",
      fmt_p(res$disp_p)
    )
  )
}
