test_that("filters remove samples and taxa", {
  toy <- make_toy()
  toy$counts[1, ] <- 0
  toy$counts[, 1] <- c(1, rep(0, 29))
  f <- apply_filters(toy$counts, toy$meta, list(
    meta_column = "group", meta_values = c("A", "B"),
    min_depth = 0, min_prevalence = 0.1, min_count = 0
  ))
  expect_equal(nrow(f$counts), 19)
  expect_false("T01" %in% colnames(f$counts))
  expect_equal(rownames(f$counts), rownames(f$meta))
  expect_length(f$log, 2)
})

test_that("robust CLR matches vegan", {
  toy <- make_toy()
  ref <- vegan::decostand(toy$counts, "rclr")
  expect_equal(as.vector(rclr(toy$counts)), as.vector(ref), tolerance = 1e-10)
  ref0 <- vegan::decostand(toy$counts, "rclr", impute = FALSE)
  ref0[is.na(ref0)] <- 0
  expect_equal(as.vector(rclr_no_impute(toy$counts)), as.vector(ref0), tolerance = 1e-10)
  degenerate <- matrix(c(0, 2, 4, 8, 1, 0, 3, 9, 5, 5, 0, 0), 3, byrow = TRUE)
  expect_equal(dim(rclr(degenerate)), c(3, 4))
})

test_that("PCoA matches cmdscale", {
  toy <- make_toy()
  d <- vegan::vegdist(vegan::decostand(toy$counts, "total"), "bray")
  ours <- ord_pcoa(d, k = 3)
  ref <- cmdscale(d, k = 3, eig = TRUE)
  expect_equal(abs(unname(ours$scores)), abs(unname(ref$points)), tolerance = 1e-8)
  expect_equal(ours$var_expl, (ref$eig / sum(ref$eig[ref$eig > 0]) * 100)[1:3], tolerance = 1e-8)
})

test_that("PCA variance matches prcomp", {
  toy <- make_toy()
  x <- clr(toy$counts)
  ref <- prcomp(x)
  ours <- ord_pca(x, k = 3)
  expect_equal(ours$var_expl, (ref$sdev^2 / sum(ref$sdev^2) * 100)[1:3], tolerance = 1e-8)
})

test_that("every registered analysis runs and gives drivers + plot", {
  toy <- make_toy()
  for (id in names(METHODS)) {
    res <- run_analysis(toy$counts, id, list(n_neighbors = 5))
    expect_equal(rownames(res$scores), rownames(toy$counts), info = id)
    expect_true(ncol(res$scores) >= 2, info = id)
    dr <- compute_drivers(res, toy$counts)
    expect_equal(nrow(dr), ncol(toy$counts), info = id)
    p <- make_ordination_plot(res, toy$meta, "Family", colour = "group", shape = "block", ellipse = TRUE)
    expect_s3_class(p, "ggplot")
    expect_match(p$labels$title, "Family level")
  }
})

test_that("highlight mode keeps chosen values and folds the rest", {
  toy <- make_toy()
  res <- run_analysis(toy$counts, "pcoa_bray")
  p <- make_ordination_plot(res, toy$meta, "Genus", colour = "group", highlight = "B")
  expect_setequal(levels(p$data$.colour), c("B", "Other"))
  p2 <- make_ordination_plot(res, toy$meta, "Genus", colour = "depth")
  expect_s3_class(p2, "ggplot")
})

test_that("PCoA on vegan's Aitchison distance equals PCA on CLR", {
  toy <- make_toy()
  pcoa <- run_analysis(toy$counts, "pcoa_aitchison", list(pseudocount = 1))
  pca <- run_analysis(toy$counts, "pca_clr", list(pseudocount = 1))
  expect_equal(abs(unname(pcoa$scores[, 1:3])), abs(unname(pca$scores[, 1:3])), tolerance = 1e-6)
  expect_equal(pcoa$var_expl[1:3], pca$var_expl[1:3], tolerance = 1e-6)
  expect_equal(c(pcoa$d), c(vegan::vegdist(toy$counts, "aitchison", pseudocount = 1)))
})

test_that("robust Aitchison PCoA uses vegan's distance, drivers on rCLR", {
  toy <- make_toy()
  res <- run_analysis(toy$counts, "pcoa_robust_aitchison")
  expect_equal(c(res$d), c(vegan::vegdist(toy$counts, "robust.aitchison")))
  expect_equal(res$log[[1]]$value, "Robust CLR")
  run_permanova(analysis_distance(res), toy$meta, "group", 19) |> expect_type("list")
})
