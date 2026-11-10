test_that("PERMANOVA matches vegan::adonis2 (categorical)", {
  toy <- make_toy()
  d <- vegan::vegdist(vegan::decostand(toy$counts, "total"), "bray")
  ours <- run_permanova(d, toy$meta, "group", permutations = 199, seed = 3)

  set.seed(3)
  perms <- permute::shuffleSet(30, control = permute::how(nperm = 199))
  ref <- vegan::adonis2(d ~ group, data = toy$meta, permutations = perms)
  expect_equal(ours$r2, ref$R2[1], tolerance = 1e-10)
  expect_equal(ours$f, ref$F[1], tolerance = 1e-10)
  expect_equal(ours$p, ref[["Pr(>F)"]][1])
})

test_that("PERMANOVA matches vegan::adonis2 (numeric, strata)", {
  toy <- make_toy()
  d <- vegan::vegdist(vegan::decostand(toy$counts, "total"), "bray")
  ours <- run_permanova(d, toy$meta, "depth", permutations = 99, strata = "block", seed = 5)

  ctrl <- permute::how(nperm = 99, blocks = factor(toy$meta$block))
  set.seed(5)
  perms <- permute::shuffleSet(30, control = ctrl)
  ref <- vegan::adonis2(d ~ depth, data = toy$meta, permutations = perms)
  expect_equal(ours$r2, ref$R2[1], tolerance = 1e-10)
  expect_equal(ours$p, ref[["Pr(>F)"]][1])
})

test_that("dispersion test matches betadisper(type = 'centroid')", {
  toy <- make_toy()
  d <- vegan::vegdist(vegan::decostand(toy$counts, "total"), "bray")
  ours <- dispersion_test(as.matrix(d)^2, factor(toy$meta$group))
  bd <- vegan::betadisper(d, factor(toy$meta$group), type = "centroid", bias.adjust = FALSE)
  expect_equal(ours, anova(bd)[["Pr(>F)"]][1], tolerance = 1e-8)
})

test_that("samples with NA in the variable are dropped", {
  toy <- make_toy()
  toy$meta$group[1:3] <- NA
  d <- vegan::vegdist(toy$counts, "bray")
  expect_equal(run_permanova(d, toy$meta, "group", 9)$n, 27)
})
