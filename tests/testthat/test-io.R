test_that("abundance file is read as samples x taxa, quotes stripped", {
  f <- tempfile(fileext = ".csv")
  writeLines(c("taxon,S1,S2,S3", "'A',1,0,2", "'B',0,5,1", "'A',1,1,1"), f)
  m <- read_abundance(f)
  expect_equal(dim(m), c(3, 2))
  expect_equal(colnames(m), c("A", "B"))
  expect_equal(unname(m["S1", "A"]), 2)      # duplicated taxa are summed
})

test_that("non-numeric matrix is rejected", {
  f <- tempfile(fileext = ".csv")
  writeLines(c("taxon,S1,S2", "A,1,x", "B,0,5"), f)
  expect_error(read_abundance(f), "not numeric")
})

test_that("ID column and rank are guessed", {
  meta <- data.frame(idx = 0:2, run = c("S1", "S2", "S9"), g = c("a", "b", "a"))
  expect_equal(guess_id_column(meta, c("S1", "S2", "S3")), "run")
  expect_equal(guess_rank("matrix_wgs_family.csv"), "Family")
  expect_equal(guess_rank("vOTU_table.tsv"), "vOTU")
  expect_equal(guess_rank("counts.csv"), "Other")
  expect_equal(rank_plural("Genus"), "genera")
})

test_that("samples are matched and aligned", {
  toy <- make_toy()
  meta <- toy$meta[rev(seq_len(nrow(toy$meta))), ][-1, ]
  m <- match_samples(toy$counts, meta, "run")
  expect_equal(rownames(m$counts), rownames(m$meta))
  expect_equal(nrow(m$counts), 29)
})

test_that("plot columns exclude IDs, index and '|' fields", {
  toy <- make_toy()
  toy$meta$V1 <- seq_len(nrow(toy$meta)) - 1
  cols <- plot_columns(toy$meta, "run")
  expect_setequal(cols$categorical, c("group", "block"))
  expect_equal(cols$numeric, "depth")
  expect_true("multi" %in% filter_columns(toy$meta, "run"))
})
