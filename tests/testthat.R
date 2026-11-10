# Run from the project root:  Rscript tests/testthat.R
for (f in list.files("R", full.names = TRUE)) source(f)
testthat::test_dir("tests/testthat", stop_on_failure = TRUE)
