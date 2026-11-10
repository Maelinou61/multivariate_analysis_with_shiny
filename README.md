# Multivariate analysis – Shiny app

Interactive ordinations (PCoA, PCA, CA, NMDS, UMAP) of an abundance matrix,
coloured by metadata, with optional PERMANOVA.

## Run

```bash
conda env create -f environment.yml      # once
conda activate shiny_mva
Rscript -e 'shiny::runApp(".", port = 8765, launch.browser = FALSE)'
```

Then open http://127.0.0.1:8765 in your browser. (The conda R has no default
browser set, so `launch.browser = TRUE` fails unless you add
`options(browser = "xdg-open");` before `shiny::runApp(...)`.)

Tests: `Rscript tests/testthat.R`

## Input files

- **Abundance matrix**: CSV/TSV, first column = taxon IDs, one column per sample, raw counts.
- **Metadata**: CSV/TSV, one row per sample. The column matching the sample names is
  detected automatically (can be changed in the app).

Columns containing `|` (multi-valued) can be used to filter samples but not to colour.

## Layout

```
app.R                 UI layout + wiring of the modules
R/fct_*.R             plain R functions (no Shiny) – usable from a notebook:
                      for (f in list.files("R", full.names = TRUE)) source(f)
R/methods_registry.R  list of analyses; add an entry to add a method
R/mod_*.R             Shiny modules (data, filter, analysis, plot, drivers, stats, summary)
www/styles.css        processing-summary strip
tests/testthat/       tests (PERMANOVA, PCoA, rCLR checked against vegan / base R)
```

## Notes on methods

- **PERMANOVA** uses the distance of the chosen analysis (never the plotted
  coordinates). It is a single-term implementation that gives the same results as
  `vegan::adonis2` but stays fast for thousands of samples. A dispersion test
  (`betadisper`, centroid) is shown next to it.
- **Robust CLR** uses `vegan::decostand(method = "rclr")` (zeros imputed by OptSpace).
- **PCoA / PCA** above 2 000 samples use truncated decompositions (RSpectra); the
  % variance of PCoA is then relative to the trace of the centred matrix.
- **NMDS** is limited to 1 500 samples (too slow beyond).
