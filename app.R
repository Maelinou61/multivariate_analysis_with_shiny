# Multivariate analysis of abundance matrices -----------------------------------------
# Run from the project root:  shiny::runApp()
# Files in R/ are sourced automatically by Shiny.

library(shiny)
library(bslib)

options(shiny.maxRequestSize = 3 * 1024^3) # uploads up to 3 GB

ui <- page_sidebar(
  title = "Multivariate analysis",
  theme = bs_theme(version = 5, primary = "#2a78d6"),
  tags$head(tags$link(rel = "stylesheet", href = "styles.css")),
  sidebar = sidebar(
    width = 340,
    accordion(
      id = "steps", open = c("Data", "Analysis"),
      mod_data_ui("data"),
      mod_filter_ui("filter"),
      mod_analysis_ui("analysis")
    ),
    actionButton("run", "Run analysis", icon = icon("play"), class = "btn-primary btn-lg")
  ),
  mod_summary_ui("summary"),
  navset_card_tab(
    id = "tabs", full_screen = TRUE,
    nav_panel("Ordination", mod_plot_ui("plot")),
    nav_panel("Scree", mod_scree_ui("scree")),
    nav_panel("Drivers", mod_drivers_ui("drivers")),
    nav_panel("PERMANOVA", mod_stats_ui("stats"))
  )
)

server <- function(input, output, session) {
  data_mod <- mod_data_server("data")
  data <- data_mod$data
  filters <- mod_filter_server("filter", data)
  analysis <- mod_analysis_server("analysis")

  # Result of the last run; cleared whenever the input data change
  result <- reactiveVal(NULL)
  observeEvent(data(), result(NULL))

  observeEvent(input$run, {
    d <- tryCatch(data(), error = function(e) NULL)
    if (is.null(d)) {
      showNotification(
        "No data loaded yet: choose your files in 'Data' and click 'Load' first.",
        type = "warning", duration = 8
      )
      return()
    }
    a <- analysis()
    label <- METHODS[[a$method_id]]$label
    showNotification(
      sprintf("Running %s on %s samples… (can take up to a minute on large data)", label, fmt_n(nrow(d$counts))),
      id = "running", type = "message", duration = NULL, closeButton = FALSE
    )
    on.exit(removeNotification("running"), add = TRUE)
    t0 <- Sys.time()

    out <- tryCatch(
      withProgress(message = METHODS[[a$method_id]]$label, value = 0.05, {
        setProgress(0.1, detail = "filtering")
        f <- apply_filters(d$counts, d$meta, filters())
        setProgress(0.3, detail = "transform, distance, ordination")
        res <- run_analysis(f$counts, a$method_id, a$params)
        setProgress(1)
        list(
          res = res,
          counts = f$counts,
          meta = f$meta,
          cols = plot_columns(f$meta, d$id_col),
          log = c(d$log, f$log, res$log)
        )
      }),
      error = function(e) {
        showNotification(conditionMessage(e), type = "error", duration = NULL)
        NULL
      }
    )
    if (!is.null(out)) {
      showNotification(
        sprintf("%s done in %.0f s.", label, as.numeric(Sys.time() - t0, units = "secs")),
        type = "message", duration = 4
      )
      nav_select("tabs", "Ordination")
    }
    result(out)
  })

  axes <- mod_plot_server("plot", result, data_mod$rank)
  mod_scree_server("scree", result)
  mod_drivers_server("drivers", result, axes)
  permanova <- mod_stats_server("stats", result)

  # Scree tab only for methods with variance per axis
  observe({
    r <- result()
    if (!is.null(r) && is.null(r$res$var_expl)) nav_hide("tabs", "Scree") else nav_show("tabs", "Scree")
  })

  mod_summary_server("summary", reactive({
    r <- result()
    if (is.null(r)) {
      d <- tryCatch(data(), error = function(e) NULL)
      return(if (is.null(d)) list() else d$log)
    }
    p <- permanova()
    c(r$log, if (!is.null(p)) list(permanova_log(p)))
  }))
}

shinyApp(ui, server)
