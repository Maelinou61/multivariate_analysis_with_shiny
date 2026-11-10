# PERMANOVA module (run on demand) -------------------------------------------------------
# Returns reactive(permanova result or NULL)

mod_stats_ui <- function(id) {
  ns <- NS(id)
  layout_sidebar(
    sidebar = sidebar(
      width = 280, open = "always",
      selectInput(ns("variable"), "Variable", choices = NULL),
      selectInput(ns("strata"), "Permute within (optional)", choices = c("None" = "")),
      selectInput(ns("perms"), "Permutations", choices = c(99, 199, 499, 999), selected = 499),
      actionButton(ns("run"), "Run PERMANOVA", icon = icon("play"), class = "btn-primary w-100"),
      uiOutput(ns("hint"))
    ),
    uiOutput(ns("result"))
  )
}

mod_stats_server <- function(id, result) {
  moduleServer(id, function(input, output, session) {
    stats_res <- reactiveVal(NULL)

    observeEvent(result(), {
      stats_res(NULL)
      cols <- result()$cols
      updateSelectInput(session, "variable",
        choices = list(Categorical = cols$categorical, Numeric = cols$numeric)
      )
      updateSelectInput(session, "strata", choices = c("None" = "", cols$categorical))
    }, ignoreNULL = FALSE)

    output$hint <- renderUI({
      r <- req(result())
      n <- nrow(r$res$scores)
      helpText(
        sprintf("Distance: %s, %s samples.", r$res$method$dist_label, fmt_n(n)),
        if (n > 3000) "Large data set: expect ~1 min per 250 permutations."
      )
    })

    observeEvent(input$run, {
      r <- req(result())
      req(input$variable)
      out <- tryCatch(
        withProgress(message = "PERMANOVA", value = 0, {
          setProgress(0.05, detail = "distance matrix")
          d <- analysis_distance(r$res)
          setProgress(0.15, detail = "permutations")
          run_permanova(d, r$meta, input$variable,
            permutations = as.integer(input$perms), strata = input$strata,
            progress = function(f) setProgress(0.15 + 0.8 * f)
          )
        }),
        error = function(e) {
          showNotification(conditionMessage(e), type = "error", duration = NULL)
          NULL
        }
      )
      stats_res(out)
    })

    output$result <- renderUI({
      s <- stats_res()
      if (is.null(s)) {
        return(div(class = "text-muted p-4", "Choose a variable and run the test."))
      }
      disp_warn <- !is.na(s$disp_p) && s$disp_p < 0.05
      tagList(
        layout_columns(
          value_box("R²", sprintf("%.3f", s$r2), showcase = icon("chart-pie")),
          value_box("Pseudo-F", sprintf("%.2f", s$f), showcase = icon("divide")),
          value_box("p-value", fmt_p(s$p), sprintf("%d permutations", s$permutations),
            showcase = icon("dice")
          ),
          value_box(
            "Dispersion", fmt_p(s$disp_p),
            if (disp_warn) "groups differ in spread" else if (is.na(s$disp_p)) "numeric variable" else "homogeneous",
            showcase = icon(if (disp_warn) "triangle-exclamation" else "circle-check"),
            theme = if (disp_warn) "warning" else NULL
          )
        ),
        tableOutput(session$ns("table"))
      )
    })

    output$table <- renderTable({
      s <- req(stats_res())
      cbind(Term = rownames(s$table), s$table)
    }, na = "", digits = 4)

    stats_res
  })
}
