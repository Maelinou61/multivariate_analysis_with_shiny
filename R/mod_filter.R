# Filter module: metadata subset + abundance thresholds --------------------------------
# Returns reactive(list(meta_column, meta_values, min_depth, min_prevalence, min_count))

mod_filter_ui <- function(id) {
  ns <- NS(id)
  accordion_panel(
    "Filters", icon = icon("filter"),
    selectInput(ns("meta_column"), "Keep samples where…", choices = c("(all samples)" = "")),
    conditionalPanel(
      "input.meta_column != ''", ns = ns,
      selectizeInput(ns("meta_values"), "…is one of", choices = NULL, multiple = TRUE)
    ),
    numericInput(ns("min_depth"), "Min. reads per sample", value = 0, min = 0, step = 100),
    sliderInput(ns("min_prevalence"), "Min. prevalence of a taxon (% samples)",
      min = 0, max = 50, value = 1, step = 0.5, post = "%"
    ),
    numericInput(ns("min_count"), "Min. total reads per taxon", value = 0, min = 0, step = 10)
  )
}

mod_filter_server <- function(id, data) {
  moduleServer(id, function(input, output, session) {
    observeEvent(data(), {
      d <- data()
      updateSelectInput(session, "meta_column",
        choices = c("(all samples)" = "", filter_columns(d$meta, d$id_col))
      )
    })

    observeEvent(input$meta_column, {
      req(input$meta_column != "")
      v <- data()$meta[[input$meta_column]]
      tab <- sort(table(v), decreasing = TRUE)
      choices <- stats::setNames(names(tab), sprintf("%s (%s)", names(tab), fmt_n(as.integer(tab))))
      updateSelectizeInput(session, "meta_values", choices = choices, selected = character(0), server = TRUE)
    })

    reactive(list(
      meta_column = input$meta_column,
      meta_values = input$meta_values,
      min_depth = input$min_depth %||% 0,
      min_prevalence = (input$min_prevalence %||% 0) / 100,
      min_count = input$min_count %||% 0
    ))
  })
}
