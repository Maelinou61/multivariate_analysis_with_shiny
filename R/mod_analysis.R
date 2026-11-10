# Analysis module: choice of method + its parameters ---------------------------------
# Returns reactive(list(method_id, params))

mod_analysis_ui <- function(id) {
  ns <- NS(id)
  accordion_panel(
    "Analysis", icon = icon("diagram-project"),
    selectInput(ns("method"), "Method", choices = method_choices()),
    conditionalPanel(
      "['pca_clr', 'umap_clr', 'pcoa_aitchison'].includes(input.method)", ns = ns,
      numericInput(ns("pseudocount"), "CLR pseudocount", value = 1, min = 1e-6, step = 0.5)
    ),
    conditionalPanel(
      "input.method.startsWith('umap')", ns = ns,
      sliderInput(ns("n_neighbors"), "Neighbours", min = 2, max = 100, value = 15),
      sliderInput(ns("min_dist"), "Min. distance", min = 0, max = 0.99, value = 0.1, step = 0.01)
    ),
    conditionalPanel(
      "input.method == 'nmds_bray'", ns = ns,
      helpText(sprintf("NMDS is limited to %s samples.", fmt_n(NMDS_MAX)))
    ),
    conditionalPanel(
      "input.method.startsWith('umap') || input.method == 'nmds_bray'", ns = ns,
      numericInput(ns("seed"), "Random seed", value = 1, min = 1, step = 1)
    )
  )
}

mod_analysis_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    reactive(list(
      method_id = input$method,
      params = list(
        pseudocount = input$pseudocount,
        n_neighbors = input$n_neighbors,
        min_dist = input$min_dist,
        seed = input$seed
      )
    ))
  })
}
