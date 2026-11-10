# Drivers module: taxa contributing most to the displayed axes ---------------------------

DRIVER_NOTES <- c(
  pca  = "Contribution = squared loading × 100 (sums to 100 per axis).",
  ca   = "Contribution = column weight × squared loading × 100 (sums to 100 per axis).",
  pcoa = "Score = covariance with the axis; contribution = share of total |covariance|.",
  nmds = "Score = Spearman correlation with the axis; contribution = share of total |correlation|.",
  umap = "Score = Spearman correlation with the axis; contribution = share of total |correlation|."
)

mod_drivers_ui <- function(id) {
  ns <- NS(id)
  tagList(
    div(
      class = "d-flex justify-content-between align-items-center mb-2",
      textOutput(ns("note"), inline = TRUE) |> tagAppendAttributes(class = "text-muted small"),
      downloadButton(ns("dl"), "CSV", class = "btn-sm btn-outline-secondary")
    ),
    DT::DTOutput(ns("table"))
  )
}

mod_drivers_server <- function(id, result, axes) {
  moduleServer(id, function(input, output, session) {
    drivers <- reactive({
      r <- req(result())
      compute_drivers(r$res, r$counts, axes())
    })

    output$note <- renderText(DRIVER_NOTES[[req(result())$res$method$engine]])

    output$table <- DT::renderDT({
      d <- drivers()
      DT::datatable(d, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE)) |>
        DT::formatSignif(2:3, 3) |>
        DT::formatRound(4:7, 2)
    })

    output$dl <- downloadHandler(
      filename = function() paste0(result()$res$method_id, "_drivers.csv"),
      content = function(file) utils::write.csv(drivers(), file, row.names = FALSE)
    )
  })
}
