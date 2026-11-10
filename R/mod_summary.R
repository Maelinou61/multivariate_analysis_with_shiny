# Processing summary: a strip of small boxes, one per pipeline step --------------------

mod_summary_ui <- function(id) {
  uiOutput(NS(id, "strip"))
}

mod_summary_server <- function(id, log) {
  moduleServer(id, function(input, output, session) {
    output$strip <- renderUI({
      entries <- log()
      if (!length(entries)) {
        return(div(class = "pipeline pipeline-empty", icon("arrow-left"),
                   "1. Choose files and click 'Load'  ·  2. Click 'Run analysis'"))
      }
      boxes <- lapply(entries, function(e) {
        div(
          class = "pipeline-step", title = e$tooltip %||% "",
          div(class = "pipeline-label", e$step),
          div(class = "pipeline-value", e$value),
          if (!is.null(e$caption)) div(class = "pipeline-caption", e$caption)
        )
      })
      arrows <- lapply(seq_along(boxes), function(i) {
        if (i < length(boxes)) list(boxes[[i]], span(class = "pipeline-arrow", "→")) else boxes[[i]]
      })
      div(class = "pipeline", arrows)
    })
  })
}
