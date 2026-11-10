# Data module: example files or uploads, sample ID column, taxonomic rank -------------
# Returns list(data = reactive(list(counts, meta, id_col, log)), rank = reactive(chr))

EXAMPLE_DIR <- "data"

example_files <- function(sub) {
  f <- list.files(file.path(EXAMPLE_DIR, sub), pattern = "\\.(csv|tsv|txt)$", full.names = TRUE)
  stats::setNames(f, basename(f))
}

mod_data_ui <- function(id) {
  ns <- NS(id)
  accordion_panel(
    "Data", icon = icon("database"),
    radioButtons(ns("source"), NULL,
      choices = c("Example files" = "example", "Upload" = "upload"),
      inline = TRUE
    ),
    conditionalPanel(
      "input.source == 'example'", ns = ns,
      selectInput(ns("ex_matrix"), "Abundance matrix", example_files("matrix")),
      selectInput(ns("ex_meta"), "Metadata", example_files("metadata"))
    ),
    conditionalPanel(
      "input.source == 'upload'", ns = ns,
      fileInput(ns("up_matrix"), "Abundance matrix (taxa × samples)", accept = c(".csv", ".tsv", ".txt")),
      fileInput(ns("up_meta"), "Metadata (one row per sample)", accept = c(".csv", ".tsv", ".txt"))
    ),
    actionButton(ns("load"), "Load", icon = icon("upload"), class = "btn-outline-primary w-100 mb-3"),
    selectInput(ns("id_col"), "Sample ID column (metadata)", choices = NULL),
    selectInput(ns("rank"), "Taxonomic rank", choices = RANKS, selected = "Other")
  )
}

mod_data_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    raw <- reactiveVal(NULL)

    observeEvent(input$load, {
      if (input$source == "example") {
        m_path <- input$ex_matrix; m_name <- basename(m_path)
        d_path <- input$ex_meta
      } else {
        if (is.null(input$up_matrix) || is.null(input$up_meta)) {
          showNotification("Upload both files first.", type = "warning")
          return()
        }
        m_path <- input$up_matrix$datapath; m_name <- input$up_matrix$name
        d_path <- input$up_meta$datapath
      }

      out <- tryCatch(
        withProgress(message = "Reading files…", value = 0.1, {
          counts <- read_abundance(m_path)
          incProgress(0.6)
          meta <- read_metadata(d_path)
          list(counts = counts, meta = meta, name = m_name)
        }),
        error = function(e) {
          showNotification(conditionMessage(e), type = "error", duration = NULL)
          NULL
        }
      )
      req(out)

      guess <- guess_id_column(out$meta, rownames(out$counts))
      if (is.na(guess)) {
        showNotification("No metadata column matches the sample names: choose the ID column.", type = "warning")
      }
      updateSelectInput(session, "id_col", choices = names(out$meta), selected = guess)
      updateSelectInput(session, "rank", selected = guess_rank(out$name))
      raw(out)
      showNotification(
        sprintf("Loaded %s: %s samples × %s taxa. Now click 'Run analysis'.",
                out$name, fmt_n(nrow(out$counts)), fmt_n(ncol(out$counts))),
        type = "message", duration = 5
      )
    })

    data <- reactive({
      r <- raw()
      req(r, input$id_col %in% names(r$meta))
      m <- tryCatch(match_samples(r$counts, r$meta, input$id_col), error = function(e) e)
      validate(need(!inherits(m, "error"), if (inherits(m, "error")) conditionMessage(m)))
      list(
        counts = m$counts,
        meta = m$meta,
        id_col = input$id_col,
        log = c(list(log_dims("Loaded", r$counts, r$name)), m$log)
      )
    })

    list(data = data, rank = reactive(input$rank))
  })
}
