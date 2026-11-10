# Ordination plot module (+ scree plot) -------------------------------------------------
# `result` : reactive(list(res, counts, meta, cols, ...)) or NULL
# `rank`   : reactive taxonomic rank (used in titles only)
# Returns reactive(c(axis_x, axis_y))

mod_plot_ui <- function(id) {
  ns <- NS(id)
  layout_sidebar(
    sidebar = sidebar(
      position = "right", width = 260, open = "desktop",
      selectInput(ns("colour"), "Colour by", choices = c("None" = "")),
      conditionalPanel(
        "output.colour_is_categorical", ns = ns,
        selectizeInput(ns("highlight"), "Highlight only", choices = NULL, multiple = TRUE,
          options = list(placeholder = "all groups")
        ),
        checkboxInput(ns("ellipse"), "95% ellipses", FALSE)
      ),
      selectInput(ns("shape"), "Shape by", choices = c("None" = "")),
      layout_columns(
        selectInput(ns("axis_x"), "X axis", choices = 1:2, selected = 1),
        selectInput(ns("axis_y"), "Y axis", choices = 1:2, selected = 2)
      ),
      sliderInput(ns("size"), "Point size", min = 0.5, max = 5, value = 1.8, step = 0.1),
      sliderInput(ns("alpha"), "Opacity", min = 0.1, max = 1, value = 0.8, step = 0.05),
      downloadButton(ns("dl_png"), "Plot (PNG)", class = "btn-sm btn-outline-secondary w-100"),
      downloadButton(ns("dl_scores"), "Coordinates (CSV)", class = "btn-sm btn-outline-secondary w-100")
    ),
    plotly::plotlyOutput(ns("plot"), height = "620px")
  )
}

mod_plot_server <- function(id, result, rank) {
  moduleServer(id, function(input, output, session) {
    observeEvent(result(), {
      r <- result()
      cols <- r$cols
      cur_colour <- input$colour %||% ""
      cur_shape <- input$shape %||% ""
      keep_colour <- if (cur_colour %in% c(cols$categorical, cols$numeric)) cur_colour else ""
      updateSelectInput(session, "colour",
        choices = c("None" = "", list(Categorical = cols$categorical, Numeric = cols$numeric)),
        selected = keep_colour
      )
      updateSelectInput(session, "shape",
        choices = c("None" = "", cols$categorical),
        selected = if (cur_shape %in% cols$categorical) cur_shape else ""
      )
      ax <- stats::setNames(seq_len(ncol(r$res$scores)), colnames(r$res$scores))
      updateSelectInput(session, "axis_x", choices = ax, selected = 1)
      updateSelectInput(session, "axis_y", choices = ax, selected = 2)
    })

    colour_is_cat <- reactive({
      r <- result()
      !is.null(r) && isTRUE((input$colour %||% "") %in% r$cols$categorical)
    })
    output$colour_is_categorical <- reactive(colour_is_cat())
    outputOptions(output, "colour_is_categorical", suspendWhenHidden = FALSE)

    observeEvent(list(input$colour, result()), {
      if (!colour_is_cat()) return()
      v <- result()$meta[[input$colour]]
      tab <- sort(table(v), decreasing = TRUE)
      updateSelectizeInput(session, "highlight",
        choices = stats::setNames(names(tab), sprintf("%s (%s)", names(tab), fmt_n(as.integer(tab)))),
        selected = character(0)
      )
    })

    axes <- reactive({
      req(result())
      a <- as.integer(c(input$axis_x, input$axis_y))
      n <- ncol(result()$res$scores)
      if (anyNA(a) || any(a > n)) a <- c(1L, 2L)
      a
    })

    gg <- reactive({
      r <- req(result())
      make_ordination_plot(
        r$res, r$meta, rank(),
        axes = axes(),
        colour = input$colour,
        shape = input$shape,
        highlight = if (colour_is_cat()) input$highlight else NULL,
        ellipse = isTRUE(input$ellipse) && colour_is_cat(),
        size = input$size, alpha = input$alpha
      )
    })

    output$plot <- plotly::renderPlotly({
      p <- gg()
      title <- sprintf("<b>%s</b><br><sup>%s</sup>", p$labels$title, p$labels$subtitle)
      p$labels$title <- NULL
      p$labels$subtitle <- NULL
      w <- plotly::ggplotly(p, tooltip = "text") |>
        plotly::layout(title = list(text = title, x = 0.02), margin = list(t = 70)) |>
        plotly::config(displaylogo = FALSE, toImageButtonOptions = list(format = "svg"))
      if (nrow(result()$res$scores) > 2000) w <- plotly::toWebGL(w)
      w
    })

    file_stub <- function() {
      r <- result()
      colour <- input$colour %||% ""
      paste(c(r$res$method_id, tolower(rank()), if (colour != "") colour), collapse = "_")
    }

    output$dl_png <- downloadHandler(
      filename = function() paste0(file_stub(), ".png"),
      content = function(file) ggplot2::ggsave(file, gg(), width = 8, height = 6, dpi = 300)
    )
    output$dl_scores <- downloadHandler(
      filename = function() paste0(file_stub(), "_coordinates.csv"),
      content = function(file) {
        r <- result()
        utils::write.csv(cbind(r$meta, r$res$scores), file, row.names = FALSE)
      }
    )

    axes
  })
}

# Scree plot -----------------------------------------------------------------------------

mod_scree_ui <- function(id) plotOutput(NS(id, "scree"), height = "420px")

mod_scree_server <- function(id, result) {
  moduleServer(id, function(input, output, session) {
    output$scree <- renderPlot({
      r <- req(result())
      validate(need(!is.null(r$res$var_expl), "No variance per axis for this method."))
      make_scree_plot(r$res)
    }, res = 110)
  })
}
