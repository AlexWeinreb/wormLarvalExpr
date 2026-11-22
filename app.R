
library(shiny)
library(dplyr)
library(tibble)
library(ggplot2)
library(stringr)
library(markdown)

dat_norm <- readRDS("data/dat_norm.rds")
gids_298 <- readRDS("data/gids_298.rds")



ui <- fluidPage(


  titlePanel("Plot gene expression during C. elegans larval development"),

  textInput("genes_str",
            label = "Gene(s)",
            value = "dpy-6, blmp-1",
            width = '50%'),
  plotOutput("distPlot"),
  uiOutput("error_msg"),
  hr(),
  includeMarkdown("README.md")

)



server <- function(input, output) {

  #~ read input ----

  genes_str <- debounce(reactive(input$genes_str), millis = 500)

  genes <- reactive({

    genes_typed <- str_split_1(genes_str(), ",") |>
      str_trim() |>
      setdiff("")

    all_gene_names <-  replace_gene_ids(genes_typed, gids_298)

    is_found <- !is.na(all_gene_names) & (all_gene_names %in% dat_norm$gene_name)

    list(gene_names = all_gene_names[is_found],
         missing = genes_typed[!is_found],
         n_found = sum(is_found))

  })

  #~ error_msg ----
  output$error_msg <- renderUI({

    g <- genes()

    msgs <- c(
      if (length(g$missing) > 0)
        paste0(
          if (length(g$gene_names) == 0) "No gene found: " else "Genes not found: ",
          paste(g$missing, collapse = ", ")
        ),
      if (length(unique(g$gene_names)) > 6)
        paste0("More than ", 6, " genes: some are plotted without markers.")
    )

    if (length(msgs) == 0) return(NULL)

    div(class = "alert alert-danger", lapply(msgs, div))
  })




  #~ plot ----
  output$distPlot <- renderPlot({

    genes_sel <- genes()$gene_names

    max_y <- filter(dat_norm, gene_name %in% genes_sel) |>
      pull(TPM_TMM) |>
      max() |>
      (\(x) 1.1 * log10(x+1) )()


    dat_norm |>
      filter(gene_name %in% genes_sel) |>
      ggplot(aes(x = time,
                 y = log10(1+TPM_TMM),
                 linetype = gene_name,
                 shape = gene_name,
                 color = gene_name)) +
      theme_classic() +
      theme(legend.position = "none") +
      scale_fill_brewer(palette = "Pastel2") +
      geom_rect(data = stage_df,
                aes(xmin = start, xmax = end,
                    ymin = -Inf, ymax = max_y - .05,
                    fill = stage),
                inherit.aes = FALSE,
                alpha = .4,
                show.legend = FALSE) +
      geom_rect(data = stage_df,
                aes(xmin = start,
                    xmax = end,
                    ymin = max_y - .08,
                    ymax = max_y + .08,
                    fill = stage),
                inherit.aes = FALSE,
                alpha = 1,
                show.legend = FALSE) +
      annotate("text",
               x = stage_df$start + (stage_df$end - stage_df$start) / 2,
               y = max_y,
               label = stage_df$stage) +
      geom_point(size = 3, alpha = .5) +
      geom_line() +
      ggtitle("Expression during larval development (Meeuse 2020 data)") +
      xlab("Developmental time (h)") +
      ylab(expression(Expression:~log[10](1 + TPM[TMM]))) +
      ggrepel::geom_text_repel(
        aes(label = gene_name),
        data = dat_norm |>
          filter(gene_name %in% genes_sel) |>
          slice_tail(n = 1, by = gene_name),
        nudge_x = 1, hjust = 0, direction = "y"
      )

  },
  res = 110)
}


shinyApp(ui = ui, server = server)



