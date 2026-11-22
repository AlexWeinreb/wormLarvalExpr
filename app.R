
library(shiny)
library(dplyr)
library(tibble)
library(ggplot2)
library(stringr)

dat_norm <- readRDS("data/dat_norm.rds")
gids_298 <- readRDS("data/gids_298.rds")



ui <- fluidPage(


  titlePanel("Plot gene expression during C. elegans larval development"),

  textInput("genes_str", "Gene(s)", value = "C26B9.7, mam-2, F45E4.5,ham-2, Y11D7A.3"),
  plotOutput("distPlot"),
  uiOutput("error_msg")

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

    missing <- genes()$missing

    if (length(missing) == 0) return(NULL)

    prefix <- if (genes()$n_found == 0) "No gene found: " else "Genes not found: "
    div(class = "alert alert-danger", paste0(prefix, paste(missing, collapse = ", ")))
  })

  #~ plot ----
  output$distPlot <- renderPlot({

    genes_sel <- genes()$gene_names

    max_y <- filter(dat_norm, gene_name %in% genes_sel) |>
      pull(TPM_TMM) |>
      max() |>
      (\(x) 1.05 * log10(x+1) )()


    dat_norm |>
      filter(gene_name %in% genes_sel) |>
      ggplot(aes(x = time,
                 y = log10(1+TPM_TMM),
                 linetype = gene_name,
                 shape = gene_name,
                 color = gene_name)) +
      theme_classic() +
      theme(legend.position = "none") +
      scale_fill_manual(values = c("L1" = "#e8dcc8",
                                   "L2" = "#c8d8e0",
                                   "L3" = "#d8e0c8",
                                   "L4" = "#e0d0c8",
                                   "adult" = "#d4c5a9")) +
      geom_rect(data = stage_df,
                aes(xmin = start,
                    xmax = end,
                    ymin = max_y - .05,
                    ymax = max_y + .05,
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
        aes(x = 1.1 * max(stage_df$end),
            y = log10(1+TPM_TMM),
            label = gene_name,
            color = gene_name),
        data = dat_norm |>
          filter(gene_name %in% genes_sel) |>
          slice_tail(n = 1, by = gene_name),
        hjust = -1, direction = "y"
      )

  })
}


shinyApp(ui = ui, server = server)



