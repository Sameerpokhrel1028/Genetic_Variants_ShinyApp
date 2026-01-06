# ============================================================
# Peanut-PanMAGIC Variant Density Shiny App
# ============================================================
# This Shiny app visualizes:
#  1) Variant density heatmap (variants per kb in 10 kb windows)
#  2) Genome clustering using a dendrogram and PCA
#
# Tabs:
#  - Variants Heatmap
#  - Genome Clusters
#  - About
#
# Inputs:
#  - 10Kbvariants.tsv : variant counts per 10 kb windows (standardized to VarPerKb)
#  - functions.R      : helper functions used by the app (required)
#
# Notes:
#  - This app assumes you run it from the project directory (the folder that
#    contains app.R, functions.R, and 10Kbvariants.tsv).
#  - Avoid using setwd() in Shiny. Instead, open the R project (or set the
#    working directory in RStudio) to the app folder before running.
#
# Install required packages (one-time):
#   install.packages("shiny")
#   install.packages("tidyverse")
#   install.packages("ape")
#   install.packages("circlize")
#   install.packages("Cairo")
#   install.packages("grid")
#   install.packages("BiocManager")
#   BiocManager::install("ComplexHeatmap")
#
# Run locally:
#   shiny::runApp()
#
# ============================================================

suppressPackageStartupMessages({
  library(shiny)
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(Cairo)
  library(grid)
  library(ape)
})

# Helper functions live in functions.R (keeps main computations out of app.R)
source("functions.R")

# Input table: variant counts per 10 kb windows (converted to variants per kb)
# Keep the file in the same directory as app.R
df <- read_tsv("10Kbvariants.tsv", show_col_types = FALSE)

genomes <- c("Bailey", "C431", "C99R", "CC41", "CC477", "CC812", "Florida07",
             "GA12Y", "GPNCWS17", "Georganic", "IAC322", "ICG1471", "Lariat",
             "Marc1", "NC94022", "NMVal.1", "NMVal.2", "TifNV", "York")

chromosomes <- paste0("chr", str_pad(1:20, width = 2, pad = "0"))

# Hypervariable windows filter (e.g., >40 variants/kb) to reduce likely artifacts
max_per_kb <- 40

# ---------------- UI ----------------

ui <- fluidPage(
  titlePanel("Variant Density Heatmap & Genome Clusters"),
  tabsetPanel(
    tabPanel(
      "Variants Heatmap",
      sidebarLayout(
        sidebarPanel(
          selectInput("chrom_AB", "Select Chromosome:",
                      choices = chromosomes, selected = "chr01"),
          selectInput("samples_AB", "Select Genomes:", choices = genomes,
                      selected = genomes[1:4], multiple = TRUE)
        ),
        mainPanel(plotOutput("heatmap_AB", height = "800px"))
      )
    ),

    tabPanel(
      "Genome Clusters",
      sidebarLayout(
        sidebarPanel(
          selectInput("chrom_div", "Select Chromosome:",
                      choices = c("All chromosomes" = "All", chromosomes),
                      selected = "All"),
          selectInput("genomes_div", "Select Genomes:",
                      choices = c("All genomes" = "All", genomes),
                      selected = "All", multiple = TRUE),
          tags$p(style = "color:black; font-size:13px; margin-top:-10px;",
                 "Select 3 or more genomes for dendrogram and PCA.")
        ),
        mainPanel(
          h4("Dendrogram (Hierarchical Clustering)"),
          plotOutput("div_dendrogram", height = "350px"),
          h4("PCA of Genome Divergence"),
          plotOutput("div_pca", height = "400px")
        )
      )
    ),

    tabPanel(
      "About",
      h3("About this app"),
      tags$p(
        style = "font-size:16px; font-weight:600;",
        "This Shiny application visualizes variant density and clusters peanut genomes using window-based variant counts."
      ),
      tags$ul(
        tags$li(strong("Variants Heatmap:"), " Variant density across selected genomes for a selected chromosome."),
        tags$li(strong("Genome Clusters:"), " Dendrogram and PCA based on similarity of variant density patterns."),
        tags$li(strong("About:"), " Summary of the app, inputs, and outputs.")
      ),
      tags$p(
        "Input: variant counts per 10 kb windows for peanut genomes (standardized to variants per kb). ",
        "Output: heatmaps, PCA, and dendrograms to visualize genome-wide divergence patterns."
      )
    )
  )
)

# ---------------- SERVER ----------------

server <- function(input, output) {

  # Apply hypervariable filter (used for heatmap, dendrogram, and PCA)
  filtered_df <- reactive({
    filter_data(df, genomes, max_per_kb)
  })

  # Variant density heatmap
  output$heatmap_AB <- renderPlot({
    chr_data_kb <- filtered_df() %>% filter(Chrom == input$chrom_AB)
    draw_heatmap(chr_data_kb, input$samples_AB)
  })

  # Divergence calculation (using filtered data)
  calc_divergence <- reactive({
    selected_genomes <- if ("All" %in% input$genomes_div || length(input$genomes_div) == 0) genomes else input$genomes_div
    chr_scope <- if (input$chrom_div == "All" || input$chrom_div == "") chromosomes else input$chrom_div

    filtered_data <- filtered_df() %>%
      filter(Chrom %in% chr_scope)

    compute_divergence(filtered_data, selected_genomes)
  })

  # Dendrogram
  output$div_dendrogram <- renderPlot({
    div <- calc_divergence()
    if (nrow(div) < 2) {
      plot.new()
      text(0.5, 0.5, "Select at least 3 genomes for clustering", cex = 1.2)
      return()
    }
    tree <- hclust(as.dist(div))
    plot(tree, main = "Genome Divergence Clustering (Dendrogram)", xlab = "", sub = "")
  })

  # PCA
  output$div_pca <- renderPlot({
    div <- calc_divergence()
    if (nrow(div) < 2) {
      plot.new()
      text(0.5, 0.5, "Select at least 3 genomes for PCA", cex = 1.2)
      return()
    }
    pca <- cmdscale(as.dist(div), eig = TRUE, k = 2)
    pca_df <- data.frame(
      Genome = rownames(div),
      PC1 = pca$points[, 1],
      PC2 = pca$points[, 2]
    )

    ggplot(pca_df, aes(x = PC1, y = PC2, label = Genome)) +
      geom_point(color = "blue", size = 3) +
      geom_text(vjust = -1, size = 4) +
      labs(title = "Genome Divergence PCA Plot",
           x = "PC1", y = "PC2") +
      theme_minimal()
  })
}

shinyApp(ui = ui, server = server)

# Deployment note:
# Keep rsconnect tokens/secrets out of public repos. If you deploy, do it locally
# and use environment variables or .Renviron (ignored by git) for credentials.
