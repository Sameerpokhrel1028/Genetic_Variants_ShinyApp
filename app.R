# ============================================================
# Variant Density + Introgression Viewer (TSV only)
# ============================================================
# Recommended workflow:
#   VCF -> scripts/parse_vcf_to_matrix.sh -> TSV -> upload TSV -> plots
# ============================================================

suppressPackageStartupMessages({
  library(shiny)
  library(dplyr)
  library(ggplot2)
  library(ape)
  library(readr)
})

source("functions.R")

# Optional demo TSV if present (your repo can include one)
default_df <- NULL
if (file.exists("10kbvariants.tsv")) {
  default_df <- readr::read_tsv("10kbvariants.tsv", show_col_types = FALSE)
  validate_variant_tsv(default_df)
}

ui <- fluidPage(
  titlePanel("Variant Density Heatmap & Genome Clusters (Introgression Viewer)"),

  tabsetPanel(
    tabPanel(
      "Data Upload",
      sidebarLayout(
        sidebarPanel(
          h4("Upload TSV (recommended workflow)"),
          fileInput("tsv_upload", "Upload variant matrix (TSV)",
                    accept = c(".tsv", ".txt", ".csv")),
          tags$p(style="font-size:12px;",
                 "TSV must have columns: Chrom, Start, End, then sample columns."),

          hr(),
          h4("Introgression-focused filtering (optional)"),
          checkboxInput("enable_filter", "Enable hypervariable-window filtering", value = FALSE),
          selectInput("filter_method", "Filtering method",
                      choices = c(
                        "Median across samples (recommended)" = "median",
                        "Max across samples (strict)" = "max",
                        "Proportion of samples above cutoff" = "prop"
                      ),
                      selected = "median"),
          sliderInput("max_per_kb", "Cutoff (variants per kb)", min = 5, max = 200, value = 40, step = 5),
          sliderInput("prop_thresh", "If proportion method: drop if ≥ this fraction exceed cutoff",
                      min = 0.05, max = 1.00, value = 0.25, step = 0.05),

          hr(),
          checkboxInput("enable_cap", "Cap extreme values for plotting (keeps windows)", value = FALSE),
          sliderInput("cap_per_kb", "Cap at (variants per kb)", min = 10, max = 300, value = 100, step = 10),

          hr(),
          h4("Export"),
          tags$p(style="font-size:12px;",
                 "Downloads the current matrix after conversion to variants/kb and after any filtering/capping."),
          downloadButton("download_processed_tsv", "Download processed TSV")
        ),

        mainPanel(
          h4("Current dataset summary"),
          verbatimTextOutput("data_summary"),
          tags$hr(),
          tags$p("Tip: Use 50–100 kb windows for smoother introgression blocks; 10 kb gives finer resolution but more noise.")
        )
      )
    ),

    tabPanel(
      "Variants Heatmap",
      sidebarLayout(
        sidebarPanel(
          selectInput("chrom_AB", "Select Chromosome:", choices = character(0)),
          selectInput("samples_AB", "Select Samples:", choices = character(0),
                      selected = NULL, multiple = TRUE)
        ),
        mainPanel(
          plotOutput("heatmap_AB", height = "800px")
        )
      )
    ),

    tabPanel(
      "Genome Clusters",
      sidebarLayout(
        sidebarPanel(
          selectInput("chrom_div", "Select Chromosome:",
                      choices = c("All chromosomes" = "All")),
          selectInput("genomes_div", "Select Samples:",
                      choices = c("All samples" = "All"),
                      selected = "All", multiple = TRUE),
          tags$p(style="color:black; font-size:13px; margin-top:-10px;",
                 "Select 3 or more samples for dendrogram and PCA.")
        ),
        mainPanel(
          h4("Dendrogram (Hierarchical Clustering)"),
          plotOutput("div_dendrogram", height = "350px"),
          h4("PCA of Sample Divergence"),
          plotOutput("div_pca", height = "400px")
        )
      )
    ),

    tabPanel(
      "About",
      h3("About this app"),
      tags$p(style="font-size:16px; font-weight:600;",
             "This Shiny app visualizes window-based, sample-wise variant density to support introgression discovery and clustering."),
      tags$ul(
        tags$li(tags$b("Input: "), "TSV matrix with columns Chrom/Start/End + sample columns (counts per window)."),
        tags$li(tags$b("Computation: "), "App converts counts to variants/kb using window length, then optionally filters/caps."),
        tags$li(tags$b("Outputs: "), "Heatmap, dendrogram clustering, and PCA.")
      ),
      tags$p("Workflow details and the VCF→TSV parser script are on GitHub:"),
      tags$p(tags$a(
        href = "https://github.com/Sameerpokhrel1028/Genetic_Variants_ShinyApp",
        target = "_blank",
        "Genetic_Variants_ShinyApp (GitHub)"
      ))
    )
  )
)

server <- function(input, output, session) {

  rv <- reactiveValues(df = default_df)

  # ---- Load TSV ----
  observeEvent(input$tsv_upload, {
    req(input$tsv_upload)
    newdf <- readr::read_tsv(input$tsv_upload$datapath, show_col_types = FALSE)
    validate_variant_tsv(newdf)
    rv$df <- newdf
    showNotification("TSV loaded successfully.", type = "message")
  })

  # ---- Derived sample names + chromosomes (universal) ----
  samples_avail <- reactive({
    req(rv$df)
    get_sample_cols(rv$df)
  })

  chroms_avail <- reactive({
    req(rv$df)
    sort(unique(rv$df$Chrom))
  })

  # Update UI dropdowns whenever dataset changes
  observe({
    req(rv$df)
    s <- samples_avail()
    cset <- chroms_avail()

    updateSelectInput(session, "chrom_AB",
                      choices = cset,
                      selected = cset[1])

    updateSelectInput(session, "samples_AB",
                      choices = s,
                      selected = head(s, 4))

    updateSelectInput(session, "chrom_div",
                      choices = c("All" = "All", cset),
                      selected = "All")

    updateSelectInput(session, "genomes_div",
                      choices = c("All" = "All", s),
                      selected = "All")
  })

  # ---- Processed data: variants/kb + optional filter/cap ----
  processed_df <- reactive({
    req(rv$df)
    df0 <- rv$df
    g <- get_sample_cols(df0)

    df_kb <- to_var_per_kb(df0, g)

    if (isTRUE(input$enable_filter)) {
      df_kb <- apply_window_filter(
        df_kb, g,
        method = input$filter_method,
        cutoff = input$max_per_kb,
        prop_thresh = input$prop_thresh
      )
    }

    if (isTRUE(input$enable_cap)) {
      df_kb <- cap_values(df_kb, g, input$cap_per_kb)
    }

    df_kb
  })

  # ---- Download processed TSV ----
  output$download_processed_tsv <- downloadHandler(
    filename = function() {
      ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
      filt <- if (isTRUE(input$enable_filter)) paste0("filter_", input$filter_method, "_", input$max_per_kb, "perkB") else "nofilter"
      cap  <- if (isTRUE(input$enable_cap)) paste0("cap_", input$cap_per_kb, "perkB") else "nocap"
      paste0("variant_density_processed_", filt, "_", cap, "_", ts, ".tsv")
    },
    content = function(file) {
      readr::write_tsv(processed_df(), file)
    }
  )

  # ---- Summary ----
  output$data_summary <- renderPrint({
    if (is.null(rv$df)) {
      cat("No dataset loaded yet.\nUpload a TSV generated by the parser script.\n")
      return()
    }
    df0 <- rv$df
    s <- get_sample_cols(df0)
    win_bp <- median(pmax(df0$End - df0$Start, 1), na.rm = TRUE)

    cat("Rows (windows):", nrow(df0), "\n")
    cat("Samples:", length(s), "\n")
    cat("Chromosomes:", length(unique(df0$Chrom)), "\n")
    cat("Median window size (bp):", win_bp, "\n")
    cat("Filtering:", ifelse(isTRUE(input$enable_filter), "ON", "OFF"), "\n")
    if (isTRUE(input$enable_filter)) {
      cat("  Method:", input$filter_method, "\n")
      cat("  Cutoff (var/kb):", input$max_per_kb, "\n")
      if (input$filter_method == "prop") cat("  Prop threshold:", input$prop_thresh, "\n")
    }
    cat("Capping:", ifelse(isTRUE(input$enable_cap), "ON", "OFF"), "\n")
    if (isTRUE(input$enable_cap)) cat("  Cap (var/kb):", input$cap_per_kb, "\n")
  })

  # ---- Heatmap ----
  output$heatmap_AB <- renderPlot({
    req(rv$df)
    req(input$chrom_AB)
    req(input$samples_AB)

    chr_data <- processed_df() %>% filter(Chrom == input$chrom_AB)
    draw_heatmap(chr_data, input$samples_AB)
  })

  # ---- Divergence + clustering ----
  calc_divergence <- reactive({
    req(rv$df)
    df_kb <- processed_df()

    s_all <- get_sample_cols(df_kb)
    selected <- if ("All" %in% input$genomes_div || length(input$genomes_div) == 0) s_all else input$genomes_div
    chr_scope <- if (input$chrom_div == "All" || input$chrom_div == "") chroms_avail() else input$chrom_div

    df_sub <- df_kb %>% filter(Chrom %in% chr_scope)
    compute_divergence(df_sub, selected)
  })

  output$div_dendrogram <- renderPlot({
    div <- calc_divergence()
    if (nrow(div) < 3) {
      plot.new()
      text(0.5, 0.5, "Select at least 3 samples for clustering", cex = 1.2)
      return()
    }
    tree <- hclust(as.dist(div))
    plot(tree, main = "Sample Divergence Clustering (Dendrogram)", xlab = "", sub = "")
  })

  output$div_pca <- renderPlot({
    div <- calc_divergence()
    if (nrow(div) < 3) {
      plot.new()
      text(0.5, 0.5, "Select at least 3 samples for PCA", cex = 1.2)
      return()
    }
    pca <- cmdscale(as.dist(div), eig = TRUE, k = 2)
    pca_df <- data.frame(
      Sample = rownames(div),
      PC1 = pca$points[, 1],
      PC2 = pca$points[, 2]
    )

    ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample)) +
      geom_point(size = 3) +
      geom_text(vjust = -0.9, size = 4) +
      labs(title = "Sample Divergence PCA", x = "PC1", y = "PC2") +
      theme_minimal()
  })
}

shinyApp(ui = ui, server = server)
