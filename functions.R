suppressPackageStartupMessages({
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# ---------- Validation ----------
validate_variant_tsv <- function(df) {
  req_cols <- c("Chrom", "Start", "End")
  if (!all(req_cols %in% colnames(df))) {
    stop("TSV must contain columns: Chrom, Start, End")
  }
  if (!is.numeric(df$Start) || !is.numeric(df$End)) {
    stop("Start and End columns must be numeric")
  }
  if (ncol(df) < 4) {
    stop("TSV must have at least one sample column after Chrom/Start/End")
  }
  invisible(TRUE)
}

get_sample_cols <- function(df) {
  setdiff(colnames(df), c("Chrom", "Start", "End"))
}

# ---------- Universal conversion: counts -> variants/kb ----------
to_var_per_kb <- function(df_in, sample_cols) {
  win_bp <- pmax(df_in$End - df_in$Start, 1)
  win_kb <- win_bp / 1000
  df_out <- df_in
  df_out[sample_cols] <- sweep(as.matrix(df_out[sample_cols]), 1, win_kb, "/")
  df_out
}

# ---------- Filtering ----------
apply_window_filter <- function(df_kb, sample_cols, method, cutoff, prop_thresh) {
  if (nrow(df_kb) == 0) return(df_kb)
  x <- as.matrix(df_kb[sample_cols])

  keep <- switch(
    method,
    "median" = apply(x, 1, median, na.rm = TRUE) <= cutoff,
    "max"    = apply(x, 1, max,    na.rm = TRUE) <= cutoff,
    "prop"   = {
      frac_above <- apply(x, 1, function(v) mean(v > cutoff, na.rm = TRUE))
      frac_above < prop_thresh
    },
    rep(TRUE, nrow(df_kb))
  )

  df_kb[keep, , drop = FALSE]
}

cap_values <- function(df_kb, sample_cols, cap_per_kb) {
  if (nrow(df_kb) == 0) return(df_kb)
  df2 <- df_kb
  df2[sample_cols] <- lapply(df2[sample_cols], function(v) pmin(v, cap_per_kb))
  df2
}

# ---------- Divergence ----------
compute_divergence <- function(df_kb, genomes) {
  mat <- df_kb %>% select(all_of(genomes)) %>% as.matrix()
  divergence <- matrix(NA_real_,
                       nrow = length(genomes),
                       ncol = length(genomes),
                       dimnames = list(genomes, genomes))
  for (i in seq_along(genomes)) {
    for (j in seq_along(genomes)) {
      divergence[i, j] <- mean(abs(mat[, i] - mat[, j]), na.rm = TRUE)
    }
  }
  divergence
}

# ---------- Heatmap ----------
draw_heatmap <- function(chr_data_kb, samples) {
  if (nrow(chr_data_kb) == 0 || length(samples) == 0) {
    plot.new()
    text(0.5, 0.5, "No data to plot", cex = 1.2)
    return(invisible(NULL))
  }

  mat <- chr_data_kb %>% select(all_of(samples)) %>% as.matrix()
  mat_log <- log10(mat + 1)

  coords_mb <- chr_data_kb$Start / 1e6
  rounded <- round(coords_mb)
  label_info <- tibble(row = seq_along(rounded), Mb = rounded) %>%
    filter(Mb %% 10 == 0) %>%
    distinct(Mb, .keep_all = TRUE)

  mb_axis <- rowAnnotation(
    Mb = anno_mark(
      at = label_info$row,
      labels = paste0(label_info$Mb, " Mb"),
      side = "left",
      labels_gp = gpar(fontsize = 9, col = "black"),
      link_width = unit(5, "mm")
    ),
    width = unit(4.0, "cm")
  )

  rng <- range(mat_log, finite = TRUE)
  mid <- mean(rng)
  col_fun <- colorRamp2(c(rng[1], mid, rng[2]), c("white", "orange", "red"))

  ht <- Heatmap(
    mat_log,
    name = "log10(var/kb + 1)",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    left_annotation = mb_axis,
    column_names_rot = 45,
    column_title = paste0("Variant density – ", unique(chr_data_kb$Chrom))
  )

  draw(ht, heatmap_legend_side = "bottom")
}
