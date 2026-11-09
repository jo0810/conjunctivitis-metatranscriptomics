# Figure 3: antimicrobial resistance gene landscape.
#   A  ARGs detected per center, stacked by primary drug class (cross-center pool excluded)
#   B  prevalence and mean abundance of the 35 high-frequency ARGs, with drug-class strip and per-center heatmap
#   C  Spearman correlation between candidate ocular pathogen species and the same ARGs; outlined cells mark
#      species-ARG pairs found on the same assembled contig; stars: two-sided P < 0.05 / 0.01 / 0.001
# Run the three Python scripts in figure3/ first, then from the repository root: Rscript figure3/figure3.R
# Output: output/figure3/Figure3.pdf (174 mm wide).

suppressPackageStartupMessages({
  library(grid)
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(patchwork)
})

# Stroke widths are kept within 0.5-1.5 pt (grid lwd 0.68-1.96).
CRM_LWD_MIN <- 0.68   # = 0.51 pt
CRM_LWD_MAX <- 1.96   # = 1.47 pt

clamp_lwd_grob <- function(g) {
  if (!is.null(g$gp) && !is.null(g$gp$lwd)) {
    g$gp$lwd <- pmin(pmax(g$gp$lwd, CRM_LWD_MIN), CRM_LWD_MAX)
  }
  if (!is.null(g$grobs)) {                      # gtable
    g$grobs <- lapply(g$grobs, clamp_lwd_grob)
  }
  if (!is.null(g$children)) {                   # gTree
    nm <- names(g$children)
    kids <- lapply(g$children, clamp_lwd_grob)
    names(kids) <- nm
    class(kids) <- "gList"
    g$children <- kids
  }
  g
}

# Inputs are read from data/ and from output/figure3/ (tables written by the Python scripts in figure3/);
# the figure is written to output/figure3/.
in_file <- function(f) file.path("data", f)
OUT <- file.path("output", "figure3"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
FIG3_IN <- OUT

# Arial is embedded through the quartz device on macOS; other platforms fall back to cairo_pdf.
open_pdf <- function(file, width, height) {
  if (capabilities("aqua")) {
    quartz(file = file, type = "pdf", width = width, height = height, family = "Arial", bg = "white")
  } else {
    cairo_pdf(file, width = width, height = height, family = "Arial", bg = "white")
  }
}

FONT_TAG           <- 6
TAG_X_LEFT         <- 3.6
FONT_AXIS_TITLE    <- 6
FONT_AXIS_TICK     <- 6
FONT_LEGEND_TITLE  <- 6
FONT_LEGEND_TEXT   <- 6
FONT_ANNOT         <- 6

# ============================================================================
# Panel A: ARGs detected per center by primary drug class
# ============================================================================
make_panel_A <- function() {

  summary_file <- file.path(FIG3_IN, "arg_feature_summary.tsv")
  dpm_file     <- file.path(FIG3_IN, "arg_dpm_matrix.tsv")
  col_file     <- in_file("library_metadata.csv")

  arg_summary <- read.delim(summary_file, check.names = FALSE, stringsAsFactors = FALSE)
  dpm         <- read.delim(dpm_file,     check.names = FALSE, stringsAsFactors = FALSE)
  col_anno    <- read.csv(col_file,       check.names = FALSE, stringsAsFactors = FALSE)

  sample_cols <- intersect(colnames(dpm), col_anno$sample)
  stopifnot(length(sample_cols) == 34)
  dpm_values  <- as.matrix(dpm[, sample_cols, drop = FALSE])
  mode(dpm_values) <- "numeric"
  rownames(dpm_values) <- dpm$feature_id

  class_order <- c(
    "fluoroquinolone", "aminoglycoside", "beta-lactam",
    "multiple_drug_classes", "multidrug_efflux_or_multiclass",
    "tetracycline", "MLS", "peptide/polymyxin", "phenicol",
    "rifamycin", "aminocoumarin", "nitroimidazole",
    "phosphonic_acid", "disinfectant/antiseptic", "diaminopyrimidine",
    "mupirocin-like", "antibacterial_free_fatty_acids"
  )
  class_palette <- c(
    "fluoroquinolone" = "#F8766D", "aminoglycoside" = "#EEC643",
    "beta-lactam" = "#F6D55C", "multiple_drug_classes" = "#F4A460",
    "multidrug_efflux_or_multiclass" = "#00A6D6", "tetracycline" = "#1F3B73",
    "MLS" = "#FFF2A8", "peptide/polymyxin" = "#FDB462", "phenicol" = "#F4C2A4",
    "rifamycin" = "#D96C9F", "aminocoumarin" = "#9E8BD9", "nitroimidazole" = "#A33B62",
    "phosphonic_acid" = "#9BD4E6", "disinfectant/antiseptic" = "#C7C7C7",
    "diaminopyrimidine" = "#8E63A9", "mupirocin-like" = "#8C5A33",
    "antibacterial_free_fatty_acids" = "#E5B7D4", "Other" = "#BDBDBD"
  )
  class_labels <- c(
    "fluoroquinolone" = "Fluoroquinolone", "aminoglycoside" = "Aminoglycoside",
    "beta-lactam" = "Beta-lactam", "multiple_drug_classes" = "Multiclass",
    "multidrug_efflux_or_multiclass" = "Efflux/multiclass", "tetracycline" = "Tetracycline",
    "MLS" = "MLS", "peptide/polymyxin" = "Peptide/polymyxin", "phenicol" = "Phenicol",
    "rifamycin" = "Rifamycin", "aminocoumarin" = "Aminocoumarin",
    "nitroimidazole" = "Nitroimidazole", "phosphonic_acid" = "Phosphonic acid",
    "disinfectant/antiseptic" = "Disinfectant", "diaminopyrimidine" = "Diaminopyrimidine",
    "mupirocin-like" = "Mupirocin-like", "antibacterial_free_fatty_acids" = "Fatty acids"
  )

  center_order  <- c("SHGH", "SY", "HG")
  center_labels <- c("SHGH" = "Center 1", "SY" = "Center 2", "HG" = "Center 3")

  center_feature_counts <- do.call(rbind, lapply(center_order, function(cc) {
    center_samples_use <- intersect(col_anno$sample[col_anno$center_code == cc], colnames(dpm_values))
    detected_features  <- rownames(dpm_values)[rowSums(dpm_values[, center_samples_use, drop = FALSE] > 0, na.rm = TRUE) > 0]
    tmp <- arg_summary[arg_summary$feature_id %in% detected_features, c("feature_id", "Primary_Drug_Class")]
    tmp <- tmp[!duplicated(tmp$feature_id), , drop = FALSE]
    tab <- table(factor(tmp$Primary_Drug_Class, levels = class_order))
    data.frame(center_code = cc, center = center_labels[[cc]],
               Primary_Drug_Class = names(tab), count = as.integer(tab),
               stringsAsFactors = FALSE)
  }))
  center_feature_counts <- center_feature_counts[center_feature_counts$count > 0, , drop = FALSE]
  center_feature_counts$class_rank  <- match(center_feature_counts$Primary_Drug_Class, class_order)
  center_feature_counts$center_rank <- match(center_feature_counts$center_code, center_order)
  center_feature_counts <- center_feature_counts[order(center_feature_counts$center_rank, center_feature_counts$class_rank), ]

  center_feature_totals <- aggregate(count ~ center_code + center, center_feature_counts, sum)
  center_feature_totals <- center_feature_totals[match(center_order, center_feature_totals$center_code), ]

  draw_panel_A <- function() {
    font_main <- "Arial"
    plot_x <- 0.155
    plot_y <- 0.155
    plot_w <- 0.315
    plot_h <- 0.590
    y_max   <- 160
    y_ticks <- c(0, 50, 100, 150)
    map_y <- function(value) plot_y + plot_h * value / y_max

    grid.rect(x = plot_x, y = plot_y, width = plot_w, height = plot_h,
              just = c("left", "bottom"), gp = gpar(fill = NA, col = "#222222", lwd = 1.15))

    tick_len <- 0.007
    for (yt in y_ticks) {
      yy <- map_y(yt)
      grid.lines(x = c(plot_x - tick_len, plot_x), y = c(yy, yy),
                 gp = gpar(col = "#222222", lwd = 0.80))
      grid.text(as.character(yt), x = plot_x - 0.012, y = yy, just = "right",
                gp = gpar(fontfamily = font_main, fontsize = FONT_AXIS_TICK, col = "#222222"))
    }

    grid.text("Number of ARG features", x = plot_x - 0.088, y = plot_y + plot_h / 2,
              rot = 90, gp = gpar(fontfamily = font_main, fontsize = FONT_AXIS_TITLE, col = "black"))

    bar_w <- plot_w * 0.235
    gap_w <- plot_w * 0.105
    total_bar_w <- 3 * bar_w + 2 * gap_w
    start_x <- plot_x + (plot_w - total_bar_w) / 2
    bar_left <- start_x + c(0, bar_w + gap_w, 2 * (bar_w + gap_w))
    bar_center <- bar_left + bar_w / 2
    stack_order <- rev(class_order)

    for (i in seq_along(center_order)) {
      cc <- center_order[i]
      bottom <- 0
      for (cls in stack_order) {
        val <- center_feature_counts$count[
          center_feature_counts$center_code == cc &
            center_feature_counts$Primary_Drug_Class == cls]
        if (length(val) == 0 || is.na(val) || val <= 0) next
        fill <- class_palette[[cls]]; if (is.null(fill)) fill <- class_palette[["Other"]]
        grid.rect(x = bar_left[i], y = map_y(bottom), width = bar_w, height = plot_h * val / y_max,
                  just = c("left", "bottom"), gp = gpar(fill = fill, col = NA))
        bottom <- bottom + val
      }
      total <- center_feature_totals$count[center_feature_totals$center_code == cc]
      grid.text(as.character(total), x = bar_center[i], y = map_y(total) + 0.020,
                gp = gpar(fontfamily = font_main, fontsize = FONT_ANNOT, col = "black"))
      grid.text(paste0("C", i), x = bar_center[i], y = plot_y - 0.052,
                gp = gpar(fontfamily = font_main, fontsize = FONT_ANNOT, lineheight = 0.88, col = "black"))
    }

    legend_classes <- class_order[class_order %in% center_feature_counts$Primary_Drug_Class]
    legend_x <- 0.540
    legend_y <- 0.720
    legend_step <- 0.030
    grid.text("Antibiotic", x = legend_x, y = legend_y + 0.045, just = "left",
              gp = gpar(fontfamily = font_main, fontsize = FONT_LEGEND_TITLE, col = "black"))
    for (k in seq_along(legend_classes)) {
      cls <- legend_classes[k]
      lab <- class_labels[[cls]]; if (is.null(lab)) lab <- cls
      fill <- class_palette[[cls]]; if (is.null(fill)) fill <- class_palette[["Other"]]
      yy <- legend_y - (k - 1) * legend_step
      grid.rect(x = legend_x, y = yy - 0.010, width = 0.022, height = 0.022,
                just = c("left", "bottom"), gp = gpar(fill = fill, col = NA))
      grid.text(lab, x = legend_x + 0.030, y = yy + 0.001, just = "left",
                gp = gpar(fontfamily = font_main, fontsize = FONT_LEGEND_TEXT, col = "#222222"))
    }
  }

  draw_panel_A
}

# ============================================================================
# Panel B: prevalence and mean abundance of the high-frequency ARGs
# ============================================================================
make_panel_B <- function() {

  summary_file <- file.path(FIG3_IN, "arg_feature_summary.tsv")
  dpm_file     <- file.path(FIG3_IN, "arg_dpm_matrix.tsv")
  col_file     <- in_file("library_metadata.csv")

  arg_summary <- read.delim(summary_file, check.names = FALSE, stringsAsFactors = FALSE)
  dpm         <- read.delim(dpm_file,     check.names = FALSE, stringsAsFactors = FALSE)
  col_anno    <- read.csv(col_file,       check.names = FALSE, stringsAsFactors = FALSE)

  sample_cols <- intersect(colnames(dpm), col_anno$sample)
  if (length(sample_cols) != 34) stop("Expected 34 shared DPM/sample annotation columns, got ", length(sample_cols))

  dpm_values  <- as.matrix(dpm[, sample_cols, drop = FALSE]); mode(dpm_values) <- "numeric"
  rownames(dpm_values) <- dpm$feature_id
  summary_idx <- match(arg_summary$feature_id, rownames(dpm_values)); stopifnot(!anyNA(summary_idx))
  arg_summary$n_passed_samples <- rowSums(dpm_values[summary_idx, , drop = FALSE] > 0, na.rm = TRUE)
  arg_summary$prevalence       <- arg_summary$n_passed_samples / length(sample_cols)
  arg_summary$prevalence_pct   <- arg_summary$prevalence * 100
  arg_summary$mean_DPM_all     <- rowMeans(dpm_values[arg_summary$feature_id, , drop = FALSE], na.rm = TRUE)

  top_n <- 35
  candidate_features <- arg_summary[arg_summary$prevalence_pct >= 50 & arg_summary$mean_DPM_all >= 3, ]
  top_features <- candidate_features[order(-candidate_features$n_passed_samples, -candidate_features$mean_DPM_all), ][seq_len(top_n), ]

  class_order <- c(
    "fluoroquinolone", "aminoglycoside", "multiple_drug_classes",
    "multidrug_efflux_or_multiclass", "tetracycline", "peptide/polymyxin",
    "aminocoumarin", "nitroimidazole"
  )
  top_features$class_rank <- match(top_features$Primary_Drug_Class, class_order)
  top_features$class_rank[is.na(top_features$class_rank)] <- length(class_order) + 1
  top_features <- top_features[order(top_features$class_rank, -top_features$mean_DPM_all, -top_features$n_passed_samples), ]
  top_features$plot_order <- seq_len(nrow(top_features))

  selected_mat <- dpm_values[top_features$feature_id, sample_cols, drop = FALSE]

  centre_samples <- split(col_anno$sample, col_anno$center_code)
  centre_samples <- centre_samples[c("SHGH", "SY", "HG")]

  mean_for <- function(samples) {
    use_samples <- intersect(samples, colnames(selected_mat))
    if (length(use_samples) == 0) rep(NA_real_, nrow(selected_mat)) else rowMeans(selected_mat[, use_samples, drop = FALSE], na.rm = TRUE)
  }
  abundance_rows <- rbind(
    Overall = rowMeans(selected_mat, na.rm = TRUE),
    SHGH = mean_for(centre_samples$SHGH),
    SY   = mean_for(centre_samples$SY),
    HG   = mean_for(centre_samples$HG)
  )
  colnames(abundance_rows) <- top_features$gene_name

  feature_table <- data.frame(
    plot_order = top_features$plot_order, feature_id = top_features$feature_id,
    gene_name = top_features$gene_name, ARO_ID = top_features$ARO_ID,
    Primary_Drug_Class = top_features$Primary_Drug_Class,
    CARD_Drug_Class = top_features$CARD_Drug_Class,
    CARD_Resistance_Mechanism = top_features$CARD_Resistance_Mechanism,
    CARD_AMR_Gene_Family = top_features$CARD_AMR_Gene_Family,
    n_passed_samples = top_features$n_passed_samples, prevalence_pct = top_features$prevalence_pct,
    mean_DPM_all = top_features$mean_DPM_all, stringsAsFactors = FALSE
  )

  overall_abundance <- as.numeric(feature_table$mean_DPM_all)
  overall_abundance_log <- log10(overall_abundance + 1)
  individual_abundance_log <- log10(selected_mat + 1)
  individual_abundance_log[!is.finite(individual_abundance_log)] <- NA_real_
  abundance_axis_max <- ceiling(max(c(overall_abundance_log, individual_abundance_log), na.rm = TRUE) * 2) / 2
  overall_abundance_scaled <- 100 * overall_abundance_log / abundance_axis_max
  individual_abundance_scaled <- 100 * individual_abundance_log / abundance_axis_max

  class_palette <- c(
    "aminoglycoside" = "#EEC643", "tetracycline" = "#1F3B73",
    "multiple_drug_classes" = "#F4A460", "multidrug_efflux_or_multiclass" = "#00A6D6",
    "aminocoumarin" = "#9E8BD9", "fluoroquinolone" = "#F8766D",
    "nitroimidazole" = "#A33B62", "peptide/polymyxin" = "#FDB462", "Other" = "#BDBDBD"
  )

  heat_fun <- colorRampPalette(c("#FFFFFF", "#FFECEC", "#FF9A9A", "#FF0000"))(256)
  heat_log <- log10(abundance_rows + 1); heat_log[!is.finite(heat_log)] <- 0
  heat_min <- 0
  heat_max <- ceiling(max(heat_log, na.rm = TRUE) * 10) / 10
  heat_scaled <- (heat_log - heat_min) / (heat_max - heat_min)
  heat_scaled[is.na(heat_scaled)] <- 0; heat_scaled[heat_scaled < 0] <- 0; heat_scaled[heat_scaled > 1] <- 1
  heat_idx <- pmax(1, pmin(256, round(heat_scaled * 255) + 1))
  heat_cols <- matrix(heat_fun[heat_idx], nrow = nrow(heat_scaled), dimnames = dimnames(heat_scaled))

  class_labels <- c(
    "aminoglycoside" = "Aminogly.", "tetracycline" = "Tetracycl.",
    "multiple_drug_classes" = "Multiclass", "multidrug_efflux_or_multiclass" = "Efflux",
    "aminocoumarin" = "Aminocoum.", "fluoroquinolone" = "Fluoroquin.",
    "nitroimidazole" = "Nitroimid.", "peptide/polymyxin" = "Peptide/poly."
  )

  draw_panel_B <- function() {
    n <- nrow(feature_table)
    plot_center <- 0.205 + 0.610 / 2
    plot_w <- 0.610 * 0.80
    plot_x <- plot_center - plot_w / 2
    cell_w <- plot_w / n
    lolli_x <- plot_x; lolli_w <- plot_w

    overall_h <- 0.024
    center_h  <- overall_h * 0.67
    gap_class_overall <- overall_h * (1 / 6)
    gap_overall_center <- overall_h * (1 / 6)
    gap_heat_lolli <- (0.020 / 3) * 1.5
    class_h <- overall_h / 3
    class_y  <- 0.825
    center1_y <- class_y - gap_class_overall - center_h
    center2_y <- center1_y - center_h
    center3_y <- center2_y - center_h
    heat_y <- center3_y

    lolli_h <- 0.232
    lolli_y <- heat_y - gap_heat_lolli - lolli_h
    gene_y  <- lolli_y - 0.011

    label_size        <- 7.3
    center_label_size <- 6.4
    heat_axis_size    <- 6.4
    axis_size         <- 6.4
    gene_size         <- 6.8

    for (j in seq_len(n)) {
      cls <- feature_table$Primary_Drug_Class[j]
      fill <- class_palette[[cls]]; if (is.null(fill)) fill <- class_palette[["Other"]]
      grid.rect(x = plot_x + (j - 1) * cell_w, y = class_y, width = cell_w, height = class_h,
                just = c("left", "bottom"), gp = gpar(fill = fill, col = NA))
    }
    grid.text("Antibiotic", x = plot_x - 0.009, y = class_y + class_h / 2,
              just = "right", gp = gpar(fontfamily = "Arial", fontsize = 6.4))

    r <- rle(feature_table$Primary_Drug_Class)

    heat_row_names <- c("Center 1", "Center 2", "Center 3")
    heat_row_y     <- c(center1_y, center2_y, center3_y)
    heat_row_heights <- c(center_h, center_h, center_h)
    for (i in seq_len(3)) {
      y <- heat_row_y[i]
      for (j in seq_len(n)) {
        grid.rect(x = plot_x + (j - 1) * cell_w, y = y, width = cell_w, height = heat_row_heights[i],
                  just = c("left", "bottom"), gp = gpar(fill = heat_cols[i + 1, j], col = NA, lwd = 0))
      }
      grid.text(heat_row_names[i], x = plot_x - 0.005, y = y + heat_row_heights[i] / 2,
                just = "right", gp = gpar(fontfamily = "Arial", fontsize = center_label_size))
    }

    scale_x <- plot_x + 24.6 * cell_w
    scale_y <- 0.415; scale_w <- 0.067; scale_h <- 0.012; n_scale <- 64
    for (k in seq_len(n_scale)) {
      grid.rect(x = scale_x + (k - 1) * scale_w / n_scale, y = scale_y,
                width = scale_w / n_scale, height = scale_h, just = c("left", "bottom"),
                gp = gpar(fill = heat_fun[round(1 + (k - 1) * 255 / (n_scale - 1))], col = NA))
    }
    grid.rect(x = scale_x, y = scale_y, width = scale_w, height = scale_h,
              just = c("left", "bottom"), gp = gpar(fill = NA, col = "#4D4D4D", lwd = 0.70))
    grid.text("Average abundance (log10[DPM + 1])", x = scale_x, y = scale_y + 0.027,
              just = "left", gp = gpar(fontfamily = "Arial", fontsize = 6.6, col = "#333333"))
    grid.lines(x = c(scale_x, scale_x), y = c(scale_y, scale_y - 0.004), gp = gpar(col = "#4D4D4D", lwd = 0.70))
    grid.lines(x = c(scale_x + scale_w, scale_x + scale_w), y = c(scale_y, scale_y - 0.004), gp = gpar(col = "#4D4D4D", lwd = 0.70))
    grid.text("0", x = scale_x, y = scale_y - 0.011, just = "center", gp = gpar(fontfamily = "Arial", fontsize = 6.4, col = "#555555"))
    grid.text(sprintf("%.1f", heat_max), x = scale_x + scale_w, y = scale_y - 0.011, just = "center", gp = gpar(fontfamily = "Arial", fontsize = 6.4, col = "#555555"))

    y_ticks <- seq(0, 100, by = 25)
    grid.rect(x = lolli_x, y = lolli_y, width = lolli_w, height = lolli_h,
              just = c("left", "bottom"), gp = gpar(fill = NA, col = "#4D4D4D", lwd = 1.05))
    inner_pad <- lolli_h * 0.070
    map_y <- function(value) lolli_y + inner_pad + (lolli_h - 2 * inner_pad) * value / 100
    tick_len_y <- 0.005
    abundance_tick_labels <- formatC(abundance_axis_max * y_ticks / 100,
                                     format = "f", digits = 1)
    for (ti in seq_along(y_ticks)) {
      yt <- y_ticks[ti]
      yy <- map_y(yt)
      grid.lines(x = c(lolli_x, lolli_x + tick_len_y), y = c(yy, yy), gp = gpar(col = "#4D4D4D", lwd = 0.85))
      grid.lines(x = c(lolli_x + lolli_w - tick_len_y, lolli_x + lolli_w), y = c(yy, yy), gp = gpar(col = "#4D4D4D", lwd = 0.85))
      grid.text(as.character(yt), x = lolli_x - 0.006, y = yy, just = "right", gp = gpar(fontfamily = "Arial", fontsize = axis_size, col = "#333333"))
      grid.text(abundance_tick_labels[ti], x = lolli_x + lolli_w + 0.006, y = yy, just = "left",
                gp = gpar(fontfamily = "Arial", fontsize = axis_size, col = "black"))
    }
    grid.text("Prevalence (%)", x = lolli_x - 0.032, y = lolli_y + lolli_h / 2,
              rot = 90, gp = gpar(fontfamily = "Arial", fontsize = label_size, col = "black"))

    dot_col <- "#3FA7E8"; stem_col <- "#BDBDBD"
    abundance_dot_col <- "#D95F02"; abundance_stem_col <- "#BDBDBD"
    xx_vec <- plot_x + (seq_len(n) - 0.5) * cell_w
    yy_ab_vec  <- map_y(overall_abundance_scaled)
    yy_prev_vec <- map_y(feature_table$prevalence_pct)
    yy_individual <- apply(individual_abundance_scaled, 2, map_y)
    sample_offsets <- seq(-0.30, 0.30, length.out = ncol(selected_mat)) * cell_w
    tick_len_x <- 0.005
    for (j in seq_len(n)) {
      grid.lines(x = c(xx_vec[j], xx_vec[j]), y = c(lolli_y, lolli_y + tick_len_x), gp = gpar(col = "#4D4D4D", lwd = 0.70))
    }
    individual_col <- grDevices::adjustcolor("#666666", alpha.f = 0.28)
    for (j in seq_len(n)) {
      for (i in seq_len(ncol(selected_mat))) {
        if (!is.finite(yy_individual[j, i])) next
        grid.circle(x = xx_vec[j] + sample_offsets[i], y = yy_individual[j, i],
                    r = unit(0.48, "pt"), gp = gpar(fill = individual_col, col = NA))
      }
    }
    for (j in seq_len(n)) {
      grid.lines(x = c(xx_vec[j], xx_vec[j]), y = c(lolli_y + inner_pad, yy_ab_vec[j]),  gp = gpar(col = abundance_stem_col, lwd = 0.90))
      grid.lines(x = c(xx_vec[j], xx_vec[j]), y = c(lolli_y + inner_pad, yy_prev_vec[j]), gp = gpar(col = stem_col, lwd = 0.78))
    }
    for (j in seq_len(n)) {
      grid.circle(x = xx_vec[j], y = yy_ab_vec[j],  r = unit(0.0082, "npc"), gp = gpar(fill = abundance_dot_col, col = NA))
      grid.circle(x = xx_vec[j], y = yy_prev_vec[j], r = unit(0.0082, "npc"), gp = gpar(fill = dot_col, col = NA))
    }

    grid.text("Abundance (log10[DPM + 1])", x = plot_x + plot_w + 0.034,
              y = lolli_y + lolli_h / 2, rot = 90,
              gp = gpar(fontfamily = "Arial", fontsize = label_size, col = "black"))

    for (j in seq_len(n)) {
      xx <- plot_x + (j - 0.5) * cell_w
      grid.text(feature_table$gene_name[j], x = xx, y = gene_y, rot = 90, just = "right",
                gp = gpar(fontfamily = "Arial", fontsize = gene_size, col = "black"))
    }
    grid.text("ARGs (prevalence ≥ 50%; mean DPM ≥ 3)", x = plot_x + plot_w / 2, y = 0.385,
              gp = gpar(fontfamily = "Arial", fontsize = 7.2, col = "black"))
  }

  draw_panel_B
}

# ============================================================================
# Panel C: species-ARG Spearman correlation
# ============================================================================
make_panel_E <- function() {

  feature_file       <- file.path(FIG3_IN, "arg_species_spearman.tsv")
  feature_order_file <- file.path(FIG3_IN, "arg_top35.csv")
  row_annot_file     <- in_file("species_annotation.csv")

  feature_pairs <- read_tsv(feature_file, show_col_types = FALSE)
  feature_order <- read_csv(feature_order_file, show_col_types = FALSE)
  row_annot     <- read_csv(row_annot_file, show_col_types = FALSE)

  species_table <- row_annot %>%
    filter(evidence_tier != "Background detection") %>%
    mutate(evidence_rank = match(evidence_tier,
      c("Conjunctivitis-associated", "Indirect ocular/infection support"))) %>%
    arrange(evidence_rank, desc(max_ra_pct), species)
  species_order <- species_table$species

  class_order <- c(
    "fluoroquinolone", "aminoglycoside", "MLS", "beta-lactam", "tetracycline",
    "multidrug_efflux_or_multiclass", "multiple_drug_classes", "peptide/polymyxin",
    "aminocoumarin", "nitroimidazole"
  )
  class_palette <- c(
    "fluoroquinolone" = "#F8766D", "aminoglycoside" = "#EEC643", "beta-lactam" = "#B79F65",
    "MLS" = "#CC79A7", "tetracycline" = "#1F3B73", "multiple_drug_classes" = "#F4A460",
    "multidrug_efflux_or_multiclass" = "#00A6D6", "peptide/polymyxin" = "#FDB462",
    "aminocoumarin" = "#9E8BD9", "nitroimidazole" = "#A33B62", "Other" = "#BDBDBD"
  )
  class_labels <- c(
    "fluoroquinolone" = "Fluoroquin.", "aminoglycoside" = "Aminogly.", "beta-lactam" = "Beta-lactam",
    "MLS" = "MLS", "tetracycline" = "Tetracycl.", "multiple_drug_classes" = "Multiclass",
    "multidrug_efflux_or_multiclass" = "Efflux/multiclass", "peptide/polymyxin" = "Peptide/poly.",
    "aminocoumarin" = "Aminocoum.", "nitroimidazole" = "Nitroimid."
  )

  feature_levels <- feature_order %>%
    filter(gene_name %in% unique(feature_pairs$ARG_gene)) %>%
    arrange(plot_order) %>% pull(gene_name)

  feature_annot_df <- feature_order %>%
    filter(gene_name %in% feature_levels) %>%
    arrange(match(gene_name, feature_levels)) %>%
    transmute(ARG_gene = factor(gene_name, levels = feature_levels),
      Primary_Drug_Class = ifelse(Primary_Drug_Class %in% names(class_palette), Primary_Drug_Class, "Other"),
      annotation = "Antibiotic", x_pos = match(gene_name, feature_levels),
      y_pos = length(species_order) + 1,
      class_fill = unname(class_palette[Primary_Drug_Class]))

  contig_file <- file.path(FIG3_IN, "arg_contig_colocated_pairs.csv")
  contig_pairs <- if (file.exists(contig_file)) read_csv(contig_file, show_col_types = FALSE)
                  else data.frame(species = character(), gene_name = character())

  feature_plot_df <- feature_pairs %>%
    filter(species %in% species_order, ARG_gene %in% feature_levels) %>%
    mutate(species = factor(species, levels = rev(species_order)),
      ARG_gene = factor(ARG_gene, levels = feature_levels),
      x_pos = match(ARG_gene, feature_levels),
      y_pos = length(species_order) - match(species, species_order) + 1,
      Primary_Drug_Class = ifelse(Primary_Drug_Class %in% names(class_palette), Primary_Drug_Class, "Other"),
      absrho = abs(spearman_rho),
      significant = !is.na(p_value) & p_value < 0.05,
      stars = case_when(p_value < 0.001 ~ "***", p_value < 0.01 ~ "**", p_value < 0.05 ~ "*", TRUE ~ ""),
      alpha_value = ifelse(significant, 0.95,
                           0.65))

  contig_ring_df <- feature_plot_df %>%
    semi_join(contig_pairs %>% transmute(species = species, ARG_gene = gene_name),
              by = c("species" = "species", "ARG_gene" = "ARG_gene"))
  message("Panel C: cells with contig-level support: ", nrow(contig_ring_df))

  theme_arg_corr <- function(base_size = 8.5) {
    theme_minimal(base_size = base_size, base_family = "Arial") +
      theme(text = element_text(family = "Arial", colour = "black"),
        panel.grid = element_blank(),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, colour = "black"),
        axis.title = element_text(colour = "black"),
        plot.title = element_text(size = base_size + 1.2, hjust = 0.5, colour = "black"),
        legend.title = element_text(size = base_size, colour = "black", margin = margin(r = 6)),
        legend.text = element_text(size = base_size - 0.6, colour = "black"),
        legend.key.height = unit(0.12, "in"), legend.key.width = unit(0.18, "in"),
        plot.margin = margin(4, 6, 4, 6))
  }

  .fl <- 0.5
  .fb <- c(-.fl, -.fl / 2, 0, .fl / 2, .fl)
  .flab <- if (abs(.fl - 1) < 1e-9) {
    c("-1.0", "-0.5", "0", "0.5", "1.0")
  } else {
    l <- vapply(.fb, function(x) sub("^-0$", "0", sprintf("%g", x)), "")
    l[1] <- paste0("≤", l[1]); l[5] <- paste0("≥", l[5]); l
  }
  common_fill <- scale_fill_gradient2(
    low = "#0072B2", mid = "#F7F7F7", high = "#D55E00", midpoint = 0,
    limits = c(-.fl, .fl), oob = scales::squish, breaks = .fb,
    labels = .flab, name = "Spearman correlation",
    guide = guide_colourbar(title.position = "left", title.hjust = 0.5,
      display = "gradient",
      barwidth = unit(0.90, "in"), barheight = unit(0.12, "in"), ticks = TRUE, order = 2))
  common_size <- scale_size_continuous(
    range = c(0.80, 3.60), limits = c(0, 1), breaks = c(0.2, 0.4, 0.6, 0.8), name = "|rho|",
    guide = guide_legend(title.position = "left",
      override.aes = list(fill = "#BDBDBD", colour = "white", alpha = 0.9), order = 3))

  p_feature <- ggplot() +
    geom_tile(data = feature_plot_df, aes(x = x_pos, y = y_pos), fill = "white", colour = "#EBEBEB", linewidth = 0.239, width = 1, height = 1) +
    geom_tile(data = feature_annot_df, aes(x = x_pos, y = y_pos), fill = feature_annot_df$class_fill, colour = NA, width = 1, height = 0.30) +
    annotate("text", x = 0.46, y = length(species_order) + 1, label = "Antibiotic", hjust = 1, vjust = 0.5, family = "Arial", size = 2.45, colour = "black") +
    geom_point(data = feature_plot_df, aes(x = x_pos, y = y_pos, size = absrho, fill = spearman_rho, alpha = alpha_value), shape = 21, colour = "white", stroke = 0.07) +
    geom_point(data = contig_ring_df, aes(x = x_pos, y = y_pos, size = absrho),
               shape = 21, fill = NA, colour = "#009E73", stroke = 0.70) +
    geom_text(data = feature_plot_df %>% filter(significant), aes(x = x_pos, y = y_pos, label = stars), size = 1.95, family = "Arial", colour = "black") +
    common_fill + common_size + scale_alpha_identity() +
    scale_x_continuous(breaks = seq_along(feature_levels), labels = feature_levels, limits = c(0.45, length(feature_levels) + 0.55), expand = c(0, 0)) +
    scale_y_continuous(breaks = seq_along(species_order), labels = rev(species_order), limits = c(0.5, length(species_order) + 0.95), expand = c(0, 0)) +
    labs(x = NULL, y = "Candidate ocular pathogen species") +
    theme_arg_corr(base_size = 7.2) + coord_fixed(ratio = 1, clip = "off") +
    theme(axis.text.y = element_text(face = "italic", size = 7.0),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7.0, margin = margin(t = 1)),
      axis.title.y = element_text(size = 7.0, margin = margin(r = 2)),
      legend.position = "bottom", legend.box = "horizontal",
      legend.spacing.x = unit(0.08, "in"), legend.box.spacing = unit(0.01, "in"),
      plot.margin = margin(0, 4, 0, 4))

  metrics_file <- tempfile(fileext = ".pdf")
  open_pdf(metrics_file, 7.01, 3.10)
  legend_table <- ggplotGrob(p_feature)
  dev.off(); unlink(metrics_file)
  legend_index <- which(legend_table$layout$name == "guide-box-bottom")
  if (length(legend_index) != 1 || inherits(legend_table$grobs[[legend_index]], "zeroGrob")) {
    stop("Unable to extract the horizontal Figure 3C legend")
  }
  legend_grob <- legend_table$grobs[[legend_index]]

  p_feature_core <- p_feature + theme(legend.position = "none")
  footer_title <- wrap_elements(full = grid::textGrob(
    "ARGs (prevalence ≥ 50%; mean DPM ≥ 3)", x = unit(1, "npc"), just = "right",
    gp = grid::gpar(fontfamily = "Arial", fontsize = 7, col = "black")))
  footer_legend <- wrap_elements(full = legend_grob)
  footer <- footer_title + footer_legend + plot_layout(widths = c(0.21, 0.79))
  p_feature_final <- p_feature_core / footer + plot_layout(heights = c(1, 0.075))

  p_feature_final
}

# ============================================================================
# Composite layout (174 mm wide)
# ============================================================================
draw_A <- make_panel_A()
draw_B <- make_panel_B()
pE     <- make_panel_E()

W_mm <- 174

A_ratio <- 3.17 / 4.35
B_ratio <- 7.37 / 4.45
E_ratio <- 7.01 / 3.10

A_w <- 55
gap_h <- 1
B_w <- W_mm - A_w - gap_h
E_w <- W_mm

A_h <- A_w / A_ratio
B_h <- B_w / B_ratio
E_h <- E_w / E_ratio

H_top <- max(A_h, B_h)
gap_v <- 4
H_tot <- H_top + gap_v + E_h

crop_A <- c(0.02, 0.07, 0.95, 0.80)
crop_B <- c(0.19, 0.30, 0.86, 0.85)
crop_E <- c(0, 0, 1, 1)

push_crop <- function(box) {
  xmin <- box[1]; ymin <- box[2]; xmax <- box[3]; ymax <- box[4]
  cw <- xmax - xmin
  ch <- ymax - ymin
  pushViewport(viewport(x = unit(-xmin / cw, "npc"), y = unit(-ymin / ch, "npc"),
                        width = unit(1 / cw, "npc"), height = unit(1 / ch, "npc"),
                        just = c("left", "bottom")))
}

draw_composite <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))

  pushViewport(viewport(x = unit(0, "mm"), y = unit(0, "mm"),
                        width = unit(E_w, "mm"), height = unit(E_h, "mm"),
                        just = c("left", "bottom"), clip = TRUE))
  grid.draw(clamp_lwd_grob(patchwork::patchworkGrob(pE)))
  upViewport()

  pushViewport(viewport(x = unit(0, "mm"),
                        y = unit(E_h + gap_v + (H_top - A_h) / 2, "mm"),
                        width = unit(A_w, "mm"), height = unit(A_h, "mm"),
                        just = c("left", "bottom"), clip = TRUE))
  push_crop(crop_A)
  draw_A()
  upViewport(2)

  pushViewport(viewport(x = unit(A_w + gap_h, "mm"),
                        y = unit(E_h + gap_v + (H_top - B_h) / 2 - 5, "mm"),
                        width = unit(B_w, "mm"), height = unit(B_h, "mm"),
                        just = c("left", "bottom"), clip = TRUE))
  push_crop(crop_B)
  draw_B()
  upViewport(2)

  tag_gp <- gpar(fontfamily = "Arial", fontsize = FONT_TAG, fontface = "bold", col = "black")
  tag_y_top <- E_h + gap_v + H_top - 3
  grid.text("A", x = unit(TAG_X_LEFT, "mm"), y = unit(tag_y_top, "mm"), just = c("left", "top"), gp = tag_gp)
  grid.text("B", x = unit(A_w + gap_h, "mm"), y = unit(tag_y_top, "mm"), just = c("left", "top"), gp = tag_gp)
  grid.text("C", x = unit(TAG_X_LEFT, "mm"), y = unit(E_h - 3, "mm"), just = c("left", "top"), gp = tag_gp)
}

W_in <- W_mm / 25.4
H_in <- H_tot / 25.4

pdf_file <- file.path(OUT, "Figure3.pdf")
open_pdf(pdf_file, W_in, H_in)
draw_composite()
dev.off()
message("Figure size: ", round(W_mm, 1), " x ", round(H_tot, 1), " mm")
message("PDF: ", pdf_file)
