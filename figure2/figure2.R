# Figure 2: microbial community and diversity.
#   A  whole-community heatmap of the displayed taxa across all libraries (relative abundance, log10 scale)
#   B  Shannon index (top) and Bray-Curtis PCoA (bottom) of the Kraken2/Bracken species profiles, per center;
#      the cross-center pooled library is not attributable to a center and is not shown in (B)
# Run from the repository root: Rscript figure2/figure2.R
# Output: output/figure2/Figure2.pdf (174 mm wide).

suppressPackageStartupMessages({
  library(grid)
  library(pheatmap)
  library(gtable)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(dplyr)
  library(vegan)
})

source("config.R")
OUT <- out_dir("figure2")

# Arial is embedded through the quartz device on macOS; other platforms fall back to cairo_pdf.
if (capabilities("aqua")) quartzFonts(Arial = quartzFont(c("Arial", "Arial Bold", "Arial Italic", "Arial Bold Italic")))
open_pdf <- function(file, width, height) {
  if (capabilities("aqua")) {
    quartz(file = file, type = "pdf", width = width, height = height, family = "Arial", bg = "white")
  } else {
    cairo_pdf(file, width = width, height = height, family = "Arial", bg = "white")
  }
}

CENTER_COLS  <- c(SHGH = "#0396D3", SY = "#8CD3CA", HG = "#DD914D", MIX = "#9F9F9F", Other = "#BDBDBD")
LIBRARY_COLS <- c(`Individual/same-patient combined` = "#00A9E0", `Cross-patient pooled` = "#B39DDB")
KINGDOM_COLS <- c(Viruses = "#E7298A", Bacteria = "#E6AB02", Fungi = "#8064A2", `Other eukaryota` = "#666666")
REL_ABUNDANCE_PAL <- colorRampPalette(c("#FFFFFF", "#FFECEC", "#FF9A9A", "#FF0000"))(256)

nudge_gt_grobs <- function(gt, layout_names, dx_pt = 0, dy_pt = 0) {
  ii <- which(gt$layout$name %in% layout_names)
  for (i in ii) {
    gt$grobs[[i]] <- grobTree(
      gt$grobs[[i]],
      vp = viewport(x = unit(0.5, "npc") + unit(dx_pt, "pt"),
                    y = unit(0.5, "npc") + unit(dy_pt, "pt"),
                    width = unit(1, "npc"), height = unit(1, "npc")))
  }
  gt
}

# ============================================================================
# Panel A: whole-community heatmap
# ============================================================================
make_panel_A <- function() {
  mat  <- as.matrix(read.csv(in_file("figure2_display_matrix.csv"), row.names = 1, check.names = FALSE)) # taxa x library (fraction)
  rowa <- read.csv(in_file("figure2_display_rowanno.csv"), row.names = 1, check.names = FALSE)             # kingdom per taxon
  ca   <- read.csv(in_file("library_metadata.csv"), check.names = FALSE)                                  # library annotation

  mat  <- mat[, ca$sample]
  rowa <- rowa[rownames(mat), , drop = FALSE]
  taxa     <- rownames(mat)
  samples  <- colnames(mat)
  kingdom  <- as.character(rowa$kingdom)         # per taxon
  center   <- as.character(ca$center_code)       # per sample
  libtype  <- as.character(ca$library_type)      # per sample
  n_row <- length(taxa); n_col <- length(samples)

  row_bounds <- head(cumsum(rle(kingdom)$lengths), -1)   # gaps between kingdoms
  col_bounds <- head(cumsum(rle(center)$lengths),  -1)   # gaps between centers

  FLOOR <- -3; CAP <- 2
  brk <- seq(FLOOR, CAP, length.out = 257)
  color_of <- function(v) {
    out <- rep("#FFFFFF", length(v))
    ok <- is.finite(v) & v > 0
    lv <- log10(v[ok] * 100)
    lv[lv < FLOOR] <- FLOOR; lv[lv > CAP] <- CAP
    idx <- findInterval(lv, brk, all.inside = TRUE)   # 1..256
    out[ok] <- REL_ABUNDANCE_PAL[idx]
    out
  }
  colmat <- matrix(color_of(as.numeric(mat)), nrow = n_row,
                   dimnames = dimnames(mat))

  draw_A <- function() {
    W <- as.numeric(convertWidth(unit(1, "npc"), "in"))
    H <- as.numeric(convertHeight(unit(1, "npc"), "in"))

    lab_w     <- 1.49
    gap_lab   <- 0.030
    king_w    <- 0.095
    gap_king  <- 0.030
    king_lab_w<- 0.150
    trackH    <- 0.090
    track_gap <- 0.014
    gap_track <- 0.040
    bp        <- 0.030
    row_gap   <- 0.022
    col_gap   <- 0.022

    heat_x0  <- king_lab_w + king_w + gap_king
    heat_w   <- W - heat_x0 - gap_lab - lab_w
    top_block<- 2 * trackH + track_gap + gap_track
    heat_top <- H - top_block
    heat_h   <- heat_top - bp
    rowH <- (heat_h - length(row_bounds) * row_gap) / n_row
    colW <- (heat_w - length(col_bounds) * col_gap) / n_col

    y_bottom <- numeric(n_row); cur <- heat_top
    for (i in seq_len(n_row)) {
      y_bottom[i] <- cur - rowH; cur <- cur - rowH
      if (i %in% row_bounds) cur <- cur - row_gap
    }
    x_left <- numeric(n_col); cur <- heat_x0
    for (j in seq_len(n_col)) {
      x_left[j] <- cur; cur <- cur + colW
      if (j %in% col_bounds) cur <- cur + col_gap
    }
    y_cen <- y_bottom + rowH / 2
    x_cen <- x_left + colW / 2

    g <- expand.grid(i = seq_len(n_row), j = seq_len(n_col))
    grid.rect(x = unit(x_cen[g$j] / W, "npc"), y = unit(y_cen[g$i] / H, "npc"),
              width = unit(colW / W, "npc"), height = unit(rowH / H, "npc"),
              gp = gpar(fill = colmat[cbind(g$i, g$j)], col = "#EEEEEE", lwd = 0.3))

    grid.text(parse(text = paste0("italic('", gsub("'", "", taxa), "')")),
              x = unit((heat_x0 + heat_w + gap_lab) / W, "npc"), y = unit(y_cen / H, "npc"),
              just = c("left", "center"),
              gp = gpar(fontfamily = "Arial", fontsize = 6))

    kx <- king_lab_w + king_w / 2
    grid.rect(x = unit(kx / W, "npc"), y = unit(y_cen / H, "npc"),
              width = unit(king_w / W, "npc"), height = unit(rowH / H, "npc"),
              gp = gpar(fill = KINGDOM_COLS[kingdom], col = NA))
    grid.text("Kingdom",
              x = unit((king_lab_w / 2) / W, "npc"),
              y = unit(mean(range(y_cen)) / H, "npc"), rot = 90, just = c("center", "center"),
              gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 6))

    lib_cen <- heat_top + gap_track + trackH / 2
    cen_cen <- heat_top + gap_track + trackH + track_gap + trackH / 2
    grid.rect(x = unit(x_cen / W, "npc"), y = unit(cen_cen / H, "npc"),
              width = unit(colW / W, "npc"), height = unit(trackH / H, "npc"),
              gp = gpar(fill = CENTER_COLS[center], col = NA))
    grid.rect(x = unit(x_cen / W, "npc"), y = unit(lib_cen / H, "npc"),
              width = unit(colW / W, "npc"), height = unit(trackH / H, "npc"),
              gp = gpar(fill = LIBRARY_COLS[libtype], col = NA))
    grid.text("Center", x = unit((heat_x0 + heat_w + gap_lab) / W, "npc"), y = unit(cen_cen / H, "npc"),
              just = c("left", "center"), gp = gpar(fontfamily = "Arial", fontsize = 6))
    grid.text("Library type", x = unit((heat_x0 + heat_w + gap_lab) / W, "npc"), y = unit(lib_cen / H, "npc"),
              just = c("left", "center"), gp = gpar(fontfamily = "Arial", fontsize = 6))
  }
  draw_A
}

# ============================================================================
# Panel B: alpha and beta diversity by center
# ============================================================================
make_panel_B <- function() {
  matrix_file     <- in_file("bracken_species_relative_abundance.csv")
  annotation_file <- in_file("library_metadata.csv")

  species_df <- read_csv(matrix_file, show_col_types = FALSE)
  col_anno <- read_csv(annotation_file, show_col_types = FALSE) %>%
    mutate(sample = as.character(sample), center_code = as.character(center_code), library_type = as.character(library_type))

  meta_cols <- c("taxonomy_id", "name", "taxonomy_lvl", "genus_name", "genus_id", "family_name", "family_id", "order_name", "order_id")
  sample_cols <- intersect(col_anno$sample, setdiff(names(species_df), meta_cols))
  stopifnot(length(sample_cols) == nrow(col_anno))

  mat <- as.matrix(species_df[, sample_cols, drop = FALSE]); mode(mat) <- "numeric"
  rownames(mat) <- species_df$name; mat[is.na(mat)] <- 0
  mat <- t(mat); mat <- mat[rowSums(mat) > 0, , drop = FALSE]; mat <- mat / rowSums(mat)
  col_anno <- col_anno[match(rownames(mat), col_anno$sample), , drop = FALSE]
  keep <- col_anno$center_code %in% c("SHGH", "SY", "HG")
  mat_plot <- mat[keep, , drop = FALSE]; anno_plot <- col_anno[keep, , drop = FALSE]

  center_order <- c("SHGH", "SY", "HG")
  center_labels <- c("SHGH" = "Center 1", "SY" = "Center 2", "HG" = "Center 3")
  center_cols <- CENTER_COLS[c("SHGH", "SY", "HG")]
  anno_plot <- anno_plot %>% mutate(
    center_code = factor(center_code, levels = center_order),
    center_label = factor(center_labels[as.character(center_code)], levels = center_labels[center_order]))

  alpha_df <- data.frame(sample = rownames(mat_plot),
                         shannon_index = vegan::diversity(mat_plot, index = "shannon"),
                         observed_species = rowSums(mat_plot > 0), stringsAsFactors = FALSE) %>%
    left_join(anno_plot, by = "sample") %>% arrange(center_code, sample)

  bray <- vegan::vegdist(mat_plot, method = "bray")
  pcoa <- cmdscale(bray, eig = TRUE, k = 2)
  positive_eig <- pcoa$eig[pcoa$eig > 0]
  var_explained <- pcoa$eig[1:2] / sum(positive_eig)
  pcoa_df <- data.frame(sample = rownames(pcoa$points), PCoA1 = pcoa$points[, 1], PCoA2 = pcoa$points[, 2],
                        stringsAsFactors = FALSE) %>% left_join(anno_plot, by = "sample") %>% arrange(center_code, sample)

  theme_figure2d <- function(base_size = 6.0) {
    theme_classic(base_size = base_size, base_family = "Arial") +
      theme(text = element_text(family = "Arial", colour = "black"),
            axis.text = element_text(size = base_size, colour = "black"),
            axis.title = element_text(size = base_size, colour = "black"),
            axis.line = element_blank(),
            axis.ticks = element_line(linewidth = 0.25, colour = "#111111"),
            axis.ticks.length = unit(2.8, "pt"),
            panel.border = element_rect(fill = NA, colour = "#111111", linewidth = 0.42),
            plot.background = element_rect(fill = NA, colour = NA),
            plot.margin = margin(5, 6, 5, 6),
            legend.title = element_blank(),
            legend.text = element_text(size = base_size, colour = "black"),
            legend.key.height = unit(0.085, "in"), legend.key.width = unit(0.125, "in"))
  }

  p_alpha <- ggplot(alpha_df, aes(x = center_label, y = shannon_index, colour = center_code, fill = center_code)) +
    geom_boxplot(width = 0.48, alpha = 0.55, outlier.shape = NA, linewidth = 0.30) +
    geom_point(position = position_jitter(width = 0.075, height = 0, seed = 20260626), size = 1.38, alpha = 0.92, stroke = 0) +
    scale_colour_manual(values = center_cols, guide = "none") + scale_fill_manual(values = center_cols, guide = "none") +
    scale_x_discrete(labels = c("Center 1" = "C1", "Center 2" = "C2", "Center 3" = "C3")) +
    labs(x = NULL, y = "Shannon index") + theme_figure2d() +
    theme(aspect.ratio = 1, axis.text.x = element_text(size = 6.0), plot.margin = margin(5, 5, 5, 5))

  p_pcoa <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, colour = center_code, fill = center_code)) +
    stat_ellipse(type = "norm", geom = "polygon", level = 0.95, alpha = 0.12, colour = NA, show.legend = FALSE) +
    stat_ellipse(type = "norm", level = 0.95, linewidth = 0.28, alpha = 0.95, show.legend = FALSE) +
    geom_point(size = 1.55, alpha = 0.94, stroke = 0) +
    scale_colour_manual(values = center_cols, breaks = center_order, labels = center_labels[center_order], guide = "none") +
    scale_fill_manual(values = center_cols, breaks = center_order, labels = center_labels[center_order], guide = "none") +
    labs(x = sprintf("PCoA1 (%.1f%%)", var_explained[1] * 100), y = sprintf("PCoA2 (%.1f%%)", var_explained[2] * 100)) +
    coord_cartesian(expand = TRUE) + theme_figure2d() +
    theme(aspect.ratio = 1, legend.background = element_rect(fill = "white", colour = NA),
          legend.box.margin = margin(0, 0, 0, 0), plot.margin = margin(3, 5, 5, 6))

  g_alpha <- ggplotGrob(p_alpha)
  g_pcoa  <- ggplotGrob(p_pcoa)
  common_widths <- do.call(grid::unit.pmax, list(g_alpha$widths, g_pcoa$widths))
  common_heights <- do.call(grid::unit.pmax, list(g_alpha$heights, g_pcoa$heights))
  g_alpha$widths <- common_widths; g_pcoa$widths <- common_widths
  g_alpha$heights <- common_heights; g_pcoa$heights <- common_heights
  b_axis_parts <- c("axis-l", "axis-b", "xlab-b", "ylab-l")
  g_alpha <- nudge_gt_grobs(g_alpha, b_axis_parts, dx_pt = -1.000, dy_pt = 0.485)
  g_pcoa  <- nudge_gt_grobs(g_pcoa,  b_axis_parts, dx_pt = -1.000, dy_pt = 0.485)

  draw_B <- function() {
    W <- as.numeric(convertWidth(unit(1, "npc"), "in"))
    H <- as.numeric(convertHeight(unit(1, "npc"), "in"))
    padx <- 0.03
    top_pad <- 0.02
    sq_gap <- 0.11
    legend_block_h <- 1.28
    avail_sq <- H - top_pad - legend_block_h - sq_gap
    side <- min(W - padx, avail_sq / 2)
    x0 <- max(padx, (W - side) / 2)

    alpha_bottom <- H - top_pad - side
    pushViewport(viewport(x = unit(x0 / W, "npc"), y = unit(alpha_bottom / H, "npc"),
                          width = unit(side / W, "npc"), height = unit(side / H, "npc"), just = c("left", "bottom")))
    grid.draw(g_alpha); upViewport()
    pcoa_top <- alpha_bottom - sq_gap; pcoa_bottom <- pcoa_top - side
    pushViewport(viewport(x = unit(x0 / W, "npc"), y = unit(pcoa_bottom / H, "npc"),
                          width = unit(side / W, "npc"), height = unit(side / H, "npc"), just = c("left", "bottom")))
    grid.draw(g_pcoa); upViewport()

    LEG_FS <- 6.0
    sw <- 5.5 / 72; pad <- 2 / 72; cg <- 2 / 72
    tw <- function(s, bold = FALSE) as.numeric(convertWidth(
      grobWidth(textGrob(s, gp = gpar(fontfamily = "Arial", fontface = if (bold) "bold" else "plain", fontsize = LEG_FS))), "in"))
    xL <- x0
    draw_row <- function(title, keys, cols, labels, y, x_start = xL) {
      grid.text(title, x = unit(x_start / W, "npc"), y = unit(y / H, "npc"), just = c("left", "center"),
                gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = LEG_FS))
      xx <- x_start + if (nzchar(title)) tw(title, TRUE) + cg else 0
      for (k in seq_along(keys)) {
        grid.rect(x = unit(xx / W, "npc"), y = unit(y / H, "npc"),
                  width = unit(sw / W, "npc"), height = unit(sw / H, "npc"), just = c("left", "center"),
                  gp = gpar(fill = cols[[keys[k]]], col = "#CCCCCC", lwd = 0.75))
        grid.text(labels[k], x = unit((xx + sw + pad) / W, "npc"), y = unit(y / H, "npc"), just = c("left", "center"),
                  gp = gpar(fontfamily = "Arial", fontsize = LEG_FS))
        xx <- xx + sw + pad + tw(labels[k]) + cg
      }
    }
    y0 <- pcoa_bottom - 0.14
    ls <- 0.165
    draw_row("Kingdom", c("Viruses", "Bacteria", "Fungi"), KINGDOM_COLS, c("Viruses", "Bacteria", "Fungi"), y0)
    draw_row("", "Other eukaryota", KINGDOM_COLS, "Other eukaryota", y0 - ls,
             x_start = xL + tw("Kingdom", TRUE) + cg)
    draw_row("Center", c("SHGH", "SY", "HG", "MIX"), CENTER_COLS,
             c("Center 1", "Center 2", "Center 3", "MIX"), y0 - 2 * ls)
    draw_row("Library type", c("Individual/same-patient combined", "Cross-patient pooled"), LIBRARY_COLS,
             c("Individual", "Pooled"), y0 - 3 * ls)

    cb_x <- xL; cb_w <- W - cb_x - 0.06; cb_h <- 0.075
    title_y <- y0 - 4 * ls - 0.02; bar_y <- title_y - 0.11
    grid.text("Relative abundance (%)", x = unit(cb_x / W, "npc"), y = unit(title_y / H, "npc"),
              just = c("left", "top"), gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = LEG_FS))
    n <- length(REL_ABUNDANCE_PAL)
    for (i in seq_len(n)) {
      grid.rect(x = unit((cb_x + (i - 1) / n * cb_w) / W, "npc"), y = unit(bar_y / H, "npc"),
                width = unit(cb_w / n / W, "npc"), height = unit(cb_h / H, "npc"),
                just = c("left", "top"), gp = gpar(fill = REL_ABUNDANCE_PAL[i], col = NA))
    }
    grid.rect(x = unit(cb_x / W, "npc"), y = unit(bar_y / H, "npc"),
              width = unit(cb_w / W, "npc"), height = unit(cb_h / H, "npc"),
              just = c("left", "top"), gp = gpar(fill = NA, col = "#999999", lwd = 0.75))
    ticks <- c("0.001", "0.01", "0.1", "1", "10", "100")
    for (t in seq_along(ticks)) {
      fx <- (t - 1) / (length(ticks) - 1)
      hj <- if (t == 1) "left" else if (t == length(ticks)) "right" else "center"
      grid.text(ticks[t], x = unit((cb_x + fx * cb_w) / W, "npc"), y = unit((bar_y - cb_h - 0.026) / H, "npc"),
                just = c(hj, "top"), gp = gpar(fontfamily = "Arial", fontsize = LEG_FS))
    }
  }
  draw_B
}

# ============================================================================
# Build both panels (a null device supplies text metrics)
# ============================================================================
opened_null_device <- FALSE
if (dev.cur() == 1L) { grDevices::pdf(file = NULL); opened_null_device <- TRUE }
draw_A <- make_panel_A()
draw_B <- make_panel_B()
if (opened_null_device) dev.off()

# ============================================================================
# Composite layout (174 mm wide)
# ============================================================================
W_mm <- 174; W_in <- W_mm / 25.4
margin_left <- 0.14; margin_right <- 0.07; margin_top <- 0.12; margin_bot <- 0.14
gap_AB <- 0.17
wA <- 4.05
xB <- margin_left + wA + gap_AB
wB <- W_in - margin_right - xB
H_in <- 5.30
hA <- H_in - margin_top - margin_bot
hB <- hA
xA <- margin_left; yA <- margin_bot; yB <- margin_bot

tag_gp <- gpar(fontfamily = "Arial", fontface = "bold", fontsize = 8)

draw_composite <- function() {
  grid.newpage(); grid.rect(gp = gpar(fill = "white", col = NA))
  pushViewport(viewport(x = unit(xA, "in"), y = unit(yA, "in"), width = unit(wA, "in"), height = unit(hA, "in"),
                        just = c("left", "bottom"), clip = "on"))
  draw_A(); upViewport()
  pushViewport(viewport(x = unit(xB, "in"), y = unit(yB, "in"), width = unit(wB, "in"), height = unit(hB, "in"),
                        just = c("left", "bottom"), clip = "on"))
  draw_B(); upViewport()
  grid.text("A", x = unit(margin_left - 0.02, "in"), y = unit(H_in - margin_top + 0.04, "in"),
            just = c("left", "top"), gp = tag_gp)
  grid.text("B", x = unit(xB - 0.02, "in"), y = unit(H_in - margin_top + 0.04, "in"),
            just = c("left", "top"), gp = tag_gp)
}

pdf_file <- file.path(OUT, "Figure2.pdf")
open_pdf(pdf_file, W_in, H_in); draw_composite(); dev.off()
message(sprintf("Canvas %.2f x %.2f mm", W_in * 25.4, H_in * 25.4))
message("PDF: ", pdf_file)
