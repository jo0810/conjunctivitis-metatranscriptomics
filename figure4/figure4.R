# Figure 4: human adenovirus (HAdV) phylogeny and HCoV-229E spike diversity.
#   A  maximum-likelihood whole-genome tree of the study HAdV genomes among HAdV-D references (fan layout);
#      study genomes colored by center and labeled with library and assigned type
#   B  maximum-likelihood tree of the study HCoV-229E genome with references (left) and spike amino-acid
#      states at polymorphic sites relative to the column consensus (right); study-specific substitutions in red
# Trees were inferred with MAFFT and IQ-TREE (see docs/upstream_processing.md).
# Run from the repository root: Rscript figure4/figure4.R
# Output: output/figure4/Figure4.pdf (174 mm wide).

suppressPackageStartupMessages({
  library(dplyr)
  library(grid)
  library(gtable)
  library(pheatmap)
  library(ape)
  library(ggplot2)
  library(ggtree)
  library(treeio)
  library(png)
  library(gridGraphics)
})

# Arial is embedded through the quartz device on macOS; other platforms fall back to cairo.
ON_MAC <- capabilities("aqua")
if (ON_MAC) quartzFonts(Arial = quartzFont(c("Arial", "Arial Bold", "Arial Italic", "Arial Bold Italic")))
PNG_TYPE <- if (ON_MAC) "quartz" else "cairo"
open_pdf8 <- function(file, width, height) {
  if (ON_MAC) quartz(type = "pdf", file = file, width = width, height = height, pointsize = 8)
  else cairo_pdf(file, width = width, height = height, pointsize = 8)
}
theme_set(theme_get() + ggplot2::theme(text = ggplot2::element_text(family = "Arial")))

source("config.R")
OUT <- out_dir("figure4")

# ============================================================================
# Layout: 174 mm wide, 6-pt text; panel A (fan) above panel B (HCoV-229E)
# ============================================================================
WIDTH_MM    <- 174
WIDTH_IN    <- WIDTH_MM / 25.4        # = 6.850 in
HEIGHT_MAX_IN <- 234 / 25.4

FONT_BC     <- 6
FONT_B      <- 6
LABEL_SIZE  <- 6
MARGIN_SIDE <- 0.14
MARGIN_BOT  <- 0.10
MARGIN_TOP  <- 0.16
GAP_ROW     <- 0.14
LABEL_PAD_Y <- 0.05

content_w   <- WIDTH_IN - 2 * MARGIN_SIDE

FAN_SIZE_IN <- 4.30
FAN_PAD_IN  <- 0.04
LIM_FAC     <- 1.55

hC <- 3.25
wC <- content_w

# ============================================================================
# Panel A: HAdV whole-genome tree (fan layout)
# ============================================================================
tree_file <- in_file("hadv_tree.treefile")
tip_file  <- in_file("hadv_tip_annotation.csv")
outgroup_tip <- "OUTGROUP__AC_000008_1"

reference_background_col <- "#DEDEDE"
study_col <- "#0396D3"
CENTER_COLS_FIG1 <- c("Center 1" = "#0396D3", "Center 2" = "#8CD3CA", "Center 3" = "#DD914D")
STUDY_OUTLINE   <- "#4D4D4D"
STUDY_LWD       <- 0.70
STUDY_CEX       <- 1.00
STUDY_LABEL_COL <- "#333333"

num <- function(x) suppressWarnings(as.numeric(x))
read_csv_base <- function(path) read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA", "NaN"))

canon_type <- function(t) {
  t <- as.character(t)
  ifelse(grepl("unresolved", t, ignore.case = TRUE), "HAdV-D (untyped)",
         sub("^HAdV-([0-9A-Za-z]+)$", "HAdV-D\\1", t))
}

base_font <- 1
draw_radial_label <- function(label, theta, r, col, cex = 1) {
  ang <- theta * 180 / pi
  left <- cos(theta) < 0
  srt <- ifelse(left, ang - 180, ang)
  adjx <- ifelse(left, 1, 0)
  text(r * cos(theta), r * sin(theta), label, srt = srt, adj = c(adjx, 0.5),
       col = col, cex = cex, font = base_font, xpd = NA)
}

tree <- read.tree(tree_file)
if (!is.null(tree$edge.length)) {
  tree$edge.length[!is.finite(tree$edge.length)] <- 0
  tree$edge.length[tree$edge.length < 0] <- 0
}
tip_tbl <- read_csv_base(tip_file)
if (!(outgroup_tip %in% tree$tip.label)) stop("Rooting outgroup missing: ", outgroup_tip)
tree <- root(tree, outgroup = outgroup_tip, resolve.root = TRUE)
tree <- ladderize(tree)
tree <- drop.tip(tree, outgroup_tip)
tip_tbl <- tip_tbl[tip_tbl$name != outgroup_tip, , drop = FALSE]

scale_factor <- 40
tree$edge.length <- tree$edge.length * scale_factor
min_terminal <- 0.006 * scale_factor
tip_edge <- tree$edge[, 2] <= Ntip(tree)
tree$edge.length[tip_edge & tree$edge.length < min_terminal] <- min_terminal

tip_idx <- match(tree$tip.label, tip_tbl$name)
if (any(is.na(tip_idx))) stop("Tree tips missing from tip table")
render <- tip_tbl[tip_idx, , drop = FALSE]
render$is_ref <- render$source_type == "reference"
render$hadv_ring <- render$hadv_major_type
render$hadv_ring[is.na(render$hadv_ring) | render$hadv_ring == ""] <- "HAdV-unknown"
# Library identifiers shown in the figure: column 1 = published library ID, column 2 = identifier used in
# the tree and tip table (Table S2 lists both).
LIBMAP_CSV <- in_file("library_id_map.csv")
.lm    <- read.csv(LIBMAP_CSV, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM")
LIBMAP <- setNames(.lm[[1]], .lm[[2]])
lib_id <- function(lab) {
  v <- unname(LIBMAP[as.character(lab)])
  if (any(is.na(v))) stop("Libraries missing from library_id_map.csv: ", paste(lab[is.na(v)], collapse = ", "))
  v
}
render$disp <- ""
.si <- which(!render$is_ref)
render$disp[.si] <- lib_id(render$label[.si])
render$center <- NA_character_
render$center[.si] <- sub("^(Center [0-9]+)-.*$", "\\1", render$disp[.si])
if (!all(render$center[.si] %in% names(CENTER_COLS_FIG1)))
  stop("Study genomes without a center color: ", paste(unique(render$center[.si]), collapse = ", "))
.mp <- unique(data.frame(old = render$label[.si], new = render$disp[.si]))
message("Panel A study tips:\n", paste(sprintf("   %-20s -> %s", .mp$old, .mp$new), collapse = "\n"))
if (anyDuplicated(.mp$new)) stop("Duplicated tip labels: ", paste(.mp$new[duplicated(.mp$new)], collapse = ", "))
render$study_label <- ifelse(render$is_ref, "",
                             paste0(render$disp, " | ", canon_type(render$hadv_ring)))

maxR <- max(node.depth.edgelength(tree)[seq_len(Ntip(tree))])
lim <- c(-1, 1) * (maxR * LIM_FAC)

cex_study <- FONT_B / 8

draw_hadv_fan <- function() {
  par(mar = c(0.2, 0.2, 0.2, 0.2), xpd = NA, family = "Arial")
  plot.phylo(tree, type = "fan", use.edge.length = TRUE, show.tip.label = FALSE,
             edge.color = "#6F6F6F", edge.width = 0.8, no.margin = TRUE,
             open.angle = 10, rotate.tree = 15, x.lim = lim, y.lim = lim)
  lp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  n_tip <- length(tree$tip.label)
  xx <- lp$xx[seq_len(n_tip)]; yy <- lp$yy[seq_len(n_tip)]
  theta <- atan2(yy, xx); r_tip <- sqrt(xx^2 + yy^2)
  study_idx <- which(!render$is_ref)
  ref_idx   <- which(render$is_ref)
  points(xx[ref_idx],   yy[ref_idx],   pch = 16, col = reference_background_col, cex = 0.80)
  points(xx[study_idx], yy[study_idx], pch = 21,
         bg = CENTER_COLS_FIG1[render$center[study_idx]],
         col = STUDY_OUTLINE, lwd = STUDY_LWD, cex = STUDY_CEX)
  label_offset <- 0.005 * scale_factor
  for (i in study_idx) {
    draw_radial_label(render$study_label[i], theta[i], r_tip[i] + label_offset,
                      STUDY_LABEL_COL, cex = cex_study)
  }
}

# ============================================================================
# Panel B: HCoV-229E tree and spike amino-acid matrix
# ============================================================================
cpt <- function(v) v / 2.8346
LW_MIN <- 0.236   # ~0.504 pt (minimum stroke 0.5 pt)
PRIV_H <- 0.95
C_tiplab   <- cpt(FONT_BC)
C_supp     <- cpt(6.0)
C_scale    <- cpt(FONT_BC)
C_domlab   <- cpt(6.0)
C_arrow    <- cpt(FONT_BC)
C_legtitle <- cpt(FONT_BC)
C_legtext  <- cpt(FONT_BC)

tr_io <- treeio::read.iqtree(in_file("hcov229e_tree.treefile"))
tr    <- as.phylo(tr_io)
tr    <- root(tr, outgroup = "NC_005831.2", resolve.root = TRUE)
tr    <- drop.tip(tr, "NC_005831.2")
supp_raw <- tr$node.label
uf <- suppressWarnings(as.numeric(sub(".*/", "", supp_raw)))

sel <- read.csv(in_file("hcov229e_reference_table.csv"), stringsAsFactors = FALSE)
STUDY <- sel$accession[sel$lineage_or_type == "Study sequence"]          # the study genome
stopifnot(length(STUDY) == 1)
sel$yr <- sapply(sel$year, function(x){ mm <- regmatches(x, regexpr("[0-9]{4}", x)); if (length(mm)) mm else "" })
sel$ctry <- trimws(sub(":.*", "", sel$country_or_region))
ann <- data.frame(name = sel$accession, ctry = sel$ctry, yr = sel$yr, stringsAsFactors = FALSE)
ann$is_study <- ann$name == STUDY
ann$attr <- ifelse(ann$ctry == "", ann$yr,
                   ifelse(ann$yr == "", ann$ctry, paste0(ann$ctry, "/", ann$yr)))
ann$lab <- ifelse(ann$is_study, paste0(lib_id(STUDY), "  HCoV-229E | China/2024"),
                  paste0(ann$name, "  HCoV-229E | ", ann$attr))
ann <- ann[ann$name %in% tr$tip.label, ]

m <- read.csv(in_file("hcov229e_spike_aa_matrix.csv"), row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
m <- m[tr$tip.label, , drop = FALSE]
cons <- apply(m, 2, function(col){ t <- table(col[col != "-"]); names(t)[which.max(t)] })
diffm <- m
for (j in seq_len(ncol(m))) diffm[[j]] <- ifelse(m[[j]] == cons[j], "·", m[[j]])
diffm <- as.data.frame(diffm); rownames(diffm) <- rownames(m)

dom_ann <- read.csv(in_file("hcov229e_spike_domains.csv"), stringsAsFactors = FALSE)
ps <- read.csv(in_file("hcov229e_private_sites.csv"), stringsAsFactors = FALSE)
ps_cols <- paste0("S", ps$ref_pos); ps_labs <- ps$label; ps_pos <- ps$ref_pos

STATUS_COL <- c(match = "#C9C9C9", gap = "#FFFFFF", sub = "#5B8AA6", private = "#D01C1F")
DOM_COL <- c(SP = "#3C3C3C", NTD = "#4E79A7", RBD = "#F28E2B", `S2-stem` = "#D9C9A3",
             FP = "#B07AA1", HR1 = "#253494", HR2 = "#A89BD4", TM = "#79706E", CT = "#A6A6A6")
DOM_RANGE <- c(SP = "1-15", NTD = "16-416", RBD = "417-547, CD13", `S2-stem` = "S2 core",
               FP = "753-773", HR1 = "774-886", HR2 = "1031-1115", TM = "1116-1135", CT = "1136-1173")
arch <- data.frame(
  dom  = c("SP","NTD","RBD","S2-stem","FP","HR1","S2-stem","HR2","TM","CT"),
  from = c(1, 16, 417, 548, 753, 774, 887, 1031, 1116, 1136),
  to   = c(15, 416, 547, 752, 773, 886, 1030, 1115, 1135, 1173))

p <- ggtree(tr, size = 0.45) %<+% ann +
  geom_tiplab(aes(label = lab, color = is_study), size = C_tiplab, align = TRUE,
              linesize = LW_MIN, offset = 0.0004, family = "Arial") +
  geom_tippoint(aes(color = is_study), size = 1.7) +
  scale_color_manual(values = c("FALSE" = "#333333", "TRUE" = "#D01C1F"), guide = "none")

d  <- p$data
ntip <- sum(d$isTip)
supp_df <- d %>% filter(!isTip) %>% arrange(node)
supp_df$uf <- uf[seq_len(nrow(supp_df))]
supp_df <- supp_df %>% filter(!is.na(uf) & uf >= 70) %>% arrange(y)

SUPP_H    <- 0.74
SUPP_MAXD <- 0.36
.tr_s     <- diff(range(d$x))
SUPP_W3   <- 0.145 * .tr_s
SUPP_W2   <- 0.097 * .tr_s
supp_df$w  <- ifelse(round(supp_df$uf) >= 100, SUPP_W3, SUPP_W2)
supp_df$yl <- supp_df$y
.ns  <- nrow(supp_df)
.xlo <- function(k) supp_df$x[k] - 1.15 * supp_df$w[k]
.xhi <- function(k) supp_df$x[k] - 0.15 * supp_df$w[k]
for (.it in seq_len(800)) {
  .moved <- FALSE
  for (i in seq_len(.ns - 1L)) for (j in (i + 1L):.ns) {
    if (.xlo(i) < .xhi(j) && .xlo(j) < .xhi(i)) {
      .dy <- supp_df$yl[j] - supp_df$yl[i]
      if (abs(.dy) < SUPP_H) {
        .need <- (SUPP_H - abs(.dy)) / 2 + 1e-9
        .sg   <- if (.dy >= 0) 1 else -1
        supp_df$yl[j] <- supp_df$yl[j] + .sg * .need
        supp_df$yl[i] <- supp_df$yl[i] - .sg * .need
        .moved <- TRUE
      }
    }
  }
  .off <- supp_df$yl - supp_df$y
  .ov  <- abs(.off) > SUPP_MAXD
  if (any(.ov)) supp_df$yl[.ov] <- supp_df$y[.ov] + sign(.off[.ov]) * SUPP_MAXD
  if (!.moved) break
}
.resid <- 0
for (i in seq_len(.ns - 1L)) for (j in (i + 1L):.ns)
  if (.xlo(i) < .xhi(j) && .xlo(j) < .xhi(i) && abs(supp_df$yl[j] - supp_df$yl[i]) < SUPP_H) .resid <- .resid + 1L
message(sprintf("Bootstrap labels: %d | max shift %.3f rows | remaining overlaps %d", .ns,
                max(abs(supp_df$yl - supp_df$y)), .resid))

p <- p + geom_text(data = supp_df, aes(x = x, y = yl, label = round(uf)),
                   inherit.aes = FALSE, size = C_supp, color = "#3A3A3A",
                   hjust = 1.15, vjust = -0.35, family = "Arial")
p <- p + geom_treescale(x = 0, y = -0.6, width = 0.01, fontsize = C_scale, linesize = 0.4, offset = 0.4, family = "Arial")

ncols     <- ncol(diffm)
treerange <- diff(range(d$x)); xmax_tip <- max(d$x[d$isTip])
message(sprintf("Panel C: treerange=%.5f  xmax_tip=%.5f  ncols=%d", treerange, xmax_tip, ncols))
OFFSET    <- 1.95 * treerange
HW_TOTAL  <- 2.30 * treerange
colw      <- HW_TOTAL / ncols
x0        <- xmax_tip + OFFSET
tipy      <- d$y[match(rownames(diffm), d$label)]
ymax      <- max(tipy, na.rm = TRUE); ymin <- min(tipy, na.rm = TRUE)
xarch <- function(pos) x0 + (pos - 1)/(1173 - 1) * HW_TOTAL

mm <- as.matrix(diffm)
long <- expand.grid(ti = seq_len(nrow(mm)), cj = seq_len(ncol(mm)))
long$aa   <- mm[cbind(long$ti, long$cj)]
long$row  <- rownames(mm)[long$ti]
long$site <- colnames(mm)[long$cj]
long$y    <- tipy[long$ti]
long$x    <- x0 + (long$cj - 0.5) * colw
is_priv <- (long$site %in% ps_cols) & (long$row == STUDY) & !(long$aa %in% c("·","-"))
long$status <- ifelse(long$aa == "·", "match",
               ifelse(long$aa == "-", "gap",
               ifelse(is_priv, "private", "sub")))
long$fill <- STATUS_COL[long$status]

p2 <- p +
  geom_tile(data = long, aes(x = x, y = y, fill = fill),
            width = colw, height = 0.62, inherit.aes = FALSE) +
  scale_fill_identity() +
  coord_cartesian(clip = "off")

priv_cells <- long[long$status == "private", ]
p2 <- p2 + geom_tile(data = priv_cells, aes(x = x, y = y), fill = STATUS_COL[["private"]],
                     width = colw, height = PRIV_H, inherit.aes = FALSE)
message(sprintf("Study-specific substitutions: %d -> %s (row %s)", nrow(priv_cells),
                paste(priv_cells$site, priv_cells$aa, sep = "=", collapse = ", "),
                paste(unique(priv_cells$row), collapse = "/")))

strip_y <- ymax + 1.0; strip_h <- 0.7
dstrip <- data.frame(cj = dom_ann$col_index, dom = dom_ann$domain)
dstrip$x <- x0 + (dstrip$cj - 0.5) * colw
dstrip$fill <- DOM_COL[dstrip$dom]
p2 <- p2 + geom_tile(data = dstrip, aes(x = x, y = strip_y, fill = fill),
                     width = colw, height = strip_h, inherit.aes = FALSE)

arch_y <- ymax + 5.7; arch_h <- 0.9
arch$xmin <- xarch(arch$from); arch$xmax <- xarch(arch$to + 1); arch$fill <- DOM_COL[arch$dom]
p2 <- p2 +
  geom_rect(data = arch, aes(xmin = xmin, xmax = xmax, ymin = arch_y - arch_h/2,
                             ymax = arch_y + arch_h/2, fill = fill),
            inherit.aes = FALSE, color = "white", linewidth = LW_MIN)
arch_lab <- arch %>% group_by(dom) %>% summarise(xc = mean(c(min(xmin), max(xmax))), .groups = "drop")
arch_in  <- arch_lab %>% filter(dom %in% c("NTD","RBD","HR1","HR2"))
arch_out <- arch_lab %>% filter(dom %in% c("FP","TM"))
arch_out$yleg <- arch_y - 2.10
p2 <- p2 +
  geom_text(data = arch_in, aes(x = xc, y = arch_y), inherit.aes = FALSE,
            size = C_domlab, color = "black", label = arch_in$dom, family = "Arial") +
  geom_segment(data = arch_out, aes(x = xc, xend = xc, y = arch_y - arch_h/2, yend = yleg + 0.45),
               inherit.aes = FALSE, color = "black", linewidth = LW_MIN) +
  geom_text(data = arch_out, aes(x = xc, y = yleg), inherit.aes = FALSE,
            size = C_domlab, color = "black", label = arch_out$dom, family = "Arial")

lab_df <- data.frame(pos = ps_pos, lab = ps_labs,
                     xcol = x0 + (match(ps_cols, colnames(mm)) - 0.5) * colw,
                     xreal = xarch(ps_pos))
arrow_y <- strip_y + 0.62
lab_df <- lab_df[order(lab_df$xcol), ]
lab_df$xlab <- lab_df$xcol
lab_df$laby <- arrow_y + c(1.30, 2.70)[seq_len(nrow(lab_df))]
p2 <- p2 +
  geom_segment(data = lab_df, aes(x = xlab, xend = xcol, y = laby - 0.50, yend = arrow_y + 0.20),
               inherit.aes = FALSE, color = "#D01C1F", linewidth = 0.25) +
  geom_point(data = lab_df, aes(x = xcol, y = arrow_y), inherit.aes = FALSE,
             shape = 25, fill = "#D01C1F", color = "#D01C1F", size = 1.00) +
  geom_text(data = lab_df, aes(x = xlab, y = laby, label = paste0(lab, " (NTD)")),
            inherit.aes = FALSE, size = C_arrow, color = "#D01C1F", fontface = "bold", family = "Arial")

lx <- x0 + HW_TOTAL + 0.25 * treerange; sw <- colw * 1.4; gap <- 0.9
leg_tile <- function(p, title, keys, cols, x, ytop) {
  p <- p + annotate("text", x = x, y = ytop, hjust = 0, size = C_legtitle, fontface = "bold", label = title, family = "Arial")
  yy <- ytop - gap
  for (k in seq_along(keys)) {
    p <- p + annotate("rect", xmin = x, xmax = x + sw, ymin = yy - 0.32, ymax = yy + 0.32,
                      fill = cols[k], color = "#888888", linewidth = LW_MIN) +
             annotate("text", x = x + sw*1.5, y = yy, hjust = 0, size = C_legtext, label = names(keys)[k], family = "Arial")
    yy <- yy - gap
  }
  p
}
stat_keys <- c("Matches consensus" = "match", "Substitution" = "sub", "Patient-private" = "private", "Gap/indel" = "gap")
p2 <- leg_tile(p2, "Residue status", stat_keys, STATUS_COL[stat_keys], lx, ymax + 2.2)
dom_keys <- setNames(names(DOM_COL), paste0(names(DOM_COL), " (", DOM_RANGE, ")"))
p2 <- leg_tile(p2, "Spike domain", dom_keys, DOM_COL[dom_keys], lx, ymax - 3.8)

p2 <- p2 +
  scale_y_continuous(expand = expansion(mult = c(0.018, 0.018))) +
  hexpand(0.05, direction =  1) +
  theme(plot.margin = margin(2, 6, 2, 2))

tmp_gc <- tempfile(fileext = ".png")
png(tmp_gc, width = 10, height = 10, units = "in", res = 72, type = PNG_TYPE)
gC <- ggplotGrob(p2)
dev.off(); unlink(tmp_gc)

# ============================================================================
# Composite: panel A cropped to its ink, panel B below, center legend
# ============================================================================

W <- WIDTH_IN

fan_replay_dev <- function(width, height) {
  open_pdf8(tempfile(fileext = ".pdf"), width, height)
  dev.control("enable")
}
measure_fan_ink <- function(size_in, dpi = 600) {
  tmp <- tempfile(fileext = ".png")
  png(tmp, width = size_in, height = size_in, units = "in", res = dpi, type = PNG_TYPE, pointsize = 8)
  draw_hadv_fan(); dev.off()
  img <- png::readPNG(tmp); unlink(tmp)
  if (length(dim(img)) == 2L) img <- array(img, dim = c(dim(img), 1L))
  dm <- dim(img)
  nz <- if (dm[3] >= 3L) pmin(img[, , 1], img[, , 2], img[, , 3]) < 0.98 else img[, , 1] < 0.98
  if (dm[3] >= 4L) nz <- nz & (img[, , 4] > 0.02)
  rr <- which(apply(nz, 1, any)); cc <- which(apply(nz, 2, any))
  if (!length(rr) || !length(cc)) return(list(left = 0, right = 0, top = 0, bottom = 0))
  list(left   = (min(cc) - 1L)/dpi, right  = (dm[2] - max(cc))/dpi,
       top    = (min(rr) - 1L)/dpi, bottom = (dm[1] - max(rr))/dpi)
}
.ink  <- measure_fan_ink(FAN_SIZE_IN)
.cutL <- max(0, .ink$left - FAN_PAD_IN); .cutR <- max(0, .ink$right  - FAN_PAD_IN)
.cutT <- max(0, .ink$top  - FAN_PAD_IN); .cutB <- max(0, .ink$bottom - FAN_PAD_IN)
fanW  <- FAN_SIZE_IN - .cutL - .cutR
fanH  <- FAN_SIZE_IN - .cutT - .cutB
message(sprintf("Panel A cropped to ink: %.2f x %.2f in", fanW, fanH))

xB_hcov <- MARGIN_SIDE;  yB_hcov <- MARGIN_BOT
hcov_top <- yB_hcov + hC
yFan <- hcov_top + GAP_ROW
xFan <- (W - fanW) / 2
fanTop <- yFan + fanH
H <- fanTop + LABEL_PAD_Y + MARGIN_TOP

gp_lab <- gpar(fontfamily = "Arial", fontface = "bold", fontsize = LABEL_SIZE)

draw_center_legend <- function(x0, ytop, dy = 0.115) {
  labs <- c(names(CENTER_COLS_FIG1), "Reference")
  for (k in seq_along(labs)) {
    y <- ytop - (k - 1) * dy
    if (k <= length(CENTER_COLS_FIG1)) {
      grid.points(unit(x0, "in"), unit(y, "in"), pch = 21, size = unit(8 * STUDY_CEX, "points"),
                  gp = gpar(fill = CENTER_COLS_FIG1[[k]], col = STUDY_OUTLINE, lwd = STUDY_LWD))
    } else {
      grid.points(unit(x0, "in"), unit(y, "in"), pch = 16, size = unit(8 * 0.80, "points"),
                  gp = gpar(col = reference_background_col))
    }
    grid.text(labs[k], x = unit(x0 + 0.075, "in"), y = unit(y, "in"), just = c("left", "centre"),
              gp = gpar(fontfamily = "Arial", fontsize = FONT_B, col = STUDY_LABEL_COL))
  }
}

draw_all <- function() {
  grid.newpage()
  pushViewport(viewport(x = unit(xFan, "in"), y = unit(yFan, "in"),
                        width = unit(fanW, "in"), height = unit(fanH, "in"),
                        just = c("left", "bottom"), clip = "on"))
  pushViewport(viewport(x = unit(-.cutL, "in"), y = unit(-.cutB, "in"),
                        width = unit(FAN_SIZE_IN, "in"), height = unit(FAN_SIZE_IN, "in"),
                        just = c("left", "bottom")))
  gridGraphics::grid.echo(draw_hadv_fan, newpage = FALSE, prefix = "fanA", device = fan_replay_dev)
  popViewport(2)
  .leg_dy <- 0.115
  .leg_r  <- 0.375 * 8 * STUDY_CEX / 72
  draw_center_legend(5.729 + .leg_r, yFan + 0.12 + 3 * .leg_dy, dy = .leg_dy)
  pushViewport(viewport(x = unit(xB_hcov, "in"), y = unit(yB_hcov, "in"),
                        width = unit(wC, "in"), height = unit(hC, "in"),
                        just = c("left", "bottom")))
  grid.draw(gC); upViewport()
  grid.text("A", x = unit(MARGIN_SIDE, "in"), y = unit(fanTop + LABEL_PAD_Y, "in"), just = c("left", "bottom"), gp = gp_lab)
  grid.text("B", x = unit(MARGIN_SIDE, "in"), y = unit(hcov_top + LABEL_PAD_Y, "in"), just = c("left", "bottom"), gp = gp_lab)
}

pdf_path <- file.path(OUT, "Figure4.pdf")
open_pdf8(pdf_path, W, H)
draw_all(); dev.off()
message(sprintf("Canvas %.2f x %.2f in (%.1f x %.1f mm)", W, H, W * 25.4, H * 25.4))
message("PDF: ", pdf_path)
