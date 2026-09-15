#!/usr/bin/env Rscript
###############################################################################
# Species that fall into more than one 98 % cluster -- the figure
#
# All numbers come from split_species_identity.py (run it first); this script
# only reads its TSVs and draws. Nothing is recomputed here, so a number in
# the figure and a number in the text cannot drift apart.
#
#   A  every species in >= 2 clusters, species on the x axis. One point per
#      cluster pair = mean identity between the two clusters' sequences,
#      sized by the sequence count of the smaller cluster; hollow grey points
#      = within-cluster identity; black diamond = species mean of the pair
#      values ("average of the averages"). Dashed line = clustering threshold.
#   B  cluster x cluster identity matrix for the first detail species
#   C  the same for the second
#      Rows and columns share one order (average-linkage clustering of
#      100 - identity) so lineages form blocks. Cells are coloured in flat
#      classes at the cut-offs the conclusions rest on (CFG$heat_bins; the
#      threshold is a class boundary). The diagonal is a cluster against
#      itself (100 %). Labels are "accession (n sequences)"; an asterisk
#      marks a centroid shorter than the length floor used in the averages.
#
# mBio CONSTRAINTS
#   width      7.5 in (two-column maximum for ASM journals)
#   resolution 300 dpi for the raster copy; a vector PDF is written too
#   type size  nothing below 8 pt (geom_text sizes are in mm: 8 pt = 2.845 mm)
#   colours    Okabe-Ito for marks; single-hue blues for the identity classes
#
# Usage:  Rscript fig_split_species_identity.R [indir]
###############################################################################

CFG <- list(
  indir = "figure_output",              # where the .py wrote its TSVs
  stem  = "fig_split_species_identity",
  threshold = 98,                       # drawn as the dashed line

  # heatmap classes: edges -> < 90, 90-95, 95-98, 98-99, >= 99
  heat_bins = c(90, 95, 98, 99),
  # which matrices become B and C; NULL = the matrix files in indir, in the
  # order the .py wrote them (most clusters first)
  detail_species = NULL,

  max_point_mm = 3.2,     # largest point in A (legend key too)
  show_within  = TRUE,    # hollow within-cluster points in A

  fig_width  = 7.5,
  # Panel sizes are PANEL areas in inches (patchwork absolute units); the
  # axis labels, titles and legends are added around them, so the canvas
  # height below is the sum of panels plus those allowances.
  panelA_in   = 2.0,      # height of A's panel (its rotated labels add ~1.7)
  panelA_w_in = 4.75,     # width of A's panel (y title + legend take the rest)
  heat_cell_in = 0.100,   # one heatmap cell; 8 pt labels are 0.111 in tall
                          # but digits only fill ~0.08 in, so rows stay legible
  dpi        = 300,
  base_pt    = 8,
  min_text_mm = 8 * 25.4 / 72,
  outdir = "figure_output"
)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) CFG$indir <- args[1]

if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
suppressPackageStartupMessages(library(ggplot2))
# ggtext italicises the binomial while leaving "(nom. nud.)" and the cluster
# count upright; without it the labels are drawn plain
has_ggtext <- requireNamespace("ggtext", quietly = TRUE)
dir.create(CFG$outdir, showWarnings = FALSE, recursive = TRUE)
qc <- function(...) cat(sprintf(...), sep = "")

PAL <- c(between = "#0072B2", within = "#999999", mean = "black",
         line = "#D55E00")

## --------------------------------------------------------------- inputs ----
rd <- function(name) {
  f <- file.path(CFG$indir, paste0(CFG$stem, "_", name, ".tsv"))
  if (!file.exists(f)) stop("missing ", f, " -- run split_species_identity.py first",
                            call. = FALSE)
  read.delim(f, sep = "\t", quote = "", check.names = FALSE,
             stringsAsFactors = FALSE, encoding = "UTF-8")
}
pairs  <- rd("cluster_pairs")
spsum  <- rd("species")
clus   <- rd("clusters")
within <- tryCatch(rd("within_cluster"), error = function(e) NULL)
if (!is.null(within) && !nrow(within)) within <- NULL

mfiles <- list.files(CFG$indir, pattern = paste0("^", CFG$stem, "_matrix_.*\\.tsv$"),
                     full.names = TRUE)
mfiles <- mfiles[order(file.info(mfiles)$mtime)]
if (!is.null(CFG$detail_species)) {
  slug <- gsub("[^A-Za-z0-9]+", "_", CFG$detail_species)
  mfiles <- file.path(CFG$indir, paste0(CFG$stem, "_matrix_", slug, ".tsv"))
}
if (length(mfiles) < 2) stop("need two matrix files for panels B and C", call. = FALSE)
mfiles <- mfiles[1:2]
read_matrix <- function(f) {
  m <- as.matrix(read.delim(f, sep = "\t", check.names = FALSE, row.names = 1))
  storage.mode(m) <- "numeric"
  m
}
mats <- lapply(mfiles, read_matrix)
# species name of a matrix: the species whose clusters are its row names
mat_species <- vapply(mats, function(m)
  clus$species[match(rownames(m)[1], clus$cluster)], character(1))
qc("panels B, C: %s\n", paste(mat_species, collapse = ", "))

strip_ver <- function(x) sub("\\.\\d+$", "", x)
label_of <- function(sp, k) {
  bin <- sub("^(\\S+\\s+\\S+).*$", "\\1", sp)
  rest <- trimws(sub("^\\S+\\s+\\S+", "", sp))
  if (has_ggtext)
    sprintf("*%s*%s (%d)", bin, if (nzchar(rest)) paste0(" ", rest) else "", k)
  else sprintf("%s (%d)", sp, k)
}

base_theme <- theme_classic(base_size = CFG$base_pt) +
  theme(text = element_text(size = CFG$base_pt),
        axis.text = element_text(size = CFG$base_pt, colour = "black"),
        axis.title = element_text(size = CFG$base_pt),
        legend.text = element_text(size = CFG$base_pt),
        legend.title = element_text(size = CFG$base_pt),
        legend.key.size = unit(8, "pt"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(3, "pt"),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.tag = element_text(size = CFG$base_pt + 2, face = "bold"),
        plot.margin = margin(4, 4, 2, 4))

## -------------------------------------------------------------- panel A ----
# species left to right from most to least similar
spsum <- spsum[order(-spsum$identity_mean, spsum$species), ]
lv <- spsum$species
spsum$label <- mapply(label_of, spsum$species, spsum$n_clusters)
pairs$species <- factor(pairs$species, levels = lv)
spsum$species <- factor(spsum$species, levels = lv)
if (!is.null(within)) within$species <- factor(within$species, levels = lv)

y_lo <- floor(min(c(pairs$identity_mean, within$identity_mean)) / 5) * 5
jit <- position_jitter(width = 0.22, height = 0, seed = 1)

pA <- ggplot(pairs, aes(x = species, y = identity_mean)) +
  geom_hline(yintercept = CFG$threshold, linetype = "22", linewidth = 0.35,
             colour = PAL[["line"]])
if (CFG$show_within && !is.null(within))
  pA <- pA + geom_point(data = within, aes(x = species, y = identity_mean),
                        shape = 1, size = 1.2, stroke = 0.35,
                        colour = PAL[["within"]], position = jit)
pA <- pA +
  geom_point(aes(size = smaller), colour = PAL[["between"]], alpha = 0.55,
             position = jit) +
  geom_point(data = spsum, aes(x = species, y = identity_mean), shape = 23,
             size = 2.2, fill = PAL[["mean"]], colour = "white", stroke = 0.3) +
  scale_size_area(name = "sequences in the\nsmaller cluster",
                  max_size = CFG$max_point_mm,
                  breaks = function(l) {
                    b <- c(1, 2, 5, 10, 20, 50, 100, 200, 500); b[b <= max(l)] }) +
  # drop = FALSE: the within layer is trained first and lacks the species
  # whose clusters are all singletons; the default would push those to the end
  scale_x_discrete(labels = setNames(spsum$label, spsum$species), drop = FALSE) +
  # the threshold is named on a secondary axis rather than by a label inside
  # the panel, where it collided with the points sitting on the line
  scale_y_continuous(limits = c(y_lo, 100),
                     breaks = seq(y_lo, 100, by = if (100 - y_lo > 30) 10 else 5),
                     labels = function(x) paste0(x, "%"),
                     sec.axis = dup_axis(name = NULL, breaks = CFG$threshold,
                                         labels = sprintf("%g%%\nthreshold",
                                                          CFG$threshold))) +
  labs(x = NULL, y = "Mean nucleotide identity between clusters") +
  base_theme +
  theme(axis.text.x = if (has_ggtext)
          ggtext::element_markdown(angle = 45, hjust = 1, vjust = 1) else
          element_text(angle = 45, hjust = 1, vjust = 1),
        axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_text(colour = PAL[["line"]], hjust = 0),
        legend.position = "right", legend.justification = c(0, 1),
        # the rotated labels run past the left edge of the panel; the margin
        # has to hold them or the first names are cut off the canvas
        plot.margin = margin(4, 4, 2, 40))

## ---------------------------------------------------------- panels B, C ----
bin_edges <- sort(CFG$heat_bins)
bin_lab <- c(sprintf("< %g%%", bin_edges[1]),
             sprintf("%g-%g%%", head(bin_edges, -1), tail(bin_edges, -1)),
             sprintf(">= %g%%", tail(bin_edges, 1)))
blues <- c("#EFF3FF", "#BDD7E7", "#6BAED6", "#2171B5", "#08306B", "#041E42")
bin_pal <- setNames(blues[seq_along(bin_lab)], bin_lab)

heat_label <- function(cl) {
  i <- match(cl, clus$cluster)
  sprintf("%s%s (%d)", strip_ver(cl),
          ifelse(!is.na(clus$centroid_short[i]) & clus$centroid_short[i] == 1,
                 "*", ""),
          clus$n_seq[i])
}

heat_panel <- function(M, sp, legend) {
  ord <- hclust(as.dist(100 - M), method = "average")$order
  cl_ord <- rownames(M)[ord]
  d <- expand.grid(a = cl_ord, b = cl_ord, stringsAsFactors = FALSE)
  d$id <- M[cbind(d$a, d$b)]
  d$a <- factor(d$a, levels = rev(cl_ord)); d$b <- factor(d$b, levels = cl_ord)
  d$bin <- factor(bin_lab[findInterval(d$id, bin_edges) + 1], levels = bin_lab)
  ggplot(d, aes(x = b, y = a, fill = bin)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    scale_fill_manual(values = bin_pal, drop = FALSE, na.value = "#E6E6E6",
                      na.translate = FALSE, name = "identity between\nclusters") +
    # the matrix is symmetric and the columns follow the row order, so the
    # column labels are dropped: rotated 8 pt accessions cost more height
    # than the cells themselves, and two label columns plus a legend cannot
    # share 7.5 in with 49 cells at 8 pt
    scale_x_discrete(labels = NULL, expand = c(0, 0)) +
    scale_y_discrete(labels = heat_label, expand = c(0, 0)) +
    # no coord_fixed(): with a fixed aspect anywhere, patchwork scales every
    # row and column to one inch-per-unit and the tall label areas shrink the
    # cells to nothing. Square cells come from the absolute panel sizes set
    # in the layout instead.
    labs(x = NULL, y = NULL, title = label_of(sp, length(cl_ord))) +
    theme_minimal(base_size = CFG$base_pt) +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = CFG$base_pt, colour = "black"),
          axis.ticks.x = element_blank(),
          legend.text = element_text(size = CFG$base_pt),
          legend.title = element_text(size = CFG$base_pt),
          legend.key.size = unit(8, "pt"),
          legend.position = legend,
          plot.title = if (has_ggtext)
            ggtext::element_markdown(size = CFG$base_pt, face = "bold", hjust = 0)
          else element_text(size = CFG$base_pt, face = "bold", hjust = 0),
          plot.tag = element_text(size = CFG$base_pt + 2, face = "bold"),
          # room on the left for the panel tag beside the row labels
          plot.margin = margin(4, 4, 2, 16))
}
pB <- heat_panel(mats[[1]], mat_species[1], legend = "none")
pC <- heat_panel(mats[[2]], mat_species[2], legend = "bottom")

## ------------------------------------------------------------- assemble ----
# patchwork was abandoned for the assembly: with a fixed-aspect panel it
# rescales every row and column to one inch-per-unit and the heatmaps shrink
# to nothing, and its absolute units apply to panels only while alignment
# then pushes the label columns off the canvas. Each panel's PANEL area is
# pinned in inches on its own gtable and the three are placed with grid.
suppressPackageStartupMessages(library(grid))
pin_panel <- function(p, w_in, h_in) {
  g <- ggplotGrob(p)
  pl <- g$layout[g$layout$name == "panel", ]
  g$widths[unique(pl$l)]  <- unit(w_in, "in")
  g$heights[unique(pl$t)] <- unit(h_in, "in")
  g
}
size_in <- function(g) c(w = convertWidth(sum(g$widths), "in", valueOnly = TRUE),
                         h = convertHeight(sum(g$heights), "in", valueOnly = TRUE))

nB <- nrow(mats[[1]]); nC <- nrow(mats[[2]])
cell <- CFG$heat_cell_in
gA <- pin_panel(pA, CFG$panelA_w_in, CFG$panelA_in)
# the class legend sits under B, the wider of the two heatmaps
gB <- pin_panel(pB + theme(legend.position = "bottom",
                           legend.direction = "horizontal"), cell * nB, cell * nB)
gC <- pin_panel(pC + theme(legend.position = "none"), cell * nC, cell * nC)
sA <- size_in(gA); sB <- size_in(gB); sC <- size_in(gC)
if (sB["w"] + sC["w"] > CFG$fig_width)
  warning(sprintf("B + C are %.2f in wide; lower heat_cell_in", sB["w"] + sC["w"]))
gap <- 0.15
fig_h <- sA["h"] + gap + max(sB["h"], sC["h"]) + 0.1
CFG$fig_height <- unname(fig_h)

draw_fig <- function() {
  grid.newpage()
  # A: top, full width, left-aligned
  pushViewport(viewport(x = unit(0, "in"), y = unit(fig_h, "in"),
                        width = unit(sA["w"], "in"), height = unit(sA["h"], "in"),
                        just = c("left", "top")))
  grid.draw(gA); upViewport()
  yB <- fig_h - sA["h"] - gap
  pushViewport(viewport(x = unit(0, "in"), y = unit(yB, "in"),
                        width = unit(sB["w"], "in"), height = unit(sB["h"], "in"),
                        just = c("left", "top")))
  grid.draw(gB); upViewport()
  pushViewport(viewport(x = unit(sB["w"] + 0.1, "in"), y = unit(yB, "in"),
                        width = unit(sC["w"], "in"), height = unit(sC["h"], "in"),
                        just = c("left", "top")))
  grid.draw(gC); upViewport()
  tag <- gpar(fontsize = CFG$base_pt + 2, fontface = "bold")
  grid.text("A", x = unit(0.05, "in"), y = unit(fig_h - 0.05, "in"),
            just = c("left", "top"), gp = tag)
  grid.text("B", x = unit(0.05, "in"), y = unit(yB - 0.02, "in"),
            just = c("left", "top"), gp = tag)
  grid.text("C", x = unit(sB["w"] + 0.15, "in"), y = unit(yB - 0.02, "in"),
            just = c("left", "top"), gp = tag)
}

## ------------------------------------------------ numbers behind the figure ----
cat("\n================ NUMBERS IN THE FIGURE ================\n")
qc("  %-40s %d\n", "species in A", nrow(spsum))
qc("  %-40s %d\n", "cluster pairs in A", nrow(pairs))
qc("  %-40s %d of %d\n", sprintf("pairs at or above %g%%", CFG$threshold),
   sum(pairs$identity_mean >= CFG$threshold), nrow(pairs))
for (i in 1:2) {
  M <- mats[[i]]; off <- M[upper.tri(M)]
  cls <- bin_lab[findInterval(off, bin_edges) + 1]
  qc("  %s (%s): %d clusters, %d short centroids (*); cluster pairs by class: %s\n",
     c("B", "C")[i], mat_species[i], nrow(M),
     sum(clus$centroid_short[match(rownames(M), clus$cluster)] == 1, na.rm = TRUE),
     paste(sprintf("%s = %d", bin_lab, table(factor(cls, levels = bin_lab))),
           collapse = ", "))
}

out_pdf <- file.path(CFG$outdir, paste0(CFG$stem, ".pdf"))
out_png <- file.path(CFG$outdir, paste0(CFG$stem, ".png"))
if (capabilities("cairo")) {
  cairo_pdf(out_pdf, width = CFG$fig_width, height = CFG$fig_height)
} else {
  pdf(out_pdf, width = CFG$fig_width, height = CFG$fig_height)
}
draw_fig(); invisible(dev.off())
png(out_png, width = CFG$fig_width, height = CFG$fig_height, units = "in",
    res = CFG$dpi, type = if (capabilities("cairo")) "cairo" else "Xlib")
draw_fig(); invisible(dev.off())
qc("\nwritten\n  %s\n  %s  (%.1f x %.1f in, %d dpi)\n", out_pdf, out_png,
   CFG$fig_width, CFG$fig_height, CFG$dpi)
