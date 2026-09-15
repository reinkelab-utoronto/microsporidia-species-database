#!/usr/bin/env Rscript
###############################################################################
# Clustering threshold sweep: three panels from threshold_sweep.tsv
#
#   A  species split across >=2 clusters   (over-splitting; worst at 100% id)
#   B  clusters holding >=2 species        (over-lumping; grows as id falls)
#   C  adjusted Rand index                 (single objective; higher is better)
#
# HOW IT WORKS
#   The analysis lives in sweep_cluster_thresholds_3.py, which writes one row
#   per (min_cov, identity) grid cell to threshold_sweep.tsv. This script only
#   draws: nothing is recomputed here, so a number in the figure cannot drift
#   from a number in the sweep. Fractions, not raw counts, are plotted --
#   the coverage filter changes how many species are scored, and only the
#   fractions are comparable across the grid.
#
#   Panel C uses the adjusted Rand index rather than the V-measure: both pick
#   out the same trade-off region, but ARI is corrected for chance agreement,
#   so it does not reward the near-trivial clusterings at the grid edges the
#   way an uncorrected score can. V, homogeneity and completeness are still
#   read from the TSV and reported in the QC output below, and the best-by-V
#   cell is compared against the best-by-ARI cell so a disagreement is loud.
#
#   The dotted guide marks the chosen threshold (CFG$chosen_id), not the ARI
#   argmax -- see the CFG comment for the plateau argument behind that.
#
# mBio CONSTRAINTS
#   width      7.5 in (two-column maximum for ASM journals)
#   height     3.0 in -- one row of three panels plus a shared legend below.
#              The legend goes UNDER the panels: collected on top it sits in
#              the same strip as the A/B/C tags and the two collide.
#   resolution 300 dpi for the raster copy; a vector PDF is written too
#   type size  nothing below 8 pt; geom_text sizes are in MILLIMETRES, where
#              8 pt = 2.845 mm, so the annotation layer is set to that
#
# Usage:  Rscript fig_cluster_sweep.R [sweep/threshold_sweep.tsv]
###############################################################################

CFG <- list(
  sweep_file = "sweep/threshold_sweep.tsv",

  fig_width  = 7.5,     # inches, mBio two-column maximum
  fig_height = 3.0,
  dpi        = 300,
  base_pt    = 8,       # smallest type anywhere in the figure
  min_text_mm = 8 * 25.4 / 72,   # 8 pt expressed in mm, for geom_text sizes

  # The dotted guide marks the CHOSEN operating threshold, which is not the
  # best-by-ARI grid cell: ARI weights lumping and splitting symmetrically,
  # and 98% is preferred as the finest threshold still on the ARI plateau --
  # it resolves more clusters and roughly halves the mixed-cluster fraction,
  # for a modest rise in split species (recoverable downstream from names),
  # while the step to 99% is where ARI collapses. The QC output below prints
  # the plateau test and the per-step costs so the results text can quote
  # them; the best-by-ARI cell is still reported there as a cross-check.
  chosen_id   = 0.98,   # NULL draws no guide
  # identities whose ARI is within this fraction of a coverage threshold's
  # maximum count as "as good as the best" for the plateau test
  plateau_tol = 0.05,

  # What the min_cov column actually thresholded. The TSV does not record
  # whether the sweep ran with --length-field length (trimmed sequence
  # length, bp) or cm_cov (covariance model coverage), so the legend title
  # is set here and must match the sweep command that made the TSV.
  filter_label = "Sequence length",
  filter_unit  = "bp",   # "nt" if the sweep thresholded model coverage

  # Okabe-Ito, in the order the coverage thresholds appear. The first two
  # match the named/provisional blue and orange used across the other figures,
  # so the palette reads as one family without implying the same categories.
  pal = c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9",
          "#D55E00", "#999999"),

  outdir = "figure_output"
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) CFG$sweep_file <- args[1]

for (p in c("ggplot2", "patchwork", "scales"))
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
dir.create(CFG$outdir, showWarnings = FALSE, recursive = TRUE)
options(width = 200)
qc <- function(...) cat(sprintf(...), sep = "")

## ------------------------------------------------------------------ data ----
if (!file.exists(CFG$sweep_file))
  stop("sweep table not found: ", CFG$sweep_file,
       "  (run sweep_cluster_thresholds_3.py first)", call. = FALSE)
d <- read.delim(CFG$sweep_file, check.names = FALSE)

need <- c("min_cov", "identity", "frac_species_split", "frac_clusters_mixed",
          "adjusted_rand", "v_measure", "homogeneity", "completeness",
          "n_species", "n_clusters_total", "split_species", "mixed_clusters",
          "n_input_seqs", "n_scored_seqs")
miss <- setdiff(need, names(d))
if (length(miss))
  stop("threshold_sweep.tsv is missing column(s): ",
       paste(miss, collapse = ", "),
       "  -- was it written by an older sweep script?", call. = FALSE)

d <- d[order(d$min_cov, d$identity), ]
covs <- sort(unique(d$min_cov))
if (length(covs) > length(CFG$pal))
  stop(length(covs), " coverage thresholds but only ", length(CFG$pal),
       " palette colours -- extend CFG$pal", call. = FALSE)
# ">=0 nt" says nothing; call the unfiltered run what it is
d$cov_lab <- factor(ifelse(d$min_cov == 0, "unfiltered",
                           sprintf("\u2265%d %s", d$min_cov, CFG$filter_unit)),
                    levels = ifelse(covs == 0, "unfiltered",
                                    sprintf("\u2265%d %s", covs, CFG$filter_unit)))
pal_use <- setNames(CFG$pal[seq_along(covs)], levels(d$cov_lab))

# best-by-ARI cell, kept as a QC cross-check against the chosen threshold;
# ties broken toward the least-filtered subset
best <- d[order(-d$adjusted_rand, d$min_cov), ][1, ]

## ---------------------------------------------------------------- panels ----
# One theme for all three panels; axis.text inherits a RELATIVE size from
# `text`, so everything that must clear 8 pt is set explicitly.
base_theme <- theme_classic(base_size = CFG$base_pt) +
  theme(
    text = element_text(size = CFG$base_pt),
    axis.text = element_text(size = CFG$base_pt, colour = "black"),
    axis.title = element_text(size = CFG$base_pt),
    legend.text = element_text(size = CFG$base_pt),
    legend.title = element_text(size = CFG$base_pt, face = "bold"),
    plot.tag = element_text(size = CFG$base_pt + 2, face = "bold"),
    plot.margin = margin(2, 6, 1, 2),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(2, "pt"),
    legend.key.size = unit(9, "pt")
  )

# identity thresholds are proportions in the TSV; percent is how they are
# spoken about ("a 97% OTU"), so the axis is drawn in percent
x_scale <- scale_x_continuous(
  labels = function(x) sprintf("%g", 100 * x),
  breaks = c(0.96, 0.97, 0.98, 0.99, 1.00))

panel_line <- function(yvar, ylab, percent = TRUE) {
  p <- ggplot(d, aes(x = identity, y = .data[[yvar]],
                     colour = cov_lab, group = cov_lab)) +
    { if (!is.null(CFG$chosen_id))
        geom_vline(xintercept = CFG$chosen_id, linetype = "dotted",
                   colour = "grey45", linewidth = 0.35) } +
    geom_line(linewidth = 0.4) +
    geom_point(size = 1.1) +
    scale_colour_manual(values = pal_use, name = CFG$filter_label) +
    x_scale +
    labs(x = "Clustering identity (%)", y = ylab) +
    base_theme
  if (percent)
    p <- p + scale_y_continuous(labels = function(x) sprintf("%g", 100 * x),
                                limits = c(0, NA),
                                expand = expansion(mult = c(0.01, 0.05)))
  p
}

pA <- panel_line("frac_species_split",
                 "Species split across \u22652 clusters (%)")
pB <- panel_line("frac_clusters_mixed",
                 "Clusters containing \u22652 species (%)")

# C: the axis starts at 0 -- ARI is a similarity, and a curve that only ever
# lives in the top tenth of a truncated axis overstates how much the
# thresholds matter. The best cell is deliberately NOT ringed or labelled:
# the figure shows the whole trade-off and the chosen guide; the argmax is a
# QC number, not the argument.
pC <- panel_line("adjusted_rand", "Adjusted Rand index", percent = FALSE) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.01, 0.05)))

## -------------------------------------------------------------- assemble ----
fig <- pA + pB + pC +
  plot_layout(nrow = 1, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom",
        legend.justification = "centre",
        legend.direction = "horizontal")

## ============================================================================
## NUMBERS BEHIND THE FIGURE
## ============================================================================
# Taken from the SAME rows the panels plotted, so a number quoted in the text
# and a number drawn in the figure cannot drift apart.

# the console may not be UTF-8 (it prints <U+2265> under a C locale), so the
# QC lines use ">=" while the figure keeps the real glyph
ascii <- function(x) gsub("\u2265", ">=", as.character(x))
fmt_cell <- function(r) sprintf(
  "id=%g%% cov=%s -> ARI=%.3f V=%.3f (hom=%.3f, com=%.3f), %d/%d species split (%.1f%%), %d/%d clusters mixed (%.1f%%)",
  100 * r$identity, ascii(r$cov_lab), r$adjusted_rand, r$v_measure,
  r$homogeneity, r$completeness,
  r$split_species, r$n_species, 100 * r$frac_species_split,
  r$mixed_clusters, r$n_clusters_with_binomials, 100 * r$frac_clusters_mixed)

best_v <- d[order(-d$v_measure, d$min_cov), ][1, ]

cat("\n================ NUMBERS BEHIND THE FIGURE ================\n")
qc("grid: %d cells = %d coverage thresholds x %d identities\n",
   nrow(d), length(covs), length(unique(d$identity)))
qc("species scored: %d-%d across the grid (coverage filter changes the set)\n",
   min(d$n_species), max(d$n_species))

# The n values the legend quotes come from the LEAST-FILTERED row: that is
# the dataset before any length cutoff, which is what "n = X sequences from
# Y species" should mean. Filtered rows describe subsets, not the dataset.
r0 <- d[d$min_cov == min(covs), ][1, ]
cat("\nfor the figure legend\n")
qc("  identity thresholds swept:   %g%% to %g%% (%d values)\n",
   100 * min(d$identity), 100 * max(d$identity), length(unique(d$identity)))
qc("  length cutoffs (colours):    %s\n",
   paste(vapply(levels(d$cov_lab), ascii, character(1)), collapse = ", "))
qc("  sequences clustered:         %d%s\n", r0$n_input_seqs,
   if (min(covs) == 0) "" else
     sprintf("  (NOTE: no unfiltered run in the grid; this is the >=%d subset)",
             min(covs)))
qc("  scored against binomials:    n = %s sequences from %s species\n",
   format(r0$n_scored_seqs, big.mark = ","),
   format(r0$n_species, big.mark = ","))
qc("  (unnamed/sp. sequences clustered but not scored: %d)\n",
   r0$n_input_seqs - r0$n_scored_seqs)
qc("\nbest adjusted-Rand cell (cross-check, not drawn):\n  %s\n",
   fmt_cell(best))
qc("best V-measure cell (cross-check, not drawn):\n  %s\n", fmt_cell(best_v))
if (best$identity == best_v$identity && best$min_cov == best_v$min_cov) {
  qc("ARI and V-measure agree on the operating point.\n")
} else {
  qc("NOTE: ARI and V-measure pick DIFFERENT cells -- look at both before\n")
  qc("quoting a threshold; V is not chance-corrected and can favour finer\n")
  qc("clusterings when many species are singletons.\n")
}

# The chosen threshold and its justification, in numbers the results text can
# quote. The claim being tested: the chosen identity is still on the ARI
# plateau (within plateau_tol of each coverage threshold's own maximum), and
# the NEXT identity step is where agreement collapses -- so the chosen
# threshold buys the finest clustering the plateau allows.
if (!is.null(CFG$chosen_id)) {
  qc("\nchosen operating threshold (the dotted guide): %g%% identity\n",
     100 * CFG$chosen_id)
  if (!CFG$chosen_id %in% d$identity)
    warning("CFG$chosen_id = ", CFG$chosen_id,
            " is not one of the identities in the sweep", call. = FALSE)

  qc("\nplateau test (ARI within %g%% of each coverage threshold's maximum)\n",
     100 * CFG$plateau_tol)
  for (cl in levels(d$cov_lab)) {
    s <- d[d$cov_lab == cl, ]
    mx <- max(s$adjusted_rand)
    edge <- max(s$identity[s$adjusted_rand >= (1 - CFG$plateau_tol) * mx])
    r <- s[s$identity == CFG$chosen_id, ]
    if (nrow(r) != 1) next
    qc("  %-12s ARI=%.3f (%3.0f%% of max %.3f)  plateau edge = %g%%  -> %s\n",
       ascii(cl), r$adjusted_rand, 100 * r$adjusted_rand / mx, mx, 100 * edge,
       if (CFG$chosen_id <= edge) "chosen id ON plateau"
       else "chosen id OFF plateau")
  }

  # what each identity step buys and costs, for the coverage threshold with
  # the best overall ARI (the curve a reader's eye lands on in panel C)
  s <- d[d$cov_lab == best$cov_lab, ]
  s <- s[order(s$identity), ]
  qc("\nstep-by-step along the best coverage threshold (%s)\n",
     ascii(best$cov_lab))
  for (i in seq_len(nrow(s)))
    qc("  %5g%%  ARI=%.3f%s  clusters=%5d  split=%5.1f%%  mixed=%5.1f%%\n",
       100 * s$identity[i], s$adjusted_rand[i],
       if (i > 1) sprintf(" (%+.3f)", s$adjusted_rand[i] - s$adjusted_rand[i - 1])
       else "         ",
       s$n_clusters_total[i], 100 * s$frac_species_split[i],
       100 * s$frac_clusters_mixed[i])
}

stats <- data.frame(
  quantity = c("grid cells", "coverage thresholds", "identity thresholds",
               "chosen identity (dotted guide)",
               "best-ARI identity", "best-ARI min_cov", "best ARI", "its V",
               "species split at best-ARI (%)", "clusters mixed at best-ARI (%)",
               "best-V identity", "best-V min_cov"),
  value = c(nrow(d), length(covs), length(unique(d$identity)),
            if (is.null(CFG$chosen_id)) NA else CFG$chosen_id,
            best$identity, best$min_cov, best$adjusted_rand, best$v_measure,
            100 * best$frac_species_split, 100 * best$frac_clusters_mixed,
            best_v$identity, best_v$min_cov))
stats <- rbind(stats, data.frame(
  quantity = c("legend n: sequences clustered (least-filtered)",
               "legend n: sequences scored (least-filtered)",
               "legend n: species scored (least-filtered)",
               "legend: lowest identity", "legend: highest identity"),
  value = c(r0$n_input_seqs, r0$n_scored_seqs, r0$n_species,
            min(d$identity), max(d$identity))))
write.table(stats, file.path(CFG$outdir, "fig_cluster_sweep_numbers.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
qc("\nwritten %s\n", file.path(CFG$outdir, "fig_cluster_sweep_numbers.tsv"))

## ----------------------------------------------------------------- output ----
out_pdf <- file.path(CFG$outdir, "fig_cluster_sweep.pdf")
out_png <- file.path(CFG$outdir, "fig_cluster_sweep.png")
ggsave(out_pdf, fig, width = CFG$fig_width, height = CFG$fig_height,
       units = "in", device = if (capabilities("cairo")) cairo_pdf else pdf)
ggsave(out_png, fig, width = CFG$fig_width, height = CFG$fig_height,
       units = "in", dpi = CFG$dpi)
qc("\nwritten\n  %s\n  %s  (%.1f x %.1f in, %d dpi = %d x %d px)\n",
   out_pdf, out_png, CFG$fig_width, CFG$fig_height, CFG$dpi,
   round(CFG$fig_width * CFG$dpi), round(CFG$fig_height * CFG$dpi))
qc("smallest type: %g pt (theme) / %.2f mm = %g pt (text layers)\n",
   CFG$base_pt, CFG$min_text_mm, CFG$min_text_mm * 72 / 25.4)
