#!/usr/bin/env Rscript
###############################################################################
# Assemble the six panels into one mBio figure
#
#   A  discovery by decade      D  microsporidia species environments (Venn)
#   B  geography (world map)    E  site of infection (tissues)
#   C  host phylum              F  spore size and shape
#
# HOW IT WORKS
#   Each panel script is run unchanged, in its own environment, with ggsave()
#   masked so nothing is written and no time is spent rendering. The ggplot
#   object is then lifted out of that environment and restyled for the
#   composite. Nothing here edits the panel scripts, so they stay runnable on
#   their own and this file cannot drift away from them.
#
# mBio CONSTRAINTS
#   width      7.5 in (two-column maximum for ASM journals)
#   height     set below; kept under ~8 in so a legend fits on the same page
#   resolution 300 dpi for the raster copy; a vector PDF is written too
#   type size  nothing below 8 pt. ggplot theme sizes are already in points, so
#              base_size = 8 is the floor; geom_text sizes are in MILLIMETRES,
#              where 8 pt = 2.845 mm, so every text layer is raised to that.
#
# Usage:  Rscript assemble_figure.R
###############################################################################

CFG <- list(
  db_path   = "Microsporidia_Characteristics_Database_merge_pro4.xlsx",
  api_file  = "host_taxonomy_output/01_species_host_long.tsv",
  pred_file = "predicted_environment/01_species_host_predicted.tsv",

  panels = c(
    A = "fig1_discovery_by_decade.R",
    B = "fig_locality_map_2.R",
    C = "fig_host_phylum_4.R",
    D = "fig_environment_venn_6.R",
    E = "fig_tissue_sites_2.R",
    F = "fig_spore_size_shape_2.R"),

  fig_width  = 7.5,     # inches, mBio two-column maximum
  fig_height = 9.0,     # taller, so no phylum has to be pooled away
  dpi        = 300,
  base_pt    = 8,       # smallest type anywhere in the figure
  # EXCEPT the Venn region counts: at 8 pt they overflow the narrow lenses.
  # Set venn_pt = base_pt to restore the 8 pt floor everywhere.
  venn_pt    = 7,
  # Everything drawn as a symbol in panel F -- legend keys AND plotted points --
  # is scaled by this factor, so the panel stays internally consistent.
  spore_key_scale = 0.70,
  # The panel script draws no-coil-data points at 1.8 mm and coil-data points at
  # 2.3 mm, which reads as two sizes of every shape. TRUE draws both at the
  # smaller size; the grey fill still tells them apart.
  spore_uniform_points = TRUE,
  min_text_mm = 8 * 25.4 / 72,   # 8 pt expressed in mm, for geom_text sizes

  # C and E are long category lists. Rotating them to vertical bars was tried
  # and abandoned: coord_flip() fights the kingdom facets in C, and 45-degree
  # category labels eat more height than the horizontal bars ever did. They are
  # kept horizontal and given a full row of their own instead.
  # Pooling keeps each list short enough for a half-width slot.
  phylum_min_n = 1,     # no pooling: every phylum keeps its own bar
  tissue_min_n = 12,

  outdir = "figure_output"
)

for (p in c("ggplot2", "patchwork"))
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
# ggtext lets one label carry several colours, so the breakdown after each
# total can be tinted to match the stack segments. Without it the same numbers
# are drawn in plain grey.
if (!requireNamespace("ggtext", quietly = TRUE))
  try(install.packages("ggtext"), silent = TRUE)
has_ggtext <- requireNamespace("ggtext", quietly = TRUE)
dir.create(CFG$outdir, showWarnings = FALSE, recursive = TRUE)
options(width = 200)
qc <- function(...) cat(sprintf(...), sep = "")

## ------------------------------------------------- canonical naming rule ----
# Named vs provisional is defined ONCE, in classify_species.R. Every panel
# script sources that file itself, and so does the flowchart, so the
# composite cannot disagree with a standalone figure or with Fig S1. The
# old machinery that injected a classifier into each panel is gone --
# there is nothing left to override. Sourced here too, so a missing file
# fails immediately rather than one panel at a time.
source("classify_species.R")

## --------------------------------------------------------------- runner ----
# Inject CFG overrides straight after the panel's own CFG list, so the panel
# script is used verbatim and only its paths and a few knobs change.
run_panel <- function(path, overrides = character(0), tag = "panel") {
  if (!file.exists(path)) stop("panel script not found: ", path, call. = FALSE)
  # Which FILE ran, and when it was last edited. The assembler executes exactly
  # what CFG$panels names -- editing a differently named copy of a panel script
  # changes nothing here, and the composite then quietly disagrees with the
  # standalone figure.
  qc("  <- %s  (edited %s)\n", path,
     format(file.info(path)$mtime, "%Y-%m-%d %H:%M"))
  # A panel writes its QC tables to CFG$outdir. Running it here with its own
  # default outdir overwrites the standalone run's tables with results computed
  # under DIFFERENT settings -- the composite passes min_n = 12 to the tissue
  # panel, so an unpooled standalone figure could end up beside a pooled counts
  # file. Each panel therefore writes into its own directory below the
  # composite's output, and standalone outputs are left alone.
  # NB: this override MUST be added before `src` is assembled -- an earlier
  # version appended it to `overrides` after the append() below, so it was
  # printed as applied but never executed, and panels kept writing into their
  # standalone output directories.
  panel_dir <- file.path(CFG$outdir, "panels", tag)
  dir.create(panel_dir, showWarnings = FALSE, recursive = TRUE)
  overrides <- c(overrides, sprintf('CFG$outdir <- "%s"', panel_dir))
  src <- readLines(path, warn = FALSE)
  i_open <- grep("^CFG <- list\\(", src)[1]
  if (is.na(i_open)) stop("no CFG block in ", path, call. = FALSE)
  i_close <- i_open - 1 + grep("^\\)\\s*$", src[i_open:length(src)])[1]
  src <- append(src, overrides, after = i_close)
  # drop the command-line override so it cannot fight with ours
  i_args <- grep("^args <- commandArgs", src)
  if (length(i_args)) src <- src[-c(i_args[1], i_args[1] + 1)]

  env <- new.env(parent = globalenv())
  # ggsave and the graphics devices are no-ops: we want the objects, not files
  env$ggsave <- function(...) invisible(NULL)
  env$png <- function(...) invisible(NULL)
  env$pdf <- function(...) invisible(NULL)
  env$cairo_pdf <- function(...) invisible(NULL)
  env$dev.off <- function(...) invisible(NULL)

  # Console output is written to a log rather than discarded, and warnings are
  # collected and re-raised rather than suppressed. Unresolved tissues, missing
  # dates and questionable taxonomy are exactly what must not go silent during
  # the final figure run.
  logfile <- file.path(panel_dir, "console.log")
  warns <- character(0)
  txt <- withCallingHandlers(
    capture.output(suppressMessages(
      eval(parse(text = paste(src, collapse = "\n")), envir = env))),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  writeLines(txt, logfile)
  attr(env, "log") <- logfile
  attr(env, "warnings") <- warns
  attr(env, "overrides") <- overrides
  if (length(warns)) {
    qc("     %d warning(s):\n", length(warns))
    for (w in unique(warns)) qc("       - %s\n", substr(w, 1, 110))
  }
  qc("     log -> %s\n", logfile)
  env
}

## ------------------------------------------------------------- restyling ----
# Panel scripts each carry their own theme, tuned for a standalone figure. In a
# composite they must share one type size and lose their outer padding.
# axis.text inherits a RELATIVE size from `text`, so setting text = 8 pt yields
# 6.4 pt axis labels. Each element that must clear 8 pt is therefore set
# explicitly -- but doing that unconditionally also un-blanks the axes of
# theme_void panels (the map and the Venn), so those get them blanked again.
compose_theme <- theme(
  text = element_text(size = CFG$base_pt),
  axis.text = element_text(size = CFG$base_pt, colour = "black"),
  axis.title = element_text(size = CFG$base_pt),
  legend.text = element_text(size = CFG$base_pt),
  legend.title = element_text(size = CFG$base_pt),
  strip.text = element_text(size = CFG$base_pt),
  plot.title = element_text(size = CFG$base_pt, face = "bold"),
  plot.subtitle = element_text(size = CFG$base_pt - 1),
  plot.margin = margin(1, 2, 1, 1),
  legend.margin = margin(0, 0, 0, 0),
  legend.box.spacing = unit(2, "pt"),
  legend.key.size = unit(7, "pt")
)

no_axes <- theme(axis.text = element_blank(), axis.title = element_blank(),
                 axis.ticks = element_blank(), panel.grid = element_blank(),
                 panel.border = element_blank())

# geom_text/geom_richtext sizes are in mm and are NOT touched by the theme, so
# they have to be raised individually or the figure fails the 8 pt rule.
bump_text_layers <- function(p, min_mm = CFG$min_text_mm) {
  if (is.null(p$layers)) return(p)
  p$layers <- lapply(p$layers, function(l) {
    g <- class(l$geom)[1]
    if (grepl("GeomText|GeomLabel|GeomRichText", g)) {
      s <- l$aes_params$size
      l$aes_params$size <- if (is.null(s)) min_mm else max(s, min_mm)
    }
    l
  })
  p
}

# Legend glyph size lives in each scale's guide as override.aes$size, not in the
# theme, so it survives every theme change. ggplot2 3.4 stores the guide as a
# plain list; 3.5 wraps it in a Guide ggproto with a params list -- both are
# handled so this does not silently do nothing after an upgrade.
scale_legend_glyphs <- function(p, factor) {
  hit <- 0
  for (sc in p$scales$scales) {
    g <- sc$guide
    if (inherits(g, "Guide") && !is.null(g$params$override.aes$size)) {
      g$params$override.aes$size <- g$params$override.aes$size * factor
      sc$guide <- g; hit <- hit + 1
    } else if (is.list(g) && !is.null(g$override.aes$size)) {
      g$override.aes$size <- g$override.aes$size * factor
      sc$guide <- g; hit <- hit + 1
    }
  }
  if (hit == 0)
    warning("no legend glyph sizes found to scale", call. = FALSE)
  p
}

# Point size lives in each layer's aes_params, entirely separate from the legend
# glyph size above. Scaling one without the other is what left panel F with
# small keys and large points.
scale_point_sizes <- function(p, factor, uniform = FALSE) {
  sz <- vapply(p$layers, function(l)
    if (grepl("GeomPoint", class(l$geom)[1]) && !is.null(l$aes_params$size))
      l$aes_params$size else NA_real_, numeric(1))
  if (!any(!is.na(sz))) {
    warning("no point layers with a fixed size in this panel", call. = FALSE)
    return(p)
  }
  target <- if (uniform) min(sz, na.rm = TRUE) else NA_real_
  for (i in seq_along(p$layers)) {
    if (is.na(sz[i])) next
    p$layers[[i]]$aes_params$size <-
      (if (uniform) target else sz[i]) * factor
  }
  p
}

restyle <- function(p, legend = NULL) {
  p <- bump_text_layers(p) + compose_theme
  if (!is.null(legend)) p <- p + theme(legend.position = legend)
  p
}

## ------------------------------------------------------------ run panels ----
cat("\n================ BUILDING PANELS ================\n")
t0 <- Sys.time()

qc("A  discovery by decade")
eA <- run_panel(CFG$panels[["A"]],
                sprintf('CFG$db_path <- "%s"', CFG$db_path), tag = "A")
pA <- eA$p


qc("B  geography")
eB <- run_panel(CFG$panels[["B"]],
                c(sprintf('CFG$db_path <- "%s"', CFG$db_path),
                  'CFG$sea_radius <- 4.6'), tag = "B")
pB <- eB$p


qc("C  host phylum")
eC <- run_panel(CFG$panels[["C"]],
                c(sprintf('CFG$api_file <- "%s"', CFG$api_file),
                  sprintf('CFG$min_n <- %d', CFG$phylum_min_n)), tag = "C")
pC <- eC$p


qc("D  environments (Venn)")
eD <- run_panel(CFG$panels[["D"]],
                c(sprintf('CFG$api_file <- "%s"', CFG$api_file),
                  sprintf('CFG$pred_file <- "%s"', CFG$pred_file),
                  'CFG$unit <- "species"'), tag = "D")
# the Venn is built inside render_unit(), so it is rebuilt here from the same
# objects rather than fished out of a closure
combo <- vapply(eD$sp_env, paste, character(1), collapse = "&")
tabD <- sort(table(combo), decreasing = TRUE)
countsD <- setNames(as.integer(tabD), names(tabD))
setsD <- eD$CFG$sets[eD$CFG$sets %in% unlist(strsplit(names(countsD), "&"))]
pD <- eD$venn_ggplot(countsD, setsD, eD$PAL)


qc("E  tissues")
eE <- run_panel(CFG$panels[["E"]],
                c(sprintf('CFG$db_path <- "%s"', CFG$db_path),
                  sprintf('CFG$min_n <- %d', CFG$tissue_min_n),
                  'CFG$run_selftest <- FALSE'), tag = "E")
pE <- eE$p


qc("F  spore size and shape")
eF <- run_panel(CFG$panels[["F"]],
                sprintf('CFG$db_path <- "%s"', CFG$db_path), tag = "F")
pF <- eF$p

qc("\npanels built in %.0f s\n",
   as.numeric(difftime(Sys.time(), t0, units = "secs")))

# What each panel actually counted. If C says "microsporidia species" when the
# standalone figure says "host species", the assembler is running an older copy
# of the panel script -- check the paths printed above.
cat("\nsettings each panel ran under (overrides applied by this script)\n")
for (nm in c("A", "B", "C", "D", "E", "F")) {
  env <- tryCatch(get(paste0("e", nm)), error = function(e) NULL)
  ov <- if (is.null(env)) NULL else attr(env, "overrides")
  if (length(ov))
    cat(sprintf("  %s  %s\n", nm,
                paste(sub("^CFG\\$", "", ov), collapse = "; ")))
}
cat("\nPanel outputs go to ", file.path(CFG$outdir, "panels"),
    "/<panel>/ so they cannot overwrite a standalone run's tables.\n", sep = "")

cat("\nwhat each panel counted\n")
for (nm in c("A", "C", "E", "F")) {
  env <- get(paste0("e", nm))
  pl <- get(paste0("p", nm))
  u <- if (!is.null(env$CFG$unit)) env$CFG$unit else "-"
  extra <- if (nm == "C" && !is.null(env$use_norm))
    sprintf("  phylum column = %s",
            if (isTRUE(env$use_norm)) "phylum_normalised" else "phylum (raw)") else ""
  cat(sprintf("  %s  unit = %-14s x axis = %-28s%s\n", nm, u,
              if (is.null(pl$labels$x)) "(none)" else pl$labels$x, extra))
}

## ------------------------------------------------- panel-specific tuning ----
# A carries the named/provisional legend for the whole figure; C and E use the
# same two colours and drop theirs, which buys three lines of height.
# A: "1850-1859" rotated is the tallest thing in the figure; "1850s" says the
# same in a third of the width. The legend goes above the panel rather than
# inside it, where it collided with the tall recent bars.
# A's per-bar counts overlap once the panel is a third of the page wide, and
# the y axis already carries the numbers. Dropped here only; the standalone
# panel keeps them. Delete the Filter() line to bring them back.
pA$layers <- Filter(function(l) !grepl("GeomText", class(l$geom)[1]), pA$layers)
pA <- restyle(pA, legend = "top") +
  labs(x = NULL) +
  # the panel writes en dashes under a UTF-8 locale and hyphens otherwise, so
  # both are matched: "1850-1859" / "1850\u20131859" -> "1850s"
  scale_x_discrete(labels = function(x) sub("[-\u2013].*$", "s", x)) +
  scale_fill_manual(values = c("Named species" = "#0072B2",
                               "Provisional species" = "#E69F00"),
                    labels = c("Named", "Provisional"), name = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.title = element_blank(),
        legend.justification = "left",
        legend.margin = margin(0, 0, 1, 0))

# The map is fixed-aspect: whatever its legend takes vertically comes straight
# out of the map. expand = FALSE and a cropped latitude range remove the empty
# ocean band left by dropping Antarctica, so the map fills its cell.
bb <- sf::st_bbox(eB$world_p)
pB <- restyle(pB, legend = "bottom") + no_axes +
  coord_sf(xlim = c(bb["xmin"], bb["xmax"]),
           ylim = c(bb["ymin"] * 1.02, bb["ymax"] * 1.02), expand = FALSE) +
  theme(legend.direction = "horizontal",
        legend.key.height = unit(4, "pt"), legend.key.width = unit(10, "pt"),
        legend.background = element_blank(),
        legend.box.spacing = unit(1, "pt"),
        legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(0, 0, 0, 0)) +
  guides(fill = guide_legend(nrow = 1, title.position = "left",
                             label.position = "bottom"))

# Both C and E carry the named+provisional breakdown after the total, in the
# stack colours. E gives up a column of width to pay for C's longer labels.
rich_breakdown <- function(total, parts, cols) {
  inner <- paste(sprintf("<span style='color:%s'>%d</span>", cols, parts),
                 collapse = "<span style='color:#7f7f7f'>+</span>")
  sprintf(paste0("<span style='color:#4d4d4d'>%d</span> ",
                 "<span style='color:#7f7f7f'>(</span>%s",
                 "<span style='color:#7f7f7f'>)</span>"), total, inner)
}

# facet_map: named vector mapping category -> facet level. A label layer with
# no facet variable is drawn in EVERY panel, which is what smeared panel C
# across all four kingdom strips the first time this was tried.
add_breakdown <- function(p, wide_counts, id_col, stack_cols, pal_use,
                          expand_mult, facet_map = NULL, facet_col = NULL) {
  p$layers <- Filter(function(l) !grepl("GeomRichText|GeomText",
                                        class(l$geom)[1]), p$layers)
  parts <- as.matrix(wide_counts[, stack_cols, drop = FALSE])
  parts[is.na(parts)] <- 0
  d <- data.frame(y = wide_counts[[id_col]], total = rowSums(parts))
  if (!is.null(facet_map) && !is.null(facet_col))
    d[[facet_col]] <- factor(unname(facet_map[as.character(d$y)]),
                             levels = levels(facet_map))
  d$lab <- vapply(seq_len(nrow(d)), function(i)
    rich_breakdown(d$total[i], as.integer(parts[i, ]), unname(pal_use[stack_cols])),
    character(1))
  d$plain <- sprintf("%d (%s)", d$total,
                     apply(parts, 1, function(r) paste(as.integer(r), collapse = "+")))
  p +
    (if (has_ggtext)
       ggtext::geom_richtext(data = d, aes(x = total, y = y, label = lab),
                             inherit.aes = FALSE, hjust = 0,
                             size = CFG$min_text_mm, fill = NA,
                             label.color = NA,
                             label.padding = unit(c(0, 0, 0, 1.2), "mm"))
     else
       geom_text(data = d, aes(x = total, y = y, label = plain),
                 inherit.aes = FALSE, hjust = -0.08, size = CFG$min_text_mm,
                 colour = "grey25")) +
    scale_x_continuous(expand = expansion(mult = c(0, expand_mult)),
                       breaks = scales::breaks_pretty(2))
}

pC <- restyle(pC, legend = "none") + labs(y = NULL)
wC <- reshape(eC$counts, idvar = "phylum_lab", timevar = "stack",
              direction = "wide")
names(wC) <- sub("^Freq\\.", "", names(wC))
wC$phylum_lab <- factor(as.character(wC$phylum_lab),
                        levels = levels(eC$counts$phylum_lab))
stackC <- levels(eC$counts$stack)
grpC <- eC$counts$group[match(levels(eC$counts$phylum_lab),
                              as.character(eC$counts$phylum_lab))]
grpC <- setNames(grpC, levels(eC$counts$phylum_lab))
pC <- add_breakdown(pC, wC, "phylum_lab", stackC,
                    setNames(c("#0072B2", "#E69F00"), stackC), 0.90,
                    facet_map = grpC, facet_col = "group") +
  theme(plot.margin = margin(9, 6, 1, 2))

pE <- restyle(pE, legend = "none") + labs(y = NULL)
wE <- reshape(as.data.frame(eE$bar), idvar = "tissue", timevar = "category",
              direction = "wide")
names(wE) <- sub("^n\\.", "", names(wE))
wE$tissue <- factor(as.character(wE$tissue), levels = levels(eE$bar$tissue))
stackE <- levels(eE$bar$category)
pE <- add_breakdown(pE, wE, "tissue", stackE,
                    setNames(c("#0072B2", "#E69F00"), stackE), 1.15) +
  # E sits against the right edge of the page, so the label headroom has to be
  # reserved in the MARGIN as well as in the scale, or the last few characters
  # fall off the canvas rather than off the panel
  theme(plot.margin = margin(9, 22, 1, 2))

# The Venn carries 0.26 of padding on each side for its outer labels, which is
# a fifth of the panel width spent on white. Cropped back in the composite.
# The region counts sit inside the ellipses and overflow the narrow lenses at
# 8 pt, so they alone are set one point smaller. The set names stay at 8 pt.
pD <- restyle(pD)
pD$layers <- lapply(pD$layers, function(l) {
  if (grepl("GeomText", class(l$geom)[1]) &&
      !identical(l$aes_params$fontface, "bold"))
    l$aes_params$size <- CFG$venn_pt * 25.4 / 72
  l
})
pD <- pD + no_axes +
  # tight enough to fill the cell, wide enough that the outer category labels
  # are not clipped by it
  coord_equal(xlim = c(-0.36, 1.36), ylim = c(0.13, 0.93), clip = "off") +
  theme(plot.margin = margin(1, 1, 1, 1))

# F's two legends stacked inside the panel covered every point above 4 um.
# Outside on the right they cost width but hide nothing.
pF <- scale_legend_glyphs(pF, CFG$spore_key_scale)
pF <- scale_point_sizes(pF, CFG$spore_key_scale,
                        uniform = CFG$spore_uniform_points)
pF <- restyle(pF, legend = "right") +
  theme(legend.justification = "top",
        legend.background = element_blank(),
        legend.key.size = unit(6, "pt"), legend.spacing.y = unit(0, "pt"),
        legend.spacing.x = unit(2, "pt"),
        legend.box = "vertical",
        legend.box.spacing = unit(2, "pt"),
        legend.title = element_text(size = CFG$base_pt, face = "bold"))

## ------------------------------------------------------------- assemble ----
# A|B across the top, the two long lists side by side in the middle, and the
# two square panels at the bottom. Keeping C and E in one row lets them share a
# height and stops either from stretching the figure on its own.
# patchwork aligns panel regions COLUMN-WISE across the whole grid, so the map
# was being pushed right to line up with the tissue panel's long category
# labels, and shrunk to match. free() releases a panel from that alignment.
# C needs the extra column: its labels are longer once the breakdown is on them
design <- "
AAABBBB
CCCCEEE
DDDFFFF
"
# Every panel is freed: none of them share an axis with another, so alignment
# only ever pushed a panel sideways to line up with a neighbour's long category
# labels -- which is what shrank the map and left a third of panel A empty.
# free() detaches a panel from patchwork's alignment, and with it from
# plot_annotation(tag_levels), so the tags are set on the panels themselves.
tag_it <- function(p, tag) p + labs(tag = tag)
pA <- tag_it(pA, "A"); pB <- tag_it(pB, "B"); pC <- tag_it(pC, "C")
pD <- tag_it(pD, "D"); pE <- tag_it(pE, "E"); pF <- tag_it(pF, "F")

fig <- free(pA) + free(pB) + free(pC) + free(pD) + free(pE) + free(pF) +
  plot_layout(design = design, heights = c(0.92, 1.20, 0.98)) &
  # "topleft" puts the tag in the plot MARGIN, so the margin has to exist --
  # with the default zero top margin the tags were drawn off-canvas
  theme(plot.tag = element_text(size = CFG$base_pt + 2, face = "bold"),
        plot.tag.position = "topleft",
        plot.margin = margin(9, 2, 1, 2))

## ============================================================================
## NUMBERS BEHIND THE FIGURE
## ============================================================================
# Every count a legend or a results sentence needs, taken from the SAME objects
# the panels plotted rather than recomputed here -- a number quoted in the text
# and a number drawn in the figure cannot then drift apart.
#
# Each panel is read defensively: if a panel script is renamed, restructured, or
# an older copy is run, the field reports NA instead of stopping the assembly.

`%||%` <- function(a, b) if (is.null(a)) b else a

n_uni <- function(x) if (is.null(x)) NA_integer_ else length(unique(x[nzchar(as.character(x)) & !is.na(x)]))
get_ <- function(env, ...) {
  for (nm in c(...)) if (!is.null(env[[nm]])) return(env[[nm]])
  NULL
}
fmt <- function(x) {
  if (length(x) == 0 || is.na(x)) return("NA")
  if (abs(x - round(x)) < 1e-9) format(round(x), big.mark = ",") else
    sprintf("%.2f", x)
}

stats_rows <- list()
add_stat <- function(panel, quantity, value, note = "") {
  v <- suppressWarnings(as.numeric(value))
  if (length(v) != 1 || is.infinite(v)) v <- NA_real_
  stats_rows[[length(stats_rows) + 1]] <<-
    data.frame(panel = panel, quantity = quantity, value = v,
               note = note, stringsAsFactors = FALSE)
}

## ---- A  discovery by decade ------------------------------------------------
plt <- get_(eA, "plt"); datA <- get_(eA, "dat")
if (!is.null(plt)) {
  add_stat("A", "microsporidia species plotted", n_uni(plt$species_name))
  # the panel labels these "Named species" / "Provisional species"; reading the
  # levels rather than hard-coding them keeps this working if they are renamed
  for (lv in levels(droplevels(factor(plt$category))))
    add_stat("A", paste0(tolower(sub(" species$", "", lv)), " species"),
             n_uni(plt$species_name[plt$category == lv]))
  add_stat("A", "decade bins drawn", n_uni(as.character(plt$bin)))
  add_stat("A", "earliest year", suppressWarnings(min(plt$year, na.rm = TRUE)))
  add_stat("A", "latest year", suppressWarnings(max(plt$year, na.rm = TRUE)))
  if (!is.null(datA))
    add_stat("A", "entries excluded, no year",
             n_uni(datA$species_name) - n_uni(plt$species_name))
}

## ---- B  geography ----------------------------------------------------------
longB <- get_(eB, "long"); cntB <- get_(eB, "counts")
seaB <- get_(eB, "sea_long"); seaC <- get_(eB, "sea_counts"); datB <- get_(eB, "dat")
if (!is.null(longB)) {
  add_stat("B", "microsporidia species plotted", n_uni(longB$species_name %||% longB$species))
  add_stat("B", "countries with >=1 species", n_uni(longB$iso3))
  if (!is.null(cntB) && nrow(cntB))
    add_stat("B", "species in the most-reported country", max(cntB$n),
             paste("in", cntB$iso3[which.max(cntB$n)]))
}
if (!is.null(seaB) && nrow(seaB)) {
  add_stat("B", "species in open water", n_uni(seaB$species_name),
           sprintf("sea_scope = %s", eB$CFG$sea_scope %||% "?"))
  add_stat("B", "named water bodies matched", n_uni(seaB$sea))
  add_stat("B", "circles drawn on the map",
           if (!is.null(seaC) && nrow(seaC)) nrow(seaC) else n_uni(seaB$basin),
           sprintf("sea_display = %s", eB$CFG$sea_display %||% "?"))
}
if (!is.null(eB$n_coastal))
  add_stat("B", "species with a coastal record also placed in a country",
           eB$n_coastal)
if (!is.null(datB) && !is.null(longB))
  add_stat("B", "entries excluded, no usable locality",
           n_uni(datB$species_name) - n_uni(longB$species_name %||% longB$species))

## ---- C  host phylum --------------------------------------------------------
longC <- get_(eC, "long"); dC <- get_(eC, "d")
if (!is.null(longC)) {
  add_stat("C", "bars count", NA,
           sprintf("unit = %s", eC$CFG$unit %||% "?"))
  add_stat("C", "host species in the panel", n_uni(dC$host_name))
  add_stat("C", "microsporidia species in the panel",
           n_uni(dC$microsporidia_species))
  add_stat("C", "entities plotted (the bar totals)", n_uni(longC$entity))
  add_stat("C", "host phyla drawn", n_uni(longC$phylum_lab))
  add_stat("C", "kingdom groups drawn", n_uni(as.character(longC$group)))
  if (!is.null(dC))
    add_stat("C", "host rows with no phylum", sum(!nzchar(dC$phylum_plot %||% dC$phylum)))
}

## ---- D  environments -------------------------------------------------------
spD <- get_(eD, "sp_env"); hoD <- get_(eD, "host_env")
if (!is.null(spD)) {
  add_stat("D", "microsporidia species plotted", length(spD))
  add_stat("D", "environments drawn", n_uni(unlist(spD)))
  add_stat("D", "species in one environment", sum(lengths(spD) == 1))
  add_stat("D", "species in >1 environment", sum(lengths(spD) > 1))
  for (e in sort(unique(unlist(spD))))
    add_stat("D", paste("species with a", e, "host"),
             sum(vapply(spD, function(v) e %in% v, logical(1))))
}
if (!is.null(hoD)) add_stat("D", "host species with an environment", length(hoD))

## ---- E  tissues ------------------------------------------------------------
longE <- get_(eE, "long"); datE <- get_(eE, "dat")
if (!is.null(longE)) {
  add_stat("E", "microsporidia species plotted", n_uni(longE$species))
  add_stat("E", "tissue categories drawn", n_uni(longE$tissue))
  per <- table(longE$species)
  add_stat("E", "species infecting one tissue", sum(per == 1))
  add_stat("E", "species infecting >1 tissue", sum(per > 1))
  if (!is.null(datE))
    add_stat("E", "entries excluded, no tissue assigned",
             n_uni(datE$species_name %||% datE$species) - n_uni(longE$species))
}

## ---- F  spore size and shape ----------------------------------------------
pltF <- get_(eF, "plt"); datF <- get_(eF, "dat")
if (!is.null(pltF)) {
  add_stat("F", "microsporidia species plotted", n_uni(pltF$species))
  add_stat("F", "spore measurements plotted", nrow(pltF))
  add_stat("F", "spore shape classes drawn", n_uni(as.character(pltF$shape)))
  if (!is.null(pltF$coil_bin)) {
    add_stat("F", "polar tubule coil bins drawn", n_uni(as.character(pltF$coil_bin)))
    add_stat("F", "species with coil data", n_uni(pltF$species[!is.na(pltF$coil_bin)]))
    add_stat("F", "species without coil data", n_uni(pltF$species[is.na(pltF$coil_bin)]))
  }
  add_stat("F", "shortest spore (um)", min(pltF$length, na.rm = TRUE))
  add_stat("F", "longest spore (um)", max(pltF$length, na.rm = TRUE))
  add_stat("F", "narrowest spore (um)", min(pltF$width, na.rm = TRUE))
  add_stat("F", "widest spore (um)", max(pltF$width, na.rm = TRUE))
  if (!is.null(datF))
    add_stat("F", "entries excluded, no length/width/shape",
             n_uni(datF$species) - n_uni(pltF$species))
}

## ---- how many species does the whole figure rest on? -----------------------
sp_sets <- list(
  A = if (!is.null(plt)) unique(plt$species_name),
  C = if (!is.null(dC)) unique(dC$microsporidia_species),
  D = if (!is.null(spD)) names(spD),
  E = if (!is.null(longE)) unique(longE$species),
  F = if (!is.null(pltF)) unique(pltF$species))
sp_sets <- sp_sets[!vapply(sp_sets, is.null, logical(1))]
if (!is.null(longB))
  sp_sets$B <- unique(longB$species_name %||% longB$species)
if (length(sp_sets)) {
  add_stat("all", "species in at least one panel",
           length(unique(unlist(sp_sets))))
  add_stat("all", "species in every panel",
           length(Reduce(intersect, sp_sets)))
}

## ---- report ----------------------------------------------------------------
stats <- do.call(rbind, stats_rows)
cat("\n================ NUMBERS BEHIND THE FIGURE ================\n")
for (pn in unique(stats$panel)) {
  s <- stats[stats$panel == pn, ]
  cat(sprintf("\n%s\n", if (pn == "all") "ACROSS PANELS" else paste("Panel", pn)))
  for (i in seq_len(nrow(s)))
    cat(sprintf("  %-42s %10s  %s\n", s$quantity[i], fmt(s$value[i]), s$note[i]))
}
cat("\nSpecies counts differ between panels because each attribute is reported\n",
    "for a different subset of entries; 'species in every panel' is the\n",
    "intersection, and is the number to quote for any cross-panel claim.\n",
    sep = "")

write.table(stats, file.path(CFG$outdir, "figure_numbers.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
qc("\nwritten %s\n", file.path(CFG$outdir, "figure_numbers.tsv"))

out_pdf <- file.path(CFG$outdir, "figure_composite.pdf")
out_png <- file.path(CFG$outdir, "figure_composite.png")
ggsave(out_pdf, fig, width = CFG$fig_width, height = CFG$fig_height,
       units = "in", device = if (capabilities("cairo")) cairo_pdf else pdf)
ggsave(out_png, fig, width = CFG$fig_width, height = CFG$fig_height,
       units = "in", dpi = CFG$dpi)

qc("\nwritten\n  %s\n  %s  (%.1f x %.1f in, %d dpi = %d x %d px)\n",
   out_pdf, out_png, CFG$fig_width, CFG$fig_height, CFG$dpi,
   round(CFG$fig_width * CFG$dpi), round(CFG$fig_height * CFG$dpi))
qc("smallest type: %g pt (theme) / %.2f mm = %g pt (text layers)\n",
   CFG$base_pt, CFG$min_text_mm, CFG$min_text_mm * 72 / 25.4)
