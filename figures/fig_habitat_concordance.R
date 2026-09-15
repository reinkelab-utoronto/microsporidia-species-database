#!/usr/bin/env Rscript
###############################################################################
# Supplemental figure -- concordance between database-derived and inferred
# host habitats
#
# The taxonomic databases (GBIF, WoRMS, Catalogue of Life) have no habitat for
# a large block of hosts, mostly terrestrial insects. Those hosts are assigned
# from a genus-level habitat table instead. This figure makes the case for
# doing so, in two panels:
#
#   A  how many hosts each method could place, and how many both could
#   B  where both placed a host, how often they agree, per habitat
#
# WHAT THIS IS NOT. Neither method is a gold standard: the databases are
# authoritative where they have data and silent where they do not, and the
# genus table always answers. So this is CONCORDANCE between two methods, not
# accuracy, and the word "accuracy" is deliberately absent from the axis
# labels. It is also measured on the hosts both methods could place, which are
# the better-known ones -- a point for the discussion, not for the figure.
#
# INPUTS   host_taxonomy_output/02_host_lookup.tsv      (databases)
#          predicted_environment/02_host_predicted.tsv  (genus table)
#          Either file may be the per-pair 01_* table instead; the host name
#          and environment columns are what matter.
#
# Usage:  Rscript fig_habitat_concordance.R
###############################################################################

CFG <- list(
  api_file  = "host_taxonomy_output/02_host_lookup.tsv",
  pred_file = "predicted_environment/02_host_predicted.tsv",

  # which habitats to show, and the order they appear in
  envs = c("freshwater", "terrestrial", "marine", "brackish"),

  outdir     = "concordance_output",
  fig_width  = 7.4,
  fig_height = 3.4
)

for (p in c("ggplot2")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages(library(ggplot2))
has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
dir.create(CFG$outdir, showWarnings = FALSE, recursive = TRUE)
options(width = 200)
qc <- function(...) cat(sprintf(...), sep = "")

## ---------------------------------------------------------------- inputs ---
read_hosts <- function(path, env_cols, label) {
  if (!file.exists(path))
    stop("required input file not found: ", path,
         "\n  produce it with host_taxonomy_environment.py / ",
         "predict_host_environment.py", call. = FALSE)
  # quote = "" because a stray double quote in a host cell makes read.delim
  # swallow every following line until it finds a partner
  n_lines <- length(readLines(path, warn = FALSE)) - 1L
  d <- read.delim(path, stringsAsFactors = FALSE, quote = "", na.strings = "",
                  comment.char = "")
  if (nrow(d) != n_lines)
    stop(sprintf("%s did not parse cleanly: %d data lines, %d rows read",
                 path, n_lines, nrow(d)), call. = FALSE)
  hit <- env_cols[env_cols %in% names(d)][1]
  if (is.na(hit))
    stop("no environment column in ", path, "\n  looked for: ",
         paste(env_cols, collapse = ", "), call. = FALSE)
  if (!"host_name" %in% names(d))
    stop("no host_name column in ", path, call. = FALSE)
  d[[hit]][is.na(d[[hit]])] <- ""
  out <- data.frame(host = trimws(d$host_name), env = d[[hit]],
                    stringsAsFactors = FALSE)
  out <- out[nzchar(out$host), ]
  # one row per host: a per-pair table lists a host once per parasite
  out <- out[!duplicated(out$host), ]
  qc("%-12s %-46s %5d host(s), env column: %s\n", label, path, nrow(out), hit)
  out
}

cat("\n============ DATABASE vs INFERRED HOST HABITAT ============\n")
api <- read_hosts(CFG$api_file, c("environment", "host_derived_environment"),
                  "databases")
prd <- read_hosts(CFG$pred_file, c("environment_all", "environment_primary"),
                  "genus table")

split_env <- function(x) lapply(strsplit(x, ";\\s*"), function(v) {
  v <- trimws(v); sort(unique(v[nzchar(v) & v %in% CFG$envs]))
})
api$set <- split_env(api$env)
prd$set <- split_env(prd$env)

hosts <- union(api$host, prd$host)
a <- setNames(api$set, api$host)[hosts]
p <- setNames(prd$set, prd$host)[hosts]
has_a <- vapply(a, function(v) length(v) > 0 && !is.null(v), logical(1))
has_p <- vapply(p, function(v) length(v) > 0 && !is.null(v), logical(1))

## ------------------------------------------------------- panel A: coverage --
cov <- data.frame(
  status = factor(c("Both methods", "Genus table only", "Databases only"),
                  levels = c("Databases only", "Genus table only",
                             "Both methods")),
  n = c(sum(has_a & has_p), sum(!has_a & has_p), sum(has_a & !has_p)))

qc("\nhost species placed by\n")
for (i in seq_len(nrow(cov)))
  qc("  %-18s %5d\n", cov$status[i], cov$n[i])
qc("  %-18s %5d\n", "neither", sum(!has_a & !has_p))

## ------------------------------------------ panel B: agreement per habitat --
both <- which(has_a & has_p)
agree_any <- mean(vapply(both, function(i) length(intersect(a[[i]], p[[i]])) > 0,
                         logical(1)))
qc("\nhosts placed by both                     : %d\n", length(both))
qc("agreeing on at least one habitat         : %.1f%%\n", 100 * agree_any)

per <- do.call(rbind, lapply(CFG$envs, function(e) {
  in_a <- both[vapply(both, function(i) e %in% a[[i]], logical(1))]
  if (!length(in_a)) return(NULL)
  same <- sum(vapply(in_a, function(i) e %in% p[[i]], logical(1)))
  data.frame(env = e, n = length(in_a), agree = same / length(in_a),
             stringsAsFactors = FALSE)
}))
qc("\nof the hosts the DATABASES place in each habitat, the share the genus\n")
qc("table places there too:\n")
for (i in seq_len(nrow(per)))
  qc("  %-12s %5d hosts   %5.1f%%\n", per$env[i], per$n[i], 100 * per$agree[i])

write.table(per, file.path(CFG$outdir, "concordance_by_habitat.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
disagree <- data.frame(
  host = hosts[both][vapply(both, function(i)
    length(intersect(a[[i]], p[[i]])) == 0, logical(1))],
  stringsAsFactors = FALSE)
if (nrow(disagree)) {
  disagree$databases <- vapply(a[disagree$host], paste, character(1), collapse = "; ")
  disagree$genus_table <- vapply(p[disagree$host], paste, character(1), collapse = "; ")
  write.table(disagree, file.path(CFG$outdir, "hosts_that_disagree.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  qc("\n%d host(s) share no habitat -> hosts_that_disagree.tsv\n", nrow(disagree))
}

## ------------------------------------------------------------------ plot ----
# Okabe-Ito, as in the main figure
pal_cov <- c("Both methods" = "#0072B2", "Genus table only" = "#009E73",
             "Databases only" = "#999999")
pal_env <- c(freshwater = "#56B4E9", brackish = "#009E73",
             marine = "#0072B2", terrestrial = "#E69F00")

base_theme <- theme_classic(base_size = 11) +
  theme(
    axis.text.y = element_text(colour = "black"),
    axis.text.x = element_text(colour = "black"),
    axis.title.x = element_text(margin = margin(t = 8), size = 10),
    legend.position = "none",
    plot.title = element_text(size = 11, face = "bold"),
    plot.margin = margin(6, 14, 4, 6)
  )

pA <- ggplot(cov, aes(x = n, y = status, fill = status)) +
  geom_col(width = 0.68, colour = "grey20", linewidth = 0.2) +
  geom_text(aes(label = n), hjust = -0.18, size = 2.9, colour = "grey20") +
  scale_fill_manual(values = pal_cov) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(x = "Host species", y = NULL) +
  base_theme

# n goes under the habitat name rather than after the percentage, where it ran
# off the panel at this width
per$label <- sprintf("%s\n(n = %d)", per$env, per$n)
per$env <- factor(per$env, levels = rev(CFG$envs))
per$label <- factor(per$label, levels = per$label[order(match(per$env,
                                                             rev(CFG$envs)))])
pB <- ggplot(per, aes(x = agree, y = label, fill = env)) +
  geom_col(width = 0.68, colour = "grey20", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * agree)),
            hjust = -0.18, size = 2.9, colour = "grey20") +
  scale_fill_manual(values = pal_env) +
  scale_x_continuous(labels = function(x) paste0(100 * x, "%"),
                     limits = c(0, 1), breaks = c(0, 0.5, 1),
                     expand = expansion(mult = c(0, 0.14))) +
  labs(x = "Host species assigned the same habitat (%)", y = NULL) +
  base_theme

if (has_patchwork) {
  library(patchwork)
  fig <- pA + pB + plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 12, face = "bold"))
} else {
  message("patchwork not installed; writing the two panels separately")
  fig <- NULL
}

save_fig <- function(plot, name, w, h) {
  ggsave(file.path(CFG$outdir, paste0(name, ".pdf")), plot, width = w, height = h,
         device = if (capabilities("cairo")) cairo_pdf else pdf)
  ggsave(file.path(CFG$outdir, paste0(name, ".png")), plot, width = w, height = h,
         dpi = 600)
}
if (!is.null(fig)) {
  save_fig(fig, "fig_habitat_concordance", CFG$fig_width, CFG$fig_height)
} else {
  save_fig(pA, "fig_habitat_coverage", 3.6, CFG$fig_height)
  save_fig(pB, "fig_habitat_agreement", 3.6, CFG$fig_height)
}

qc("\nwritten to %s/\n  fig_habitat_concordance.pdf / .png\n  concordance_by_habitat.tsv\n  hosts_that_disagree.tsv\n",
   CFG$outdir)
qc("\nFor the legend: %.0f%% of the %d host species placed by both methods were\n",
   100 * agree_any, length(both))
qc("assigned at least one habitat in common.\n")
