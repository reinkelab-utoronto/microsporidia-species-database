#!/usr/bin/env Rscript
###############################################################################
# Prune the microsporidia 18S ML tree to the tips shown in Fig 2.
#
# The maximum-likelihood tree is inferred from every centroid that entered the
# alignment. Figure 2 shows only clusters with at least one reported host, so
# the tree is pruned to those tips here. This is a DISPLAY step, not part of the
# inference -- kept as its own script (rather than a throwaway snippet) so the
# 680 -> N tip accounting is reproducible and auditable, which is exactly the
# step a reviewer will ask about.
#
# The tips to keep come from 03_CLUSTER_HOST_ATTRIBUTES.tsv -- the SAME file the
# tree figure reads for its rings and symbols -- so the pruning and the
# annotations can never disagree about which clusters have a host.
#
# Three things this does that the original snippet did not:
#   * coerces the host count to a NUMBER (the old "n > 0" compared a character
#     column, a lexicographic accident that only happened to work)
#   * matches tree tips to accessions the same way plot_18S_tree.R does
#     (vsearch "NNN|" prefix and version suffix tolerated), so a spelling
#     difference cannot silently drop a tip
#   * ALWAYS retains the outgroup, even with no host record -- otherwise
#     plot_18S_tree.R's root(outgroup = ...) stops with "OUTGROUP not found"
#
# Usage:
#   Rscript prune_tree.R
#   Rscript prune_tree.R <full.treefile> <attributes.tsv> <out.treefile>
###############################################################################

CFG <- list(
  in_tree   = "mito_muscle01.treefile",
  attr_file = "cluster_output/03_CLUSTER_HOST_ATTRIBUTES.tsv",
  # No version number in the output name on purpose: point plot_18S_tree.R's
  # TREE_FILE at THIS file. Numbered copies (…_2_pruned, …_8_pruned) are how the
  # tree and the pruning drifted apart in the first place.
  out_tree  = "CENTROIDS_WITH_HOSTS_pruned.treefile",

  acc_col        = "centroid_accession",  # tip labels are these accessions
  host_count_col = "n_distinct_hosts",    # keep a cluster when this is > 0
  # tips ALWAYS kept even with no host: the outgroup must survive pruning
  always_keep = "MF278562.1"              # Mitosporidium daphniae
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) CFG$in_tree   <- args[1]
if (length(args) >= 2) CFG$attr_file <- args[2]
if (length(args) >= 3) CFG$out_tree  <- args[3]

if (!requireNamespace("ape", quietly = TRUE)) install.packages("ape")
suppressPackageStartupMessages(library(ape))
qc <- function(...) cat(sprintf(...), sep = "")

## ------------------------------------------------------------------ inputs ---
if (!file.exists(CFG$in_tree))
  stop("tree not found: ", CFG$in_tree, call. = FALSE)
if (!file.exists(CFG$attr_file))
  stop("attribute table not found: ", CFG$attr_file, call. = FALSE)

full <- read.tree(CFG$in_tree)
qc("full tree      : %d tips  (%s)\n", Ntip(full), CFG$in_tree)

attrs <- read.delim(CFG$attr_file, stringsAsFactors = FALSE, na.strings = "",
                    quote = "\"")
for (col in c(CFG$acc_col, CFG$host_count_col))
  if (!col %in% names(attrs))
    stop("column '", col, "' not in ", CFG$attr_file,
         "\n  found: ", paste(names(attrs), collapse = ", "), call. = FALSE)

## ------------------------------------------------------- which tips to keep ---
# host count as an integer, blanks / NA treated as zero
n_host <- suppressWarnings(as.integer(attrs[[CFG$host_count_col]]))
n_host[is.na(n_host)] <- 0L
keep_acc  <- unique(attrs[[CFG$acc_col]][n_host > 0])
qc("attribute table: %d rows, %d cluster(s) with >=1 host\n",
   nrow(attrs), length(keep_acc))

# Tolerant accession matching, as in plot_18S_tree.R: strip a vsearch "NNN|"
# prefix and the version suffix before comparing.
strip_pre <- function(x) sub("^[0-9]+\\|", "", x)
base_acc  <- function(x) sub("\\.[0-9]+$", "", strip_pre(x))
keep_base <- base_acc(keep_acc)

tip      <- full$tip.label
has_host <- tip %in% keep_acc | base_acc(tip) %in% keep_base
protect  <- tip %in% CFG$always_keep | base_acc(tip) %in% base_acc(CFG$always_keep)
keep_tip <- has_host | protect

if (any(protect & !has_host))
  qc("outgroup kept despite no host record: %s\n",
     paste(tip[protect & !has_host], collapse = ", "))
if (!any(protect))
  qc("!! WARNING: none of the always-keep tips (%s) is in the tree; ",
     paste(CFG$always_keep, collapse = ", "),
     "plot_18S_tree.R will fail to root.\n")
if (!sum(keep_tip))
  stop("no tips left to keep -- do the accession columns match the tip labels?",
       call. = FALSE)

# 680-vs-N accounting, made explicit: host-bearing accessions not found in the
# tree (and therefore silently absent from the figure) are reported, not hidden.
absent <- keep_acc[!(keep_base %in% base_acc(tip))]
if (length(absent))
  qc("!! %d host-bearing accession(s) in the table are NOT in the tree: %s%s\n",
     length(absent), paste(head(absent, 6), collapse = ", "),
     if (length(absent) > 6) ", ..." else "")

## ------------------------------------------------------------------- prune ---
pruned <- drop.tip(full, tip[!keep_tip])
qc("pruned         : %d -> %d tips (%d dropped: no reported host)\n",
   Ntip(full), Ntip(pruned), Ntip(full) - Ntip(pruned))

if (!all(CFG$always_keep %in% pruned$tip.label |
         base_acc(CFG$always_keep) %in% base_acc(pruned$tip.label)))
  qc("!! WARNING: an always-keep tip is missing from the pruned tree\n")

write.tree(pruned, CFG$out_tree)
qc("written        : %s\n", CFG$out_tree)
qc("  -> set TREE_FILE in plot_18S_tree.R to this exact file\n")
