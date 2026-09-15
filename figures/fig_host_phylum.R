#!/usr/bin/env Rscript
###############################################################################
# Figure -- microsporidia species by host phylum
#
# Stacked bars of the number of microsporidia species recorded from each host
# phylum, split by named vs provisional, matching the discovery-by-decade and
# tissue panels.
#
# INPUT  host_taxonomy_output/01_species_host_long.tsv, written by
#        host_taxonomy_environment.py. That is the only table carrying a host
#        lineage; the offline predictor gives a coarse host_group but no
#        phylum, so it cannot substitute here. A missing file is a hard error.
#
# THREE THINGS THE PREVIOUS VERSION GOT WRONG, WORTH KNOWING ABOUT
#
#   1. A LOG AXIS AND STACKED BARS ARE INCOMPATIBLE. On a log scale the two
#      segments of a stack do not add up to the total -- the visual length of
#      "named + provisional" is not the length of their sum. The old figure got
#      away with a log axis because its bars were DODGED. This panel is stacked,
#      so the axis is linear. Set CFG$log_panel = TRUE to get a companion
#      dodged log panel for the long tail, which is legitimate.
#
#   2. Phyla were hardcoded in two lists and anything outside them was silently
#      dropped. Here the phyla come from the data, the kingdom column decides
#      the Metazoa/Protista grouping, and anything unplaced is REPORTED.
#
#   3. n_distinct(Host_formatted) counted HOSTS. Named vs provisional is a
#      property of the microsporidian, so the unit here is microsporidia
#      species (CFG$unit = "microsporidia"). CFG$unit = "host" counts host
#      species instead, stacked by whether their parasites are named,
#      provisional, or both -- a real partition, so nothing is double counted.
#
# Usage:  Rscript fig_host_phylum.R [path/to/01_species_host_long.tsv]
###############################################################################

CFG <- list(
  api_file = "host_taxonomy_output/01_species_host_long.tsv",

  # WHICH HOST RECORDS TO COUNT.
  # host_taxonomy_environment.py harvests both "Natural Host(s)" and
  # "Experimental Host(s)" by default, so 01_species_host_long.tsv contains
  # laboratory infections as well as natural ones. An experimental infection
  # shows what a parasite CAN infect, not what it does infect, so counting it
  # overstates realised host range -- and the phylogeny figure already excludes
  # them, so including them here would make the two figures disagree.
  #   "natural"      Natural Host(s) only (default)
  #   "experimental" Experimental Host(s) only
  #   "all"          both, i.e. the previous behaviour
  host_role_filter = "natural",

  # "host"          one bar segment per HOST species (default)
  # "microsporidia"  one per microsporidian instead
  unit = "host",

  # WHAT ONE BAR COUNTS.
  #   "species" one microsporidian counted ONCE per phylum, however many hosts
  #             in that phylum it infects (default; the y label says species)
  #   "records" every species-host row counted separately
  # The two differ a lot: a parasite with three nematode hosts is three rows
  # and one species, so the row count in the input table is always the larger
  # number. Both are printed per phylum either way.
  count_unit = "species",
  subset = "all",             # "all" | "named" | "provisional" (unit=microsporidia)

  min_n = 1,                  # phyla with fewer than this are pooled as "Other"
  log_panel = FALSE,          # add a dodged log-scale companion panel
  show_extant = FALSE,        # add a panel comparing with described diversity
  # TRUE removes the unresolved bar AND the unresolved kingdom strip. The
  # counts are still printed and written to the QC files.
  drop_unresolved = TRUE,
  # host_taxonomy_environment.py marks every row it is unsure about in a
  # needs_review column. "report" counts them in the console, "exclude" also
  # drops them from the bars so nothing unreviewed reaches the figure.
  flagged = "report",         # "report" | "exclude" | "ignore"
  # TRUE drops the "(... hyperhost)" rows of a hyperparasite chain, so one
  # microsporidian is not counted in two phyla. Works on an existing table --
  # the annotation is preserved in host_raw, so the lookup scripts need not be
  # re-run.
  drop_hyperhosts = TRUE,

  # WHICH PHYLUM COLUMN TO PLOT.
  #   "auto"   phylum_normalised where the table has it, else phylum (default)
  #   "raw"    the phylum exactly as the winning source gave it
  # The sources do not agree on phylum names: GBIF files gregarines under
  # Chromista / MYZOZOA while NCBI and COL call them APICOMPLEXA, so plotting
  # the raw column splits one real group across two bars and undercounts both.
  # host_taxonomy_environment.py reconciles them into phylum_normalised.
  phylum_col = "auto",    # TRUE hides the "(no phylum)" bar; it is still reported

  outdir = "phylum_output",
  fig_width = 7.6,
  fig_height = 6.0
)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) CFG$api_file <- args[1]

for (p in c("ggplot2")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages(library(ggplot2))
# ggtext lets one label carry several colours, so the breakdown after the total
# can be tinted to match the stack segments. Without it the same label is drawn
# in plain grey -- the numbers are identical either way.
if (!requireNamespace("ggtext", quietly = TRUE))
  try(install.packages("ggtext"), silent = TRUE)
has_ggtext <- requireNamespace("ggtext", quietly = TRUE)
has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
dir.create(CFG$outdir, showWarnings = FALSE, recursive = TRUE)
options(width = 200)
qc <- function(...) cat(sprintf(...), sep = "")
`%||%` <- function(a, b) if (length(a) == 0 || is.na(a[1])) b else a

## ---------------------------------------------------------------- inputs ----
if (!file.exists(CFG$api_file))
  stop("required input file not found: ", CFG$api_file,
       "\n\nProduce it with:\n",
       "  python host_taxonomy_environment.py --db <database>.xlsx --email <you>\n",
       "\nThe offline predictor's table cannot be used here: it has no phylum.",
       call. = FALSE)

cat("\n================ MICROSPORIDIA BY HOST PHYLUM ================\n")
qc("INPUT: %s  (modified %s)\n", CFG$api_file,
   format(file.info(CFG$api_file)$mtime, "%Y-%m-%d %H:%M"))

# TSVs are read with quote = "" -- a single stray double quote in a host cell
# (the database has several, e.g. 'Lecudina (Ophioidina) sp." in ...') makes
# read.delim swallow every following line until it finds a closing quote, which
# silently merges hundreds of rows into one and drops most phyla from the
# figure. The row count is checked against the file's own line count so this
# can never happen quietly again.
n_lines <- length(readLines(CFG$api_file, warn = FALSE)) - 1L
d <- read.delim(CFG$api_file, stringsAsFactors = FALSE, quote = "",
                na.strings = "", comment.char = "")
if (nrow(d) != n_lines) {
  d2 <- read.delim(CFG$api_file, stringsAsFactors = FALSE, quote = "\"",
                   na.strings = "")
  stop(sprintf(paste0("the table did not parse cleanly: %d data lines in the ",
                      "file but %d rows read (%d with quoting on).\n",
                      "  A field almost certainly contains an unbalanced quote ",
                      "or an embedded tab.\n",
                      "  Re-run host_taxonomy_environment.py to rewrite it."),
               n_lines, nrow(d), nrow(d2)), call. = FALSE)
}
need_cols <- c("microsporidia_species", "microsporidia_category", "host_name",
               "phylum", "kingdom")
missing_cols <- setdiff(need_cols, names(d))
if (length(missing_cols))
  stop("columns missing from ", CFG$api_file, ": ",
       paste(missing_cols, collapse = ", "), call. = FALSE)
for (v in need_cols) d[[v]][is.na(d[[v]])] <- ""
d <- d[nzchar(d$microsporidia_species), ]

# Fragments that were descriptions rather than names parse to an empty
# host_name (e.g. 'A gregarine, "form en comete," in ...'). They are real
# species-host records but they are not host SPECIES, so they are dropped when
# the unit is hosts, and reported either way.
n_blank_host <- sum(!nzchar(trimws(d$host_name)))
if (n_blank_host)
  qc("rows with no host name       : %d%s\n", n_blank_host,
     if (identical(CFG$unit, "host"))
       "  (excluded: not a host species)" else "")
if (identical(CFG$unit, "host")) d <- d[nzchar(trimws(d$host_name)), ]
qc("rows (species x host pairs) : %d\n", nrow(d))


## ------------------------------------------- hyperparasite chain rows -------
# The database records a whole parasite chain in one cell, with the roles
# annotated:
#   "Selenidium pygospionis (direct host); Pygospio elegans (polychaete hyperhost)"
#   "Enterocystis rhithrogenae (hyperparasitic); Rhithrogena semicolorata"
# A semicolon means separate hosts, so both fragments become rows and one
# microsporidian is counted in two phyla at once. The annotation says which
# fragment is the microsporidian's host, and it survives in host_raw -- so the
# figure can drop the wrong rows without re-running the lookup scripts.
HYPERHOST_RX <- "hyper[- ]?host|host of the host|carrier host"
PARASITE_RX  <- "hyperparasit|archigregarine|\\bgregarine\\b"
DIRECT_RX    <- "direct host|primary host|true host|actual host"

drop_hyperhost_rows <- function(d, label = "rows") {
  if (!"host_raw" %in% names(d)) {
    qc("note: no host_raw column; hyperparasite rows cannot be identified\n")
    return(d)
  }
  raw <- ifelse(is.na(d$host_raw), "", as.character(d$host_raw))
  cellkey <- if ("cell" %in% names(d) && any(nzchar(as.character(d$cell))))
    paste(d$microsporidia_species, d$cell) else
    paste(d$microsporidia_species, d$host_column)

  is_hyper  <- grepl(HYPERHOST_RX, raw, ignore.case = TRUE)
  is_parasi <- grepl(PARASITE_RX, raw, ignore.case = TRUE)
  is_direct <- grepl(DIRECT_RX, raw, ignore.case = TRUE)

  # rule 1: an explicit hyperhost, where the same cell names a direct host
  cells_with_direct <- unique(cellkey[is_direct | is_parasi])
  drop1 <- is_hyper & cellkey %in% cells_with_direct
  # rule 2: an unannotated companion in a cell whose other fragment is marked
  # as the parasite itself
  cells_with_parasite <- unique(cellkey[is_parasi])
  drop2 <- !is_parasi & !is_direct & !is_hyper & cellkey %in% cells_with_parasite
  drop <- drop1 | drop2

  if (any(drop)) {
    qc("\nhyperparasite chains: %d %s dropped (the host's host, not the "
       , sum(drop), label)
    qc("microsporidian's)\n")
    show <- unique(data.frame(species = d$microsporidia_species[drop],
                              dropped = raw[drop],
                              phylum = if ("phylum_plot" %in% names(d))
                                d$phylum_plot[drop] else d$phylum[drop],
                              stringsAsFactors = FALSE))
    print(utils::head(show, 12), row.names = FALSE, right = FALSE)
    write.table(show, file.path(CFG$outdir, "QC_hyperhost_rows_dropped.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
  } else {
    qc("hyperparasite chains: none found to drop\n")
  }
  d[!drop, , drop = FALSE]
}

## ------------------------------------------- natural vs experimental -------
## Applied before any counting, so every downstream number -- per phylum, per
## species, the QC files -- reflects the same set of records.
if (!identical(CFG$host_role_filter, "all")) {
  if (!"host_column" %in% names(d)) {
    qc("note: no host_column in the input; cannot separate natural from\n")
    qc("      experimental hosts, so ALL rows are counted\n")
  } else {
    hc <- ifelse(is.na(d$host_column), "", as.character(d$host_column))
    want <- if (identical(CFG$host_role_filter, "natural")) "natural" else "experimental"
    is_nat <- grepl("natural", hc, ignore.case = TRUE)
    is_exp <- grepl("experimental|lab", hc, ignore.case = TRUE)
    keep <- if (want == "natural") is_nat else is_exp
    ## a row from an input file rather than a database column has no role;
    ## keep it, since dropping it would silently discard user-supplied hosts
    keep <- keep | (!is_nat & !is_exp)
    qc("\nhost records by source column:\n")
    for (k in sort(unique(hc))) qc("  %-28s %d\n", k, sum(hc == k))
    qc("host_role_filter = %s -> keeping %d of %d row(s)\n",
       CFG$host_role_filter, sum(keep), nrow(d))
    if (!any(keep))
      stop("host_role_filter removed every row; check the host_column values")
    d <- d[keep, , drop = FALSE]
  }
}

## ------------------------------------------------- flagged rows -------------
if (CFG$drop_hyperhosts) d <- drop_hyperhost_rows(d, "row(s)")

if ("needs_review" %in% names(d)) {
  d$needs_review <- ifelse(is.na(d$needs_review), "", as.character(d$needs_review))
  fl <- d[nzchar(d$needs_review), ]
  cat("rows flagged needs_review    : ", nrow(fl), " (",
      sprintf("%.1f", 100 * nrow(fl) / nrow(d)), "%) [flagged = ",
      CFG$flagged, "]\n", sep = "")
  if (nrow(fl) > 0 && !identical(CFG$flagged, "ignore")) {
    rs <- unlist(strsplit(fl$needs_review, "; "))
    tb <- sort(table(rs[nzchar(rs)]), decreasing = TRUE)
    for (i in seq_along(tb)) qc("  %6d  %s\n", as.integer(tb)[i], names(tb)[i])
    write.table(unique(fl[, c("microsporidia_species", "host_name", "phylum",
                              "phylum_agreement", "phylum_dispute",
                              "needs_review")]),
                file.path(CFG$outdir, "QC_flagged_rows.tsv"), sep = "\t",
                row.names = FALSE, quote = FALSE)
  }
  if (identical(CFG$flagged, "exclude")) {
    n0 <- nrow(d)
    d <- d[!nzchar(d$needs_review), ]
    qc("  -> %d row(s) excluded from the figure\n", n0 - nrow(d))
  }
} else {
  qc("note: no needs_review column; this table predates the cross-source check\n")
}

## -------------------------------------------------------------- kingdoms ----
# GBIF says Animalia, NCBI says Metazoa, COL and WoRMS say other things again.
# The grouping is derived from whatever the lineage gives rather than from a
# hardcoded phylum list, so a phylum new to the database still lands somewhere.
# When a source gives a phylum but no kingdom -- COL frequently does for
# protists -- the kingdom-based grouping sent a perfectly well resolved host to
# the "(unresolved)" strip. The phylum decides the group whenever the kingdom
# cannot.
PHYLUM_GROUP <- c(
  apicomplexa = "Protista", ciliophora = "Protista", amoebozoa = "Protista",
  cercozoa = "Protista", myzozoa = "Protista", miozoa = "Protista",
  euglenozoa = "Protista", percolozoa = "Protista", foraminifera = "Protista",
  heliozoa = "Protista", bigyra = "Protista", ochrophyta = "Protista",
  arthropoda = "Metazoa", chordata = "Metazoa", nematoda = "Metazoa",
  annelida = "Metazoa", mollusca = "Metazoa", platyhelminthes = "Metazoa",
  rotifera = "Metazoa", bryozoa = "Metazoa", cnidaria = "Metazoa",
  nemertea = "Metazoa", acanthocephala = "Metazoa", porifera = "Metazoa",
  gastrotricha = "Metazoa", phoronida = "Metazoa", echinodermata = "Metazoa",
  tardigrada = "Metazoa", nematomorpha = "Metazoa", dicyemida = "Metazoa",
  cephalorhyncha = "Metazoa", entoprocta = "Metazoa", chaetognatha = "Metazoa",
  ascomycota = "Fungi", basidiomycota = "Fungi", microsporidia = "Fungi",
  chytridiomycota = "Fungi", mucoromycota = "Fungi")

norm_kingdom <- function(k) {
  kl <- tolower(trimws(k))
  ifelse(kl %in% c("animalia", "metazoa"), "Metazoa",
  ifelse(kl %in% c("chromista", "protozoa", "protista", "harosa", "discoba",
                   "eukaryota", "sar", "excavata", "rhizaria", "alveolata",
                   "amoebozoa"), "Protista",
  ifelse(kl == "fungi", "Fungi",
  ifelse(kl %in% c("plantae", "viridiplantae"), "Plants & algae",
  ifelse(nzchar(kl), "Other", "(unresolved)")))))
}
d$group <- norm_kingdom(d$kingdom)
## ------------------------------------------------ which phylum column -------
use_norm <- identical(CFG$phylum_col, "auto") &&
  "phylum_normalised" %in% names(d)
if (use_norm) {
  d$phylum_normalised <- ifelse(is.na(d$phylum_normalised), "",
                                as.character(d$phylum_normalised))
  d$phylum_plot <- ifelse(nzchar(d$phylum_normalised), d$phylum_normalised,
                          d$phylum)
  changed <- d[nzchar(d$phylum) & nzchar(d$phylum_normalised) &
               d$phylum != d$phylum_normalised, ]
  cat("phylum column plotted        : phylum_normalised\n")
  if (nrow(changed)) {
    cat("  reconciled across sources   : ", nrow(changed), " row(s)\n", sep = "")
    tb <- sort(table(paste(changed$phylum, "->", changed$phylum_normalised)),
               decreasing = TRUE)
    for (i in seq_along(tb))
      cat(sprintf("  %6d  %s\n", as.integer(tb)[i], names(tb)[i]))
  }
} else {
  d$phylum_plot <- d$phylum
  cat("phylum column plotted        : phylum (raw)",
      if (!"phylum_normalised" %in% names(d))
        "  << no phylum_normalised column in this table" else "", "\n", sep = "")
}
from_phylum <- PHYLUM_GROUP[tolower(trimws(ifelse(nzchar(d$phylum_plot),
                                                  d$phylum_plot, d$phylum)))]
rescued <- d$group %in% c("(unresolved)", "Other") & !is.na(from_phylum)
if (any(rescued)) {
  qc("kingdom missing, grouped by phylum instead: %d row(s) (%s)\n",
     sum(rescued), paste(unique(d$phylum_plot[rescued]), collapse = ", "))
  d$group[rescued] <- unname(from_phylum[rescued])
}

d$phylum_lab <- ifelse(nzchar(d$phylum_plot), d$phylum_plot,
                       "(no phylum resolved)")

# A phylum name that is numeric, one character, or otherwise not name-shaped
# means the columns have shifted -- the classic symptom is a bar labelled "0"
# picked up from an is_marine column.
bad_name <- unique(d$phylum_plot[nzchar(d$phylum_plot) &
                                 !grepl("^[A-Z][A-Za-z .-]{2,}$", d$phylum_plot)])
if (length(bad_name)) {
  qc("\n!! %d phylum value(s) do not look like taxon names: %s\n",
     length(bad_name), paste(utils::head(bad_name, 8), collapse = ", "))
  qc("   That usually means the columns are misaligned in %s.\n", CFG$api_file)
  qc("   Column names read: %s\n", paste(names(d), collapse = ", "))
}
qc("distinct phyla to plot       : %d\n",
   length(unique(d$phylum_lab[d$phylum_lab != "(no phylum resolved)"])))

## ------------------------------------------------------- unresolved hosts ---
no_phy <- d[!nzchar(d$phylum_plot), ]
qc("\nhost rows with no phylum    : %d (%.1f%%)\n", nrow(no_phy),
   100 * nrow(no_phy) / nrow(d))
if (nrow(no_phy)) {
  aff <- length(unique(no_phy$microsporidia_species))
  only <- setdiff(unique(no_phy$microsporidia_species),
                  unique(d$microsporidia_species[nzchar(d$phylum_plot)]))
  qc("  microsporidia species affected            : %d\n", aff)
  qc("  ... of which have NO placed host at all   : %d  <- lost from the figure\n",
     length(only))
  tb <- sort(table(no_phy$host_name), decreasing = TRUE)
  qc("  most frequent unplaced hosts: %s\n",
     paste(sprintf("%s (%d)", names(tb)[seq_len(min(6, length(tb)))],
                   as.integer(tb)[seq_len(min(6, length(tb)))]), collapse = ", "))
  write.table(unique(no_phy[, c("host_name", "microsporidia_species")]),
              file.path(CFG$outdir, "QC_hosts_without_phylum.tsv"), sep = "\t",
              row.names = FALSE, quote = FALSE)
}

## ------------------------------------------------------------- long table ---
if (CFG$unit == "microsporidia") {
  if (CFG$subset != "all") {
    n0 <- length(unique(d$microsporidia_species))
    d <- d[d$microsporidia_category == CFG$subset, ]
    qc("\nsubset = %s: %d of %d species kept\n", CFG$subset,
       length(unique(d$microsporidia_species)), n0)
  }
  long <- d[, c("microsporidia_species", "microsporidia_category",
                "phylum_lab", "group")]
  # a species with several hosts in one phylum is ONE species but several rows
  if (identical(CFG$count_unit, "species")) long <- unique(long)
  names(long) <- c("entity", "stack", "phylum_lab", "group")
  long$stack <- factor(long$stack, levels = c("named", "provisional"),
                       labels = c("Named species", "Provisional species"))
  unit_label <- if (identical(CFG$count_unit, "records"))
    "species-host records" else "microsporidia species"
  stack_title <- NULL
} else {
  # ONE ROW PER HOST SPECIES. The stack is decided by the best parasite the host
  # carries: a host infected by at least one NAMED microsporidian counts as a
  # host of a named species, whatever else also infects it; only hosts whose
  # every parasite is provisional count as provisional. The two categories are
  # therefore mutually exclusive and no host is counted twice.
  any_named <- tapply(d$microsporidia_category, d$host_name,
                      function(x) any(x == "named"))
  key <- unique(d[, c("host_name", "phylum_lab", "group")])
  key$stack <- factor(ifelse(unname(any_named[key$host_name]),
                             "Host of a named species",
                             "Host of provisional species only"),
                      levels = c("Host of a named species",
                                 "Host of provisional species only"))
  long <- data.frame(entity = key$host_name, stack = key$stack,
                     phylum_lab = key$phylum_lab, group = key$group,
                     stringsAsFactors = FALSE)
  # a host species belongs to one phylum, so unique() here is already one row
  # per host and count_unit does not apply
  unit_label <- "host species"
  stack_title <- NULL
}
long <- long[!is.na(long$stack), ]

## ------------------------------------------------- multi-phylum entities ----
n_phy <- tapply(long$phylum_lab, long$entity, function(x) length(unique(x)))
qc("\n%s in the figure : %d\n", unit_label, length(n_phy))
qc("  in 1 phylum   : %d (%.1f%%)\n", sum(n_phy == 1), 100 * mean(n_phy == 1))
if (any(n_phy > 1)) {
  qc("  in >1 phylum  : %d  (counted once in EACH, so column totals exceed the\n",
     sum(n_phy > 1))
  qc("                  number of %s)\n", unit_label)
  multi <- names(n_phy)[n_phy > 1]
  mt <- do.call(rbind, lapply(utils::head(multi[order(-n_phy[multi])], 8),
    function(s) data.frame(entity = s,
                           phyla = paste(sort(unique(long$phylum_lab[long$entity == s])),
                                         collapse = "; "))))
  print(mt, row.names = FALSE, right = FALSE)
  write.table(data.frame(entity = multi, n_phyla = as.integer(n_phy[multi])),
              file.path(CFG$outdir, "QC_multi_phylum.tsv"), sep = "\t",
              row.names = FALSE, quote = FALSE)
}

if (CFG$drop_unresolved) {
  n0 <- nrow(long)
  long <- long[long$phylum_lab != "(no phylum resolved)" &
               long$group != "(unresolved)", ]
  qc("\ndrop_unresolved = TRUE: %d %s hidden from the figure\n",
     n0 - nrow(long), unit_label)
  qc("  (still counted in the console table and the QC files above)\n")
}

## -------------------------------------------------------------- counting ----
counts <- as.data.frame(table(phylum_lab = long$phylum_lab, stack = long$stack),
                        stringsAsFactors = FALSE)
tot <- aggregate(Freq ~ phylum_lab, counts, sum)
names(tot)[2] <- "total"
grp <- unique(long[, c("phylum_lab", "group")])
grp <- grp[!duplicated(grp$phylum_lab), ]
tot <- merge(tot, grp, by = "phylum_lab")
tot <- tot[order(-tot$total), ]

pooled <- tot$phylum_lab[tot$total < CFG$min_n]
if (length(pooled)) {
  qc("\npooled into 'Other' (fewer than %d %s): %s\n", CFG$min_n, unit_label,
     paste(pooled, collapse = ", "))
  long$phylum_lab[long$phylum_lab %in% pooled] <- "Other"
  long$group[long$phylum_lab == "Other"] <- "Other"
  counts <- as.data.frame(table(phylum_lab = long$phylum_lab, stack = long$stack),
                          stringsAsFactors = FALSE)
  tot <- aggregate(Freq ~ phylum_lab, counts, sum); names(tot)[2] <- "total"
  grp <- unique(long[, c("phylum_lab", "group")])
  grp <- grp[!duplicated(grp$phylum_lab), ]
  tot <- merge(tot, grp, by = "phylum_lab"); tot <- tot[order(-tot$total), ]
}

wide <- reshape(counts, idvar = "phylum_lab", timevar = "stack",
                direction = "wide")
names(wide) <- sub("^Freq\\.", "", names(wide))
wide <- merge(wide, tot[, c("phylum_lab", "group", "total")], by = "phylum_lab")
wide <- wide[order(-wide$total), ]
wide$`% of all` <- sprintf("%.1f", 100 * wide$total / length(unique(long$entity)))

# Both numbers, side by side, so a bar can be reconciled against the input
# table without arithmetic: the table has ROWS, the figure counts SPECIES.
rec <- as.data.frame(table(phylum_lab = d$phylum_lab), stringsAsFactors = FALSE)
names(rec)[2] <- "host_records"
rec$unique_hosts <- as.integer(tapply(d$host_name, d$phylum_lab,
                                      function(x) length(unique(x)))[rec$phylum_lab])
rec$unique_microsporidia <- as.integer(
  tapply(d$microsporidia_species, d$phylum_lab,
         function(x) length(unique(x)))[rec$phylum_lab])
wide <- merge(wide, rec, by = "phylum_lab", all.x = TRUE)
wide <- wide[order(-wide$total), ]

cat(sprintf("\n================ %s PER HOST PHYLUM ================\n",
            toupper(unit_label)))
print(wide, row.names = FALSE, right = FALSE)
cat("\nhost_records = rows in the input table. unique_hosts and",
    " unique_microsporidia\nare the two things a bar could count; this figure",
    " counts ", unit_label, ".\n", sep = "")
write.table(wide, file.path(CFG$outdir, "counts_by_phylum.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(long, file.path(CFG$outdir, "per_entity_phylum.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)

qc("\nby kingdom-level group:\n")
gg <- aggregate(Freq ~ group, merge(counts, grp, by = "phylum_lab"), sum)
print(gg[order(-gg$Freq), ], row.names = FALSE, right = FALSE)

## ------------------------------------------------------------------ plot ----
# Okabe-Ito, same blue/orange as the discovery and tissue panels
pal <- if (CFG$unit == "microsporidia")
  c("Named species" = "#0072B2", "Provisional species" = "#E69F00") else
  c("Host of a named species" = "#0072B2",
    "Host of provisional species only" = "#E69F00")

lvl <- rev(tot$phylum_lab)
counts$phylum_lab <- factor(counts$phylum_lab, levels = lvl)
counts$stack <- factor(counts$stack, levels = names(pal))
lab <- aggregate(Freq ~ phylum_lab, counts, sum)

# Bar labels: the total, then the stack breakdown in parentheses, each part
# tinted to match its segment. For a phylum with 8 species it is otherwise
# impossible to read the split off a 3-pixel bar.
brk <- reshape(counts, idvar = "phylum_lab", timevar = "stack",
               direction = "wide")
names(brk) <- sub("^Freq\\.", "", names(brk))
stack_levels <- names(pal)[names(pal) %in% as.character(counts$stack)]
lab$parts_plain <- vapply(as.character(lab$phylum_lab), function(ph) {
  v <- vapply(stack_levels, function(s)
    as.integer(brk[[s]][brk$phylum_lab == ph] %||% 0), integer(1))
  paste(v, collapse = "+")
}, character(1))
lab$parts_rich <- vapply(as.character(lab$phylum_lab), function(ph) {
  v <- vapply(stack_levels, function(s)
    as.integer(brk[[s]][brk$phylum_lab == ph] %||% 0), integer(1))
  paste(sprintf("<span style='color:%s'>%d</span>", pal[stack_levels], v),
        collapse = "<span style='color:#7f7f7f'>+</span>")
}, character(1))
lab$txt_plain <- sprintf("%d (%s)", lab$Freq, lab$parts_plain)
lab$txt_rich <- sprintf("<span style='color:#4d4d4d'>%d</span> "
                        , lab$Freq)
lab$txt_rich <- paste0(lab$txt_rich,
                       "<span style='color:#7f7f7f'>(</span>", lab$parts_rich,
                       "<span style='color:#7f7f7f'>)</span>")

# The kingdom group must be attached to EVERY layer's data before the plot is
# built. ggplot captures a layer's data frame at the time the layer is added,
# so adding the column afterwards leaves the label layer without it -- and a
# layer with no facetting variable is drawn in every panel, which silently
# replicated all the labels and all the y-axis levels across the facets.
tot$group <- factor(tot$group, levels = c("Metazoa", "Protista", "Fungi",
                                          "Plants & algae", "Other",
                                          "(unresolved)"))
grp_of <- setNames(as.character(tot$group), tot$phylum_lab)
counts$group <- factor(grp_of[as.character(counts$phylum_lab)],
                       levels = levels(tot$group))
lab$group <- factor(grp_of[as.character(lab$phylum_lab)],
                    levels = levels(tot$group))
facet_it <- length(unique(na.omit(counts$group))) > 1

p <- ggplot(counts, aes(x = Freq, y = phylum_lab, fill = stack)) +
  geom_col(width = 0.74, colour = "grey20", linewidth = 0.2,
           position = position_stack(reverse = TRUE)) +
  {if (has_ggtext)
     ggtext::geom_richtext(data = lab, aes(x = Freq, y = phylum_lab,
                                           label = txt_rich),
                           inherit.aes = FALSE, hjust = 0, size = 2.7,
                           fill = NA, label.color = NA,
                           label.padding = unit(c(0, 0, 0, 1.6), "mm"))
   else
     geom_text(data = lab, aes(x = Freq, y = phylum_lab, label = txt_plain),
               inherit.aes = FALSE, hjust = -0.10, size = 2.7,
               colour = "grey25")} +
  scale_fill_manual(values = pal, name = stack_title) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.30))) +
  scale_y_discrete(drop = TRUE) +
  labs(x = sprintf("Number of %s", unit_label), y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    axis.text = element_text(colour = "black"),
    axis.title.x = element_text(margin = margin(t = 8)),
    legend.position = if (facet_it) "top" else c(0.98, 0.06),
    legend.justification = if (facet_it) "right" else c(1, 0),
    legend.background = element_blank(),
    legend.key.size = unit(0.9, "lines"),
    plot.margin = margin(6, 16, 4, 6)
  )

if (facet_it)
  p <- p + facet_grid(rows = vars(group), scales = "free_y", space = "free_y",
                      switch = "y", drop = TRUE) +
    theme(strip.placement = "outside",
          strip.background = element_rect(fill = "grey95", colour = NA),
          strip.text.y.left = element_text(angle = 0, face = "bold", size = 9),
          panel.spacing.y = unit(0.35, "lines"))

ggsave(file.path(CFG$outdir, "fig_host_phylum.pdf"), p,
       width = CFG$fig_width, height = CFG$fig_height,
       device = if (capabilities("cairo")) cairo_pdf else pdf)
ggsave(file.path(CFG$outdir, "fig_host_phylum.png"), p,
       width = CFG$fig_width, height = CFG$fig_height, dpi = 600)
qc("\nwrote fig_host_phylum.pdf / .png\n")

## ------------------------------------------------- optional log companion ---
if (CFG$log_panel) {
  # DODGED, not stacked: segment lengths on a log axis do not add, so a stacked
  # log bar is a lie. Dodged bars on a log axis are fine.
  cl <- counts[counts$Freq > 0, ]
  pl <- ggplot(cl, aes(x = Freq, y = phylum_lab, fill = stack)) +
    geom_col(width = 0.7, position = position_dodge(width = 0.75),
             colour = "grey20", linewidth = 0.2) +
    scale_fill_manual(values = pal, name = stack_title) +
    scale_x_log10(expand = expansion(mult = c(0, 0.08))) +
    annotation_logticks(sides = "b", outside = FALSE,
                        short = unit(0.05, "cm"), mid = unit(0.08, "cm"),
                        long = unit(0.12, "cm")) +
    labs(x = sprintf("Number of %s (log scale)", unit_label), y = NULL) +
    theme_classic(base_size = 11) +
    theme(axis.text = element_text(colour = "black"),
          legend.position = "top", legend.key.size = unit(0.9, "lines"))
  ggsave(file.path(CFG$outdir, "fig_host_phylum_log.pdf"), pl,
         width = CFG$fig_width, height = CFG$fig_height,
         device = if (capabilities("cairo")) cairo_pdf else pdf)
  ggsave(file.path(CFG$outdir, "fig_host_phylum_log.png"), pl,
         width = CFG$fig_width, height = CFG$fig_height, dpi = 600)
  qc("wrote fig_host_phylum_log.pdf / .png (dodged, not stacked -- see comment)\n")
}

## ------------------------------------- optional comparison with diversity ---
# Described extant species per phylum. These are literature values, not
# something this pipeline can derive, so they live in an editable table and
# every phylum without one is reported rather than silently dropped.
EXTANT <- data.frame(
  phylum = c("Arthropoda", "Chordata", "Platyhelminthes", "Annelida", "Mollusca",
             "Nematoda", "Bryozoa", "Cnidaria", "Rotifera", "Nemertea",
             "Acanthocephala", "Porifera", "Gastrotricha", "Cephalorhyncha",
             "Dicyemida", "Phoronida", "Apicomplexa", "Ciliophora", "Amoebozoa",
             "Cercozoa"),
  extant = c(1082297, 69913, 18616, 14399, 65442, 3455, 5434, 11151, 2014, 1373,
             1330, 9092, 852, 237, 122, 19, 5000, 8613, 48, 100),
  stringsAsFactors = FALSE)

if (CFG$show_extant) {
  m <- merge(tot, EXTANT, by.x = "phylum_lab", by.y = "phylum", all.x = TRUE)
  no_ref <- m$phylum_lab[is.na(m$extant)]
  if (length(no_ref))
    qc("\nno extant-species reference for: %s\n", paste(no_ref, collapse = ", "))
  m <- m[!is.na(m$extant), ]
  m$per_1000 <- 1000 * m$total / m$extant
  m <- m[order(-m$per_1000), ]
  cat("\nhosts per 1000 described species in the phylum:\n")
  print(within(m[, c("phylum_lab", "total", "extant", "per_1000")],
               per_1000 <- sprintf("%.2f", per_1000)),
        row.names = FALSE, right = FALSE)
  write.table(m, file.path(CFG$outdir, "counts_vs_extant.tsv"), sep = "\t",
              row.names = FALSE, quote = FALSE)

  pe <- ggplot(m, aes(x = per_1000, y = reorder(phylum_lab, per_1000))) +
    geom_col(width = 0.7, fill = "grey45", colour = "grey20", linewidth = 0.2) +
    geom_text(aes(label = sprintf("%.2f", per_1000)), hjust = -0.18, size = 2.7,
              colour = "grey20") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.16))) +
    labs(x = sprintf("%s per 1000 described species", unit_label), y = NULL,
         caption = "described-species counts are literature values, see EXTANT in this script") +
    theme_classic(base_size = 11) +
    theme(axis.text = element_text(colour = "black"),
          plot.caption = element_text(size = 7, colour = "grey40"))
  ggsave(file.path(CFG$outdir, "fig_host_phylum_per_diversity.pdf"), pe,
         width = CFG$fig_width, height = CFG$fig_height * 0.8,
         device = if (capabilities("cairo")) cairo_pdf else pdf)
  ggsave(file.path(CFG$outdir, "fig_host_phylum_per_diversity.png"), pe,
         width = CFG$fig_width, height = CFG$fig_height * 0.8, dpi = 600)
  qc("wrote fig_host_phylum_per_diversity.pdf / .png\n")
}

qc("\nwritten to %s/\n  fig_host_phylum.pdf / .png\n  counts_by_phylum.tsv\n  per_entity_phylum.tsv\n  QC_*.tsv\n",
   CFG$outdir)
