#!/usr/bin/env Rscript
###############################################################################
# Figure 1A -- Microsporidia species descriptions over time
#
# Stacked bar chart of species per decade, split by NAMED (valid binomial)
# vs PROVISIONAL (no specific epithet / uncultured / unnamed etc.).
#
# Tested under R 4.3.3 with readxl 1.4.x + ggplot2 3.4.x.
#
# Design decisions, all explicit rather than silent:
#   * Column lookup is on a normalised key (NFKD -> ASCII -> [a-z0-9] only),
#     so the micro sign (U+00B5) vs Greek mu (U+03BC) trap cannot bite. A
#     missing column is a warning, never a silent zero.
#   * Date rule: every 4-digit year in the cell is extracted and the MINIMUM is
#     kept (date_rule = "min_all"). Set date_rule = "min_publication" to ignore
#     years inside a "(Samples Collected ...)" / "*collected" clause.
#     The script reports how many entries the two rules disagree on.
#   * Every cell that yields no year is counted, listed, and excluded --
#     nothing is silently coerced to NA. Missing dates are broken down by
#     category so a batch import with no dates cannot hide inside the total.
#   * Repeated header rows pasted into the data are detected and dropped.
#   * Classification uses the shared canonical rule in classify_species.R,
#     the single implementation used by every figure script and the flowchart.
#
# Usage:  Rscript fig1_discovery_by_decade.R [path/to/database.xlsx]
###############################################################################

## ---------------------------------------------------------------- config ----
CFG <- list(
  db_path    = "Microsporidia_Characteristics_Database_merge_pro4.xlsx",
  sheet      = "Actively Updated Masterlist",
  name_col   = "Species Name",
  date_col   = "Date Identified (year)",
  remarks_col = "Important Remarks",

  date_rule  = "min_all",      # "min_all" | "min_publication"
  year_min   = 1600,           # anything outside [year_min, year_max] is rejected
  year_max   = as.integer(format(Sys.Date(), "%Y")),

  # Dates written by fill_dates_from_ncbi.py carry a provenance tag in
  # Important Remarks ("date basis: GenBank submission year 2021 (OL117026)").
  # A submission year is not a description year, so those entries are counted
  # separately and can be dropped outright.
  #   "include" -- plot them, report how many per bin (default)
  #   "exclude" -- drop them from the figure
  derived_dates = "include",

  # opt-in second-chance parse: "collection year: YYYY" in Important Remarks,
  # for rows the fill script could not resolve. Also a proxy, also reported.
  remarks_fallback = FALSE,

  bin_width  = 10,
  bin_start  = 1850,


  max_print  = 20,             # rows shown per QC table on the console
  outdir     = "fig1_output",
  fig_width  = 7.5,
  fig_height = 4.6
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) CFG$db_path <- args[1]

## -------------------------------------------------------------- packages ----
need <- c("readxl", "ggplot2")
for (p in need) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({ library(readxl); library(ggplot2) })
has_stringi <- requireNamespace("stringi", quietly = TRUE)

dir.create(CFG$outdir, showWarnings = FALSE, recursive = TRUE)
old_width <- getOption("width"); options(width = 200)   # keep QC tables on one line
on.exit(options(width = old_width), add = TRUE)

# En dash in bin labels breaks graphics devices under a non-UTF-8 locale
# (R silently substitutes dots). Fall back to a hyphen when that is the case.
utf8_ok <- grepl("UTF-?8", Sys.getlocale("LC_CTYPE"), ignore.case = TRUE)
DASH <- if (utf8_ok) "\u2013" else "-"
if (!utf8_ok)
  message("note: locale is not UTF-8 (", Sys.getlocale("LC_CTYPE"),
          "); using '-' in axis labels")

## --------------------------------------------------------------- helpers ----
deaccent <- function(x) {
  if (has_stringi) stringi::stri_trans_general(x, "Latin-ASCII")
  else iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT", sub = "")
}

# normalised column key: kills the micro-sign / mu ambiguity and all punctuation
norm_key <- function(x) gsub("[^a-z0-9]", "", tolower(deaccent(as.character(x))))

get_col <- function(df, wanted, required = TRUE) {
  hit <- which(norm_key(names(df)) == norm_key(wanted))
  if (length(hit) == 0) {
    msg <- sprintf("column not found: %s\n  available: %s",
                   wanted, paste(names(df), collapse = " | "))
    if (required) stop(msg, call. = FALSE) else { warning(msg, call. = FALSE); return(NULL) }
  }
  if (length(hit) > 1)
    warning(sprintf("column %s matched %d headers; using the first", wanted, length(hit)),
            call. = FALSE)
  df[[hit[1]]]
}

# console-safe printing of QC tables: truncate long cells, cap rows
show_tbl <- function(df, n = CFG$max_print, width = 46, file = NULL) {
  if (nrow(df) == 0) return(invisible(NULL))
  d <- head(df, n)
  d[] <- lapply(d, function(col) {
    s <- as.character(col); s[is.na(s)] <- ""
    ifelse(nchar(s) > width, paste0(substr(s, 1, width - 3), "..."), s)
  })
  print(d, row.names = FALSE, right = FALSE)
  if (nrow(df) > n)
    cat(sprintf("  ... %d more row(s)%s\n", nrow(df) - n,
                if (is.null(file)) "" else paste0(", see ", file)))
}

write_qc <- function(df, file) {
  if (nrow(df) == 0) return(invisible(NULL))
  write.table(df, file.path(CFG$outdir, file), sep = "\t",
              row.names = FALSE, quote = TRUE, fileEncoding = "UTF-8")
}

## ----------------------------------------------------- classifier ----------
# Named vs provisional comes from classify_species.R -- the ONE canonical
# implementation shared by every figure script, the assembler and the
# flowchart. RX_NOM_1 / RX_NOM_2 (used by norm_binomial below) are defined
# there too. Do not redefine classify_species() here.
source("classify_species.R")

# genus + epithet key, for duplicate detection only
norm_binomial <- function(name) {
  vapply(name, function(nm) {
    s <- tolower(deaccent(nm))
    s <- gsub("\\(.*?\\)", " ", s, perl = TRUE)
    s <- gsub(RX_NOM_1, " ", s, perl = TRUE, ignore.case = TRUE)
    s <- gsub(RX_NOM_2, " ", s, perl = TRUE, ignore.case = TRUE)
    s <- gsub("[^a-z0-9\\s]", " ", s, perl = TRUE)
    tok <- setdiff(strsplit(trimws(s), "\\s+")[[1]], c("nov", "gen", "et", "and", ""))
    paste(head(tok, 2), collapse = " ")
  }, character(1), USE.NAMES = FALSE)
}

## ---------------------------------------------------------- date parsing ----
# Strip "(Samples Collected 2014-2016)" / "*collected 2019" style clauses so the
# collection year cannot masquerade as the description year.
COLLECT_RX <- "collect|sampl|receiv|source fish"

strip_collection_clauses <- function(s) {
  s <- gsub(sprintf("\\([^)]*(?:%s)[^)]*\\)", COLLECT_RX), " ", s, perl = TRUE, ignore.case = TRUE)
  gsub(sprintf("\\*[^;]*(?:%s)[^;]*", COLLECT_RX), " ", s, perl = TRUE, ignore.case = TRUE)
}

extract_years <- function(s) {
  m <- regmatches(s, gregexpr("(?<!\\d)\\d{4}(?!\\d)", s, perl = TRUE))[[1]]
  y <- suppressWarnings(as.integer(m))
  y[!is.na(y) & y >= CFG$year_min & y <= CFG$year_max]
}

earliest_year <- function(cells, rule = CFG$date_rule) {
  vapply(cells, function(x) {
    s <- if (is.na(x)) "" else as.character(x)
    if (rule == "min_publication") s <- strip_collection_clauses(s)
    y <- extract_years(s)
    if (length(y) == 0) NA_integer_ else min(y)
  }, integer(1), USE.NAMES = FALSE)
}

# NB: regmatches(x, regexpr(...)) returns ONLY the elements that matched, so
# indexing it positionally silently misaligns names and years. Match positions
# are carried explicitly instead.
match_text <- function(s, rx) {
  out <- rep(NA_character_, length(s))
  m <- regexpr(rx, s, perl = TRUE, ignore.case = TRUE)
  hit <- which(m > 0L)
  if (length(hit))
    out[hit] <- substring(s[hit], m[hit], m[hit] + attr(m, "match.length")[hit] - 1L)
  out
}

year_from_text <- function(txt) {
  y <- suppressWarnings(as.integer(sub(".*?((?:1[6-9]|20)\\d{2}).*", "\\1", txt, perl = TRUE)))
  y[is.na(y) | y < CFG$year_min | y > CFG$year_max] <- NA_integer_
  y
}

# "date basis: GenBank submission year 2021 (OL117026)" written by
# fill_dates_from_ncbi.py -- returns the basis label, or NA where absent.
date_basis_from_remarks <- function(s) {
  txt <- match_text(s, "date basis:[^;]*")
  basis <- rep(NA_character_, length(s))
  hit <- which(!is.na(txt))
  if (length(hit)) {
    b <- tolower(txt[hit])
    basis[hit] <- ifelse(grepl("submission", b), "GenBank submission year",
                  ifelse(grepl("publication", b), "publication year",
                  ifelse(grepl("collection", b), "sample collection year",
                  ifelse(grepl("record date", b), "GenBank record date", "other"))))
  }
  basis
}

collection_year_from_remarks <- function(s) {
  year_from_text(match_text(s, "collection year:\\s*\\d{4}"))
}

## ------------------------------------------------------------------ load ----
if (!file.exists(CFG$db_path)) stop("database not found: ", CFG$db_path, call. = FALSE)
raw <- read_excel(CFG$db_path, sheet = CFG$sheet, col_types = "text",
                  .name_repair = "minimal")
message(sprintf("read %d rows x %d columns from sheet '%s'",
                nrow(raw), ncol(raw), CFG$sheet))

remarks_raw <- get_col(raw, CFG$remarks_col, required = FALSE)
dat <- data.frame(
  row_xlsx     = seq_len(nrow(raw)) + 1L,          # +1 for the header row
  species_name = trimws(get_col(raw, CFG$name_col)),
  date_raw     = get_col(raw, CFG$date_col),
  remarks      = if (is.null(remarks_raw)) "" else remarks_raw,
  stringsAsFactors = FALSE
)
dat$date_raw[is.na(dat$date_raw)] <- ""
dat$remarks[is.na(dat$remarks)]   <- ""

qc <- function(...) cat(sprintf(...), sep = "")
cat("\n================ PARSING REPORT ================\n")

# repeated header rows pasted in when sheets are appended
hdr <- which(!is.na(dat$species_name) &
             norm_key(dat$species_name) == norm_key(CFG$name_col))
if (length(hdr)) {
  qc("!! %d repeated HEADER row(s) found inside the data (xlsx row %s) -- dropped\n",
     length(hdr), paste(dat$row_xlsx[hdr], collapse = ", "))
  dat <- dat[-hdr, ]
}

n_blank <- sum(is.na(dat$species_name) | !nzchar(dat$species_name))
if (n_blank) qc("dropping %d row(s) with a blank Species Name\n", n_blank)
dat <- dat[!is.na(dat$species_name) & nzchar(dat$species_name), ]

dat$category <- factor(classify_species(dat$species_name),
                       levels = c("named", "provisional"),
                       labels = c("Named species", "Provisional species"))
dat$year      <- earliest_year(dat$date_raw)
dat$year_alt  <- earliest_year(dat$date_raw,
                               rule = if (CFG$date_rule == "min_all") "min_publication" else "min_all")

## ------------------------------------------------- date provenance ----------
# A year in the date column is a described-species date UNLESS the row carries a
# "date basis:" tag, in which case it was back-filled from GenBank.
dat$basis     <- date_basis_from_remarks(dat$remarks)
dat$derived   <- !is.na(dat$year) & !is.na(dat$basis)
dat$year_src  <- ifelse(is.na(dat$year), NA_character_,
                        ifelse(dat$derived, dat$basis, "date column"))

## ------------------------------------------------- optional remarks rescue ---
rescue <- data.frame()
if (any(is.na(dat$year))) {
  i <- which(is.na(dat$year))
  cy <- collection_year_from_remarks(dat$remarks[i])
  hit <- i[!is.na(cy)]
  if (length(hit)) {
    rescue <- data.frame(species_name = dat$species_name[hit],
                         collection_year = cy[!is.na(cy)],
                         stringsAsFactors = FALSE)
    if (CFG$remarks_fallback) {
      dat$year[hit]     <- cy[!is.na(cy)]
      dat$year_alt[hit] <- cy[!is.na(cy)]
      dat$year_src[hit] <- "remarks collection year"
      dat$derived[hit]  <- TRUE
    }
  }
}

## ---------------------------------------------------------- QC / warnings ----
qc("date rule            : %s%s\n", CFG$date_rule,
   if (CFG$remarks_fallback) " (+ remarks collection-year fallback)" else "")
qc("entries with a name  : %d  (named %d, provisional %d)\n",
   nrow(dat), sum(dat$category == "Named species"),
   sum(dat$category == "Provisional species"))

if (any(dat$derived)) {
  qc("\ndates back-filled from sequence metadata: %d\n", sum(dat$derived))
  bs <- as.data.frame(table(basis = dat$year_src[dat$derived],
                            category = dat$category[dat$derived]),
                      responseName = "n")
  bs <- bs[bs$n > 0, ]
  print(bs, row.names = FALSE, right = FALSE)
  qc("   a submission year is an upper bound on the description date, not the\n")
  qc("   description date itself; CFG$derived_dates = \"%s\"\n", CFG$derived_dates)
  if (identical(CFG$derived_dates, "exclude")) {
    qc("   -> excluded from the figure\n")
    dat$year[dat$derived] <- NA_integer_
  }
}

no_year <- dat[is.na(dat$year), ]
qc("\ncells yielding no year: %d  (excluded from the figure)\n", nrow(no_year))
if (nrow(no_year)) {
  # break the loss down, so a dateless batch import cannot hide in the total
  brk <- as.data.frame(table(category = no_year$category), responseName = "missing")
  brk$total <- as.integer(table(dat$category)[as.character(brk$category)])
  brk$`% of category` <- sprintf("%.1f", 100 * brk$missing / brk$total)
  print(brk, row.names = FALSE, right = FALSE)
  if (max(brk$missing / brk$total) > 0.10) {
    qc("!! a category is missing >10%% of its dates: the stack under-counts that\n")
    qc("   category in the affected bins. Interpret those bars with care.\n")
  }
  show_tbl(data.frame(xlsx_row = no_year$row_xlsx,
                      species  = no_year$species_name,
                      category = as.character(no_year$category),
                      cell     = ifelse(nzchar(no_year$date_raw), no_year$date_raw, "<empty>")),
           file = "QC_no_year_parsed.tsv")
  write_qc(no_year[, c("row_xlsx", "species_name", "category", "date_raw")],
           "QC_no_year_parsed.tsv")
}

if (nrow(rescue)) {
  qc("\nof those, %d carry 'collection year: YYYY' in %s (%s)\n",
     nrow(rescue), CFG$remarks_col,
     if (CFG$remarks_fallback) "USED" else "not used; set CFG$remarks_fallback = TRUE")
  show_tbl(rescue, n = 5, file = "QC_remarks_collection_years.tsv")
  write_qc(rescue, "QC_remarks_collection_years.tsv")
}

disagree <- dat[!is.na(dat$year) & !is.na(dat$year_alt) & dat$year != dat$year_alt, ]
qc("\nentries where 'min_all' and 'min_publication' disagree: %d\n", nrow(disagree))
if (nrow(disagree)) {
  show_tbl(data.frame(species = disagree$species_name, cell = disagree$date_raw,
                      used = disagree$year, other_rule = disagree$year_alt),
           file = "QC_date_rule_disagreement.tsv")
  qc("  -> collection/sampling years, not description years. Switch CFG$date_rule to compare.\n")
  write_qc(disagree[, c("species_name", "date_raw", "year", "year_alt")],
           "QC_date_rule_disagreement.tsv")
}

n_multi <- sum(vapply(dat$date_raw, function(x) length(unique(extract_years(x))) > 1,
                      logical(1)))
qc("\ncells listing more than one distinct year: %d (minimum taken)\n", n_multi)

## duplicates -----------------------------------------------------------------
dat$key <- norm_binomial(dat$species_name)

exact <- table(tolower(dat$species_name))
exact <- names(exact)[exact > 1]
qc("exact duplicate Species Name strings: %d name(s), %d row(s)\n",
   length(exact), sum(tolower(dat$species_name) %in% exact))
if (length(exact)) {
  d <- dat[tolower(dat$species_name) %in% exact, c("row_xlsx", "species_name", "date_raw")]
  d <- d[order(tolower(d$species_name)), ]
  show_tbl(d, file = "QC_exact_duplicate_names.tsv")
  write_qc(d, "QC_exact_duplicate_names.tsv")
}

dup_keys <- names(which(table(dat$key[dat$category == "Named species"]) > 1))
qc("\nduplicate genus+epithet keys among named species: %d key(s), %d row(s)\n",
   length(dup_keys), sum(dat$key %in% dup_keys & dat$category == "Named species"))
if (length(dup_keys)) {
  d <- dat[dat$key %in% dup_keys & dat$category == "Named species",
           c("key", "species_name", "date_raw")]
  d <- d[order(d$key), ]
  show_tbl(d, file = "QC_duplicate_named_keys.tsv")
  write_qc(d, "QC_duplicate_named_keys.tsv")
}

## name-format problems -------------------------------------------------------
odd <- dat[!grepl("^[A-Z][a-z]", dat$species_name), c("row_xlsx", "species_name")]
if (nrow(odd)) {
  qc("\nspecies names not starting with a capitalised genus: %d\n", nrow(odd))
  show_tbl(odd, n = 10, file = "QC_odd_species_names.tsv")
  write_qc(odd, "QC_odd_species_names.tsv")
}

## near-duplicate genus spellings (Plistophora vs Pleistophora etc.) ----------
genus <- vapply(strsplit(dat$key, " "), function(z) if (length(z)) z[1] else "", character(1))
genus <- sort(unique(genus[nzchar(genus)]))
dm <- utils::adist(genus)
pairs <- which(dm == 1 & upper.tri(dm), arr.ind = TRUE)
if (nrow(pairs)) {
  qc("\ngenus names differing by one character (possible misspellings):\n")
  gp <- data.frame(genus_a = genus[pairs[, 1]], genus_b = genus[pairs[, 2]],
                   n_a = as.integer(table(genus)[genus[pairs[, 1]]]),
                   n_b = as.integer(table(genus)[genus[pairs[, 2]]]))
  print(gp, row.names = FALSE, right = FALSE)
  write_qc(gp, "QC_similar_genus_names.tsv")
}

## ---------------------------------------------------------------- binning ----
plt <- dat[!is.na(dat$year), ]
if (nrow(plt) == 0) stop("no entries have a parseable year", call. = FALSE)
breaks <- seq(CFG$bin_start, CFG$year_max + CFG$bin_width, by = CFG$bin_width)
if (min(plt$year) < CFG$bin_start)
  stop(sprintf("%d entries fall before bin_start (%d); lower CFG$bin_start",
               sum(plt$year < CFG$bin_start), CFG$bin_start), call. = FALSE)

plt$bin_start <- breaks[findInterval(plt$year, breaks)]
# the final bin is truncated at year_max rather than running to the end of the
# decade, so the label states the years the data actually cover
lab_of <- function(b) sprintf("%d%s%d", b, DASH,
                              pmin(b + CFG$bin_width - 1, CFG$year_max))

used <- seq(min(plt$bin_start), max(plt$bin_start), by = CFG$bin_width)  # keep internal zeros
plt$bin <- factor(lab_of(plt$bin_start), levels = lab_of(used))

tab <- as.data.frame(table(bin = plt$bin, category = plt$category),
                     responseName = "n")
tab$bin <- factor(tab$bin, levels = lab_of(used))

wide <- reshape(tab, idvar = "bin", timevar = "category", direction = "wide")
names(wide) <- sub("^n\\.", "", names(wide))
wide$Total <- rowSums(wide[, -1, drop = FALSE])
wide$`% provisional` <- sprintf("%.1f", 100 * wide$`Provisional species` / wide$Total)
wide$`% provisional`[wide$Total == 0] <- "-"
if (any(plt$derived)) {
  d <- table(factor(plt$bin[plt$derived], levels = lab_of(used)))
  wide$`of which back-filled` <- as.integer(d[as.character(wide$bin)])
}

cat("\n================ COUNTS BEHIND THE FIGURE ================\n")
print(wide, row.names = FALSE, right = FALSE)
qc("\nplotted: %d entries (%d named, %d provisional); %d excluded for no year\n",
   nrow(plt), sum(plt$category == "Named species"),
   sum(plt$category == "Provisional species"), nrow(no_year))
last_bin <- max(used)
last_span <- CFG$year_max - last_bin + 1
if (last_span < CFG$bin_width)
  qc("NOTE: final bin %s spans %d years, not %d -- its height is not directly\n      comparable with the full bins.\n",
     lab_of(last_bin), last_span, CFG$bin_width)

write_qc(wide, "fig1_counts_by_decade.tsv")
write_qc(plt[, c("row_xlsx", "species_name", "category", "date_raw",
                 "year", "year_src", "derived", "bin")], "fig1_per_species_year.tsv")

## ------------------------------------------------------------------ plot ----
# Okabe-Ito: blue + orange, safe for deuteranopia/protanopia and in greyscale
pal <- c("Named species" = "#0072B2", "Provisional species" = "#E69F00")

totals <- aggregate(n ~ bin, tab, sum)

p <- ggplot(tab, aes(x = bin, y = n, fill = category)) +
  # reverse = TRUE puts the first factor level (Named) at the BOTTOM of the stack
  geom_col(width = 0.78, colour = "grey20", linewidth = 0.2,
           position = position_stack(reverse = TRUE)) +
  geom_text(data = totals[totals$n > 0, ], aes(x = bin, y = n, label = n),
            inherit.aes = FALSE, vjust = -0.5, size = 2.9, colour = "grey20") +
  scale_fill_manual(values = pal, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  labs(x = "Decade of description", y = "Number of species reported") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    axis.title.x     = element_text(margin = margin(t = 8)),
    axis.title.y     = element_text(margin = margin(r = 8)),
    legend.position  = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.key.size  = unit(0.9, "lines"),
    plot.margin      = margin(6, 10, 4, 6)
  )

pdf_dev <- if (capabilities("cairo")) cairo_pdf else pdf
ggsave(file.path(CFG$outdir, "fig1_discovery_by_decade.pdf"), p,
       width = CFG$fig_width, height = CFG$fig_height, device = pdf_dev)
ggsave(file.path(CFG$outdir, "fig1_discovery_by_decade.png"), p,
       width = CFG$fig_width, height = CFG$fig_height, dpi = 600)

qc("\nwritten to %s/\n  fig1_discovery_by_decade.pdf / .png\n  fig1_counts_by_decade.tsv\n  fig1_per_species_year.tsv\n  QC_*.tsv\n",
   CFG$outdir)
