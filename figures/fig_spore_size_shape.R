#!/usr/bin/env Rscript
###############################################################################
# Figure -- spore size and shape
#
# Scatter of spore length vs width, coloured by shape class, parsed out of the
# free-text "Spore Shape (Class; Condition)", "Spore Length Average (um)" and
# "Spore Width Average (um)" columns.
#
# WHAT MAKES THIS HARD
#   Shape cells hold several descriptors with conditions:
#     "oval; pyriform"
#     "pyriform (c1); oval (c2)"
#     "ovoid (meiospore, fresh); collapsed (meiospore, stained)"
#     "rod (fresh); ovocylindrical (fixed)"
#   Measurement cells hold several values, errors, ranges and conditions:
#     "3.15(fresh); 3.1(stained)"
#     "4.03 +/- 0.28 (3.71 - 4.33 um) (fresh)"
#     "3.7 - 4.9 um (avg: 4.4 um)"
#     "4 (normal spores); 6 (macrospores)"
#
# EVERY RULE IS STATED AND EVERY DROPPED VALUE IS REPORTED. Nothing is silently
# coerced to NA and nothing is silently clipped off the axes.
#
# Usage:  Rscript fig_spore_size_shape.R [path/to/database.xlsx]
###############################################################################

## ---------------------------------------------------------------- config ----
CFG <- list(
  db_path   = "Microsporidia_Characteristics_Database_merge_pro4.xlsx",
  sheet     = "Actively Updated Masterlist",
  name_col  = "Species Name",
  len_col   = "Spore Length Average (\u00b5m)",   # MICRO SIGN U+00B5
  wid_col   = "Spore Width Average (\u00b5m)",
  shape_col = "Spore Shape (Class; Condition)",
  coil_avg_col   = "Polar Tubule Coils Average",
  coil_range_col = "Polar Tubule Coils Range",

  # which measurement to take when a cell offers several, best first
  condition_priority = c("fresh", "live", "unspecified", "fixed", "stained",
                         "preserved", "other"),
  # spore types to avoid when an ordinary spore measurement is also present
  deprioritise_types = c("macrospore", "microspore", "meiospore", "aberrant",
                         "binucleate", "transovarial", "octospore"),

  # a species whose shape cell names more than one class
  #   "first"  use the first descriptor (as entered = the primary shape)
  #   "drop"   exclude the species from the figure
  #   "all"    plot it once per class (points are then not independent)
  multi_shape = "first",

  # how the panel encodes shape and coils
  #   "colour_circles"  every point a circle, colour = shape class (as before)
  #   "colour_symbols"  colour AND symbol = shape class
  #   "coil_heat"       symbol = shape class, colour = polar tubule coils,
  #                     entries with no coil data still drawn, in grey
  style = "coil_heat",
  coil_breaks = c(0, 4, 6, 8, 10, 13, 16, Inf),
  # recover coil ranges that Excel silently turned into dates (see COIL NOTES)
  repair_excel_dates = TRUE,
  # TRUE excludes entries whose shape and dimensions describe different spores
  # and cannot be reconciled; FALSE keeps and flags them.
  require_same_condition = TRUE,

  # "all" | "named" | "provisional"
  subset = "all",
  # NULL = fit to the data. Setting limits reports, not hides, what falls out.
  xlim = NULL,
  ylim = NULL,
  facet_by_category = FALSE,

  outdir     = "spore_output",
  fig_width  = 7.2,
  fig_height = 5.4
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) CFG$db_path <- args[1]

## -------------------------------------------------------------- packages ----
for (p in c("readxl", "ggplot2")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({ library(readxl); library(ggplot2) })
has_stringi <- requireNamespace("stringi", quietly = TRUE)

dir.create(CFG$outdir, showWarnings = FALSE, recursive = TRUE)
old_width <- getOption("width"); options(width = 200)
on.exit(options(width = old_width), add = TRUE)
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a[1])) b else a

qc <- function(...) cat(sprintf(...), sep = "")

## --------------------------------------------------------------- helpers ----
deaccent <- function(x) {
  if (has_stringi) stringi::stri_trans_general(x, "Latin-ASCII")
  else iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT", sub = "")
}
norm_key <- function(x) gsub("[^a-z0-9]", "", tolower(deaccent(as.character(x))))

# The spore columns use the MICRO SIGN (U+00B5) and the polar tubule columns use
# GREEK SMALL LETTER MU (U+03BC). They render identically and are not equal as
# strings, so headers are matched on a normalised key, never literally.
get_col <- function(df, wanted, required = TRUE) {
  hit <- which(norm_key(names(df)) == norm_key(wanted))
  if (length(hit) == 0) {
    msg <- sprintf("column not found: %s\n  available: %s",
                   wanted, paste(names(df), collapse = " | "))
    if (required) stop(msg, call. = FALSE) else { warning(msg, call. = FALSE); return(NULL) }
  }
  df[[hit[1]]]
}

show_tbl <- function(df, n = 25, width = 62, file = NULL) {
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  d <- head(df, n)
  d[] <- lapply(d, function(col) {
    s <- as.character(col); s[is.na(s)] <- ""
    ifelse(nchar(s) > width, paste0(substr(s, 1, width - 3), "..."), s)
  })
  print(d, row.names = FALSE, right = FALSE)
  if (nrow(df) > n)
    cat(sprintf("  ... %d more%s\n", nrow(df) - n,
                if (is.null(file)) "" else paste0(", see ", file)))
}

write_qc <- function(df, file) {
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  write.table(df, file.path(CFG$outdir, file), sep = "\t",
              row.names = FALSE, quote = TRUE, fileEncoding = "UTF-8")
}

## ------------------------------------------------------- species classifier --
# classify_species() comes from classify_species.R -- the ONE canonical
# named/provisional implementation, shared by every figure script, the
# assembler and the flowchart. Do not define a local copy here.
source("classify_species.R")

## ============================================================================
## SHAPE CLASSIFICATION
## ============================================================================
# Applied to one descriptor at a time, in this order. The first pattern that
# matches wins, so the specific compounds must come before the generic words:
# "oval-cylindrical" is Ovocylindrical, not Oval; "elongate-pyriform" is
# Pyriform, not Rod-shaped. This ordering replaces the old scheme of collecting
# every match and then applying a fixed priority, which classified
# "short-cylindrical" as Rod-shaped and "elongate-ovoid" as Rod-shaped too.
SHAPE_RULES <- list(
  Ovocylindrical = paste0("ovo[- ]?cylindric|oval[- ]?cylindric|barrel|biconic",
                          "|rounded[- ]?cylindric|short[- ]?cylindric",
                          "|shortly cylindric|hemispherical end"),
  Pyriform       = paste0("pyriform|pyrifrom|pyrifom|pear|flask|lageniform",
                          "|bottle[- ]?shape|drop[- ]?shape|tear[- ]?drop"),
  # \\bround\\b deliberately does NOT match "rounded", which is a modifier
  # ("broadly rounded posterior"), only the bare descriptor "round"
  Spherical      = paste0("spheric|spheroid|globular|coccoid|diplococc",
                          "|roundish|\\bround(?:ed)?\\b(?![^,]*(?:posterior|anterior|end|pole))",
                          "|almost uniform in diameter"),
  # typos are in the data and must be matched: ellispoidal, ellpitical
  Elliptical     = paste0("ellips|elliptic|episoid|ellispoid|ellpitic",
                          "|elipsoid|lozenge"),
  Rod            = paste0("rod|bacill|bascill|cylindric|tubular|fusiform",
                          "|elongate narrow|elongate-narrow|narrow|slender",
                          "|bacculiform|needle|acicular|long thin|straight"),
  Oval           = paste0("oval|ovoid|ovate|oviform|ovoidal|egg[- ]?shape",
                          "|oblong|conical|coniform|cucumiform|elongat")
)
SHAPE_LABEL <- c(Elliptical = "Elliptical", Oval = "Oval",
                 Ovocylindrical = "Ovocylindrical", Pyriform = "Pyriform",
                 Rod = "Rod-shaped", Spherical = "Spherical")

# Real shapes that are none of the six classes. Counted and reported as "Other"
# rather than dropped, because they are informative morphology, not noise.
OTHER_SHAPES <- paste0("lanceolate|horseshoe|crescent|c[- ]shaped|bean[- ]?shape",
                       "|kidney|clavate|club[- ]?shape|cordiform|sphenoidal",
                       "|rice[- ]?grain|reniform|spindle|comma|spiral",
                       "|spirilliform|dumbbell|cuneiform|wedge|triangular",
                       "|rectangular|semicircular|bow[- ]?shape|irregular")

# Descriptors that qualify a shape rather than being one. Ignored silently,
# because reporting them as unresolved would bury the real misses.
MODIFIERS <- paste0("attenuat|collaps|curv|bent|incurvat|truncat|blunt|point",
                    "|aberrant|abberant|deformed|approx|twisted|sharpen",
                    "|^broad$|^broadly|^long$|^short$|^wide$|wider in",
                    "|uniform in|uniform shape|variable in shape|^uniform$",
                    "|spores in pairs|nucleate|mononuclear|sporophorous",
                    "|double-walled|flagellated|polar filament|rounded,",
                    "|rounded (?:anterior|posterior)|posterior|anterior|pole",
                    "|concave|invaginat|depress|taper|equally|ends|end ",
                    "|similar to|smaller|larger|narrowed|inflated|constrict",
                    "|slightly more narrow|original description|not stated",
                    "|unknown|translation|from figures|dominant|sequence",
                    "|diameter|section|refractile|vacuole|nucleus|nuclei")

classify_shape_term <- function(term) {
  t <- tolower(deaccent(term))
  t <- gsub("\\(.*?\\)", " ", t, perl = TRUE)     # drop the condition
  t <- trimws(gsub("\\s+", " ", t))
  if (!nzchar(t)) return(c(class = NA, kind = "empty"))
  for (nm in names(SHAPE_RULES))
    if (grepl(SHAPE_RULES[[nm]], t, perl = TRUE))
      return(c(class = unname(SHAPE_LABEL[nm]), kind = "class"))
  if (grepl(OTHER_SHAPES, t, perl = TRUE)) return(c(class = "Other", kind = "other"))
  if (grepl(MODIFIERS, t, perl = TRUE)) return(c(class = NA, kind = "modifier"))
  c(class = NA, kind = "unresolved")
}

# Condition text in the parentheses of a shape descriptor, e.g. "(meiospore,
# stained)". Used to prefer a fresh, ordinary spore over a stained or aberrant
# one when a cell describes several.
term_condition <- function(term) {
  m <- regmatches(term, regexpr("\\(([^)]*)\\)", term, perl = TRUE))
  if (!length(m)) "" else tolower(gsub("[()]", "", m))
}

parse_shape_cell <- function(cell) {
  if (is.na(cell) || !nzchar(trimws(cell)))
    return(list(classes = character(0), primary = NA_character_,
                unresolved = character(0), kinds = character(0)))
  terms <- trimws(strsplit(as.character(cell), ";")[[1]])
  terms <- terms[nzchar(terms)]
  res <- lapply(terms, classify_shape_term)
  cls <- vapply(res, function(r) r[["class"]], character(1))
  knd <- vapply(res, function(r) r[["kind"]], character(1))
  cond <- vapply(terms, term_condition, character(1), USE.NAMES = FALSE)

  ok <- which(!is.na(cls) & cls != "Other")
  # prefer an unconditioned or fresh/live descriptor over a stained or
  # aberrant one: "rod (fresh); ovocylindrical (fixed)" is a rod
  if (length(ok)) {
    score <- rep(1L, length(ok))
    score[grepl("fresh|live", cond[ok])] <- 0L
    score[!nzchar(cond[ok])] <- 0L
    score[grepl("stain|fixed|collaps|preserv|aberrant", cond[ok])] <- 2L
    score[grepl(paste(CFG$deprioritise_types, collapse = "|"), cond[ok])] <-
      score[grepl(paste(CFG$deprioritise_types, collapse = "|"), cond[ok])] + 1L
    ok <- ok[order(score, ok)]
  }
  # the condition of the descriptor that was chosen, so the shape can be checked
  # against the condition the MEASUREMENTS came from
  norm_cond <- function(x) trimws(gsub("\\s+", " ", gsub("[^a-z ]", " ", tolower(x))))
  list(classes = unique(cls[ok]),
       primary = if (length(ok)) cls[ok[1]] else
                 if (any(cls == "Other", na.rm = TRUE)) "Other" else NA_character_,
       primary_cond = if (length(ok)) norm_cond(cond[ok[1]]) else NA_character_,
       n_terms = as.integer(length(terms)),
       unresolved = terms[knd == "unresolved"],
       kinds = knd)
}

## ============================================================================
## MEASUREMENT PARSING
## ============================================================================
# One cell -> one number, by this rule:
#   1. split on ";" into terms; each term is one reported measurement
#   2. rank terms by condition (fresh/live > unstated > fixed/stained) and
#      deprioritise macrospores, meiospores and other non-ordinary spore types
#   3. within the winning term take the FIRST number, which is the mean --
#      "4.03 +/- 0.28 (3.71 - 4.33 um)" is a mean of 4.03, not 0.28 or 3.71
#   4. if the term states an average explicitly ("avg: 4.4"), use that instead,
#      because "3.7 - 4.9 um (avg: 4.4 um)" leads with a range minimum
#   5. if the term is only a range ("0.5-0.8"), use the midpoint and flag it
NUM <- "[0-9]+(?:\\.[0-9]+)?"

parse_measure_term <- function(term) {
  t <- tolower(deaccent(term))
  t <- gsub("\u00b5|\u03bc", "u", t)
  # explicit average wins over anything else in the term
  m <- regmatches(t, regexpr(paste0("(?:avg|average|mean)\\s*[:=]?\\s*", NUM),
                             t, perl = TRUE))
  if (length(m))
    return(list(value = as.numeric(regmatches(m, regexpr(NUM, m, perl = TRUE))),
                basis = "stated average"))
  # a bare range with no leading mean
  stripped <- gsub("\\(.*?\\)", " ", t, perl = TRUE)
  rng <- regmatches(stripped, regexpr(paste0(NUM, "\\s*[-\u2013]\\s*", NUM),
                                      stripped, perl = TRUE))
  lead <- regmatches(stripped, regexpr(NUM, stripped, perl = TRUE))
  if (length(rng) && length(lead) && startsWith(trimws(stripped), sub("^\\s+", "", rng))) {
    n <- as.numeric(regmatches(rng, gregexpr(NUM, rng, perl = TRUE))[[1]])
    return(list(value = mean(n), basis = "range midpoint"))
  }
  if (length(lead)) return(list(value = as.numeric(lead), basis = "first value"))
  # nothing outside the parentheses; try inside, e.g. "(0.5-0.8 um diameter)"
  inside <- regmatches(t, gregexpr(NUM, t, perl = TRUE))[[1]]
  if (length(inside) >= 2) return(list(value = mean(as.numeric(inside[1:2])),
                                       basis = "range midpoint (parenthetical)"))
  if (length(inside) == 1) return(list(value = as.numeric(inside),
                                       basis = "first value (parenthetical)"))
  list(value = NA_real_, basis = "no number found")
}

cond_rank <- function(cond) {
  cond <- tolower(cond)
  base <- match("other", CFG$condition_priority)
  if (!nzchar(trimws(cond))) base <- match("unspecified", CFG$condition_priority)
  for (i in seq_along(CFG$condition_priority)) {
    k <- CFG$condition_priority[i]
    if (k == "unspecified") next
    if (grepl(k, cond)) { base <- i; break }
  }
  if (grepl(paste(CFG$deprioritise_types, collapse = "|"), cond)) base <- base + 10L
  base
}

parse_measure_cell <- function(cell) {
  if (is.na(cell) || !nzchar(trimws(cell)))
    return(list(value = NA_real_, basis = "empty", condition = "", n_terms = 0L))
  terms <- trimws(strsplit(as.character(cell), ";")[[1]])
  terms <- terms[nzchar(terms)]
  if (!length(terms)) return(list(value = NA_real_, basis = "empty",
                                  condition = "", n_terms = 0L))
  cond <- vapply(terms, function(t) {
    p <- regmatches(t, gregexpr("\\(([^)]*)\\)", t, perl = TRUE))[[1]]
    paste(tolower(gsub("[()]", "", p)), collapse = " ")
  }, character(1), USE.NAMES = FALSE)
  # condition words also appear bare: "Type I 5.10 fixed/5.75 fresh"
  cond <- paste(cond, ifelse(grepl("fresh|live|fixed|stained", tolower(terms)),
                             tolower(terms), ""))
  ord <- order(vapply(cond, cond_rank, integer(1)), seq_along(terms))
  for (i in ord) {
    p <- parse_measure_term(terms[i])
    if (!is.na(p$value))
      return(list(value = p$value, basis = p$basis,
                  condition = trimws(cond[i]), n_terms = length(terms),
                  index = i))
  }
  list(value = NA_real_, basis = "no number found", condition = "",
       n_terms = as.integer(length(terms)), index = NA_integer_)
}

# Every term in a cell, with its value and its condition label. Used to pair the
# length and width of the SAME reported spore.
# "5.10 fixed/5.75 fresh" is two measurements written with a slash instead of a
# semicolon. Left as one term the parser takes the first number -- the FIXED
# value -- while the condition scan sees the word "fresh" further along, so the
# entry is ranked as fresh but carries the fixed number. A slash is therefore
# treated as a separator when BOTH sides carry a number and at least one side
# names a condition; ranges and ratios ("1.5-2.0/0.8-1.0") are left alone.
COND_WORD <- "fresh|live|fixed|stained|preserv|collaps"
split_slash <- function(term) {
  # "+/-" and "+ / -" are the plus-minus sign, not a separator
  term <- gsub("\\+\\s*/\\s*-", "+-", term)
  if (!grepl("/", term)) return(term)
  parts <- trimws(strsplit(term, "/")[[1]])
  if (length(parts) < 2) return(term)
  has_num <- vapply(parts, function(x) grepl("[0-9]", x), logical(1))
  if (!all(has_num)) return(term)
  if (!any(grepl(COND_WORD, tolower(parts)))) return(term)
  parts
}

measure_terms <- function(cell) {
  if (is.na(cell) || !nzchar(trimws(cell))) return(NULL)
  terms <- trimws(strsplit(as.character(cell), ";")[[1]])
  terms <- unlist(lapply(terms, split_slash), use.names = FALSE)
  terms <- terms[nzchar(terms)]
  if (!length(terms)) return(NULL)
  cond <- vapply(terms, function(t) {
    p <- regmatches(t, gregexpr("\\(([^)]*)\\)", t, perl = TRUE))[[1]]
    paste(tolower(gsub("[()]", "", p)), collapse = " ")
  }, character(1), USE.NAMES = FALSE)
  cond <- paste(cond, ifelse(grepl("fresh|live|fixed|stained", tolower(terms)),
                             tolower(terms), ""))
  vals <- vapply(terms, function(t) parse_measure_term(t)$value, numeric(1),
                 USE.NAMES = FALSE)
  data.frame(idx = seq_along(terms), value = vals,
             cond = trimws(cond),
             cond_key = trimws(gsub("\\s+", " ", gsub("[^a-z ]", " ", cond))),
             rank = vapply(cond, cond_rank, integer(1), USE.NAMES = FALSE),
             stringsAsFactors = FALSE)
}

## ============================================================================
## POLAR TUBULE COILS
## ============================================================================
# COIL NOTES -- Excel date corruption
# Many "Polar Tubule Coils Range" cells were typed as "5-6", "8-9" and so on,
# and Excel autocorrected them into DATES. They now sit in the sheet as 5-digit
# serial numbers (44322) or ISO datetimes (2021-05-06), and the range is
# recoverable because Excel read the two numbers as month and day:
#   44322 -> 2021-05-06 -> "5 - 6"
# Every repair is reported so the source cells can be fixed and reformatted as
# TEXT. Ranges whose first number exceeds 12 ("13-14") were never valid dates
# and survived untouched, which is why only some of the column is affected.
EXCEL_EPOCH <- as.Date("1899-12-30")

repair_excel_date <- function(cell) {
  x <- trimws(as.character(cell))
  d <- NA
  if (grepl("^[0-9]{5}$", x)) {
    d <- EXCEL_EPOCH + as.integer(x)
  } else if (grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}", x)) {
    d <- as.Date(substr(x, 1, 10))
  }
  if (is.na(d)) return(list(text = x, repaired = FALSE))
  list(text = sprintf("%d - %d", as.integer(format(d, "%m")),
                      as.integer(format(d, "%d"))), repaired = TRUE)
}

# One number per cell: the stated average, else the midpoint of the range.
# Multi-value cells are ranked like the size columns, so "7.3 (meiospores);
# 7 (macrospores)" prefers the ordinary spore.
parse_coils <- function(avg_cell, range_cell) {
  a <- parse_measure_cell(avg_cell)
  if (!is.na(a$value)) return(list(value = a$value, basis = "stated average",
                                   repaired = FALSE))
  rep <- repair_excel_date(range_cell)
  if (!nzchar(rep$text)) return(list(value = NA_real_, basis = "empty",
                                     repaired = FALSE))
  if (grepl("uncoiled", rep$text, ignore.case = TRUE))
    return(list(value = 0, basis = "uncoiled", repaired = rep$repaired))
  terms <- trimws(strsplit(rep$text, ";")[[1]])
  terms <- terms[nzchar(terms)]
  for (t in terms) {
    stripped <- gsub("\\(.*?\\)", " ", t, perl = TRUE)
    n <- as.numeric(regmatches(stripped,
                               gregexpr(NUM, stripped, perl = TRUE))[[1]])
    n <- n[!is.na(n)]
    if (length(n) >= 2) return(list(value = mean(n[1:2]), basis = "range midpoint",
                                    repaired = rep$repaired))
    if (length(n) == 1) return(list(value = n[1], basis = "single value",
                                    repaired = rep$repaired))
  }
  list(value = NA_real_, basis = "no number found", repaired = rep$repaired)
}

## ------------------------------------------------------------------ load ----
if (!file.exists(CFG$db_path)) stop("database not found: ", CFG$db_path, call. = FALSE)
raw <- read_excel(CFG$db_path, sheet = CFG$sheet, col_types = "text",
                  .name_repair = "minimal")
message(sprintf("read %d rows x %d columns from sheet '%s'",
                nrow(raw), ncol(raw), CFG$sheet))

dat <- data.frame(
  row_xlsx  = seq_len(nrow(raw)) + 1L,
  species   = trimws(get_col(raw, CFG$name_col)),
  len_raw   = get_col(raw, CFG$len_col),
  wid_raw   = get_col(raw, CFG$wid_col),
  shape_raw = get_col(raw, CFG$shape_col),
  coil_avg_raw   = get_col(raw, CFG$coil_avg_col, required = FALSE),
  coil_range_raw = get_col(raw, CFG$coil_range_col, required = FALSE),
  stringsAsFactors = FALSE
)
for (v in c("len_raw", "wid_raw", "shape_raw", "coil_avg_raw", "coil_range_raw"))
  if (!is.null(dat[[v]])) dat[[v]][is.na(dat[[v]])] <- "" else dat[[v]] <- ""

cat("\n================ SPORE PARSING REPORT ================\n")
hdr <- which(norm_key(dat$species) == norm_key(CFG$name_col))
if (length(hdr)) {
  qc("!! %d repeated HEADER row(s) inside the data (xlsx row %s) -- dropped\n",
     length(hdr), paste(dat$row_xlsx[hdr], collapse = ", "))
  dat <- dat[-hdr, ]
}
dat <- dat[nzchar(dat$species), ]
dat$category <- classify_species(dat$species)
if (CFG$subset != "all") {
  n0 <- nrow(dat); dat <- dat[dat$category == CFG$subset, ]
  qc("subset = %s: %d of %d entries kept\n", CFG$subset, nrow(dat), n0)
}

## ------------------------------------------------------------- parse all ----
sh <- lapply(dat$shape_raw, parse_shape_cell)
dat$shape_primary <- vapply(sh, function(s) s$primary, character(1))
dat$shape_cond <- vapply(sh, function(s) s$primary_cond %||% NA_character_,
                         character(1))
dat$shape_nterms <- vapply(sh, function(s) as.integer(s$n_terms %||% 0L),
                           integer(1))
dat$shape_all     <- vapply(sh, function(s) paste(s$classes, collapse = "; "), character(1))
dat$n_shape_class <- vapply(sh, function(s) length(s$classes), integer(1))

lp <- lapply(dat$len_raw, parse_measure_cell)
wp <- lapply(dat$wid_raw, parse_measure_cell)
dat$length      <- vapply(lp, function(x) x$value, numeric(1))
dat$len_basis   <- vapply(lp, function(x) x$basis, character(1))
dat$len_cond    <- vapply(lp, function(x) x$condition, character(1))
dat$width       <- vapply(wp, function(x) x$value, numeric(1))
dat$wid_basis   <- vapply(wp, function(x) x$basis, character(1))
dat$wid_cond    <- vapply(wp, function(x) x$condition, character(1))

# Length and width are ranked independently, which can pull the two numbers
# from DIFFERENT reported spores: "5.5; 4.2 (fresh)" and "3.3; 7.0 (fresh)"
# gave a 4.2 x 7.0 spore that exists in neither source. When the two cells list
# the same number of terms they describe the same spores in the same order, so
# the width is re-read at the index chosen for the length.
# Return the shape class of the descriptor whose condition matches `want`,
# or NA when the shape cell has no such descriptor.
shape_terms_matching <- function(cell, want) {
  # `want` is the measurement condition; a descriptor matches when it does not
  # CONFLICT with it, using the same token comparison
  if (is.na(cell) || !nzchar(trimws(cell))) return(NA_character_)
  terms <- trimws(strsplit(as.character(cell), ";")[[1]])
  terms <- terms[nzchar(terms)]
  if (length(terms) < 2) return(NA_character_)
  for (t in terms) {
    if (conditions_conflict(term_condition(t), want)) next
    cl <- classify_shape_term(t)[["class"]]
    if (!is.na(cl) && cl != "Other") return(cl)
  }
  NA_character_
}

repaired <- list()
for (i in seq_len(nrow(dat))) {
  lt <- measure_terms(dat$len_raw[i]); wt <- measure_terms(dat$wid_raw[i])
  # Pairing used to be attempted only when BOTH cells had two or more terms, so
  # "3 (fresh); 2.1 (fixed)" against a single unlabelled width was never checked.
  if (is.null(lt) || is.null(wt) || (nrow(lt) < 2 && nrow(wt) < 2)) next
  lt <- lt[!is.na(lt$value), ]; wt <- wt[!is.na(wt$value), ]
  if (!nrow(lt) || !nrow(wt)) next
  li <- lt[order(lt$rank, lt$idx), ][1, ]
  wi_own <- wt[order(wt$rank, wt$idx), ][1, ]

  # 1. match on the condition label, which survives the two cells being
  #    written in different orders ("3 (stained); 1.7" vs "1.45; 3 (stained)")
  # an EMPTY condition is itself a label: it matches the unlabelled width term,
  # which is what makes "3 (stained); 1.7" pair with "1.45; 3 (stained)"
  hit <- wt[wt$cond_key == li$cond_key, ]
  # 2. otherwise, equal term counts mean the same spores in the same order
  if (!nrow(hit) && nrow(lt) == nrow(wt)) hit <- wt[wt$idx == li$idx, ]
  if (!nrow(hit) || hit$idx[1] == wi_own$idx) next

  repaired[[length(repaired) + 1]] <- data.frame(
    row_xlsx = dat$row_xlsx[i], species = dat$species[i],
    len_raw = dat$len_raw[i], wid_raw = dat$wid_raw[i],
    length = li$value, condition = li$cond,
    width_before = wi_own$value, width_after = hit$value[1],
    stringsAsFactors = FALSE)
  dat$width[i] <- hit$value[1]
  dat$wid_basis[i] <- paste(dat$wid_basis[i], "(paired to length term)")
  dat$wid_cond[i] <- hit$cond[1]
}
## ------------------------------- shape vs measurement condition ------------
# Shape, length and width were each chosen independently, so a point could
# combine a fresh shape with fixed dimensions. Where BOTH the shape cell and the
# measurement cells name conditions, they are now compared, and the shape term
# matching the measurement condition is used instead where one exists.
# Comparing the condition strings verbatim is useless: the measurement label
# carries the number and units too ("meiospore, fresh 7.2 +/- 0.3 (meiospore,
# fresh)"), so it never equals the shape label "meiospore fresh". The two are
# reduced to the tokens that actually identify a spore, and PREPARATION is
# compared separately from SPORE TYPE -- a mismatch only counts when both sides
# name the same kind of thing and disagree.
PREP_TOKENS <- c("fresh", "live", "fixed", "stained", "preserved", "collapsed",
                 "dried", "smear")
TYPE_TOKENS <- c("meiospore", "macrospore", "microspore", "megaspore",
                 "octospore", "monospore", "binucleate", "uninucleate",
                 "diplokaryotic", "primary", "secondary", "lanceolate")
cond_tokens <- function(x, vocab) {
  w <- unlist(strsplit(tolower(gsub("[^a-z ]", " ", x)), " +"))
  sort(unique(intersect(w, vocab)))
}
conditions_conflict <- function(a, b) {
  pa <- cond_tokens(a, PREP_TOKENS); pb <- cond_tokens(b, PREP_TOKENS)
  ta <- cond_tokens(a, TYPE_TOKENS); tb <- cond_tokens(b, TYPE_TOKENS)
  (length(pa) && length(pb) && !identical(pa, pb)) ||
    (length(ta) && length(tb) && !identical(ta, tb))
}

shape_fix <- list(); shape_mismatch <- list()
for (i in seq_len(nrow(dat))) {
  mc <- dat$len_cond[i] %||% ""
  sc <- dat$shape_cond[i] %||% ""
  if (!nzchar(trimws(mc)) || !nzchar(trimws(sc))) next
  if (!conditions_conflict(mc, sc)) next
  alt <- shape_terms_matching(dat$shape_raw[i], mc)
  if (!is.na(alt)) {
    shape_fix[[length(shape_fix) + 1]] <- data.frame(
      row_xlsx = dat$row_xlsx[i], species = dat$species[i],
      shape_raw = dat$shape_raw[i], measurement_condition = mc,
      shape_before = dat$shape_primary[i], shape_after = alt,
      stringsAsFactors = FALSE)
    dat$shape_primary[i] <- alt
    dat$shape_cond[i] <- mc
  } else {
    shape_mismatch[[length(shape_mismatch) + 1]] <- data.frame(
      row_xlsx = dat$row_xlsx[i], species = dat$species[i],
      shape_raw = dat$shape_raw[i], shape = dat$shape_primary[i],
      shape_condition = sc, measurement_condition = mc,
      stringsAsFactors = FALSE)
  }
}
if (length(shape_fix)) {
  sf <- do.call(rbind, shape_fix)
  qc("\n%d entr%s had the shape re-taken from the spore the measurements came from:\n",
     nrow(sf), ifelse(nrow(sf) == 1, "y", "ies"))
  show_tbl(sf[, c("species", "measurement_condition", "shape_before",
                  "shape_after")], n = 8, file = "QC_shape_repaired.tsv")
  write_qc(sf, "QC_shape_repaired.tsv")
}
if (length(shape_mismatch)) {
  sm <- do.call(rbind, shape_mismatch)
  qc("\n%d entr%s report the shape and the measurements for DIFFERENT spores,\n",
     nrow(sm), ifelse(nrow(sm) == 1, "y", "ies"))
  qc("  with no matching shape descriptor available")
  qc(if (identical(CFG$require_same_condition, TRUE))
       " -- excluded (CFG$require_same_condition = TRUE)\n" else
       " -- kept, flagged (set CFG$require_same_condition = TRUE to exclude)\n")
  show_tbl(sm[, c("species", "shape", "shape_condition",
                  "measurement_condition")], n = 8,
           file = "QC_shape_condition_mismatch.tsv")
  write_qc(sm, "QC_shape_condition_mismatch.tsv")
  if (identical(CFG$require_same_condition, TRUE))
    dat$shape_primary[dat$row_xlsx %in% sm$row_xlsx] <- NA_character_
}

if (length(repaired)) {
  rp <- do.call(rbind, repaired)
  qc("\n%d entr%s had length and width taken from different reported spores;\n",
     nrow(rp), if (nrow(rp) == 1) "y" else "ies")
  qc("  width re-read from the same term as the length:\n")
  show_tbl(rp[, c("species", "len_raw", "length", "width_before", "width_after")],
           n = 10, file = "QC_paired_measurements.tsv")
  write_qc(rp, "QC_paired_measurements.tsv")
}

cp <- Map(parse_coils, dat$coil_avg_raw,
          if (CFG$repair_excel_dates) dat$coil_range_raw
          else rep("", nrow(dat)))
dat$coils       <- vapply(cp, function(x) x$value, numeric(1))
dat$coil_basis  <- vapply(cp, function(x) x$basis, character(1))
dat$coil_fixed  <- vapply(cp, function(x) isTRUE(x$repaired), logical(1))

qc("entries                       : %d\n", nrow(dat))
qc("  shape recorded              : %d\n", sum(nzchar(trimws(dat$shape_raw))))
qc("  length recorded             : %d   parsed: %d\n",
   sum(nzchar(trimws(dat$len_raw))), sum(!is.na(dat$length)))
qc("  width recorded              : %d   parsed: %d\n",
   sum(nzchar(trimws(dat$wid_raw))), sum(!is.na(dat$width)))

has_coil_txt <- nzchar(trimws(dat$coil_avg_raw)) | nzchar(trimws(dat$coil_range_raw))
qc("  polar tubule coils recorded : %d   parsed: %d\n",
   sum(has_coil_txt), sum(!is.na(dat$coils)))

## --------------------------------------------- Excel date corruption --------
fixed_rows <- dat[dat$coil_fixed, ]
n_corrupt <- sum(grepl("^[0-9]{5}$|^[0-9]{4}-[0-9]{2}-[0-9]{2}",
                       trimws(dat$coil_range_raw)))
if (n_corrupt) {
  qc("\n!! %d '%s' cell(s) hold an Excel DATE, not a coil range.\n",
     n_corrupt, CFG$coil_range_col)
  qc("   Excel read \"5-6\" as 6 May and stored the serial number. The month and\n")
  qc("   day recover the original range; %d of these were used here (the rest\n",
     nrow(fixed_rows))
  qc("   already had a stated average). FIX AT SOURCE: format that column as TEXT.\n")
  cr <- dat[grepl("^[0-9]{5}$|^[0-9]{4}-[0-9]{2}-[0-9]{2}", trimws(dat$coil_range_raw)), ]
  cr$recovered <- vapply(cr$coil_range_raw,
                         function(x) repair_excel_date(x)$text, character(1))
  show_tbl(data.frame(species = cr$species, stored = cr$coil_range_raw,
                      recovered_range = cr$recovered,
                      stated_average = cr$coil_avg_raw), n = 10,
           file = "QC_excel_date_coils.tsv")
  write_qc(cr[, c("row_xlsx", "species", "coil_range_raw", "recovered",
                  "coil_avg_raw", "coils")], "QC_excel_date_coils.tsv")
}
coil_bad <- dat[has_coil_txt & is.na(dat$coils), ]
if (nrow(coil_bad)) {
  qc("\ncoil text present but NO number parsed: %d\n", nrow(coil_bad))
  show_tbl(data.frame(species = coil_bad$species, avg = coil_bad$coil_avg_raw,
                      range = coil_bad$coil_range_raw), n = 10,
           file = "QC_coils_unparsed.tsv")
  write_qc(coil_bad[, c("row_xlsx", "species", "coil_avg_raw", "coil_range_raw")],
           "QC_coils_unparsed.tsv")
}

## ------------------------------------------- unresolved shape descriptors ---
unres <- unlist(lapply(seq_along(sh), function(i) sh[[i]]$unresolved))
if (length(unres)) {
  u <- as.data.frame(table(descriptor = tolower(trimws(unres))),
                     stringsAsFactors = FALSE)
  u <- u[order(-u$Freq), ]; names(u)[2] <- "n"
  qc("\nUNRESOLVED shape descriptors: %d occurrences, %d distinct\n",
     sum(u$n), nrow(u))
  qc("  (matched no shape class, no 'other' shape and no modifier)\n")
  show_tbl(u, n = 25, file = "QC_unresolved_shape_terms.tsv")
  write_qc(u, "QC_unresolved_shape_terms.tsv")
  ur <- dat[vapply(sh, function(s) length(s$unresolved) > 0, logical(1)),
            c("row_xlsx", "species", "shape_raw", "shape_primary")]
  write_qc(ur, "QC_unresolved_shape_rows.tsv")
}

other_rows <- dat[!is.na(dat$shape_primary) & dat$shape_primary == "Other", ]
if (nrow(other_rows)) {
  qc("\nshapes outside the six classes ('Other'): %d entries\n", nrow(other_rows))
  show_tbl(other_rows[, c("species", "shape_raw")], n = 12,
           file = "QC_other_shapes.tsv")
  write_qc(other_rows[, c("row_xlsx", "species", "shape_raw")], "QC_other_shapes.tsv")
}

no_class <- dat[nzchar(trimws(dat$shape_raw)) & is.na(dat$shape_primary), ]
if (nrow(no_class)) {
  qc("\nshape text present but NO class assigned: %d entries\n", nrow(no_class))
  show_tbl(no_class[, c("species", "shape_raw")], n = 10,
           file = "QC_shape_no_class.tsv")
  write_qc(no_class[, c("row_xlsx", "species", "shape_raw")], "QC_shape_no_class.tsv")
}

multi <- dat[dat$n_shape_class > 1, ]
qc("\nentries naming more than one shape class: %d (multi_shape = \"%s\")\n",
   nrow(multi), CFG$multi_shape)
if (nrow(multi)) {
  show_tbl(multi[, c("species", "shape_raw", "shape_all", "shape_primary")],
           n = 10, file = "QC_multi_shape.tsv")
  write_qc(multi[, c("row_xlsx", "species", "shape_raw", "shape_all",
                     "shape_primary")], "QC_multi_shape.tsv")
}

## -------------------------------------------- measurement parse failures ----
for (side in c("len", "wid")) {
  raws <- dat[[paste0(side, "_raw")]]
  vals <- dat[[if (side == "len") "length" else "width"]]
  bad <- dat[nzchar(trimws(raws)) & is.na(vals), ]
  if (nrow(bad)) {
    qc("\n%s: %d cell(s) with text but NO parseable number\n",
       if (side == "len") "LENGTH" else "WIDTH", nrow(bad))
    show_tbl(data.frame(species = bad$species, cell = bad[[paste0(side, "_raw")]]),
             n = 10, file = sprintf("QC_%s_unparsed.tsv", side))
    write_qc(bad[, c("row_xlsx", "species", paste0(side, "_raw"))],
             sprintf("QC_%s_unparsed.tsv", side))
  }
}
derived <- dat[!is.na(dat$length) &
               (grepl("midpoint", dat$len_basis) | grepl("midpoint", dat$wid_basis)), ]
if (nrow(derived)) {
  qc("\n%d measurement(s) taken as a RANGE MIDPOINT (no mean was stated):\n",
     nrow(derived))
  show_tbl(derived[, c("species", "len_raw", "length", "wid_raw", "width")],
           n = 8, file = "QC_range_midpoints.tsv")
  write_qc(derived[, c("row_xlsx", "species", "len_raw", "length", "len_basis",
                       "wid_raw", "width", "wid_basis")], "QC_range_midpoints.tsv")
}
qc("\nmeasurement basis (length / width):\n")
print(table(length = dat$len_basis[nzchar(trimws(dat$len_raw))]), right = FALSE)
print(table(width = dat$wid_basis[nzchar(trimws(dat$wid_raw))]), right = FALSE)

## ------------------------------------------------------------ build plot ----
plt <- dat[!is.na(dat$length) & !is.na(dat$width) &
           dat$length > 0 & dat$width > 0 & !is.na(dat$shape_primary), ]

if (CFG$multi_shape == "drop") {
  n0 <- nrow(plt); plt <- plt[plt$n_shape_class <= 1, ]
  qc("\nmulti_shape = drop: %d multi-class entries removed\n", n0 - nrow(plt))
} else if (CFG$multi_shape == "all") {
  plt <- do.call(rbind, lapply(seq_len(nrow(plt)), function(i) {
    cl <- strsplit(plt$shape_all[i], "; ")[[1]]
    cl <- cl[nzchar(cl)]
    if (!length(cl)) cl <- plt$shape_primary[i]
    r <- plt[rep(i, length(cl)), ]; r$shape_primary <- cl; r
  }))
  qc("\nmulti_shape = all: expanded to %d points\n", nrow(plt))
}

qc("\nplottable (shape + length + width): %d entries\n", nrow(plt))
qc("  lost for want of a shape class : %d\n",
   sum(!is.na(dat$length) & !is.na(dat$width) & is.na(dat$shape_primary)))
qc("  lost for want of a measurement : %d\n",
   sum(!is.na(dat$shape_primary) & (is.na(dat$length) | is.na(dat$width))))

lev <- c("Elliptical", "Oval", "Ovocylindrical", "Pyriform", "Rod-shaped",
         "Spherical", "Other")
plt$shape <- factor(plt$shape_primary, levels = lev)
plt <- plt[!is.na(plt$shape), ]

cat("\n================ SHAPE CLASS SUMMARY ================\n")
summ <- do.call(rbind, lapply(levels(droplevels(plt$shape)), function(s) {
  d <- plt[plt$shape == s, ]
  data.frame(shape = s, n = nrow(d),
             length_median = round(median(d$length), 2),
             length_range = sprintf("%.1f-%.1f", min(d$length), max(d$length)),
             width_median = round(median(d$width), 2),
             width_range = sprintf("%.1f-%.1f", min(d$width), max(d$width)),
             ratio_median = round(median(d$length / d$width), 2),
             stringsAsFactors = FALSE)
}))
print(summ, row.names = FALSE, right = FALSE)
write_qc(summ, "shape_class_summary.tsv")
if (CFG$facet_by_category)
  print(table(shape = plt$shape, category = plt$category), right = FALSE)

## ------------------------------------------- length shorter than width ------
# By convention "length" is the long axis. A point above the 1:1 line is either
# a genuinely wider-than-long spore, a swapped pair of columns, or two values
# read from different conditions -- worth eyeballing before publication.
flip <- plt[plt$width > plt$length * 1.02, ]
if (nrow(flip)) {
  qc("\n!! %d entries have width > length (above the 1:1 line):\n", nrow(flip))
  show_tbl(data.frame(species = flip$species, shape = as.character(flip$shape),
                      length = flip$length, width = flip$width,
                      len_raw = flip$len_raw, wid_raw = flip$wid_raw),
           n = 12, file = "QC_width_exceeds_length.tsv")
  write_qc(flip[, c("row_xlsx", "species", "shape_primary", "length", "width",
                    "len_raw", "wid_raw", "len_cond", "wid_cond")],
           "QC_width_exceeds_length.tsv")
}

## ------------------------------------------------------- axis limits ------
xr <- range(plt$length); yr <- range(plt$width)
qc("\nobserved length %.1f-%.1f um, width %.1f-%.1f um\n", xr[1], xr[2], yr[1], yr[2])
if (!is.null(CFG$xlim) || !is.null(CFG$ylim)) {
  xl <- if (is.null(CFG$xlim)) xr else CFG$xlim
  yl <- if (is.null(CFG$ylim)) yr else CFG$ylim
  out <- plt[plt$length < xl[1] | plt$length > xl[2] |
             plt$width  < yl[1] | plt$width  > yl[2], ]
  if (nrow(out)) {
    qc("!! %d point(s) fall OUTSIDE the requested axis limits and are not drawn:\n",
       nrow(out))
    show_tbl(out[, c("species", "shape_primary", "length", "width")],
             n = 12, file = "QC_outside_axis_limits.tsv")
    write_qc(out[, c("row_xlsx", "species", "shape_primary", "length", "width",
                     "len_raw", "wid_raw")], "QC_outside_axis_limits.tsv")
  }
}

## ------------------------------------------------------------------ plot ----
# Okabe-Ito: distinguishable under deuteranopia and protanopia
pal <- c("Elliptical" = "#0072B2", "Oval" = "#009E73",
         "Ovocylindrical" = "#D55E00", "Pyriform" = "#E69F00",
         "Rod-shaped" = "#CC79A7", "Spherical" = "#56B4E9",
         "Other" = "#999999")
SYMBOLS <- c(16, 17, 15, 18, 8, 4, 3)   # one per level of `lev`
# Filled symbols for the coil overlay: 21-25 take a fill, so the marker reads
# as a solid block of colour rather than a thin outline. Only five fillable
# symbols exist, so the two least common classes fall back to 8 and 13, which
# take `colour` instead -- both aesthetics are mapped to the coil bin so all
# seven look the same. Order follows `lev`.
SYMBOLS_FILL <- c(Elliptical = 23, Oval = 24, Ovocylindrical = 8,
                  Pyriform = 25, `Rod-shaped` = 22, Spherical = 21, Other = 13)

base_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey93", linewidth = 0.3),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    axis.text = element_text(colour = "black"),
    legend.background = element_rect(fill = "white", colour = "grey40",
                                     linewidth = 0.3),
    legend.key = element_blank(),
    legend.key.size = unit(0.85, "lines"),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    plot.margin = margin(8, 10, 6, 8)
  )

p <- ggplot(plt, aes(x = length, y = width)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted",
              colour = "grey70", linewidth = 0.3) +
  labs(x = "Spore length (\u00b5m)", y = "Spore width (\u00b5m)") +
  base_theme

if (CFG$style == "colour_circles") {
  # every point a circle; shape class carried by colour alone
  p <- p + geom_point(aes(colour = shape), shape = 16, size = 2.1, alpha = 0.85) +
    scale_colour_manual(values = pal, name = "Spore shape", drop = TRUE) +
    theme(legend.position = c(0.985, 0.985), legend.justification = c(1, 1))

} else if (CFG$style == "colour_symbols") {
  p <- p + geom_point(aes(colour = shape, shape = shape), size = 1.9,
                      alpha = 0.85, stroke = 0.5) +
    scale_colour_manual(values = pal, name = "Spore shape", drop = TRUE) +
    scale_shape_manual(values = SYMBOLS, name = "Spore shape", drop = TRUE) +
    theme(legend.position = c(0.985, 0.985), legend.justification = c(1, 1))

} else if (CFG$style == "coil_heat") {
  # Shape class -> symbol, coil count -> colour. Two variables on one panel:
  # the symbol says what the spore looks like, the colour says how many coils
  # its polar tubule has. Entries with no coil data are NOT dropped -- they are
  # drawn first, in grey, so the extent of the missing data is visible rather
  # than hidden by a smaller point count.
  cl <- CFG$coil_breaks
  clab <- vapply(seq_len(length(cl) - 1), function(i) {
    lo <- cl[i]; hi <- cl[i + 1]
    if (is.infinite(hi)) sprintf("%g+", lo) else sprintf("%g\u2013%g", lo, hi)
  }, character(1))
  plt$coil_bin <- cut(plt$coils, breaks = cl, labels = clab,
                      right = TRUE, include.lowest = TRUE)
  known <- plt[!is.na(plt$coil_bin), ]
  unknown <- plt[is.na(plt$coil_bin), ]
  coil_cols <- setNames(grDevices::hcl.colors(length(clab), "viridis",
                                              rev = TRUE), clab)

  p <- p +
    geom_point(data = unknown, aes(shape = shape), colour = "grey80",
               fill = "grey85", size = 1.8, alpha = 0.8, stroke = 0.3) +
    geom_point(data = known, aes(shape = shape, fill = coil_bin,
                                 colour = coil_bin),
               size = 2.3, alpha = 0.95, stroke = 0.35) +
    scale_fill_manual(values = coil_cols, name = "Polar tubule coils",
                      drop = FALSE, na.translate = FALSE,
                      guide = guide_legend(order = 1,
                                           override.aes = list(shape = 21,
                                                               colour = "grey30",
                                                               size = 2.6))) +
    scale_colour_manual(values = coil_cols, guide = "none",
                        na.translate = FALSE) +
    scale_shape_manual(values = SYMBOLS_FILL[lev], name = "Spore shape",
                       drop = TRUE,
                       guide = guide_legend(order = 2,
                                            override.aes = list(fill = "grey60",
                                                                colour = "grey30",
                                                                size = 2.4))) +
    theme(legend.position = "right", legend.justification = "center",
          legend.box.background = element_blank(),
          legend.background = element_blank())
  qc("\ncoil overlay: %d point(s) with coil data, %d drawn in grey (no coil data)\n",
     nrow(known), nrow(unknown))
  print(table(`coils` = plt$coil_bin, useNA = "ifany"), right = FALSE)
  cat("\ncoils by shape class (median, n with coil data):\n")
  cs <- do.call(rbind, lapply(levels(droplevels(plt$shape)), function(sh) {
    d <- plt[plt$shape == sh & !is.na(plt$coils), ]
    if (!nrow(d)) return(NULL)
    data.frame(shape = sh, n_with_coils = nrow(d),
               coils_median = round(median(d$coils), 1),
               coils_range = sprintf("%g-%g", min(d$coils), max(d$coils)),
               stringsAsFactors = FALSE)
  }))
  print(cs, row.names = FALSE, right = FALSE)
  write_qc(cs, "coils_by_shape.tsv")

} else {
  stop("unknown CFG$style: ", CFG$style, call. = FALSE)
}

if (!is.null(CFG$xlim) || !is.null(CFG$ylim))
  p <- p + coord_cartesian(xlim = CFG$xlim, ylim = CFG$ylim)
if (CFG$facet_by_category) {
  p <- p + facet_wrap(~ category) +
    theme(legend.position = "right", legend.justification = "center")
}

ggsave(file.path(CFG$outdir, "fig_spore_size_shape.pdf"), p,
       width = CFG$fig_width, height = CFG$fig_height,
       device = if (capabilities("cairo")) cairo_pdf else pdf)
ggsave(file.path(CFG$outdir, "fig_spore_size_shape.png"), p,
       width = CFG$fig_width, height = CFG$fig_height, dpi = 600)

write_qc(dat[, c("row_xlsx", "species", "category", "shape_raw", "shape_primary",
                 "shape_all", "n_shape_class", "len_raw", "length", "len_basis",
                 "len_cond", "wid_raw", "width", "wid_basis", "wid_cond")],
         "per_entry_spore_data.tsv")

qc("\nwritten to %s/\n  fig_spore_size_shape.pdf / .png\n  shape_class_summary.tsv\n  per_entry_spore_data.tsv\n  QC_*.tsv\n",
   CFG$outdir)
