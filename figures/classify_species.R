###############################################################################
# classify_species.R -- THE canonical named/provisional rule.
#
# This file is the single source of truth for the R side. Every R script that
# classifies microsporidia names must source() this file and must NOT define
# its own classify_species() / is_provisional_name(). The Python twin is
# classify_species.py; check_classifier_parity.py verifies the two agree on
# every Species Name in the database.
#
# The rule (identical to the CANON_* block formerly in assemble_figure_8.R):
#   provisional if:  empty name; uninformative prefix (uncultured/unnamed/
#   unknown/undescribed/mutated/...); fewer than 2 informative tokens;
#   sp./spp. among the first three tokens; purely numeric second token.
#   Otherwise named. "Microsporidium" stays a valid collective-group genus,
#   and nomenclatural acts ("sp. n.", "n. gen. sp.") are stripped BEFORE the
#   sp. test so "Ameson earli sp. n." is named.
#
# ONE DELIBERATE FIX vs the old assembler block: tokens were previously
# passed through setdiff(), which also DE-DUPLICATED them as a side effect,
# so a name with a repeated token could classify differently in R than in
# Python (which kept duplicates). Duplicates are now kept in both languages.
###############################################################################

CANON_deaccent <- function(x) {
  if (requireNamespace("stringi", quietly = TRUE))
    stringi::stri_trans_general(x, "Latin-ASCII")
  else iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT", sub = "")
}

RX_NOM_1 <- "\\b(?:nov|n)\\.(?:\\s*,?\\s*(?:n\\.\\s*)?(?:gen|sp|spp|comb)\\.?)+"
RX_NOM_2 <- paste0("\\b(?:gen|sp|spp|ssp|subsp|var|comb|stat)\\.?\\s*,?\\s*",
                   "(?:et\\s*|and\\s*|&\\s*)?(?:nov\\.?(?![a-z])|n\\.)")
UNINF <- paste0("^(?:uncultured|unidentified|unclassified|environmental|",
                "microsporidian?|fungal|eukaryot\\w*|",
                "unnamed|unknown|undescribed|mutated)\\b")

classify_species <- function(name, ...) {
  vapply(name, function(nm) {
    if (is.na(nm) || !nzchar(trimws(nm))) return("provisional")
    s <- tolower(CANON_deaccent(gsub("\u00a0", " ", nm, fixed = TRUE)))
    s <- gsub("\\(.*?\\)", " ", s, perl = TRUE)
    s <- gsub(RX_NOM_1, " ", s, perl = TRUE, ignore.case = TRUE)
    s <- gsub(RX_NOM_2, " ", s, perl = TRUE, ignore.case = TRUE)
    s <- gsub("\\b(?:cf|aff|nr)\\.?\\s+", " ", s, perl = TRUE)
    s <- gsub("[^a-z0-9\\s]", " ", s, perl = TRUE)
    tok <- strsplit(trimws(s), "\\s+")[[1]]
    # drop stopword tokens but KEEP duplicates and order (see header note)
    tok <- tok[!tok %in% c("nov", "gen", "et", "and", "comb", "stat", "")]
    if (length(tok) == 0) return("provisional")
    if (grepl(UNINF, paste(tok, collapse = " "), perl = TRUE))
      return("provisional")
    if (length(tok) < 2) return("provisional")
    if (any(tok[seq_len(min(3, length(tok)))] %in% c("sp", "spp")))
      return("provisional")
    if (grepl("^[0-9]+$", tok[2])) return("provisional")
    "named"
  }, character(1), USE.NAMES = FALSE)
}

is_provisional_name <- function(x) classify_species(x) == "provisional"
