#!/usr/bin/env Rscript
###############################################################################
# Figure -- site of infection (tissue tropism)
#
# Number of species reported from each tissue system, parsed out of the
# free-text "Site of Infection" column and split by named vs provisional.
#
# WHY THE OLD APPROACH NEEDED SO MANY MANUAL OVERRIDES
#   The previous script collected EVERY tissue keyword found anywhere in a
#   string, so any phrase that located one tissue relative to another produced
#   two categories, and each case had to be hardcoded:
#     "cells of hindgut near Malpighian tubules"  -> gut + Malpighian
#     "longitudinal muscles of body wall"         -> muscle + integument
#     "subcutaneous layer of muscles"             -> integument + muscle
#     "connective tissue of the coelomic cavity"  -> connective + body cavity
#   Two rules replace all of those overrides:
#     1. HEAD-NOUN RULE. Within one term the infected tissue is the FIRST
#        keyword; anything after it is describing where that tissue sits.
#     2. Positional adjectives (subcutaneous, longitudinal, skeletal, ...) are
#        not category keywords on their own, only in fixed phrases such as
#        "subcutaneous connective tissue", so "subcutaneous layer of muscles"
#        resolves to muscle rather than to integument.
#   The twelve overrides in the old script are checked as test cases at the end
#   of this file: set CFG$run_selftest = TRUE.
#
# A species may genuinely infect several tissues. Terms separated by ";" (and
# by "," where each fragment names a tissue) are parsed independently and the
# species is counted once in each category, exactly as for the country map.
#
# Usage:  Rscript fig_tissue_sites.R [path/to/database.xlsx]
###############################################################################

## ---------------------------------------------------------------- config ----
CFG <- list(
  db_path   = "Microsporidia_Characteristics_Database_merge_pro4.xlsx",
  sheet     = "Actively Updated Masterlist",
  name_col  = "Species Name",
  site_col  = "Site of Infection",

  subset       = "all",     # "all" | "named" | "provisional"
  min_n        = 1,         # categories with fewer species are pooled as "Other"
  # "category"  keep Systemic / general as its own bar (default)
  # "exclude"   drop the systemic term; species naming only a systemic site are
  #             then lost from the panel entirely
  # A species naming BOTH a systemic and a specific site keeps the specific one
  # either way -- this switch controls only the systemic term itself.
  systemic     = "category",
  show_percent = FALSE,     # label bars with % of entries as well as n
  run_selftest = TRUE,      # check the head-noun rule against known phrasings

  outdir     = "tissue_output",
  fig_width  = 7.6,
  fig_height = 5.8
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
qc <- function(...) cat(sprintf(...), sep = "")

## --------------------------------------------------------------- helpers ----
deaccent <- function(x) {
  if (has_stringi) stringi::stri_trans_general(x, "Latin-ASCII")
  else iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT", sub = "")
}
norm_key <- function(x) gsub("[^a-z0-9]", "", tolower(deaccent(as.character(x))))

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
## TISSUE GAZETTEER
## ============================================================================
# alias -> category. Aliases are matched as whole words; the EARLIEST match in
# a term wins (head-noun rule), with the longest alias winning at equal start.
# Misspellings present in the data are included deliberately: Malphigian,
# hyopdermal, conective, mesentary, hematocoel, oenocystes, peritenoum,
# fatbody, silkglands, traceae, havily, gult.
TISSUE <- c(
  # --- fat body / adipose -------------------------------------------------
  "fat body" = "Fat body / adipose", "fatbody" = "Fat body / adipose",
  "fat bodies" = "Fat body / adipose", "fat cells" = "Fat body / adipose",
  "fat tissue" = "Fat body / adipose", "adipose" = "Fat body / adipose",
  "adipose body" = "Fat body / adipose", "adipocyte" = "Fat body / adipose",
  "oenocyte" = "Fat body / adipose", "oenocystes" = "Fat body / adipose",
  "storage cells" = "Fat body / adipose", "storages cells" = "Fat body / adipose",
  "chloragocyte" = "Fat body / adipose",

  # --- muscle --------------------------------------------------------------
  "muscle" = "Muscle", "muscles" = "Muscle", "muscular" = "Muscle",
  "musculature" = "Muscle", "myofibril" = "Muscle", "myofiber" = "Muscle",
  "myocyte" = "Muscle", "sarcoplasm" = "Muscle", "sarcolemma" = "Muscle",
  "sacrolemma" = "Muscle", "skeletal fibers" = "Muscle",
  "muscle fibers" = "Muscle", "myoblast" = "Muscle",

  # --- gut / digestive -----------------------------------------------------
  "gut" = "Gut / digestive tract", "gut wall" = "Gut / digestive tract",
  "gut-wall" = "Gut / digestive tract", "gult" = "Gut / digestive tract",
  "midgut" = "Gut / digestive tract", "mid-gut" = "Gut / digestive tract",
  "hindgut" = "Gut / digestive tract", "foregut" = "Gut / digestive tract",
  "intestin" = "Gut / digestive tract", "enterocyte" = "Gut / digestive tract",
  "digestive tract" = "Gut / digestive tract",
  "digestive tube" = "Gut / digestive tract",
  "digestive tissue" = "Gut / digestive tract",
  "digestive epithelium" = "Gut / digestive tract",
  "gastrointestinal" = "Gut / digestive tract",
  "gastric caeca" = "Gut / digestive tract", "gastric ceca" = "Gut / digestive tract",
  "gastric serosa" = "Gut / digestive tract", "caecal mass" = "Gut / digestive tract",
  "pyloric" = "Gut / digestive tract", "cecae" = "Gut / digestive tract",
  "stomach" = "Gut / digestive tract", "ventriculus" = "Gut / digestive tract",
  "colon" = "Gut / digestive tract", "cloaca" = "Gut / digestive tract",
  "goblet cells" = "Gut / digestive tract", "villus" = "Gut / digestive tract",
  "alimentary" = "Gut / digestive tract", "peritrophic" = "Gut / digestive tract",

  # --- Malpighian tubules --------------------------------------------------
  "malpighian" = "Malpighian tubules", "malphigian" = "Malpighian tubules",
  "malpighan" = "Malpighian tubules",

  # --- reproductive --------------------------------------------------------
  "ovary" = "Reproductive / gonad", "ovaries" = "Reproductive / gonad",
  "ovarian" = "Reproductive / gonad", "oocyte" = "Reproductive / gonad",
  "ova" = "Reproductive / gonad", "egg" = "Reproductive / gonad",
  "eggs" = "Reproductive / gonad", "egg sacs" = "Reproductive / gonad",
  "gonad" = "Reproductive / gonad", "gonadal" = "Reproductive / gonad",
  "testis" = "Reproductive / gonad", "testes" = "Reproductive / gonad",
  "testicular" = "Reproductive / gonad", "germ cells" = "Reproductive / gonad",
  "germline" = "Reproductive / gonad", "nurse cells" = "Reproductive / gonad",
  "follicular" = "Reproductive / gonad", "uterus" = "Reproductive / gonad",
  "vas deferens" = "Reproductive / gonad", "genital duct" = "Reproductive / gonad",
  "embryo" = "Reproductive / gonad", "reproductive" = "Reproductive / gonad",

  # --- nervous -------------------------------------------------------------
  "nervous system" = "Nervous system", "nerve" = "Nervous system",
  "nerves" = "Nervous system", "ganglia" = "Nervous system",
  "ganglion" = "Nervous system", "brain" = "Nervous system",
  "spinal cord" = "Nervous system", "meninges" = "Nervous system",
  "cerebrospinal" = "Nervous system", "neural" = "Nervous system",

  # --- integument ----------------------------------------------------------
  "hypodermis" = "Integument / epidermis", "hypoderm" = "Integument / epidermis",
  "hypodermal" = "Integument / epidermis", "hyopdermal" = "Integument / epidermis",
  "epidermis" = "Integument / epidermis", "epidermal" = "Integument / epidermis",
  "integument" = "Integument / epidermis", "cuticle" = "Integument / epidermis",
  "carapace" = "Integument / epidermis", "derma cells" = "Integument / epidermis",
  "skin" = "Integument / epidermis", "body wall" = "Integument / epidermis",
  "body-wall" = "Integument / epidermis",

  # --- connective ----------------------------------------------------------
  "connective" = "Connective tissue", "conective" = "Connective tissue",
  "mesenchyme" = "Connective tissue", "mesenter" = "Connective tissue",
  "mesentary" = "Connective tissue", "mensenteries" = "Connective tissue",
  "fibroblast" = "Connective tissue", "stroma" = "Connective tissue",
  "lamina propria" = "Connective tissue", "submucosa" = "Connective tissue",
  "peritoneal" = "Connective tissue",
  "subcutaneous tissue" = "Connective tissue",
  "subcuticular connective" = "Connective tissue",
  "subcutaneous connective" = "Connective tissue",
  "bone marrow" = "Connective tissue",

  # --- body fluids / cavity ------------------------------------------------
  "haemolymph" = "Body fluid / cavity", "hemolymph" = "Body fluid / cavity",
  "haemocoel" = "Body fluid / cavity", "hemocoel" = "Body fluid / cavity",
  "haemocele" = "Body fluid / cavity", "hemocoele" = "Body fluid / cavity",
  "hematocoel" = "Body fluid / cavity", "blood" = "Body fluid / cavity",
  "blood sinus" = "Body fluid / cavity", "hemal" = "Body fluid / cavity",
  "body cavity" = "Body fluid / cavity", "coelom" = "Body fluid / cavity",
  "coelomic" = "Body fluid / cavity", "abdominal cavity" = "Body fluid / cavity",
  "visceral cavity" = "Body fluid / cavity", "general cavity" = "Body fluid / cavity",
  "haemal" = "Body fluid / cavity",

  # --- immune cells --------------------------------------------------------
  "haemocyte" = "Immune cells", "hemocyte" = "Immune cells",
  "coelomocyte" = "Immune cells", "lymphocyte" = "Immune cells",
  "macrophage" = "Immune cells", "phagocyte" = "Immune cells",
  "leukocyte" = "Immune cells", "hematopoietic" = "Immune cells",
  "haematopoietic" = "Immune cells", "spleen" = "Immune cells",
  "lymph node" = "Immune cells", "podocyte" = "Immune cells",

  # --- hepatic -------------------------------------------------------------
  "liver" = "Liver / hepatopancreas", "hepatic" = "Liver / hepatopancreas",
  "hepatopancrea" = "Liver / hepatopancreas",
  "digestive gland" = "Liver / hepatopancreas",
  "bile duct" = "Liver / hepatopancreas", "bilary" = "Liver / hepatopancreas",
  "gall bladder" = "Liver / hepatopancreas", "gallbladder" = "Liver / hepatopancreas",

  # --- excretory -----------------------------------------------------------
  "kidney" = "Excretory / renal", "renal" = "Excretory / renal",
  "urinary" = "Excretory / renal", "nephridia" = "Excretory / renal",
  "glomeruli" = "Excretory / renal", "antennal gland" = "Excretory / renal",
  "excretory" = "Excretory / renal", "excetory" = "Excretory / renal",

  # --- respiratory ---------------------------------------------------------
  "gill" = "Respiratory / gills", "gills" = "Respiratory / gills",
  "trachea" = "Respiratory / gills", "tracheal" = "Respiratory / gills",
  "traceae" = "Respiratory / gills", "lung" = "Respiratory / gills",
  "pseudobranch" = "Respiratory / gills",

  # --- secretory glands ----------------------------------------------------
  "silk gland" = "Silk / salivary glands", "silkgland" = "Silk / salivary glands",
  "salivary" = "Silk / salivary glands", "gland cells" = "Silk / salivary glands",
  "pancreas" = "Silk / salivary glands", "adrenal" = "Silk / salivary glands",

  # --- circulatory ---------------------------------------------------------
  "heart" = "Heart / circulatory", "cardiac" = "Heart / circulatory",
  "blood vessel" = "Heart / circulatory", "endothelial" = "Heart / circulatory",
  "endothelium" = "Heart / circulatory",

  # --- eye -----------------------------------------------------------------
  "eye" = "Eye", "cornea" = "Eye", "conjunctiva" = "Eye", "iris" = "Eye",
  "choroid" = "Eye",

  # --- systemic ------------------------------------------------------------
  "systemic" = "Systemic / general", "general infection" = "Systemic / general",
  "generalized" = "Systemic / general", "all organs" = "Systemic / general",
  "most organs" = "Systemic / general", "visceral organs" = "Systemic / general",
  "internal organs" = "Systemic / general", "other organs" = "Systemic / general",
  "host tissues" = "Systemic / general", "disseminated" = "Systemic / general",
  "throughout" = "Systemic / general", "body completely filled" = "Systemic / general",
  "somatic tissue" = "Systemic / general", "whole body" = "Systemic / general",

  # --- xenoma --------------------------------------------------------------
  "xenoma" = "Xenoma / cyst", "cyst" = "Xenoma / cyst",
  "metacercaria" = "Xenoma / cyst", "sporocyst" = "Xenoma / cyst",

  # --- added after reviewing the unresolved-term report ---------------------
  "fatty body" = "Fat body / adipose",
  "rodlet cells" = "Immune cells",
  "mensetery" = "Connective tissue",
  "verson" = "Silk / salivary glands",
  "germline tissues" = "Reproductive / gonad",
  "ileum" = "Gut / digestive tract",
  "jejunal" = "Gut / digestive tract",
  "caecum" = "Gut / digestive tract",
  "caeca" = "Gut / digestive tract",
  "ceca" = "Gut / digestive tract",
  "cecal" = "Gut / digestive tract",
  "rectum" = "Gut / digestive tract",
  "oesophagus" = "Gut / digestive tract",
  "esophagus" = "Gut / digestive tract",
  "proventriculus" = "Gut / digestive tract",
  "gastrodermis" = "Gut / digestive tract",
  "diverticulum" = "Gut / digestive tract",
  "mucosa" = "Gut / digestive tract",
  "gut epithelium" = "Gut / digestive tract",
  "nervous" = "Nervous system",
  "neuron" = "Nervous system",
  "neuronal" = "Nervous system",
  "motor neurons" = "Nervous system",
  "hanglia" = "Nervous system",
  "subesophageal" = "Nervous system",
  "photoreceptor" = "Nervous system",
  "sensory bulb" = "Nervous system",
  "receptors" = "Nervous system",
  "cranial" = "Nervous system",
  "oviduct" = "Reproductive / gonad",
  "vitellaria" = "Reproductive / gonad",
  "vitelline" = "Reproductive / gonad",
  "ovocyte" = "Reproductive / gonad",
  "ovotestis" = "Reproductive / gonad",
  "ovules" = "Reproductive / gonad",
  "ovum" = "Reproductive / gonad",
  "gonods" = "Reproductive / gonad",
  "spermatogonia" = "Reproductive / gonad",
  "seminal vesicle" = "Reproductive / gonad",
  "germaria" = "Reproductive / gonad",
  "germ balls" = "Reproductive / gonad",
  "sex cells" = "Reproductive / gonad",
  "sex organs" = "Reproductive / gonad",
  "genital organs" = "Reproductive / gonad",
  "tropharia" = "Reproductive / gonad",
  "uterine" = "Reproductive / gonad",
  "copulation sac" = "Reproductive / gonad",
  "female glands" = "Reproductive / gonad",
  "myocardium" = "Heart / circulatory",
  "pericardium" = "Heart / circulatory",
  "pericardial" = "Heart / circulatory",
  "epicardium" = "Heart / circulatory",
  "vessel walls" = "Heart / circulatory",
  "dermal" = "Integument / epidermis",
  "epidermia" = "Integument / epidermis",
  "tegument" = "Integument / epidermis",
  "fins" = "Integument / epidermis",
  "fin" = "Integument / epidermis",
  "skeleton" = "Integument / epidermis",
  "snail mantle" = "Integument / epidermis",
  "head capsule" = "Integument / epidermis",
  "seam cells" = "Integument / epidermis",
  "parenchym" = "Connective tissue",
  "paranchyma" = "Connective tissue",
  "parenchimatous" = "Connective tissue",
  "peritoneum" = "Connective tissue",
  "peritenoum" = "Connective tissue",
  "peritoneal" = "Connective tissue",
  "septa" = "Connective tissue",
  "inter-membrane" = "Connective tissue",
  "diaphram" = "Connective tissue",
  "mesentery" = "Connective tissue",
  "hemoctye" = "Immune cells",
  "nephrocytes" = "Immune cells",
  "archaeocytes" = "Immune cells",
  "amphidiskoblasts" = "Immune cells",
  "spongioblasts" = "Immune cells",
  "adipo-phagocytic" = "Immune cells",
  "trachael" = "Respiratory / gills",
  "peritracheal" = "Respiratory / gills",
  "tracheole" = "Respiratory / gills",
  "air sac" = "Respiratory / gills",
  "nasal epithelium" = "Respiratory / gills",
  "nasal mucosa" = "Respiratory / gills",
  "nasal sinus" = "Respiratory / gills",
  "labial glands" = "Silk / salivary glands",
  "pharyngeal glands" = "Silk / salivary glands",
  "verson's gland" = "Silk / salivary glands",
  "versons gland" = "Silk / salivary glands",
  "thyroid" = "Silk / salivary glands",
  "parathyroid" = "Silk / salivary glands",
  "protonephridial" = "Excretory / renal",
  "maxillary gland" = "Excretory / renal",
  "body fluid" = "Body fluid / cavity",
  "peritoneal cavity" = "Body fluid / cavity",
  "mouth cavity" = "Body fluid / cavity",
  "conjuctiva" = "Eye",
  "other tissues" = "Systemic / general",
  "general in adults" = "Systemic / general",
  "infection general" = "Systemic / general",
  "all segments" = "Systemic / general",
  "swimming organs" = "Systemic / general"
)

# Coverings and linings. They name a real tissue, but when an organ is also
# named in the same term the organ is the site: "under peritoneum of liver" is
# a liver infection, not a connective-tissue one. Applied only when nothing
# else in the term matched.
WEAK_ALIASES <- c("peritoneum", "peritenoum", "peritoneal", "mucosa",
                  "submucosa", "serosa", "stroma", "septa", "diaphram",
                  "inter-membrane", "somatic tissue", "other tissues")

# Adjectives of position or quality. They qualify a tissue, they are not one,
# so "subcutaneous layer of muscles" is muscle and "skeletal muscle of the
# abdominal cavity" is muscle. Fixed phrases such as "subcutaneous connective
# tissue" are in TISSUE above and win because they are longer.
POSITIONAL <- paste0("subcutaneous|subcuticular|longitudinal|skeletal|striated",
                     "|smooth|abdominal|thoracic|lateral|dorsal|ventral",
                     "|anterior|posterior|superficial|deep|proximal|distal",
                     "|primary|secondary|main|rarely|occassionally|severe",
                     "|upper|lower|mid|large|small|locomotor|flight|trunk")

# Body regions and non-tissue words: ignored silently, so they do not crowd out
# the genuinely unrecognised terms in the report.
IGNORE_WORDS <- c(
  "abdomen", "thorax", "limbs", "appendages", "head", "body", "trunk",
  "larva", "larvae", "adult", "adults", "imago", "queens", "male", "female",
  "type 1", "type 2", "type i", "type ii", "c1", "c2", "c3",
  "epithelium", "epithelia", "epithelial", "epithelial cells", "cells",
  "tissue", "tissues", "site", "sites", "organ", "organs", "wall", "walls",
  "layer", "surface", "cavity", "region", "parenchyma", "cytoplasm",
  "nucleus", "intranuclear", "primary cells", "host", "hosts", "unknown",
  "not resolved", "not recovered", "not reported", "uncertain", "various",
  "several", "other", "others", "rare", "common", "first", "second",
  "and", "or", "of", "in", "at", "the", "with", "within", "including",
  "wasp host", "archigregarine host", "gregarine cytoplasm", "flagellum",
  "around centre of animal", "tissues around basal bulb",
  "tissues in protocephalic area", "binucleate cells",
  "binucleate cells (pre-pansporocyst stage) of the myxosporean",
  "host-cell cytoplasm", "transovarian transmission was mentioned",
  "site/tissue not recovered from accessible sources", "endodermal",
  "tissue site not resolved in the description", "2)", "posterior body",
  "excluding", "except", "etc", "vessels", "lesion", "fresh", "stained"
)
# Matched as exact strings rather than compiled into one alternation: entries
# like "2)" are literal text and would break a regex with an unmatched paren.
is_ignorable <- function(x) x %in% IGNORE_WORDS

## -------------------------------------------------------------- matching ----
# Aliases sorted longest first so "subcutaneous connective tissue" beats
# "connective", and "gastric caeca" beats "gut".
# An alias whose only distinguishing feature is whitespace is a trap: aliases
# are matched against a trimmed, space-collapsed string, so "maxillary gland "
# matched only when other text followed, and the same organ was classified two
# different ways depending on the rest of the sentence. Aliases are therefore
# trimmed, and any alias mapping to two categories is a hard error rather than
# a silent coin toss.
names(TISSUE) <- trimws(names(TISSUE))
dup <- names(TISSUE)[duplicated(names(TISSUE))]
if (length(dup)) {
  conflict <- vapply(unique(dup), function(a)
    length(unique(unname(TISSUE[names(TISSUE) == a]))) > 1, logical(1))
  if (any(conflict))
    stop("alias(es) mapped to more than one tissue category: ",
         paste(sprintf("%s -> {%s}", names(conflict)[conflict],
                       vapply(names(conflict)[conflict], function(a)
                         paste(unique(unname(TISSUE[names(TISSUE) == a])),
                               collapse = ", "), character(1))),
               collapse = "; "), call. = FALSE)
  TISSUE <- TISSUE[!duplicated(names(TISSUE))]
}

ALIASES <- names(TISSUE)[order(-nchar(names(TISSUE)))]

# Find the FIRST tissue keyword in a term. Everything after it locates that
# tissue rather than naming a second one -- the head-noun rule.
first_tissue <- function(term) {
  s <- tolower(deaccent(term))
  s <- gsub("\\(.*?\\)", " ", s, perl = TRUE)      # host names, "(primary site)"
  s <- gsub("[^a-z0-9 /-]", " ", s, perl = TRUE)
  s <- gsub("\\s+", " ", trimws(s))
  if (!nzchar(s)) return(list(category = NA_character_, alias = NA_character_,
                              kind = "empty", text = s))
  best_pos <- Inf; best_alias <- NA_character_
  weak_pos <- Inf; weak_alias <- NA_character_
  for (a in ALIASES) {
    m <- regexpr(paste0("\\b", gsub("([.|()\\^{}+$*?])", "\\\\\\1", a)), s, perl = TRUE)
    if (m <= 0) next
    if (a %in% WEAK_ALIASES) {
      if (m < weak_pos) { weak_pos <- m; weak_alias <- a }
    } else if (m < best_pos || (m == best_pos && nchar(a) > nchar(best_alias))) {
      best_pos <- m; best_alias <- a
    }
  }
  if (is.na(best_alias) && !is.na(weak_alias)) best_alias <- weak_alias
  if (!is.na(best_alias))
    return(list(category = unname(TISSUE[[best_alias]]), alias = best_alias,
                kind = "matched", text = s))
  stripped <- gsub(paste0("\\b(?:", POSITIONAL, ")\\b"), " ", s, perl = TRUE)
  stripped <- trimws(gsub("\\s+", " ", stripped))
  if (!nzchar(stripped) || is_ignorable(stripped))
    return(list(category = NA_character_, alias = NA_character_,
                kind = "non-specific", text = s))
  list(category = NA_character_, alias = NA_character_, kind = "unresolved",
       text = s)
}

# Split a cell into terms. ";" always separates; "," only when the pieces each
# look like a tissue, so "Muscles in abdomen, thorax, limbs" stays one term.
split_terms <- function(cell) {
  parts <- trimws(strsplit(as.character(cell), ";")[[1]])
  parts <- parts[nzchar(parts)]
  out <- character(0)
  for (p in parts) {
    if (!grepl(",", p)) { out <- c(out, p); next }
    frag <- trimws(strsplit(p, ",")[[1]])
    frag <- frag[nzchar(frag)]
    hits <- vapply(frag, function(f) first_tissue(f)$kind == "matched",
                   logical(1), USE.NAMES = FALSE)
    out <- c(out, if (sum(hits) >= 2) frag else p)
  }
  out
}

parse_site <- function(cell) {
  if (is.na(cell) || !nzchar(trimws(cell)))
    return(list(categories = character(0), unresolved = character(0),
                kinds = character(0)))
  terms <- split_terms(cell)
  res <- lapply(terms, first_tissue)
  list(categories = unique(na.omit(vapply(res, function(r) r$category, character(1)))),
       unresolved = terms[vapply(res, function(r) r$kind, character(1)) == "unresolved"],
       kinds = vapply(res, function(r) r$kind, character(1)))
}

## ------------------------------------------------------------- self-test ----
# Every phrase the old script needed a hardcoded override for.
selftest <- function() {
  cases <- list(
    c("cells of hindgut near Malpighian tubules", "Gut / digestive tract"),
    c("longitudinal muscles of body wall", "Muscle"),
    c("connective tissue of the coelomic cavity", "Connective tissue"),
    c("under peritoneum of liver", "Liver / hepatopancreas"),
    c("skeletal muscle of the abdominal cavity", "Muscle"),
    c("subcutaneous connective tissue", "Connective tissue"),
    c("lymphocytes lying in body cavity", "Immune cells"),
    c("Skeletal musculature of abdomen and thorax", "Muscle"),
    c("subcutaneous connective tissue at base of fins", "Connective tissue"),
    c("muscles of the lateral wall of the abdomen near the anus", "Muscle"),
    c("subcutaneous layer of muscles", "Muscle"),
    c("connective tissue of muscle of dorsal fin", "Connective tissue"),
    c("Muscle surrounding ovaries", "Muscle"),
    c("oenocystes at surface of fat body", "Fat body / adipose"),
    c("connective tissue around oesophagus", "Connective tissue"),
    c("Gut epithelium of colon", "Gut / digestive tract"),
    c("mucosa of small intestine", "Gut / digestive tract"),
    c("Epithelium of the hepatopancreas", "Liver / hepatopancreas"),
    c("liver; kidney (tubule epithelium)", "Liver / hepatopancreas")
  )
  fails <- 0
  for (cs in cases) {
    got <- parse_site(cs[1])$categories
    ok <- length(got) && got[1] == cs[2]
    if (!ok) {
      fails <- fails + 1
      qc("  FAIL  %-58s expected %-24s got %s\n", cs[1], cs[2],
         if (length(got)) paste(got, collapse = "+") else "<none>")
    }
  }
  qc("head-noun self-test: %d/%d phrasings resolve to the intended tissue%s\n",
     length(cases) - fails, length(cases),
     if (fails == 0) " (all 12 old manual overrides now unnecessary)" else "")
}

## ------------------------------------------------------------------ load ----
if (!file.exists(CFG$db_path)) stop("database not found: ", CFG$db_path, call. = FALSE)
raw <- read_excel(CFG$db_path, sheet = CFG$sheet, col_types = "text",
                  .name_repair = "minimal")
message(sprintf("read %d rows x %d columns from sheet '%s'",
                nrow(raw), ncol(raw), CFG$sheet))

dat <- data.frame(
  row_xlsx = seq_len(nrow(raw)) + 1L,
  species  = trimws(get_col(raw, CFG$name_col)),
  site_raw = get_col(raw, CFG$site_col),
  stringsAsFactors = FALSE
)
dat$site_raw[is.na(dat$site_raw)] <- ""

cat("\n================ SITE OF INFECTION PARSING REPORT ================\n")
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
if (CFG$run_selftest) selftest()

parsed <- lapply(dat$site_raw, parse_site)
dat$n_cat     <- vapply(parsed, function(p) length(p$categories), integer(1))
dat$cat_list  <- vapply(parsed, function(p) paste(p$categories, collapse = "; "),
                        character(1))

has_site <- nzchar(trimws(dat$site_raw))
qc("\nentries                     : %d\n", nrow(dat))
qc("  no Site of Infection      : %d (%.1f%%)\n",
   sum(!has_site), 100 * sum(!has_site) / nrow(dat))
qc("  site recorded             : %d\n", sum(has_site))
qc("    resolved to >=1 tissue  : %d\n", sum(dat$n_cat > 0))
qc("    tissues per entry       : median %g, max %d\n",
   median(dat$n_cat[dat$n_cat > 0]), max(dat$n_cat))

## ---------------------------------------------------- problematic entries ---
unres <- unlist(lapply(parsed, function(p) p$unresolved))
if (length(unres)) {
  u <- as.data.frame(table(term = tolower(trimws(unres))), stringsAsFactors = FALSE)
  u <- u[order(-u$Freq), ]; names(u)[2] <- "n"
  qc("\nUNRESOLVED terms: %d occurrences, %d distinct\n", sum(u$n), nrow(u))
  qc("  (no tissue keyword, and not a body region or generic word)\n")
  show_tbl(u, n = 30, file = "QC_unresolved_terms.tsv")
  write_qc(u, "QC_unresolved_terms.tsv")
}

nocat <- dat[has_site & dat$n_cat == 0, ]
if (nrow(nocat)) {
  qc("\nsite text present but NO tissue assigned: %d entries\n", nrow(nocat))
  show_tbl(data.frame(species = nocat$species, site = nocat$site_raw),
           n = 20, file = "QC_no_tissue_assigned.tsv")
  write_qc(nocat[, c("row_xlsx", "species", "site_raw")],
           "QC_no_tissue_assigned.tsv")
}

## -------------------------------------------------------------- counting ----
long <- do.call(rbind, lapply(which(dat$n_cat > 0), function(i) {
  data.frame(species = dat$species[i], category = dat$category[i],
             tissue = strsplit(dat$cat_list[i], "; ")[[1]],
             stringsAsFactors = FALSE)
}))
long <- long[!duplicated(long[, c("species", "tissue")]), ]

# ---- systemic policy -------------------------------------------------------
# Reported here in full, because the count depends on it and the Methods must
# state which was used.
sys_sp <- unique(long$species[long$tissue == "Systemic / general"])
sys_only <- setdiff(sys_sp, unique(long$species[long$tissue != "Systemic / general"]))
qc("\nsystemic / general infections\n")
qc("  species with a systemic term            : %d\n", length(sys_sp))
qc("  ... of which systemic is the ONLY site  : %d\n", length(sys_only))
qc("  ... which also name a specific tissue   : %d  (the specific tissue is kept)\n",
   length(sys_sp) - length(sys_only))
if (identical(CFG$systemic, "exclude")) {
  n0 <- nrow(long)
  long <- long[long$tissue != "Systemic / general", ]
  qc("  CFG$systemic = \"exclude\": %d row(s) removed, %d species lost entirely\n",
     n0 - nrow(long), length(sys_only))
} else {
  qc("  CFG$systemic = \"category\": plotted as its own bar\n")
}

counts <- as.data.frame(table(tissue = long$tissue, category = long$category),
                        stringsAsFactors = FALSE)
tot <- as.data.frame(table(tissue = long$tissue), stringsAsFactors = FALSE)
names(tot)[2] <- "n"
tot <- tot[order(-tot$n), ]

pooled <- tot$tissue[tot$n < CFG$min_n]
if (length(pooled)) {
  qc("\npooled into 'Other' (fewer than %d species): %s\n",
     CFG$min_n, paste(pooled, collapse = ", "))
  long$tissue[long$tissue %in% pooled] <- "Other"
  counts <- as.data.frame(table(tissue = long$tissue, category = long$category),
                          stringsAsFactors = FALSE)
  tot <- as.data.frame(table(tissue = long$tissue), stringsAsFactors = FALSE)
  names(tot)[2] <- "n"; tot <- tot[order(-tot$n), ]
}

wide <- reshape(counts, idvar = "tissue", timevar = "category", direction = "wide")
names(wide) <- sub("^Freq\\.", "", names(wide))
wide$Total <- rowSums(wide[, -1, drop = FALSE])
wide <- wide[order(-wide$Total), ]
wide$`% of entries with a site` <-
  sprintf("%.1f", 100 * wide$Total / sum(dat$n_cat > 0))

cat("\n================ SPECIES PER TISSUE ================\n")
qc("%d tissue categories, %d species-tissue pairs, %d species placed\n\n",
   nrow(wide), nrow(long), length(unique(long$species)))
print(wide, row.names = FALSE, right = FALSE)
write_qc(wide, "counts_by_tissue.tsv")

multi <- dat[dat$n_cat >= 5, ]
if (nrow(multi)) {
  qc("\n%d entries name 5 or more tissues (broad or systemic infections):\n",
     nrow(multi))
  show_tbl(data.frame(species = multi$species, n = multi$n_cat,
                      tissues = multi$cat_list), n = 8,
           file = "QC_many_tissue_entries.tsv")
  write_qc(multi[, c("row_xlsx", "species", "site_raw", "n_cat", "cat_list")],
           "QC_many_tissue_entries.tsv")
}

## ------------------------------------------------------------------ plot ----
lvl <- rev(wide$tissue)
long$tissue <- factor(long$tissue, levels = lvl)
long$category <- factor(long$category, levels = c("named", "provisional"),
                        labels = c("Named species", "Provisional species"))
pal <- c("Named species" = "#0072B2", "Provisional species" = "#E69F00")

bar <- as.data.frame(table(tissue = long$tissue, category = long$category),
                     responseName = "n")
lab <- aggregate(n ~ tissue, bar, sum)
lab$txt <- if (CFG$show_percent)
  sprintf("%d (%.0f%%)", lab$n, 100 * lab$n / sum(dat$n_cat > 0)) else
  as.character(lab$n)

p <- ggplot(bar, aes(x = n, y = tissue, fill = category)) +
  geom_col(width = 0.74, colour = "grey20", linewidth = 0.2,
           position = position_stack(reverse = TRUE)) +
  geom_text(data = lab, aes(x = n, y = tissue, label = txt), inherit.aes = FALSE,
            hjust = -0.18, size = 2.9, colour = "grey20") +
  scale_fill_manual(values = pal, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Number of species reported", y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y = element_text(colour = "black"),
    axis.text.x = element_text(colour = "black"),
    axis.title.x = element_text(margin = margin(t = 8)),
    legend.position = c(0.98, 0.06),
    legend.justification = c(1, 0),
    legend.background = element_blank(),
    legend.key.size = unit(0.9, "lines"),
    plot.margin = margin(6, 14, 4, 6)
  )

ggsave(file.path(CFG$outdir, "fig_tissue_sites.pdf"), p,
       width = CFG$fig_width, height = CFG$fig_height,
       device = if (capabilities("cairo")) cairo_pdf else pdf)
ggsave(file.path(CFG$outdir, "fig_tissue_sites.png"), p,
       width = CFG$fig_width, height = CFG$fig_height, dpi = 600)

write_qc(dat[, c("row_xlsx", "species", "category", "site_raw", "n_cat",
                 "cat_list")], "per_entry_tissues.tsv")

qc("\nwritten to %s/\n  fig_tissue_sites.pdf / .png\n  counts_by_tissue.tsv\n  per_entry_tissues.tsv\n  QC_*.tsv\n",
   CFG$outdir)
