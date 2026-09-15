# Microsporidia database assembly flowchart (v3: exact step accounting)
#
# Run with the two spreadsheets beside this script:
#   Rscript microsporidia_database_flowchart_3.R
#
# Or supply paths explicitly:
#   Rscript microsporidia_database_flowchart_3.R database.xlsx \
#     "Actively Updated Masterlist" output/microsporidia_database_flowchart \
#     mbio.0149021st001.xlsx
#
# The script reads the current workbook AND the original published database
# (mBio 2021 Table S1, 1,440 rows), reconstructs what happened at each
# assembly step, and creates:
#   <stem>.png          high-resolution talk figure
#   <stem>.pdf          vector supplemental figure
#   <stem>_counts.csv   metric audit (calculated / override / used)
#   <stem>_steps.csv    step-by-step ledger with running totals
#   <stem>_rows.csv     per-row provenance assignment (for checking)
#
# How the counts are derived (all from "Important Remarks" + the original
# species list; nothing is hard-coded except the documented alias table):
#
#   cluster     "PROVISIONAL SPECIES: defined by 18S sequence cluster"
#   llm         any "[ADDED <date> ...]" marker (2026-08-06 audits and the
#               2026-08-09 accession/paper reconciliation)
#   original    row carries a species from the original database: exact name
#               match, a name listed in a "combined records:" merge note, or
#               a documented rename (original_name_aliases below)
#   publication new species from 2021-2026 submissions ("[TRANSFERRED/ADDED ...]")
#   llm (unmarked) any remaining row with no marker and no original match is an
#               addition made during the LLM-assisted audit without a marker
#               (currently one row: Unikaryon sp. 2); counted with llm_search
#
#   Merges are parsed from the two annotation styles:
#     "[MERGED 2026-08-06: original rows 258, 1483]"
#         rows <= original_row_limit are original entries, larger rows are
#         2021-2026 submissions
#     "[MERGED 2026-08-10 - ...; combined records: A | B | C; ...]"
#         n combined names; the number of "[TRANSFERRED/ADDED" markers on the
#         row tells how many of them were submissions, the rest are originals

library(grid)

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop(
    "Package 'readxl' is required. Install it with: ",
    "install.packages('readxl')"
  )
}

# ============================ INPUTS ================================
all_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", all_args, value = TRUE)
script_dir <- if (length(file_arg) == 1) {
  dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  getwd()
}

args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[1] else file.path(
  script_dir, "Microsporidia_Characteristics_Database_merge_pro4.xlsx"
)
sheet_name <- if (length(args) >= 2) args[2] else "Actively Updated Masterlist"
output_stem <- if (length(args) >= 3) args[3] else file.path(
  script_dir, "microsporidia_database_flowchart"
)
# Original published database (mBio 2021, Table S1). If the file is missing
# the script falls back to published_original_count and marker-only logic,
# with a warning (see below).
original_file <- if (length(args) >= 4) args[4] else file.path(
  script_dir, "mbio.0149021st001.xlsx"
)
original_sheet <- "Masterlist"
published_original_count <- 1440L

# Set any value below to an integer to replace the spreadsheet calculation.
# Leave it as NA_integer_ to use the calculated value.
overrides <- list(
  original_entries = NA_integer_,
  merged_redundant_records = NA_integer_,
  publication_additions = NA_integer_,
  llm_search_additions = NA_integer_,
  genbank_cluster_additions = NA_integer_,
  named_species = NA_integer_,
  provisional_species = NA_integer_,
  species_with_18s = NA_integer_,
  total_species = NA_integer_
)

# Provenance markers used in the current database. Edit these regular
# expressions if the wording in Important Remarks changes in a future version.
cluster_pattern <- "PROVISIONAL SPECIES:\\s*defined by 18S sequence cluster"
llm_pattern <- "\\[ADDED \\d{4}-\\d{2}-\\d{2}"
transfer_pattern <- "\\[TRANSFERRED/ADDED \\d{4}-\\d{2}-\\d{2}"
merge_rows_pattern <- "\\[MERGED \\d{4}-\\d{2}-\\d{2}:\\s*original rows ([0-9, ]+)\\]"
merge_names_pattern <- "combined records:\\s*([^;\\]]+)"
# In the "original rows" merge annotations, row numbers up to this value are
# entries of the original database; larger numbers are 2021-2026 submissions.
original_row_limit <- 1440L

# Original species whose current row cannot be found by name or merge note.
# original name = current "Species Name". Two originals pointing at the same
# current row means a duplicate was consolidated before the 2026 curation.
original_name_aliases <- c(
  "Microspordium sp. C. aurata" = "Microsporidium sp. isolate TR",
  "Microsporidium sp. 24" = "Microsporidium prosopium",
  "Pleistophora manningeri" = "Pleistophora manningeri (Plistophora manningeri)",
  "Plistophora manningeri" = "Pleistophora manningeri (Plistophora manningeri)",
  "Episeptuni inversum" = "Episeptum inversum"
)

# Optional wording changes
figure_title <- "Building a comprehensive microsporidia species database"
llm_note <- "Broader LLM-assisted publication audit"
genbank_note <- "Sequence clusters treated as provisional species"
# publication_note is built from the counts below; set a string to replace it.
publication_note <- NULL
# ====================================================================

if (!file.exists(input_file)) {
  stop("Spreadsheet not found: ", input_file)
}

database <- readxl::read_excel(
  input_file,
  sheet = sheet_name,
  col_types = "text",
  .name_repair = "minimal"
)
database <- as.data.frame(database, check.names = FALSE)

normalize_header <- function(x) {
  tolower(trimws(gsub("\\s+", " ", as.character(x), perl = TRUE)))
}

find_column <- function(data, candidates) {
  headers <- normalize_header(names(data))
  matches <- match(normalize_header(candidates), headers, nomatch = 0L)
  matches <- matches[matches > 0L]
  if (length(matches) == 0L) {
    stop(
      "Could not find required column. Expected one of: ",
      paste(candidates, collapse = ", ")
    )
  }
  matches[1]
}

species_col <- find_column(database, c("Species Name"))
remarks_col <- find_column(database, c("Important Remarks", "Remarks"))
accession_col <- find_column(
  database,
  c("18S Accession #", "18S Accession", "18S accession number")
)

as_clean_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("\u00a0", " ", x, fixed = TRUE)
  trimws(gsub("\\s+", " ", x, perl = TRUE))
}
name_key <- function(x) tolower(as_clean_text(x))

# Remove empty spreadsheet rows before calculating anything.
species_names <- as_clean_text(database[[species_col]])
database <- database[nzchar(species_names), , drop = FALSE]
species_names <- as_clean_text(database[[species_col]])
remarks <- as_clean_text(database[[remarks_col]])
accessions <- as_clean_text(database[[accession_col]])
n_rows <- nrow(database)

# ---------------------- Marker flags per row ------------------------
cluster_flag <- grepl(cluster_pattern, remarks, ignore.case = TRUE, perl = TRUE)
llm_flag <- grepl(llm_pattern, remarks, perl = TRUE)
transfer_count <- vapply(
  regmatches(remarks, gregexpr(transfer_pattern, remarks, perl = TRUE)),
  length, integer(1)
)

# "[MERGED <date>: original rows a, b, c]" -> list of integer vectors
merge_rows <- lapply(remarks, function(x) {
  m <- regmatches(x, regexec(merge_rows_pattern, x, perl = TRUE))[[1]]
  if (length(m) < 2L) return(integer(0))
  as.integer(trimws(strsplit(m[2], ",", fixed = TRUE)[[1]]))
})
# "combined records: A | B | C" -> list of character vectors
merge_names <- lapply(remarks, function(x) {
  m <- regmatches(x, regexec(merge_names_pattern, x, perl = TRUE))[[1]]
  if (length(m) < 2L) return(character(0))
  as_clean_text(strsplit(m[2], "|", fixed = TRUE)[[1]])
})

# ---------------------- Original-database matching ------------------
have_original <- file.exists(original_file)
if (have_original) {
  original <- readxl::read_excel(
    original_file, sheet = original_sheet, col_types = "text",
    .name_repair = "minimal"
  )
  original <- as.data.frame(original, check.names = FALSE)
  orig_species_col <- find_column(original, c("Species", "Species Name"))
  original_names <- as_clean_text(original[[orig_species_col]])
  original_names <- original_names[nzchar(original_names)]
  original_entries_calculated <- length(original_names)

  # For every original species, find the current row that carries it.
  current_keys <- name_key(species_names)
  alias_keys <- name_key(original_name_aliases)
  names(alias_keys) <- name_key(names(original_name_aliases))
  carrier_row <- vapply(name_key(original_names), function(key) {
    hit <- which(current_keys == key)
    if (length(hit) == 0L) {
      hit <- which(vapply(merge_names, function(nm) key %in% name_key(nm),
                          logical(1)))
    }
    if (length(hit) == 0L && key %in% names(alias_keys)) {
      hit <- which(current_keys == alias_keys[[key]])
    }
    if (length(hit) == 0L) NA_integer_ else hit[1]
  }, integer(1))

  unaccounted <- original_names[is.na(carrier_row)]
  if (length(unaccounted) > 0L) {
    warning(
      length(unaccounted), " original species could not be located in the ",
      "current database (add them to original_name_aliases): ",
      paste(head(unaccounted, 10), collapse = "; ")
    )
  }
  original_flag <- seq_len(n_rows) %in% carrier_row
} else {
  warning(
    "Original database not found (", original_file, "). Falling back to ",
    "marker-only logic: rows without addition/cluster markers are treated ",
    "as original entries, and original_entries = ", published_original_count,
    ". Direct curated additions without a marker cannot be separated from ",
    "original entries in this mode."
  )
  original_entries_calculated <- published_original_count
  # Transferred rows whose 2026-08-10 merge note combined them with at least
  # one original record are original entries as well.
  originals_in_group <- pmax(lengths(merge_names) - transfer_count, 0L)
  original_flag <- !cluster_flag & !llm_flag &
    (transfer_count == 0L | originals_in_group > 0L)
  unaccounted <- character(0)
}

# Rows that are clusters or LLM additions can never be originals.
conflict <- original_flag & (cluster_flag | llm_flag)
if (any(conflict)) {
  warning(sum(conflict), " rows match an original species but carry an ",
          "LLM/cluster marker; they are counted as original entries.")
  cluster_flag[conflict] <- FALSE
  llm_flag[conflict] <- FALSE
}

# Mutually exclusive provenance per row
provenance <- ifelse(
  cluster_flag, "genbank_cluster",
  ifelse(llm_flag, "llm_search",
    ifelse(original_flag, "original",
      ifelse(transfer_count > 0L, "publication_submission",
             "llm_search_unmarked")))
)

# ---------------------- Merge accounting ----------------------------
# 2026-08-06 merges ("original rows a, b"): a <= limit is an original entry.
m06_orig <- vapply(merge_rows, function(r) sum(r <= original_row_limit),
                   integer(1))
m06_subm <- vapply(merge_rows, function(r) sum(r > original_row_limit),
                   integer(1))
# The retained row keeps its own TRANSFERRED marker when it is a submission.
m06_subm_absorbed <- pmax(m06_subm - as.integer(m06_orig == 0L), 0L)
m06_orig_redundant <- pmax(m06_orig - 1L, 0L)
m06_subm_into_original <- ifelse(m06_orig > 0L, m06_subm_absorbed, 0L)
m06_subm_duplicates <- ifelse(m06_orig == 0L, m06_subm_absorbed, 0L)

# 2026-08-10 identity merges ("combined records: A | B | C")
m10_n <- lengths(merge_names)
m10_subm <- ifelse(m10_n > 0L, pmin(transfer_count, m10_n), 0L)
m10_orig <- pmax(m10_n - m10_subm, 0L)
m10_orig_redundant <- pmax(m10_orig - 1L, 0L)
m10_subm_into_original <- ifelse(m10_orig > 0L, m10_subm, 0L)
m10_subm_duplicates <- ifelse(m10_orig == 0L, pmax(m10_subm - 1L, 0L), 0L)

identity_merge_groups <- sum(m10_orig >= 2L)
identity_merged_originals <- sum(m10_orig_redundant)
original_rows_retained <- sum(original_flag)
# Originals that disappeared without an identity-merge note were consolidated
# as duplicates before the 2026 curation (e.g. spelling-variant rows).
pre2026_consolidated <- original_entries_calculated -
  original_rows_retained - identity_merged_originals
merged_redundant_records_calculated <- identity_merged_originals +
  pre2026_consolidated

submissions_new <- sum(provenance == "publication_submission")
submissions_into_original <- sum(m06_subm_into_original) +
  sum(m10_subm_into_original)
submissions_duplicates <- sum(m06_subm_duplicates) + sum(m10_subm_duplicates)
submissions_total <- submissions_new + submissions_into_original +
  submissions_duplicates
llm_unmarked_additions <- sum(provenance == "llm_search_unmarked")

publication_additions_calculated <- submissions_new
llm_search_additions_calculated <- sum(provenance == "llm_search") +
  llm_unmarked_additions
llm_dates <- table(regmatches(
  remarks[provenance == "llm_search"],
  regexpr("(?<=\\[ADDED )\\d{4}-\\d{2}-\\d{2}",
          remarks[provenance == "llm_search"], perl = TRUE)
))
genbank_cluster_additions_calculated <- sum(provenance == "genbank_cluster")
total_species_calculated <- n_rows

# --------------------- Named/provisional species --------------------
# Named vs provisional comes from classify_species.R -- the SAME canonical
# rule used by every figure panel, so this flowchart's totals and Fig 1
# cannot disagree. Sourced from the script's own directory because this
# script supports being run from elsewhere.
source(file.path(script_dir, "classify_species.R"))

provisional_flag <- vapply(
  species_names, is_provisional_name, logical(1)
)
provisional_species_calculated <- sum(provisional_flag)
named_species_calculated <- sum(!provisional_flag)
species_with_18s_calculated <- sum(nzchar(accessions))

# ----------------------------- Overrides ----------------------------
use_override <- function(name, calculated) {
  value <- overrides[[name]]
  if (is.null(value) || length(value) == 0L || is.na(value)) {
    return(as.integer(calculated))
  }
  as.integer(value)
}

original_entries <- use_override(
  "original_entries", original_entries_calculated
)
merged_redundant_records <- use_override(
  "merged_redundant_records", merged_redundant_records_calculated
)
publication_additions <- use_override(
  "publication_additions", publication_additions_calculated
)
llm_search_additions <- use_override(
  "llm_search_additions", llm_search_additions_calculated
)
genbank_cluster_additions <- use_override(
  "genbank_cluster_additions", genbank_cluster_additions_calculated
)
named_species <- use_override("named_species", named_species_calculated)
provisional_species <- use_override(
  "provisional_species", provisional_species_calculated
)
species_with_18s <- use_override(
  "species_with_18s", species_with_18s_calculated
)
total_species <- use_override("total_species", total_species_calculated)

post_merge_entries <- original_entries - merged_redundant_records
workflow_total <- post_merge_entries + publication_additions +
  llm_search_additions + genbank_cluster_additions

if (workflow_total != total_species) {
  warning(
    "Workflow components sum to ", workflow_total,
    ", but total_species is ", total_species,
    ". Check provenance patterns or overrides."
  )
}

if (named_species + provisional_species != total_species) {
  warning("Named + provisional species does not equal total species.")
}

fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

if (is.null(publication_note)) {
  publication_note <- paste0(
    fmt(submissions_into_original + submissions_duplicates),
    " merged into existing entries"
  )
}

# ----------------------------- Audit files --------------------------
metric_names <- c(
  "original_entries", "merged_redundant_records", "publication_additions",
  "llm_search_additions", "genbank_cluster_additions", "named_species",
  "provisional_species", "species_with_18s", "total_species"
)
calculated_values <- c(
  original_entries_calculated, merged_redundant_records_calculated,
  publication_additions_calculated, llm_search_additions_calculated,
  genbank_cluster_additions_calculated, named_species_calculated,
  provisional_species_calculated, species_with_18s_calculated,
  total_species_calculated
)
used_values <- c(
  original_entries, merged_redundant_records, publication_additions,
  llm_search_additions, genbank_cluster_additions, named_species,
  provisional_species, species_with_18s, total_species
)
override_values <- vapply(
  metric_names,
  function(name) {
    value <- overrides[[name]]
    if (is.null(value) || length(value) == 0L || is.na(value)) NA_integer_
    else as.integer(value)
  },
  integer(1)
)
audit_table <- data.frame(
  metric = metric_names,
  calculated_from_spreadsheet = calculated_values,
  override = override_values,
  value_used_in_figure = used_values,
  stringsAsFactors = FALSE
)
write.csv(
  audit_table,
  paste0(output_stem, "_counts.csv"),
  row.names = FALSE,
  na = ""
)

# Step ledger: what happened at each stage, with running totals
step_change <- c(
  original_entries_calculated,
  -pre2026_consolidated,
  -identity_merged_originals,
  0L,
  submissions_new,
  llm_search_additions_calculated,
  genbank_cluster_additions_calculated
)
step_table <- data.frame(
  step = seq_along(step_change),
  description = c(
    "Original published database (mBio 2021 Table S1)",
    "Duplicate rows consolidated before the 2026 curation",
    paste0("Identity-based merges among original entries (",
           identity_merge_groups, " groups, 2026-08-10)"),
    paste0("2021-2026 submissions merged into existing entries (",
           submissions_into_original, " into original entries; ",
           submissions_duplicates, " duplicate submissions merged; ",
           submissions_total, " submissions in total)"),
    "New species from 2021-2026 submissions",
    paste0("LLM-assisted literature audit additions (",
           paste(names(llm_dates), llm_dates, sep = ": ", collapse = "; "),
           "; unmarked: ", llm_unmarked_additions, ")"),
    "GenBank 18S sequence clusters added as provisional species"
  ),
  change = step_change,
  running_total = cumsum(step_change),
  stringsAsFactors = FALSE
)
write.csv(step_table, paste0(output_stem, "_steps.csv"), row.names = FALSE)

row_table <- data.frame(
  sheet_row = which(nzchar(as_clean_text(
    readxl::read_excel(input_file, sheet = sheet_name, col_types = "text",
                       .name_repair = "minimal")[[species_col]]
  ))) + 1L,
  species_name = species_names,
  provenance = provenance,
  provisional_name = provisional_flag,
  has_18s = nzchar(accessions),
  transferred_markers = transfer_count,
  merge_0806_original_rows = vapply(merge_rows, paste, character(1),
                                    collapse = " "),
  merge_0810_combined_records = vapply(merge_names, paste, character(1),
                                       collapse = " | "),
  stringsAsFactors = FALSE
)
write.csv(row_table, paste0(output_stem, "_rows.csv"), row.names = FALSE)

# ----------------------------- Figure -------------------------------
colors <- list(
  background = "#F7F9FC",
  ink = "#183042",
  muted = "#5C6B76",
  arrow = "#8293A1",
  source = "#1F4E79",
  source_fill = "#EAF2F8",
  publication = "#2A7F9E",
  publication_fill = "#E7F5F8",
  llm = "#6F58A6",
  llm_fill = "#F1EDF8",
  genbank = "#C77819",
  genbank_fill = "#FFF3E3",
  final = "#18766E",
  final_fill = "#E5F5F2",
  white = "#FFFFFF"
)

draw_round_box <- function(x, y, w, h, border, fill, title, number,
                           line1 = NULL, line2 = NULL) {
  grid.roundrect(
    x = x, y = y, width = w, height = h,
    r = unit(0.018, "npc"),
    gp = gpar(fill = fill, col = border, lwd = 2.2)
  )
  grid.text(
    title, x = x, y = y + h * 0.29,
    gp = gpar(col = border, fontsize = 12.5, fontface = "bold")
  )
  grid.text(
    number, x = x, y = y + h * 0.01,
    gp = gpar(col = colors$ink, fontsize = 22, fontface = "bold")
  )
  if (!is.null(line1)) {
    grid.text(
      line1, x = x, y = y - h * 0.22,
      gp = gpar(col = colors$ink, fontsize = 9.5, fontface = "bold")
    )
  }
  if (!is.null(line2)) {
    grid.text(
      line2, x = x, y = y - h * 0.36,
      gp = gpar(col = colors$muted, fontsize = 8.7)
    )
  }
}

draw_arrow <- function(x0, y0, x1, y1) {
  grid.lines(
    x = unit(c(x0, x1), "npc"), y = unit(c(y0, y1), "npc"),
    arrow = arrow(type = "closed", length = unit(0.13, "inches")),
    gp = gpar(col = colors$arrow, lwd = 2.2, lineend = "round")
  )
}

draw_stat_chip <- function(x, y, w, h, value, label) {
  grid.roundrect(
    x = x, y = y, width = w, height = h,
    r = unit(0.012, "npc"),
    gp = gpar(fill = colors$white, col = NA)
  )
  grid.text(
    fmt(value), x = x, y = y + h * 0.12,
    gp = gpar(col = colors$final, fontsize = 19, fontface = "bold")
  )
  grid.text(
    label, x = x, y = y - h * 0.22,
    gp = gpar(col = colors$ink, fontsize = 10.2, fontface = "bold")
  )
}

draw_flowchart <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = colors$background, col = NA))

  grid.text(
    figure_title, x = 0.5, y = 0.92,
    gp = gpar(col = colors$ink, fontsize = 24, fontface = "bold")
  )
  grid.text(
    "Database assembly and taxon provenance",
    x = 0.5, y = 0.865,
    gp = gpar(col = colors$muted, fontsize = 12)
  )

  box_y <- 0.635
  box_w <- 0.205
  box_h <- 0.29
  xs <- c(0.13, 0.375, 0.625, 0.87)

  draw_round_box(
    xs[1], box_y, box_w, box_h,
    colors$source, colors$source_fill,
    "1  Existing\ndatabase", paste0(fmt(original_entries), " entries"),
    paste0(fmt(merged_redundant_records), " redundant records merged"),
    paste0(fmt(post_merge_entries), " retained")
  )
  draw_round_box(
    xs[2], box_y, box_w, box_h,
    colors$publication, colors$publication_fill,
    "2  Publications\n2021–2026", paste0("+", fmt(publication_additions)),
    "species added", publication_note
  )
  draw_round_box(
    xs[3], box_y, box_w, box_h,
    colors$llm, colors$llm_fill,
    "3  LLM-assisted\nsearches", paste0("+", fmt(llm_search_additions)),
    "species added", llm_note
  )
  draw_round_box(
    xs[4], box_y, box_w, box_h,
    colors$genbank, colors$genbank_fill,
    "4  GenBank 18S\nclustering", paste0("+", fmt(genbank_cluster_additions)),
    "provisional species", genbank_note
  )

  for (i in 1:3) {
    draw_arrow(xs[i] + box_w / 2 + 0.006, box_y,
               xs[i + 1] - box_w / 2 - 0.006, box_y)
  }
  draw_arrow(xs[4], box_y - box_h / 2 - 0.015, 0.77, 0.395)

  final_x <- 0.5
  final_y <- 0.255
  final_w <- 0.74
  final_h <- 0.255
  grid.roundrect(
    x = final_x, y = final_y, width = final_w, height = final_h,
    r = unit(0.02, "npc"),
    gp = gpar(fill = colors$final_fill, col = colors$final, lwd = 2.6)
  )
  grid.text(
    paste0(fmt(total_species), " total species"),
    x = final_x, y = final_y + final_h * 0.31,
    gp = gpar(col = colors$final, fontsize = 26, fontface = "bold")
  )

  chip_y <- final_y - final_h * 0.12
  chip_w <- 0.19
  chip_h <- 0.105
  draw_stat_chip(0.27, chip_y, chip_w, chip_h, named_species, "Named species")
  draw_stat_chip(0.50, chip_y, chip_w, chip_h, provisional_species, "Provisional species")
  draw_stat_chip(0.73, chip_y, chip_w, chip_h, species_with_18s, "With 18S accessions")

  grid.text(
    "Named = species-level name; provisional = unnamed or “sp.” designation",
    x = 0.5, y = 0.055,
    gp = gpar(col = colors$muted, fontsize = 9.5)
  )
}

# High-resolution PNG for talks
png(
  paste0(output_stem, ".png"),
  width = 3200, height = 1800, res = 240, bg = colors$background
)
draw_flowchart()
invisible(dev.off())

# Vector PDF for publication/supplemental use
pdf(
  paste0(output_stem, ".pdf"),
  width = 13.333, height = 7.5, useDingbats = FALSE, bg = colors$background
)
draw_flowchart()
invisible(dev.off())

# ----------------------------- Console summary ----------------------
message("Created ", output_stem, ".png and ", output_stem, ".pdf")
message("Spreadsheet: ", normalizePath(input_file))
if (have_original) message("Original database: ", normalizePath(original_file))
message("Audit files: ", output_stem, "_counts.csv, _steps.csv, _rows.csv")
message("")
message("Step ledger:")
for (i in seq_len(nrow(step_table))) {
  message(sprintf("  %d. %-90s %+6d  -> %s", step_table$step[i],
                  step_table$description[i], step_table$change[i],
                  fmt(step_table$running_total[i])))
}
message("")
message("Existing database: ", fmt(original_entries), " entries; ",
        fmt(merged_redundant_records), " redundant records merged (",
        identity_merged_originals, " identity merges + ",
        pre2026_consolidated, " earlier duplicates); ",
        fmt(post_merge_entries), " retained")
message("Publications 2021-2026: +", fmt(publication_additions), " (",
        submissions_new, " new from ", submissions_total, " submissions; ",
        submissions_into_original, " merged into existing entries, ",
        submissions_duplicates, " duplicate submissions)")
message("LLM-assisted searches: +", fmt(llm_search_additions), " (",
        llm_unmarked_additions, " without an [ADDED] marker)")
message("GenBank 18S clusters: +", fmt(genbank_cluster_additions))
message("Workflow total: ", fmt(workflow_total))
message("Database rows: ", fmt(total_species))
