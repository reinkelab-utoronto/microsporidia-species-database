#!/usr/bin/env Rscript
###############################################################################
# Environmental range of each microsporidia species
#
# A species infects one or more hosts; each host lives in one or more
# environments. This asks: how many species are confined to one environment,
# and how many span several? The answer is drawn as an AREA-PROPORTIONAL Euler
# diagram (a Venn where the regions are sized by their counts) with the counts
# printed inside.
#
# TWO WAYS TO DEFINE "THE ENVIRONMENTS OF A SPECIES", AND THEY GIVE VERY
# DIFFERENT FIGURES. Choose deliberately:
#
#   env_field = "environment_all"       (DEFAULT -- matches the databases)
#     every habitat each host occupies. This is what GBIF and WoRMS themselves
#     return: WoRMS gives four independent flags (isMarine, isBrackish,
#     isFreshwater, isTerrestrial) and a euryhaline fish or an insect with an
#     aquatic larva has several set at once. A Venn built from database
#     annotations is therefore an "all habitats" Venn by construction, and this
#     setting reproduces it. Roughly a third of species come out
#     multi-environment.
#
#   env_field = "environment_primary"
#     one habitat per host: where the host mainly lives. A species then spans
#     two environments only if its HOSTS differ -- one host in freshwater,
#     another terrestrial. About 99% of species come out single-environment.
#     This is a stricter and arguably more biological claim ("does the parasite
#     bridge habitats?"), but it is NOT comparable to a database-annotated Venn
#     and it will look nothing like the published panel.
#
#   Note the two answer different questions and neither is wrong. What is wrong
#   is mixing them, which is why the merge below forces both sources onto the
#   same footing (see PREDICTION_FIELD).
#
# The diagram is drawn with eulerr, which fits circles or ellipses whose areas
# approximate the counts. FOUR-SET EULER DIAGRAMS CANNOT ALWAYS BE DRAWN
# EXACTLY, so the fit error is computed and printed; if it is large the exact
# UpSet-style panel written alongside is the honest figure.
#
# Usage:  Rscript fig_environment_euler.R
###############################################################################

CFG <- list(
  # tables produced earlier. Either or both; hosts present in both are resolved
  # by source_rule.
  api_file  = "host_taxonomy_output/01_species_host_long.tsv",
  pred_file = "predicted_environment/01_species_host_predicted.tsv",
  # Which sources to use. EVERY source listed here MUST exist -- a missing file
  # is a hard error, never a silent fall-back to whatever else is on disk.
  # To use only one deliberately, remove the other from this vector.
  sources = c("api", "prediction"),
  # drop the "(... hyperhost)" rows of a hyperparasite chain, as panel C does
  drop_hyperhosts = TRUE,

  # WHICH HOST RECORDS TO COUNT -- matches fig_host_phylum.R's host_role_filter.
  # The lookup tables carry both Natural Host(s) and Experimental Host(s) rows
  # (host_taxonomy_environment.py harvests both), and an experimental infection
  # shows what a parasite CAN infect in the lab, not where it lives -- counting
  # it adds the lab host's habitat to the species and makes this panel disagree
  # with panel C and the phylogeny, which both exclude experimental hosts.
  #   "natural"      Natural Host(s) only (default; the host_column contains
  #                  "natural")
  #   "experimental" Experimental Host(s) only
  #   "all"          both, i.e. the previous behaviour
  # Rows with no host_column (e.g. a hand-made input file) carry no role and are
  # always kept, so this never silently drops user-supplied hosts.
  host_role_filter = "natural",

  # "api_first"        use the API habitat where it exists, prediction to fill
  # "prediction_first" the reverse
  # "union"            take both (widest, most generous)
  source_rule = "api_first",

  env_field = "environment_all",       # or "environment_primary"; see above
  # environments to draw. Dropping "brackish" (small, and almost always
  # alongside marine) gives a 3-set diagram, which can be drawn exactly and
  # scaled by every package. Species are then classified on the rest.
  # drawing order matters for a four-ellipse Venn: the two outermost ellipses
  # meet in only one small lens, so put the least-overlapping sets there.
  sets = c("terrestrial", "brackish", "freshwater", "marine"),
  subset    = "all",                   # "all" | "named" | "provisional"

  # "auto" = scaled Euler for 2-3 sets, classic four-ellipse Venn for 4.
  # "classic" forces the fixed Venn shape, "scaled" forces area-proportional.
  shape_rule = "auto",
  # "host"    one row per host species, parasites ignored
  # "species" one row per microsporidian, environments unioned over its hosts
  # "both"    draw both (default)
  unit = "both",
  shape      = "ellipse",              # eulerr shape when scaling (2-3 sets)
  outdir     = "environment_euler",
  fig_width  = 6.6,
  fig_height = 5.4
)

for (p in c("ggplot2")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages(library(ggplot2))
has_eulerr <- requireNamespace("eulerr", quietly = TRUE)
has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
dir.create(CFG$outdir, showWarnings = FALSE, recursive = TRUE)
options(width = 200)
qc <- function(...) cat(sprintf(...), sep = "")
cap1 <- function(x) sub("^(\\w)", "\\U\\1", as.character(x), perl = TRUE)

## ---------------------------------------------------------------- palette ---
# Okabe-Ito, chosen so the four fills stay distinct under deuteranopia and
# protanopia AND differ in lightness, so the overlaps remain readable in
# greyscale. Blue-family for the wet habitats, orange for dry, which also makes
# the diagram legible at a glance.
PAL <- c(
  freshwater  = "#56B4E9",   # sky blue
  brackish    = "#009E73",   # bluish green
  marine      = "#0072B2",   # blue
  terrestrial = "#7F4F24"    # orange
)
ENVS <- names(PAL)

## --------------------------------------------------------------- helpers ----
read_long <- function(path, env_col_candidates, label = path) {
  if (!file.exists(path))
    stop("input file not found: ", path, call. = FALSE)
  # quote = "" and a row-count check, as in the other readers: a stray double
  # quote in a host cell otherwise makes read.delim swallow following lines,
  # silently merging rows and dropping most of the table without any error
  n_lines <- length(readLines(path, warn = FALSE)) - 1L
  d <- read.delim(path, stringsAsFactors = FALSE, quote = "", na.strings = "",
                  comment.char = "")
  if (nrow(d) != n_lines)
    stop(sprintf("%s did not parse cleanly: %d data lines, %d rows read",
                 path, n_lines, nrow(d)), call. = FALSE)
  hit <- env_col_candidates[env_col_candidates %in% names(d)][1]
  if (is.na(hit))
    stop("no environment column in ", path, "\n  looked for: ",
         paste(env_col_candidates, collapse = ", "), "\n  found: ",
         paste(names(d), collapse = ", "), call. = FALSE)
  sp <- c("microsporidia_species", "species")[c("microsporidia_species", "species")
                                              %in% names(d)][1]
  cat <- if ("microsporidia_category" %in% names(d)) d$microsporidia_category else ""
  out <- data.frame(species = d[[sp]], category = cat,
                    host = d$host_name, env = d[[hit]],
                    host_raw = if ("host_raw" %in% names(d)) d$host_raw else "",
                    host_column = if ("host_column" %in% names(d))
                      d$host_column else "",
                    cell = if ("cell" %in% names(d)) d$cell else "",
                    stringsAsFactors = FALSE)
  out$env[is.na(out$env)] <- ""
  out <- out[nzchar(out$species), ]

  # NATURAL vs EXPERIMENTAL, applied here so every count downstream -- the merge,
  # both units, the trace table -- reflects the same set of rows. Same rule as
  # fig_host_phylum.R: keep rows whose host_column names the wanted role, plus
  # any row that carries no role at all (so a hand-made input file is not lost).
  if (!identical(CFG$host_role_filter, "all")) {
    if (!"host_column" %in% names(d)) {
      qc("  %-10s no host_column; cannot separate natural from experimental, "
         , label); qc("all rows kept\n")
    } else {
      hc <- ifelse(is.na(out$host_column), "", tolower(out$host_column))
      is_nat <- grepl("natural", hc)
      is_exp <- grepl("experimental|lab", hc)
      want   <- identical(CFG$host_role_filter, "natural")
      keep <- (if (want) is_nat else is_exp) | (!is_nat & !is_exp)
      qc("  %-10s host_role_filter = %s: keeping %d of %d row(s)\n",
         label, CFG$host_role_filter, sum(keep), nrow(out))
      if (!any(keep))
        stop("host_role_filter removed every row from ", path,
             "; check the host_column values", call. = FALSE)
      out <- out[keep, , drop = FALSE]
    }
  }
  out
}

split_env <- function(x) strsplit(x, ";\\s*")

## ------------------------------------------------------------------ load ----
# --- required inputs, checked before anything else --------------------------
# Proceeding with one source when the other is absent produces a figure that
# looks fine and answers a different question. That is always an error here;
# opting out has to be explicit via CFG$sources.
want <- CFG$sources
if (!length(want) || !all(want %in% c("api", "prediction")))
  stop("CFG$sources must contain \"api\" and/or \"prediction\"", call. = FALSE)
need <- c(api = CFG$api_file, prediction = CFG$pred_file)[want]
absent <- need[!file.exists(need)]
if (length(absent))
  stop(paste0(
    "required input file(s) not found:\n",
    paste0("  ", names(absent), ": ", absent, collapse = "\n"),
    "\n\nProduce them with:\n",
    if ("api" %in% names(absent))
      "  python host_taxonomy_environment.py --db <database>.xlsx --email <you>\n" else "",
    if ("prediction" %in% names(absent))
      "  python predict_host_environment.py --db <database>.xlsx\n" else "",
    "\nOr, to deliberately use only the source you have, set\n",
    sprintf("  CFG$sources <- c(%s)\n",
            paste0("\"", setdiff(want, names(absent)), "\"", collapse = ", "))),
    call. = FALSE)
if (length(want) == 1)
  message("NOTE: using the ", want, " source only, because CFG$sources says so.")

cat("\n---- host-role filter ----\n")
api <- if ("api" %in% want)
  read_long(CFG$api_file, c("environment", "host_derived_environment"),
            "API") else NULL
# The API column is inherently multi-label (a union of WoRMS/GBIF flags), so
# pairing it with the prediction's single-label "primary" column would let
# API-covered hosts be multi-environment while prediction-only hosts could
# never be -- inflating the intersections exactly where the API has coverage.
# The prediction column is therefore chosen to match the API's semantics.
PREDICTION_FIELD <- if (identical(CFG$env_field, "environment_primary"))
  "environment_primary" else "environment_all"
prd <- if ("prediction" %in% want)
  read_long(CFG$pred_file, c(PREDICTION_FIELD, "environment_all",
                             "environment_primary"), "prediction") else NULL
if (identical(CFG$env_field, "environment_primary"))
  qc("")  # nothing; the warning below covers it
cat("\n================ ENVIRONMENTAL RANGE ================\n")
qc("INPUT api  : %s%s\n", CFG$api_file,
   if (file.exists(CFG$api_file))
     sprintf("  (modified %s)",
             format(file.info(CFG$api_file)$mtime, "%Y-%m-%d %H:%M")) else
     "  << NOT FOUND")
qc("INPUT pred : %s%s\n", CFG$pred_file,
   if (file.exists(CFG$pred_file))
     sprintf("  (modified %s)",
             format(file.info(CFG$pred_file)$mtime, "%Y-%m-%d %H:%M")) else
     "  << NOT FOUND")
qc("(the 02_* per-host tables are NOT read; counts taken from them will not\n")
qc(" match unless both files came from the same run)\n")
qc("API table       : %s\n", if (is.null(api)) "not found" else
   sprintf("%d rows", nrow(api)))
qc("prediction table: %s  (field: %s)\n", if (is.null(prd)) "not found" else
   sprintf("%d rows", nrow(prd)), PREDICTION_FIELD)
if (identical(CFG$env_field, "environment_primary") && !is.null(api)) {
  qc("\n!! env_field = \"environment_primary\" but the API column is multi-label.\n")
  qc("   API-covered hosts can still carry several habitats, so the two sources\n")
  qc("   are not on the same footing. Use \"environment_all\", or set\n")
  qc("   source_rule = \"prediction_first\" to make the definition uniform.\n")
}
qc("source rule     : %s\n", CFG$source_rule)

# merge host by host, not species by species: the choice of source is a
# property of the HOST lookup
key <- function(d) paste(d$species, d$host, sep = "\r")
merged <- NULL
if (!is.null(api) && !is.null(prd)) {
  ak <- key(api); pk <- key(prd)
  all_k <- union(ak, pk)
  ae <- setNames(api$env, ak)[all_k]
  pe <- setNames(prd$env, pk)[all_k]
  ae[is.na(ae)] <- ""; pe[is.na(pe)] <- ""
  env <- switch(CFG$source_rule,
    api_first        = ifelse(nzchar(ae), ae, pe),
    prediction_first = ifelse(nzchar(pe), pe, ae),
    union            = trimws(paste(ae, pe, sep = "; ")),
    stop("unknown source_rule: ", CFG$source_rule, call. = FALSE))
  # carry host_raw / host_column / cell through the merge, not just
  # species/category/host -- drop_hyperhost_rows() needs host_raw, and without
  # it the merged (default, two-source) frame silently skipped hyperhost removal
  meta_cols <- c("species", "category", "host", "host_raw", "host_column", "cell")
  meta <- rbind(api[, meta_cols], prd[, meta_cols])
  meta <- meta[!duplicated(key(meta)), ]
  merged <- data.frame(meta[match(all_k, key(meta)), ], env = unname(env),
                       stringsAsFactors = FALSE)
  qc("host rows: API only %d, prediction only %d, in both %d\n",
     length(setdiff(ak, pk)), length(setdiff(pk, ak)),
     length(intersect(ak, pk)))
  n_filled <- sum(!nzchar(ae) & nzchar(pe))
  qc("  habitat taken from the prediction because the API had none: %d\n",
     n_filled)
} else {
  merged <- if (is.null(api)) prd else api
}
merged$env <- gsub("^;\\s*|;\\s*$", "", merged$env)


## ------------------------------------------- hyperparasite chain rows -------
# See fig_host_phylum.R for the reasoning: in "X (direct host); Y (hyperhost)"
# only X is the microsporidian's host, and counting Y as well adds its habitat
# to the species. The annotation survives in host_raw, so no re-run is needed.
# Same two rules as fig_host_phylum.R, so panels C and D exclude the same
# hyperparasite-chain rows. Scoped to the CELL (species + cell, or species +
# host_column when no cell column is present), not to the whole species: a
# genus that is a genuine host in one cell must not be dropped because it is a
# hyperhost in another.
#   rule 1  an explicit hyperhost in a cell that also names a direct host/parasite
#   rule 2  the unannotated companion in a cell whose other fragment is the parasite
drop_hyperhost_rows <- function(d) {
  if (!all(c("host_raw", "species") %in% names(d))) return(d)
  raw <- ifelse(is.na(d$host_raw), "", as.character(d$host_raw))
  cellsrc <- if ("cell" %in% names(d) && any(nzchar(as.character(d$cell))))
    d$cell else if ("host_column" %in% names(d)) d$host_column else d$host
  cellkey <- paste(d$species, ifelse(is.na(cellsrc), "", as.character(cellsrc)))
  is_hyper  <- grepl("hyper[- ]?host|host of the host|carrier host", raw,
                     ignore.case = TRUE)
  is_parasi <- grepl("hyperparasit|archigregarine|\\bgregarine\\b", raw,
                     ignore.case = TRUE)
  is_direct <- grepl("direct host|primary host|true host|actual host", raw,
                     ignore.case = TRUE)
  cells_with_direct   <- unique(cellkey[is_direct | is_parasi])
  cells_with_parasite <- unique(cellkey[is_parasi])
  drop1 <- is_hyper & cellkey %in% cells_with_direct
  drop2 <- !is_parasi & !is_direct & !is_hyper & cellkey %in% cells_with_parasite
  drop  <- drop1 | drop2
  if (any(drop))
    qc("hyperparasite chains: %d host row(s) dropped (%d explicit hyperhost, "
       , sum(drop), sum(drop1))
  if (any(drop)) qc("%d unannotated companion)\n", sum(drop2 & !drop1))
  d[!drop, , drop = FALSE]
}

if (isTRUE(CFG$drop_hyperhosts)) merged <- drop_hyperhost_rows(merged)

## ---------------------------------------------------- provenance tracing ----
# Where a host count goes missing between the input tables and the diagram.
# Three things can drop it, and they are easy to confuse:
#   1. the host name did not join between the two tables
#   2. source_rule discarded one source's habitats for that host
#   3. the habitat is not in CFG$sets
cat("\n---- how each environment survives the merge ----\n")
count_env <- function(vec, e) sum(vapply(split_env(vec), function(v)
  e %in% trimws(v), logical(1)))
trace <- do.call(rbind, lapply(ENVS, function(e) {
  data.frame(environment = e,
             hosts_api = if (is.null(api)) NA_integer_ else count_env(api$env, e),
             hosts_pred = if (is.null(prd)) NA_integer_ else count_env(prd$env, e),
             hosts_after_merge = count_env(merged$env, e),
             in_CFG_sets = e %in% CFG$sets,
             stringsAsFactors = FALSE)
}))
# the most either source alone would have supported, vs what survived
trace$best_single_source <- pmax(trace$hosts_api, trace$hosts_pred, na.rm = TRUE)
trace$lost_vs_best <- pmax(0, trace$best_single_source - trace$hosts_after_merge)
print(trace, row.names = FALSE, right = FALSE)
qc("\n(counts are HOST ROWS, i.e. species-host pairs, not species)\n")
if (any(trace$lost_vs_best > 0, na.rm = TRUE)) {
  qc("\n!! source_rule = \"%s\" discards the other source's habitats for a host\n",
     CFG$source_rule)
  qc("   whenever the preferred source has ANY value. If WoRMS returns only\n")
  qc("   \"marine\" for a host the prediction calls \"marine; brackish\", the\n")
  qc("   brackish flag is dropped. Set source_rule = \"union\" to keep both.\n")
  qc("   Worst affected: %s (%d host rows fewer than the better source alone).\n",
     trace$environment[which.max(trace$lost_vs_best)], max(trace$lost_vs_best))
}
if (!all(ENVS %in% CFG$sets))
  qc("\nNOT DRAWN because they are absent from CFG$sets: %s\n",
     paste(setdiff(ENVS, CFG$sets), collapse = ", "))

if (!is.null(api) && !is.null(prd)) {
  only_api <- setdiff(api$host, prd$host)
  only_prd <- setdiff(prd$host, api$host)
  if (length(only_api) || length(only_prd)) {
    qc("\nhost names present in only one table: API %d, prediction %d\n",
       length(only_api), length(only_prd))
    if (length(only_api))
      qc("  e.g. %s\n", paste(utils::head(only_api, 4), collapse = "; "))
    qc("  these are joined on the host name, so a difference in name parsing\n")
    qc("  between the two runs will silently split a host into two rows.\n")
  }
}

if (CFG$subset != "all" && any(nzchar(merged$category))) {
  n0 <- nrow(merged)
  merged <- merged[merged$category == CFG$subset, ]
  qc("subset = %s: %d of %d host rows kept\n", CFG$subset, nrow(merged), n0)
}

## ============================================================================
## TWO UNITS OF ANALYSIS
## ============================================================================
# 1. HOSTS   -- what environments do the host species themselves occupy?
#              One row per unique host, parasites ignored. This is the figure
#              that answers "where do microsporidian hosts live", and it is what
#              a host-annotation Venn built straight from GBIF/WoRMS shows.
# 2. SPECIES -- what environments does each MICROSPORIDIAN occupy, inferred as
#              the union over all of its hosts? One row per parasite. A species
#              lands in an intersection when its hosts differ, or when one host
#              is itself flagged for several habitats.
#
# The two answer different questions and their totals are not comparable: a
# parasite with twelve brackish hosts is twelve rows in (1) and one in (2).
# CFG$unit picks which to draw; "both" writes both, which is the default.

env_sets_by <- function(id) {
  lst <- split_env(merged$env)
  out <- tapply(seq_len(nrow(merged)), id, function(i) {
    e <- unique(trimws(unlist(lst[i])))
    sort(e[nzchar(e) & e %in% CFG$sets])
  }, simplify = FALSE)
  out[lengths(out) > 0]
}

# --- self-contained four-ellipse Venn ---------------------------------------
ellipse_pts <- function(x0, y0, a, b, rot, n = 200) {
  t <- seq(0, 2 * pi, length.out = n)
  r <- rot * pi / 180
  data.frame(x = x0 + a * cos(t) * cos(r) - b * sin(t) * sin(r),
             y = y0 + a * cos(t) * sin(r) + b * sin(t) * cos(r))
}
in_ellipse <- function(px, py, x0, y0, a, b, rot) {
  r <- -rot * pi / 180
  dx <- px - x0; dy <- py - y0
  u <- dx * cos(r) - dy * sin(r)
  v <- dx * sin(r) + dy * cos(r)
  (u / a)^2 + (v / b)^2 <= 1
}
# classic layout: sets 1 and 4 outermost, 2 and 3 inner
VENN4 <- data.frame(x0 = c(0.35, 0.45, 0.55, 0.65),
                    y0 = c(0.47, 0.57, 0.57, 0.47),
                    a  = rep(0.35, 4), b = rep(0.20, 4),
                    rot = c(-35, -35, 35, 35))
VENN3 <- data.frame(x0 = c(0.40, 0.60, 0.50), y0 = c(0.60, 0.60, 0.40),
                    a = rep(0.28, 3), b = rep(0.28, 3), rot = c(0, 0, 0))
VENN2 <- data.frame(x0 = c(0.40, 0.60), y0 = c(0.50, 0.50),
                    a = rep(0.28, 2), b = rep(0.28, 2), rot = c(0, 0))

# Label positions are computed from the geometry, not hardcoded: rasterise the
# ellipses, find the pixels belonging to each exact region, and place the label
# at the point of that region furthest from its own boundary.
region_anchors <- function(par, res = 320) {
  g <- expand.grid(x = seq(0, 1, length.out = res), y = seq(0, 1, length.out = res))
  memb <- vapply(seq_len(nrow(par)), function(i)
    in_ellipse(g$x, g$y, par$x0[i], par$y0[i], par$a[i], par$b[i], par$rot[i]),
    logical(nrow(g)))
  key <- apply(memb, 1, function(r) paste(which(r), collapse = "&"))
  out <- list()
  for (k in unique(key[key != ""])) {
    idx <- which(key == k)
    cx <- mean(g$x[idx]); cy <- mean(g$y[idx])
    # furthest-from-edge point, approximated on a thinned grid
    sub <- idx[seq(1, length(idx), length.out = min(400, length(idx)))]
    edge <- which(key != k)
    edge <- edge[seq(1, length(edge), length.out = min(3000, length(edge)))]
    d <- vapply(sub, function(i)
      min((g$x[i] - g$x[edge])^2 + (g$y[i] - g$y[edge])^2), numeric(1))
    best <- sub[which.max(d)]
    out[[k]] <- c(x = g$x[best], y = g$y[best], cx = cx, cy = cy)
  }
  out
}

venn_ggplot <- function(counts, sets, pal) {
  n <- length(sets)
  par_ <- switch(as.character(n), "4" = VENN4, "3" = VENN3, "2" = VENN2,
                 stop("the built-in Venn handles 2-4 sets; got ", n,
                      call. = FALSE))
  anchors <- region_anchors(par_)

  shapes <- do.call(rbind, lapply(seq_len(n), function(i) {
    d <- ellipse_pts(par_$x0[i], par_$y0[i], par_$a[i], par_$b[i], par_$rot[i])
    d$set <- sets[i]; d
  }))
  shapes$set <- factor(shapes$set, levels = sets)

  # map "freshwater&marine" -> the numeric key "1&3" used by the anchors
  lab <- do.call(rbind, lapply(names(counts), function(k) {
    idx <- sort(match(strsplit(k, "&")[[1]], sets))
    if (anyNA(idx)) return(NULL)
    key <- paste(idx, collapse = "&")
    a <- anchors[[key]]
    if (is.null(a)) return(NULL)
    data.frame(x = a[["x"]], y = a[["y"]], n = counts[[k]], combo = k)
  }))
  missing_regions <- setdiff(names(counts), lab$combo)
  if (length(missing_regions))
    warning("no place on the diagram for: ",
            paste(missing_regions, collapse = ", "), call. = FALSE)

  # Category labels are placed geometrically: the point on each outline
  # furthest from the centre of the whole diagram, nudged a little further out.
  # Hardcoded offsets break as soon as the layout or the set order changes.
  cat_lab <- do.call(rbind, lapply(seq_len(n), function(i) {
    d <- ellipse_pts(par_$x0[i], par_$y0[i], par_$a[i], par_$b[i], par_$rot[i])
    r <- sqrt((d$x - 0.5)^2 + (d$y - 0.5)^2)
    j <- which.max(r)
    ux <- (d$x[j] - 0.5) / r[j]; uy <- (d$y[j] - 0.5) / r[j]
    off <- if (n >= 4) 0.045 else 0.075   # circles need more clearance
    data.frame(set = sets[i], x = d$x[j] + off * ux, y = d$y[j] + off * uy,
               hj = 0.5 - 0.5 * ux)
  }))
  cat_lab$set <- factor(cat_lab$set, levels = sets)

  # leave room for the labels themselves, which sit outside the outlines
  xr <- range(c(shapes$x, cat_lab$x)); yr <- range(c(shapes$y, cat_lab$y))
  # x needs more padding than y: the outer labels are left/right-justified and
  # extend well past their anchor point
  pad_x <- 0.26; pad_y <- 0.10

  ggplot() +
    geom_polygon(data = shapes, aes(x, y, group = set, fill = set),
                 alpha = 0.55, colour = "grey25", linewidth = 0.5) +
    geom_text(data = lab, aes(x, y, label = n), size = 3.1, colour = "grey10") +
    geom_text(data = cat_lab, aes(x, y, label = cap1(set), colour = set,
                              hjust = hj),
                size = 3.9, fontface = "bold") +
    scale_fill_manual(values = pal[sets]) +
    scale_colour_manual(values = pal[sets]) +
    coord_equal(xlim = xr + c(-pad_x, pad_x), ylim = yr + c(-pad_y, pad_y),
                clip = "off") +
    theme_void() +
    theme(legend.position = "none",
          plot.background = element_rect(fill = "white", colour = NA))
}


upset_ggplot <- function(counts, sets, pal, title) {
  ord <- order(-counts)
  bars <- data.frame(combo = names(counts)[ord], n = as.integer(counts)[ord])
  bars$x <- seq_len(nrow(bars))
  bars$solo <- ifelse(grepl("&", bars$combo), NA_character_, bars$combo)
  grid_df <- do.call(rbind, lapply(seq_len(nrow(bars)), function(i) {
    members <- strsplit(bars$combo[i], "&")[[1]]
    data.frame(x = i, env = sets, on = sets %in% members)
  }))
  grid_df$env <- factor(grid_df$env, levels = rev(sets))
  seg <- do.call(rbind, lapply(split(grid_df[grid_df$on, ], grid_df$x[grid_df$on]),
                               function(d) if (nrow(d) < 2) NULL else
                                 data.frame(x = d$x[1], ymin = min(as.integer(d$env)),
                                            ymax = max(as.integer(d$env)))))
  th <- theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(), legend.position = "none",
          axis.text = element_text(colour = "black"),
          plot.title = element_text(size = 10, face = "bold"))
  p1 <- ggplot(bars, aes(x = x, y = n, fill = solo)) +
    geom_col(width = 0.68) +
    scale_fill_manual(values = pal, na.value = "grey45") +
    geom_text(aes(label = n), vjust = -0.35, size = 2.9, colour = "grey20") +
    scale_x_continuous(limits = c(0.4, nrow(bars) + 0.6), expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
    labs(title = title, x = NULL, y = NULL) +
    th + theme(axis.text.x = element_blank(),
               panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3))
  p2 <- ggplot(grid_df, aes(x = x, y = env)) +
    geom_point(colour = "grey88", size = 3) +
    {if (!is.null(seg)) geom_segment(data = seg,
       aes(x = x, xend = x, y = ymin, yend = ymax), inherit.aes = FALSE,
       colour = "grey30", linewidth = 0.5)} +
    geom_point(data = grid_df[grid_df$on, ], aes(colour = env), size = 3) +
    scale_colour_manual(values = pal) +
    scale_x_continuous(limits = c(0.4, nrow(bars) + 0.6), expand = c(0, 0)) +
    labs(x = NULL, y = NULL) + th +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  list(top = p1, dots = p2, n_combo = nrow(bars))
}

render_unit <- function(entity_env, tag, unit_label) {
  if (!length(entity_env)) {
    qc("\n[%s] nothing to draw\n", tag); return(invisible(NULL))
  }
  combo <- vapply(entity_env, paste, character(1), collapse = "&")
  tabl <- sort(table(combo), decreasing = TRUE)
  counts <- setNames(as.integer(tabl), names(tabl))
  sets_present <- CFG$sets[CFG$sets %in% unlist(strsplit(names(counts), "&"))]
  n_env <- lengths(entity_env)

  cat(sprintf("\n================ %s ================\n", toupper(unit_label)))
  qc("%s with at least one environment : %d\n", unit_label, length(entity_env))
  for (k in sort(unique(n_env)))
    qc("  in %d environment%s : %4d (%.1f%%)\n", k, ifelse(k == 1, " ", "s"),
       sum(n_env == k), 100 * mean(n_env == k))
  per_env <- vapply(CFG$sets, function(e)
    sum(vapply(entity_env, function(v) e %in% v, logical(1))), integer(1))
  qc("\nper environment: %s\n",
     paste(sprintf("%s %d", names(per_env), per_env), collapse = ", "))
  cat("\nexact combinations:\n")
  print(data.frame(combination = names(counts), n = as.integer(counts),
                   pct = sprintf("%.1f%%", 100 * counts / length(entity_env))),
        row.names = FALSE, right = FALSE)

  write.table(data.frame(combination = names(counts), n = as.integer(counts)),
              file.path(CFG$outdir, sprintf("combination_counts_%s.tsv", tag)),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(data.frame(entity = names(entity_env),
                         environments = vapply(entity_env, paste, character(1),
                                               collapse = "; "),
                         n_environments = lengths(entity_env)),
              file.path(CFG$outdir, sprintf("%s_environments.tsv", tag)),
              sep = "\t", row.names = FALSE, quote = FALSE)

  pv <- venn_ggplot(counts, sets_present, PAL)
  ggsave(file.path(CFG$outdir, sprintf("fig_venn_%s.pdf", tag)), pv,
         width = CFG$fig_width, height = CFG$fig_height,
         device = if (capabilities("cairo")) cairo_pdf else pdf)
  ggsave(file.path(CFG$outdir, sprintf("fig_venn_%s.png", tag)), pv,
         width = CFG$fig_width, height = CFG$fig_height, dpi = 600)
  qc("\nwrote fig_venn_%s.pdf / .png (%d sets%s)\n", tag, length(sets_present),
     if (length(sets_present) == 4) ", four-ellipse layout, not area-proportional"
     else "")

  if (has_patchwork) {
    library(patchwork)
    u <- upset_ggplot(counts, sets_present, PAL,
                      sprintf("%s by environmental range", unit_label))
    up <- u$top / u$dots + plot_layout(heights = c(2.6, 1))
    ggsave(file.path(CFG$outdir, sprintf("fig_upset_%s.pdf", tag)), up,
           width = max(5, 0.55 * u$n_combo + 2.6), height = 4.6,
           device = if (capabilities("cairo")) cairo_pdf else pdf)
    ggsave(file.path(CFG$outdir, sprintf("fig_upset_%s.png", tag)), up,
           width = max(5, 0.55 * u$n_combo + 2.6), height = 4.6, dpi = 600)
    qc("wrote fig_upset_%s.pdf / .png (exact, no geometric distortion)\n", tag)
  }
  invisible(counts)
}

## ------------------------------------------------------------------ run ----
host_env <- env_sets_by(merged$host)
sp_env   <- env_sets_by(merged$species)

n_host_total <- length(unique(merged$host))
n_sp_total   <- length(unique(merged$species))
qc("\nunique hosts %d (with an environment: %d)\n", n_host_total, length(host_env))
qc("microsporidia species %d (with an environment: %d)\n", n_sp_total,
   length(sp_env))

if (CFG$unit %in% c("host", "both"))
  render_unit(host_env, "hosts", "host species")
if (CFG$unit %in% c("species", "both"))
  render_unit(sp_env, "microsporidia", "microsporidia species")

qc("\nwritten to %s/\n", CFG$outdir)
qc("  fig_venn_hosts.*          environments of the HOST species\n")
qc("  fig_venn_microsporidia.*  environments of each MICROSPORIDIAN\n")
qc("  fig_upset_*.*             the same counts without geometric distortion\n")
qc("  combination_counts_*.tsv, *_environments.tsv\n")
