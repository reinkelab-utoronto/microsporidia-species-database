## =============================================================================
## Circular 18S tree with clade highlighting that tolerates FRAGMENTED clades.
##
## Two kinds of clade definition:
##
##   MAJOR clades  - given as a pair of boundary tips, highlighted as the single
##                   MRCA clade (the original behaviour, unchanged).
##
##   MINOR clades  - given as a LIST OF MEMBER ACCESSIONS. These are not
##                   monophyletic in an 18S tree, so instead of one MRCA the
##                   script finds every MAXIMAL PURE subclade (a node whose
##                   descendants are all members of that clade and nothing
##                   else) and highlights each piece in the same colour.
##                   Single-accession clades fall out of this automatically as
##                   one-tip pieces.
##
## It also reports, per minor clade, how many pieces it is broken into and how
## many foreign tips sit inside its overall MRCA -- i.e. how fragmented it is.
## =============================================================================

library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
library(tidytree)
library(treeio)
library(grid)

## capitalise the first letter, display only ("brackish" -> "Brackish")
cap1 <- function(x) sub("^(\\w)", "\\U\\1", as.character(x), perl = TRUE)

## ---- settings ---------------------------------------------------------------
TREE_FILE <- "CENTROIDS_WITH_HOSTS_8_pruned.treefile"
OUTGROUP  <- "MF278562.1"
OUT_STEM  <- "18S_tree_circular_with_legends"  # a timestamp is appended, so
                                               # each run writes new files and
                                               # nothing is silently overwritten
OUT_DIR   <- "."          # where the figures go
OUT_PNG   <- TRUE         # also write a PNG alongside the PDF
PNG_DPI   <- 600          # ASM allows up to 600
## mBio / ASM: recommended figure dimensions are 7 x 9 inches or less at
## 300 dpi. 7 inches is the maximum printed width, so the figure is built at
## exactly that and everything else is scaled to it.
## LAYOUT, derived rather than guessed (see the arithmetic below).
## ASM allows 7 x 9 in. The tree is a circle, so its diameter is the SHORTER
## side of its panel. Giving it the full 7 in width means it also needs 7 in of
## height, which leaves 9 - 7 = 2 in for a legend band underneath. That is the
## largest possible tree inside the ASM box: 7.00 in across, versus 4.62 in
## when the legend sat in a side column.
##
## Legend capacity at 6 pt with 1.55x leading: 0.129 in per row, so 15 rows in
## the 2 in band. A worst-case legend is ~56 rows (3 status + 10 clades +
## 6 environment + 10 phyla + 16 orders + headers and notes), so 4 columns
## (60 slots) suffice; 5 leaves headroom. The longest label,
## "Eupercaria incertae sedis", measures 0.60 in at 6 pt and fits a 1.40 in
## column (5 columns) with room for the swatch.
FIG_W       <- 7          # inches -- ASM maximum width
FIG_H       <- 9          # inches -- ASM maximum height
## The tree diameter is min(FIG_W, FIG_H - LEGEND_IN). At LEGEND_IN = 2.6 the
## tree was HEIGHT-limited (6.4 in available against 7 in of width), so packing
## the legend into a shorter band directly buys tree size. 2.0 in makes the two
## limits equal at 7.0 in, which is the most the ASM box allows.
LEGEND_IN   <- 2.0        # height of the legend band beneath the tree
LEGEND_COLS <- 6          # columns across that band (7/6 = 1.17 in each; the
                          # longest label is 0.60 in at 6 pt, so it fits)
ORDER_GROUP_NEW_COL <- TRUE
                          # start each host-order group (Chordata, Arthropoda)
                          # in its own column, so a group is never split across
                          # a column break
LEGEND_FS   <- 6          # legend body font size (pt) -- ASM floor
TREE_MARGIN <- -0.80      # inches. ggtree's circular layout leaves the circle
                          # filling only ~80% of its panel, and neither
                          # scale_*_continuous(expand=) nor coord clipping
                          # changes that. Negative plot margins push the panel
                          # outside the viewport, scaling the circle up.
                          # Measured by rendering and counting non-white pixels
                          # on the border rows/columns, then bisecting:
                          #   -0.70 -> 96.4% fill, no clipping
                          #   -0.80 -> 98.7% fill, no clipping  <- used
                          #   -0.83 -> 99.3% fill, 15 px clipped
                          #   -0.90 -> 99.9% fill, 204 px clipped
                          # Raise it only if you also shrink HOSTSYM_MAX or
                          # HOSTSYM_GAP, since the symbol band sets the limit.

## A second, separate PDF with the tip NAMES readable. All ring annotations are
## dropped except the clade ring, because at a legible label size the figure has
## to be very large and the extra tracks only add clutter -- this version is for
## looking up which taxon sits where, not for presenting the data.
OUT_LABELS   <- TRUE
LABEL_SIZE   <- 1.2       # tip label text size
LABEL_FIG_W  <- 26        # inches; labelled version needs far more room
LABEL_FIG_H  <- 28  # 2 in taller than wide: the legend strip
                          # sits below a full-width circle
LABEL_RING_W <- 0.02      # clade ring thickness in the labelled version
LABEL_TREE_MARGIN <- -4.4   # inches, negative: inflates the circle past
                          # CoordPolar's built-in ~20% rim, as TREE_MARGIN
                          # does for the main figure. Tuned by rendering and
                          # measuring; raise toward 0 if labels clip.
LABEL_LEGEND_IN   <- 1.4  # height of the legend band under the tree
LABEL_LEGEND_COLW <- 2.6  # inches per legend column (5 columns x 2 rows)
LABEL_LEGEND_FS   <- 14   # legend font size (pt)
LABEL_XLIM   <- 1.25      # radial room for the labels. 1.60 left a third of
                          # the panel radius empty, shrinking the circle and
                          # pushing the legend ~5 in below the tips; 1.25 just
                          # clears "Species name (ACC123.1)" labels at size 1.2.
                          # Raise it if the longest names clip at the rim.
CLADOGRAM <- TRUE       # TRUE = branch.length "none" (equidistant tips).
                        # NOTE: the rings assume equidistant tips. With real
                        # branch lengths, tips sit at different radii and the
                        # rings detach from them.
TIPLAB    <- FALSE      # tip labels collide with the rings; off by default
TIPLAB_SZ <- 0.5

## rings (fractions of the tree's x range)
RING_WIDTH <- 0.022     # radial thickness of each ring
RING_GAP   <- 0.012     # gap between the tips and the first ring
RING_SPACE <- 0.006     # gap between successive rings

## environment rings, drawn inside-out in this order
## Environments are resolved PER HOST and then unioned over every host in the
## cluster -- not looked up by parasite name. This is what lets provisional
## clusters ("Enterocytozoon sp. MSP-C0042") carry environments: their hosts
## are known even though the parasite name appears in no upstream table.
## TWO ROUTES, and the union of both is used:
##   route A (hosts)   cluster -> all_hosts -> per-host habitat. Works for
##                     provisional clusters, whose parasite name is in no table.
##   route B (species)  cluster_label -> microsporidia_environments.tsv. Works
##                     for named species and is the route the euler figure uses.
##                     APPLIED ONLY to the cluster that owns the curated database
##                     record (owns_database_record / attribute_basis in the
##                     attribute table): the lookup is by species name, so a
##                     species split across clusters would otherwise stamp its
##                     full habitat union onto every piece. Non-owning clusters
##                     use route A alone -- the habitats of their own GenBank
##                     hosts -- matching how make_cluster_table.py credits the
##                     curated record to a single cluster.
## Using both means a path failing cannot blank the rings on its own.
ATTR_FILE  <- "cluster_output/03_CLUSTER_HOST_ATTRIBUTES.tsv"  # cluster -> hosts
SPECIES_ENV_FILE <- "environment_euler/microsporidia_environments.tsv"
HOST_ENV_FILES <- c("host_taxonomy_output/01_species_host_long.tsv",
                    "predicted_environment/01_species_host_predicted.tsv")
HOST_ENV_RULE  <- "api_first"   # "api_first" | "prediction_first" | "union"
                                # matches source_rule in fig_environment_euler.R
ENV_ORDER  <- c("terrestrial", "freshwater", "marine", "brackish")
ENV_COLORS <- c(terrestrial = "#7F4F24",   # Okabe-Ito, matching
                freshwater  = "#56B4E9",   # fig_environment_euler.R
                marine      = "#0072B2",
                brackish    = "#009E73")
## named vs provisional (branch colour)
ENTITY_LABELS <- c("Named species", "Provisional species")
ENTITY_COLORS <- c("Named species" = "#0072B2",
                   "Provisional species" = "#E69F00")
## Some Natural Host(s) entries flag an organism as a hyperhost -- the host of
## the host, not of the microsporidian, e.g.
##   "Pygospio elegans (polychaete hyperhost)"
## Plotting these as host taxa would inflate host breadth and put the wrong
## phylum symbol on a tip, so they are dropped before any host is used.
DROP_HYPERHOSTS <- TRUE
HYPERHOST_RE    <- "hyper[ -]?(host|parasit)"

## Some hyperhosts are not annotated as such in the source data. Rather than
## guessing from the host name -- which risks deleting genuine records for a
## species that is a real host elsewhere -- name the exceptions explicitly:
## for these microsporidian species, hosts in the listed phyla are not plotted.
## Metchnikovellids are hyperparasites of gregarines living in polychaetes, so
## the polychaete is the hyperhost, not the host.
EXEMPT_HOST_PHYLUM <- list(
  "Metchnikovella incurvata" = "Annelida",
  "Metchnikovella spiralis"  = "Annelida"
)

ENTITY_NA_COL <- "grey88"   # tips with no entity_type
BRANCH_COL    <- "grey20"   # branches carry no colour at all; every category
                            # in the figure is encoded by a ring

## ---- host taxonomy symbols ---------------------------------------------------
## Plotted at PHYLUM level, except Chordata and Arthropoda which are resolved to
## ORDER. A cluster may infect hosts in several groups, so the symbols for one
## tip are STACKED radially outward rather than overplotted.
HOST_TAX_FILE <- "host_taxonomy_output/01_species_host_long.tsv"
SPLIT_TO_ORDER <- c("Chordata", "Arthropoda")   # these go to order level
HOSTSYM_MAX    <- 10        # max symbols stacked per tip; extras are counted
HOSTSYM_SIZE   <- 0.94     # symbol size (a further 10% down from 1.04). Below ~1.5 only the first five
                           # PHYLUM_SHAPES (circle, square, triangle, asterisk,
                           # plus) stay clearly distinguishable -- the compound
                           # glyphs later in the pool turn into small blobs. If
                           # you have more than five phyla, raise this or raise
                           # HOSTTAX_MIN_N to pool the tail.
HOSTSYM_GAP    <- 0.032    # radial spacing between stacked symbols
HOSTSYM_START  <- 0.038    # gap between the last ring and the first symbol
HOSTTAX_MIN_N  <- 5        # phyla AND orders found on fewer than this many
                           # clusters are pooled into "Other phylum" /
                           # "other order". Counted in CLUSTERS (tips), not in
                           # symbols, so one widespread host does not inflate a
                           # group and a cluster with two Diptera counts once.
COLOUR_BY_PHYLUM_IF_NO_ORDERS <- TRUE
                           # if no order can be resolved at all, colour the
                           # symbols by phylum instead of leaving them all grey

## PHYLUM sets the SHAPE, ORDER sets the COLOUR (Chordata / Arthropoda only).
## Solid glyphs are used so that mapping order to `colour` fills the whole
## symbol -- the `fill` aesthetic is already taken by the six rings, and two
## fill scales cannot coexist in one ggplot without ggnewscale.
## SHAPE SET -- two constraints, both learned by rendering and looking:
##
## 1. NO ROTATION-CONJUGATE PAIRS. geom_point glyphs are drawn aligned to the
##    PAGE, not to the ring, and ggplot2 cannot rotate them per point. Around a
##    circular tree the ring runs at every angle, so a square appears tilted
##    where the ring is diagonal. That is only a PROBLEM if the set also
##    contains a diamond -- the two become confusable. So a square is fine on
##    its own; a square AND a diamond are not. The conjugate pairs are
##    square/diamond (15,0 vs 18,5), plus/cross (3 vs 4), and circle-plus /
##    circle-cross (10 vs 13). Only one member of each pair is ever used.
##
## 2. DISTINCT BY OUTLINE AT SMALL SIZE. Compound glyphs (circle-X, square-X,
##    hexagram, square-triangle) turn into indistinguishable blobs below about
##    size 1.5, and pairs that differ only in FILL (open vs filled triangle) or
##    in SIZE (pch 16 vs 19 vs 20) are not distinguishable either. Those are
##    all excluded, which leaves a small set -- five shapes that genuinely read
##    apart at this density. Colour carries the rest of the information, and
##    HOSTTAX_MIN_N pools the tail so the set is rarely exhausted.
PHYLUM_SHAPES <- c(Arthropoda      = 16,   # filled circle
                   Chordata        = 15,   # filled square
                   Nematoda        = 17,   # filled triangle
                   Mollusca        = 8,    # asterisk
                   Annelida        = 3,    # plus
                   Platyhelminthes = 13,   # circle-cross (needs a larger size)
                   Rotifera        = 7,    # square-cross
                   Bryozoa         = 11,   # hexagram
                   Cnidaria        = 14,   # square-triangle
                   Other           = 1)    # open circle
## first five read clearly at any size; the rest only above ~1.5
PHYLUM_POOL   <- c(16, 15, 17, 8, 3, 13, 7, 11, 14, 1)

NOORDER_COL   <- "grey25"  # colour for anything unassigned

## A qualitative palette of maximally distinct colours (Kelly / Tol lineage).
## Sequential ramps like hcl.colors("Reds 3") are unusable here: with 9
## Arthropoda and 5 Chordata orders the neighbouring steps differ by a few
## percent of lightness and cannot be told apart at symbol size.
DISTINCT_PAL <- c(
  "#E6194B", "#3CB44B", "#4363D8", "#F58231", "#911EB4", "#42D4F4",
  "#F032E6", "#BFEF45", "#FABED4", "#469990", "#DCBEFF", "#9A6324",
  "#800000", "#AAFFC3", "#808000", "#FFD8B1", "#000075", "#A9A9A9",
  "#FFE119", "#00A5A5", "#B03060", "#7CFC00", "#6A5ACD", "#FF8C69")

ABSENT_COL <- "grey88"  # habitat recorded, species not in it
NODATA_COL <- "grey96"  # no habitat record for this species

## ---- read and root ----------------------------------------------------------
tree <- read.tree(TREE_FILE)
cat("Tree has", Ntip(tree), "tips\n")

if (!OUTGROUP %in% tree$tip.label) {
  cat("\n*** OUTGROUP", OUTGROUP, "not found. Candidates present:\n")
  print(grep("Metchnikovella|Amphiamblys|Rozella|Paramicrosporidium|Mitosporidium",
             tree$tip.label, value = TRUE, ignore.case = TRUE))
  stop("Fix OUTGROUP and re-run.")
}
tree <- root(tree, outgroup = OUTGROUP, resolve.root = TRUE)

## ---- cluster -> hosts -> environments ---------------------------------------
## Two joins, both keyed on things that exist for every cluster:
##   tree tip (centroid accession) -> all_hosts   (03_CLUSTER_HOST_ATTRIBUTES)
##   host name -> habitat flags                   (01_* per-host tables)
## The parasite's own name is never used, so named and provisional clusters are
## treated identically.
norm_host <- function(x) tolower(trimws(gsub("\\s+", " ", x)))

## Split a "; "-joined host list and drop any entry flagged as a hyperhost.
## Applied everywhere a host list is consumed -- environment rings and host
## taxonomy symbols both read the same lists, so filtering in one place keeps
## them consistent.
## Split a "; "-joined host list and drop any entry annotated as a hyperhost.
split_hosts <- function(x) {
  v <- trimws(unlist(strsplit(as.character(x), ";")))
  v <- v[nzchar(v)]
  if (DROP_HYPERHOSTS && length(v)) {
    drop <- grepl(HYPERHOST_RE, v, ignore.case = TRUE)
    if (any(drop)) {
      hyper_dropped <<- unique(c(hyper_dropped, v[drop]))
      v <- v[!drop]
    }
  }
  v
}
hyper_dropped <- character(0)

## --- cluster -> hosts (RESTORED) ---------------------------------------------
## attrs, attr_idx, tip_hosts and tip_label_sp are all consumed further down;
## an earlier edit removed the block that creates them, which surfaced as a
## length mismatch rather than "object not found" because a stale copy from a
## previous session was still in the workspace.
attrs <- NULL
attr_idx     <- rep(NA_integer_, Ntip(tree))
tip_hosts    <- rep(NA_character_, Ntip(tree))
tip_label_sp <- rep(NA_character_, Ntip(tree))
## Does this tip's cluster OWN the curated database record for its species?
## Set from 03_CLUSTER_HOST_ATTRIBUTES.tsv below. Only the owning cluster may
## pull curated (species-name) environments via route B; the others get GenBank
## hosts only. Default TRUE so an old attribute table lacking the column keeps
## the previous behaviour rather than blanking rings (a NOTE is printed then).
tip_owns     <- rep(TRUE, Ntip(tree))
acc_col <- NA

if (file.exists(ATTR_FILE)) {
  attrs <- read.delim(ATTR_FILE, stringsAsFactors = FALSE, na.strings = "",
                      quote = "\"")
  acc_col <- intersect(c("centroid_accession", "centroid"), names(attrs))[1]
  h_col   <- intersect(c("all_hosts", "hosts"), names(attrs))[1]
  if (is.na(acc_col) || is.na(h_col)) {
    cat("*** ", ATTR_FILE, " lacks a centroid or host column ***\n", sep = "")
    attrs <- NULL
  } else {
    ## match on the accession as written, then with any vsearch "NNN|" prefix
    ## removed, then without the version suffix
    strip <- function(x) sub("^[0-9]+\\|", "", x)
    base  <- function(x) sub("\\.[0-9]+$", "", strip(x))
    i <- match(strip(tree$tip.label), attrs[[acc_col]])
    j <- match(base(tree$tip.label),  base(attrs[[acc_col]]))
    i[is.na(i)] <- j[is.na(i)]
    attr_idx  <- i
    tip_hosts <- attrs[[h_col]][i]
    lab_col <- intersect(c("cluster_label", "species_name"), names(attrs))[1]
    if (!is.na(lab_col)) tip_label_sp <- attrs[[lab_col]][i]
    ## Ownership: cluster_to_database.py assigns each database row's curated
    ## attributes to ONE cluster (the one holding most of its accessions) and
    ## records that in owns_database_record (species names it owns) / the
    ## attribute_basis string. A species split across clusters -- Vittaforma
    ## corneae spans many -- therefore has ONE owning cluster; the rest carry
    ## only their own GenBank member hosts. Route B below honours this so the
    ## curated species-level environment is not copied onto every same-named tip.
    own_col <- intersect(c("owns_database_record", "attribute_basis"),
                         names(attrs))[1]
    if (!is.na(own_col)) {
      ov <- attrs[[own_col]][i]; ov <- ifelse(is.na(ov), "", ov)
      tip_owns <- if (own_col == "owns_database_record") nzchar(ov)
                  else grepl("owns the database record", ov)
    } else {
      cat("NOTE: ", ATTR_FILE, " has neither owns_database_record nor ",
          "attribute_basis;\n  route B cannot respect record ownership, so a ",
          "curated species-level habitat\n  may appear on every cluster that ",
          "shares the species name. Re-run\n  cluster_to_database.py to add the ",
          "column.\n", sep = "")
    }
    cat("Cluster table: ", sum(!is.na(i)), " of ", Ntip(tree),
        " tips matched in ", ATTR_FILE, "\n", sep = "")
    if (any(is.na(i)))
      cat("  unmatched e.g.: ",
          paste(head(tree$tip.label[is.na(i)], 3), collapse = ", "),
          "\n    (an outgroup added to the tree but absent from the cluster\n",
          "     table is expected here -- it will show as unclassified)\n",
          sep = "")
  }
} else {
  cat("NOTE: ", ATTR_FILE, " not found -- rings will be blank\n", sep = "")
}

## ---- host taxonomy per cluster ----------------------------------------------
## host name -> taxon label, then unioned over every host in the cluster. Same
## host-level join as the environment rings, so provisional clusters work too.
host_taxon <- character(0)
if (file.exists(HOST_TAX_FILE)) {
  ht <- read.delim(HOST_TAX_FILE, stringsAsFactors = FALSE, na.strings = "",
                   quote = "")
  hcol <- intersect(c("host_name", "host"), names(ht))[1]
  pcol <- intersect(c("phylum_normalised", "phylum"), names(ht))[1]
  ocol <- intersect(c("order", "order_normalised", "class"), names(ht))[1]
  if (is.na(hcol) || is.na(pcol)) {
    cat("*** ", HOST_TAX_FILE, " lacks host/phylum columns (found: ",
        paste(names(ht), collapse = ", "), ") ***\n", sep = "")
  } else {
    ph <- ifelse(is.na(ht[[pcol]]) | !nzchar(ht[[pcol]]), NA, ht[[pcol]])
    if (pcol == "phylum_normalised" && "phylum" %in% names(ht))
      ph <- ifelse(is.na(ph), ht$phylum, ph)
    ord <- if (!is.na(ocol)) ht[[ocol]] else rep(NA_character_, nrow(ht))
    ## Order is recorded only for the phyla being split; everything else keeps
    ## a blank order and is drawn in one neutral colour.
    ord <- ifelse(ph %in% SPLIT_TO_ORDER & !is.na(ord) & nzchar(ord), ord, "")
    keyv <- norm_host(ht[[hcol]])
    ok <- nzchar(keyv) & !is.na(ph) & nzchar(ph)
    host_phylum <- setNames(ph[ok], keyv[ok])
    host_order  <- setNames(ord[ok], keyv[ok])
    host_phylum <- host_phylum[!duplicated(names(host_phylum))]
    host_order  <- host_order[!duplicated(names(host_order))]
    host_taxon  <- host_phylum          # presence test used elsewhere
    ## build the message first: a multi-line if/else inside cat() is not
    ## parseable, because R closes the call at the end of the first branch
    ord_msg <- if (is.na(ocol)) {
      " (no order column -- all one colour)"
    } else {
      paste0("; ", sum(nzchar(host_order)), " with an order (", ocol, ")")
    }
    cat("Host taxonomy: ", length(host_phylum), " hosts with a phylum",
        ord_msg, "\n", sep = "")
  }
} else {
  cat("NOTE: ", HOST_TAX_FILE, " not found -- host symbols omitted\n", sep = "")
}

host_phylum <- if (exists("host_phylum")) host_phylum else character(0)
host_order  <- if (exists("host_order"))  host_order  else character(0)

exempt_log <- NULL
tip_taxa <- vector("list", Ntip(tree))   # each element: data.frame(phylum, order)
for (i in seq_len(Ntip(tree))) {
  hs <- tip_hosts[i]
  tip_taxa[[i]] <- data.frame(phylum = character(0), order = character(0))
  if (is.na(hs) || !nzchar(hs) || !length(host_phylum)) next
  hv <- norm_host(split_hosts(hs))
  hv <- hv[hv %in% names(host_phylum)]
  if (!length(hv)) next
  df <- data.frame(host   = hv,
                   phylum = unname(host_phylum[hv]),
                   order  = unname(host_order[hv]), stringsAsFactors = FALSE)

  ## named exemptions: for these microsporidian species, hosts in the listed
  ## phyla are hyperhosts and are not plotted
  sp <- tip_label_sp[i]
  if (!is.na(sp) && nzchar(sp) && length(EXEMPT_HOST_PHYLUM)) {
    ex <- EXEMPT_HOST_PHYLUM[[sp]]
    if (is.null(ex)) {                     # tolerate label variants
      k <- which(tolower(trimws(names(EXEMPT_HOST_PHYLUM))) ==
                 tolower(trimws(sp)))
      if (length(k)) ex <- EXEMPT_HOST_PHYLUM[[k[1]]]
    }
    if (!is.null(ex)) {
      hit <- df$phylum %in% ex
      if (any(hit)) {
        exempt_log <<- rbind(exempt_log, data.frame(
          species = sp, dropped_host = df$host[hit],
          phylum = df$phylum[hit], stringsAsFactors = FALSE))
        df <- df[!hit, , drop = FALSE]
      }
    }
  }
  if (!nrow(df)) next

  df <- df[!duplicated(paste(df$phylum, df$order)), , drop = FALSE]
  tip_taxa[[i]] <- df[order(df$phylum, df$order),
                      c("phylum", "order"), drop = FALSE]
}
n_tx <- vapply(tip_taxa, nrow, integer(1))
if (!is.null(exempt_log) && nrow(exempt_log)) {
  cat("Host(s) excluded by EXEMPT_HOST_PHYLUM (hyperhosts):\n")
  print(unique(exempt_log), row.names = FALSE)
} else if (length(EXEMPT_HOST_PHYLUM)) {
  cat("EXEMPT_HOST_PHYLUM matched nothing. The names must equal the\n",
      "  cluster_label exactly; present labels include:\n   ",
      paste(head(unique(tip_label_sp[!is.na(tip_label_sp)]), 4),
            collapse = "; "), "\n", sep = "")
}
cat("Tips with >=1 host taxon: ", sum(n_tx > 0), " of ", Ntip(tree), "\n", sep = "")
if (DROP_HYPERHOSTS && length(hyper_dropped)) {
  cat("Hyperhosts excluded (", length(hyper_dropped), "):\n", sep = "")
  for (h in hyper_dropped) cat("   ", h, "\n", sep = "")
} else if (DROP_HYPERHOSTS) {
  cat("Hyperhosts excluded: none matched /", HYPERHOST_RE, "/\n", sep = "")
}
if (any(n_tx > HOSTSYM_MAX))
  cat("  ", sum(n_tx > HOSTSYM_MAX), " tip(s) have more than ", HOSTSYM_MAX,
      " taxa; only the first ", HOSTSYM_MAX,
      " are drawn (all kept in clade_assignments.csv)\n", sep = "")

## ---- tip classification: named vs provisional species -----------------------
## entity_type comes from 03_CLUSTER_HOST_ATTRIBUTES.tsv, so the figure uses the
## same named/provisional call as the database outputs.
##
## Shown as the innermost RING, not as branch colour: propagating the label up
## the tree would imply that an ancestral lineage is "named" or "provisional",
## which is meaningless -- the distinction is a property of the terminal taxa.
tip_entity <- rep(NA_character_, Ntip(tree))
if (!is.null(attrs) && "entity_type" %in% names(attrs)) {
  tip_entity <- attrs$entity_type[attr_idx]
} else {
  cat("*** entity_type not available -- tip points will be grey ***\n")
}
tip_entity <- ifelse(tip_entity == "named", ENTITY_LABELS[1],
              ifelse(tip_entity == "provisional", ENTITY_LABELS[2], NA))

cat("\nTip classification: ",
    sum(tip_entity == ENTITY_LABELS[1], na.rm = TRUE), " named, ",
    sum(tip_entity == ENTITY_LABELS[2], na.rm = TRUE), " provisional, ",
    sum(is.na(tip_entity)), " unassigned\n", sep = "")

## ---- MAJOR clades: pairs of boundary tips -----------------------------------
major_boundaries <- list(
  Glugeida          = c("MT462180.1", "KU302782.1"),
  Enterocytozoonida = c("MZ079008.1", "FJ756174.1"),
  Ovavesiculida     = c("KX360153.1", "EF564602.1"),
  Nosematida        = c("PV877882.1", "FJ794881.1"),
  Amblyosporida     = c("AY090051.1", "AY090068.1"),  
  Neopereziida      = c("ON054957.1", "OK112607.1"),
  Metchnikovella    = c("KX214672.1", "PP057786.1")
)

## ---- MINOR clades: explicit member lists, may be fragmented -----------------
minor_members <- list(
  Astathelohaniida = c("FJ755981.1", "AJ438963.1", "AF492593.1", "OM630069.1",
                       "OM630067.1", "AY183664.1", "AF294779.1", "AJ302319.1",
                       "OP859151.1", "PQ670044.1", "GQ206147.1", "KC111802.1",
                       "PV577712.1"),
  Caudosporida     = c("AM411639.1", "OL583677.1", "KR704917.1", "AY973624.1",
                       "AF132544.1", "KC855552.1", "AY090069.1"),
  Areosporida      = c("PP852508.1")
)

## Clade colours deliberately AVOID the four environment colours below
## (#E69F00, #56B4E9, #0072B2, #009E73). Reusing them makes the clade ring
## indistinguishable from the terrestrial / marine rings sitting next to it.
clade_colors <- c(
  Glugeida          = "#882255",   # wine
  Astathelohaniida  = "#88CCEE",   # pale cyan
  Caudosporida      = "#44AA99",   # teal
  Areosporida       = "#AA4499",   # purple
  Enterocytozoonida = "#DDCC77",   # sand
  Ovavesiculida     = "#332288",   # indigo
  Nosematida        = "#D55E00",   # vermilion
  Amblyosporida     = "#CC79A7",   # pink
  Neopereziida      = "#000000",   # black
  Metchnikovella    = "#999999"    # grey
)

## ---- helpers ----------------------------------------------------------------
## Descendant tip indices for every node, computed once from the edge matrix.
## Avoids depending on offspring() semantics across tidytree versions.
build_descendants <- function(tree) {
  n_tip  <- Ntip(tree)
  n_node <- n_tip + tree$Nnode
  kids   <- split(tree$edge[, 2], tree$edge[, 1])
  desc   <- vector("list", n_node)
  for (i in seq_len(n_tip)) desc[[i]] <- i
  ## post-order: process nodes in decreasing depth
  ord <- rev(unique(c(tree$edge[, 1])))
  ord <- ord[order(node.depth(tree)[ord])]
  for (nd in ord) {
    ch <- kids[[as.character(nd)]]
    desc[[nd]] <- sort(unlist(desc[ch], use.names = FALSE))
  }
  desc
}

parent_of <- function(tree) {
  pv <- integer(Ntip(tree) + tree$Nnode)
  pv[tree$edge[, 2]] <- tree$edge[, 1]
  pv
}

## Maximal nodes whose descendant tips are ALL in target_tips.
maximal_pure_nodes <- function(tree, target_tips, desc, par) {
  tset <- target_tips
  pure <- which(vapply(seq_along(desc), function(i) {
    d <- desc[[i]]
    length(d) > 0 && all(d %in% tset)
  }, logical(1)))
  ## drop any node whose parent is also pure -> keeps only the maximal ones
  pure[!(par[pure] %in% pure)]
}

## ---- resolve MAJOR clades ---------------------------------------------------
desc <- build_descendants(tree)
par  <- parent_of(tree)
n_tip <- Ntip(tree)

clade_pieces <- list()   # clade -> integer vector of nodes to highlight
clade_report <- list()

cat("\n=== MAJOR clades (MRCA of boundary pair) ===\n")
for (nm in names(major_boundaries)) {
  pr <- major_boundaries[[nm]]
  miss <- pr[!pr %in% tree$tip.label]
  if (length(miss)) {
    cat("SKIPPED", nm, "- missing tip(s):", paste(miss, collapse = ", "), "\n")
    next
  }
  nd <- getMRCA(tree, pr)
  clade_pieces[[nm]] <- nd
  n_in <- length(desc[[nd]])
  cat(sprintf("%-18s node %-5s %4d tips\n", nm, nd, n_in))
  clade_report[[nm]] <- data.frame(clade = nm, type = "major", pieces = 1L,
                                   members_given = 2L, tips_covered = n_in,
                                   foreign_in_mrca = NA_integer_)
}

## ---- resolve MINOR clades (fragmented) --------------------------------------
cat("\n=== MINOR clades (maximal pure subclades) ===\n")
for (nm in names(minor_members)) {
  accs <- minor_members[[nm]]
  miss <- accs[!accs %in% tree$tip.label]
  present <- accs[accs %in% tree$tip.label]
  if (length(miss))
    cat("  NOTE", nm, "- not in tree:", paste(miss, collapse = ", "), "\n")
  if (!length(present)) {
    cat("SKIPPED", nm, "- no members present\n")
    next
  }
  idx <- match(present, tree$tip.label)
  nodes <- maximal_pure_nodes(tree, idx, desc, par)
  clade_pieces[[nm]] <- nodes

  ## how contaminated would a single MRCA be?
  foreign <- NA_integer_
  if (length(idx) > 1) {
    mrca <- getMRCA(tree, idx)
    foreign <- length(setdiff(desc[[mrca]], idx))
  }
  n_internal <- sum(nodes > n_tip)
  n_single   <- sum(nodes <= n_tip)
  cat(sprintf("%-18s %2d members -> %2d piece(s) (%d subclade, %d single tip); ",
              nm, length(present), length(nodes), n_internal, n_single))
  cat(if (is.na(foreign)) "single member\n"
      else sprintf("MRCA would include %d foreign tips\n", foreign))

  clade_report[[nm]] <- data.frame(
    clade = nm, type = "minor", pieces = length(nodes),
    members_given = length(accs),
    tips_covered = length(unlist(desc[nodes], use.names = FALSE)),
    foreign_in_mrca = foreign)
}

report <- bind_rows(clade_report)
cat("\n")
print(report, row.names = FALSE)
write.csv(report, "clade_fragmentation_report.csv", row.names = FALSE)

## ---- plot -------------------------------------------------------------------
p <- ggtree(tree, layout = "circular",
            branch.length = if (CLADOGRAM) "none" else "branch.length",
            size = 0.25, colour = BRANCH_COL) +
  theme_tree() +
  theme(legend.position = "none",
        plot.margin = margin(TREE_MARGIN, TREE_MARGIN, TREE_MARGIN,
                             TREE_MARGIN, "in"))

if (TIPLAB) p <- p + geom_tiplab(size = TIPLAB_SZ)

## ---- rings ------------------------------------------------------------------
## Ring 1 (innermost)  named vs provisional species
## Ring 2              clade, including singletons
## Rings 3-6           terrestrial, freshwater, marine, brackish
##
## Branches are uncoloured and no clade shading is drawn: every category in the
## figure is read from the rings, so nothing is encoded twice and no ancestral
## branch is implied to have a species status or a habitat.
tip_clade <- rep(NA_character_, n_tip)
for (nm in names(clade_pieces)) {
  tips <- unlist(desc[clade_pieces[[nm]]], use.names = FALSE)
  tips <- tips[tips <= n_tip]
  tip_clade[tips] <- nm
}

pd <- p$data
xmax <- max(pd$x, na.rm = TRUE)
tip_y <- pd$y[match(seq_len(n_tip), pd$node)]

## ---- how many ARCS each clade forms on the ring, in THIS tree ---------------
## The legend used to report length(clade_pieces[[nm]]) -- the number of
## maximal pure subclades. That is computed per tree, but it can exceed what a
## reader can count: two pure subclades that happen to sit side by side in the
## tip order fuse into ONE contiguous coloured arc. The legend now reports
## contiguous arcs in plotting order (the ring is closed, so a clade spanning
## the seam is one arc, not two). Both counts go to the console and the
## fragmentation report so they can be reconciled.
clade_arcs <- local({
  ord <- order(tip_y)
  cl  <- ifelse(is.na(tip_clade[ord]), "__none__", tip_clade[ord])
  r   <- rle(cl)
  v   <- r$values
  if (length(v) > 1 && v[1] != "__none__" && v[1] == v[length(v)])
    v <- v[-1]                              # merge the run across the seam
  tab <- table(v[v != "__none__"])
  setNames(as.integer(tab), names(tab))
})
cmp <- data.frame(clade = names(clade_pieces),
                  highlighted_pieces = as.integer(lengths(clade_pieces)),
                  arcs_on_ring = as.integer(clade_arcs[names(clade_pieces)]),
                  stringsAsFactors = FALSE)
cmp$arcs_on_ring[is.na(cmp$arcs_on_ring)] <- 0L
cat("\n---- clade arcs on the ring (the count the legend reports) ----\n")
print(cmp, row.names = FALSE)
if (any(cmp$arcs_on_ring != cmp$highlighted_pieces))
  cat("  (pieces > arcs: adjacent pure subclades fused into one visible arc)\n")
report$arcs_on_ring <- cmp$arcs_on_ring[match(report$clade, cmp$clade)]
write.csv(report, "clade_fragmentation_report.csv", row.names = FALSE)

## ---- environments: union of route A (hosts) and route B (species name) ------
## Route A resolves each host in the cluster and unions the results, so a
## provisional cluster still gets habitats. Route B looks the cluster's species
## name up in the euler script's output -- but ONLY for the cluster that owns the
## curated database record (tip_owns), because that lookup is keyed by species
## name and would otherwise copy the whole species' habitat union onto every
## same-named cluster (all the Vittaforma corneae pieces, say). Non-owning
## clusters keep route A alone, i.e. the habitats of their own GenBank hosts.
## A tip is "no data" only when both applicable routes come back empty.
sp_env_lookup <- character(0)
if (file.exists(SPECIES_ENV_FILE)) {
  et <- read.delim(SPECIES_ENV_FILE, stringsAsFactors = FALSE, na.strings = "",
                   quote = "\"")
  sc <- intersect(c("entity", "species", "microsporidia_species"), names(et))[1]
  ec <- intersect(c("environments", "environment", "environment_all"), names(et))[1]
  if (!is.na(sc) && !is.na(ec)) {
    sp_env_lookup <- setNames(ifelse(is.na(et[[ec]]), "", et[[ec]]),
                              norm_host(et[[sc]]))
    sp_env_lookup <- sp_env_lookup[!duplicated(names(sp_env_lookup))]
    cat("Route B: ", sum(nzchar(sp_env_lookup)), " species with habitats\n", sep = "")
  }
} else {
  cat("Route B: ", SPECIES_ENV_FILE, " not found\n", sep = "")
}

## --- route A lookup: host name -> habitat (RESTORED) -------------------------
## host_env is consumed in the loop below. The block that builds it was lost in
## an earlier edit, so the script only ran where a stale host_env survived in
## the workspace from a previous session; in a clean session it failed with
## "object 'host_env' not found".
host_env <- character(0)
for (f in HOST_ENV_FILES) {
  if (!file.exists(f)) { cat("NOTE: host env file not found: ", f, "\n", sep = ""); next }
  d <- read.delim(f, stringsAsFactors = FALSE, na.strings = "", quote = "\"")
  ecol <- intersect(c("environment", "host_derived_environment",
                      "environment_all", "environment_primary"), names(d))[1]
  hcol <- intersect(c("host_name", "host"), names(d))[1]
  if (is.na(ecol) || is.na(hcol)) {
    cat("NOTE: no host/environment columns in ", f, "\n", sep = ""); next
  }
  v <- setNames(ifelse(is.na(d[[ecol]]), "", d[[ecol]]), norm_host(d[[hcol]]))
  v <- v[nzchar(names(v))]
  v <- v[!duplicated(names(v))]
  if (!length(host_env)) {
    host_env <- v
  } else if (identical(HOST_ENV_RULE, "union")) {
    k <- union(names(host_env), names(v))
    a <- host_env[k]; b <- v[k]; a[is.na(a)] <- ""; b[is.na(b)] <- ""
    host_env <- setNames(trimws(paste(a, b, sep = "; ")), k)
  } else {
    keep <- if (identical(HOST_ENV_RULE, "prediction_first")) v else host_env
    fill <- if (identical(HOST_ENV_RULE, "prediction_first")) host_env else v
    k <- union(names(keep), names(fill))
    a <- keep[k]; b <- fill[k]; a[is.na(a)] <- ""; b[is.na(b)] <- ""
    host_env <- setNames(ifelse(nzchar(a), a, b), k)
  }
  cat("Loaded ", length(v), " host habitat record(s) from ", f, "\n", sep = "")
}
host_env <- gsub("^;\\s*|;\\s*$", "", host_env)
cat("Route A: ", sum(nzchar(host_env)), " host name(s) with a habitat\n", sep = "")
if (!length(host_env))
  cat("*** no host habitat data loaded -- the environment rings will be blank;",
      " check HOST_ENV_FILES ***\n", sep = "")

env_mat <- matrix(NA, nrow = n_tip, ncol = length(ENV_ORDER),
                  dimnames = list(tree$tip.label, ENV_ORDER))
n_A <- n_B <- n_any <- n_B_skip <- 0L
host_hits <- host_miss <- character(0)

for (i in seq_len(n_tip)) {
  envs <- character(0)
  hs <- tip_hosts[i]
  if (!is.na(hs) && nzchar(hs)) {
    hv <- norm_host(split_hosts(hs))
    found <- hv %in% names(host_env)
    host_hits <- c(host_hits, hv[found]); host_miss <- c(host_miss, hv[!found])
    if (any(found)) {
      a <- unique(trimws(unlist(strsplit(
        paste(host_env[hv[found]], collapse = "; "), ";"))))
      a <- a[nzchar(a)]
      if (length(a)) { envs <- union(envs, a); n_A <- n_A + 1L }
    }
  }
  sp <- tip_label_sp[i]
  if (!is.na(sp) && nzchar(sp) && length(sp_env_lookup)) {
    v <- sp_env_lookup[norm_host(sp)]
    if (!is.na(v) && nzchar(v)) {
      # route B applies only to the cluster that owns the curated record; a
      # same-named non-owning cluster would otherwise inherit the species-level
      # habitat union and every V. corneae piece would look identical
      if (isTRUE(tip_owns[i])) {
        b <- unique(trimws(unlist(strsplit(v, ";")))); b <- b[nzchar(b)]
        if (length(b)) { envs <- union(envs, b); n_B <- n_B + 1L }
      } else {
        n_B_skip <- n_B_skip + 1L
      }
    }
  }
  if (!length(envs)) next
  n_any <- n_any + 1L
  env_mat[i, ] <- ENV_ORDER %in% envs
}

cat("\n---- environment rings ----\n")
cat("tips resolved via hosts (A):", n_A, "\n")
cat("tips resolved via name  (B):", n_B, "\n")
cat("route B suppressed (non-owning cluster shares a species name):",
    n_B_skip, "\n")
cat("tips with >=1 habitat      :", n_any, "of", n_tip, "\n")
cat("host names with no habitat :", length(unique(host_miss)), "\n")
if (length(unique(host_miss)))
  cat("  e.g.", paste(head(unique(host_miss), 4), collapse = "; "), "\n")
for (j in seq_along(ENV_ORDER))
  cat(sprintf("  %-12s %d tips\n", ENV_ORDER[j], sum(env_mat[, j], na.rm = TRUE)))

## --- assemble every ring into ONE data frame so a single fill scale works ----
rings <- list()
r <- 1 + RING_GAP

## ring 1: named vs provisional
rings[[1]] <- data.frame(
  y = tip_y, x = xmax * r,
  fill_key = ifelse(is.na(tip_entity), "__nostatus__", tip_entity),
  ring = "status", stringsAsFactors = FALSE)

## ring 2: clade
r <- r + RING_WIDTH + RING_SPACE
rings[[2]] <- data.frame(
  y = tip_y, x = xmax * r,
  fill_key = tip_clade, ring = "clade", stringsAsFactors = FALSE)

## rings 3-6: environments
for (j in seq_along(ENV_ORDER)) {
  r <- r + RING_WIDTH + RING_SPACE
  v <- env_mat[, j]
  rings[[j + 2]] <- data.frame(
    y = tip_y, x = xmax * r,
    fill_key = ifelse(is.na(v), "__nodata__",
                      ifelse(v, ENV_ORDER[j], "__absent__")),
    ring = ENV_ORDER[j], stringsAsFactors = FALSE)
}
ring_df <- do.call(rbind, rings)
ring_df <- ring_df[!is.na(ring_df$fill_key), ]

fill_values <- c(ENTITY_COLORS,
                 `__nostatus__` = ENTITY_NA_COL,
                 clade_colors[names(clade_pieces)],
                 ENV_COLORS,
                 `__absent__` = ABSENT_COL,
                 `__nodata__` = NODATA_COL)
ring_df$fill_key <- factor(ring_df$fill_key, levels = names(fill_values))

p <- p +
  geom_tile(data = ring_df, aes(x = x, y = y, fill = fill_key),
            width = xmax * RING_WIDTH, height = 1, inherit.aes = FALSE) +
  scale_fill_manual(values = fill_values, na.value = NA, guide = "none")

cat("Rings drawn: status, clade, +", length(ENV_ORDER),
    "environment tracks\n")

## ---- host taxonomy symbols, stacked outside the rings -----------------------
## PHYLUM sets the SHAPE. COLOUR is the host ORDER within Chordata / Arthropoda,
## and the PHYLUM itself for every other phylum -- so every symbol is coloured
## and nothing falls back to grey. Symbols for one tip stack radially outward.
all_ph  <- unlist(lapply(tip_taxa, function(d) d$phylum), use.names = FALSE)
all_ord <- unlist(lapply(tip_taxa, function(d) d$order),  use.names = FALSE)

## which phylum does each order belong to? the legend draws an order using its
## parent phylum's shape
ord_parent <- character(0)
if (length(all_ord)) {
  ok <- !is.na(all_ord) & nzchar(all_ord) & !is.na(all_ph)
  if (any(ok)) {
    sp <- split(all_ph[ok], all_ord[ok])
    ord_parent <- vapply(sp, function(v)
      names(sort(table(v), decreasing = TRUE))[1], character(1))
    ## each pooled bucket maps to the phylum named in its own label
    pooled_ph_names <- sort(unique(all_ph[ok]))
    if (length(pooled_ph_names))
      ord_parent <- c(ord_parent,
        setNames(pooled_ph_names, paste0("other ", pooled_ph_names, " order")))
  }
}

## defaults so the legend can run even when nothing resolves
ph_levels <- character(0); colour_levels <- character(0)
shape_map <- integer(0);   colour_map <- character(0)
ord_group <- character(0)
sym <- NULL

if (length(all_ph)) {
  ## --- counts are per CLUSTER, not per symbol ---
  ph_by_tip  <- lapply(tip_taxa, function(d) unique(d$phylum))
  ord_by_tip <- lapply(tip_taxa, function(d) unique(d$order[nzchar(d$order)]))
  ph_n  <- table(unlist(ph_by_tip,  use.names = FALSE))
  ord_n <- table(unlist(ord_by_tip, use.names = FALSE))

  keep_ph    <- names(ph_n)[ph_n  >= HOSTTAX_MIN_N]
  pooled_ph  <- setdiff(names(ph_n), keep_ph)
  keep_ord   <- names(ord_n)[ord_n >= HOSTTAX_MIN_N]
  pooled_ord <- setdiff(names(ord_n), keep_ord)

  if (length(pooled_ph))
    cat("Phyla pooled as 'Other phylum' (<", HOSTTAX_MIN_N, " clusters): ",
        paste(sprintf("%s(%d)", pooled_ph, ph_n[pooled_ph]), collapse = ", "),
        "\n", sep = "")
  other_lab <- function(phy) paste0("other ", phy, " order")

  ## which phyla actually contribute a pooled order?
  pooled_parent <- character(0)
  if (length(pooled_ord)) {
    okp <- !is.na(all_ord) & all_ord %in% pooled_ord & !is.na(all_ph)
    pooled_parent <- sort(unique(all_ph[okp]))
    cat("Orders pooled (<", HOSTTAX_MIN_N, " clusters): ",
        paste(sprintf("%s(%d)", pooled_ord, ord_n[pooled_ord]), collapse = ", "),
        "\n  -> buckets: ",
        paste(other_lab(pooled_parent), collapse = ", "), "\n", sep = "")
  }

  ## audit table: every phylum and order with its cluster count and whether it
  ## survived the threshold. Without this a taxon that quietly falls into
  ## "Other phylum" looks like a bug rather than a threshold decision.
  tax_audit <- rbind(
    data.frame(rank = "phylum", taxon = names(ph_n),
               n_clusters = as.integer(ph_n),
               shown = ifelse(names(ph_n) %in% keep_ph, "yes",
                              "no - pooled as 'Other phylum'"),
               stringsAsFactors = FALSE),
    if (length(ord_n)) data.frame(rank = "order", taxon = names(ord_n),
               n_clusters = as.integer(ord_n),
               shown = ifelse(names(ord_n) %in% keep_ord, "yes",
                              "no - pooled as 'other order'"),
               stringsAsFactors = FALSE))
  tax_audit <- tax_audit[order(tax_audit$rank, -tax_audit$n_clusters), ]
  write.csv(tax_audit, "host_taxon_counts.csv", row.names = FALSE)
  cat("\n---- host taxa: cluster counts and threshold (>= ", HOSTTAX_MIN_N,
      ") ----\n", sep = "")
  print(tax_audit, row.names = FALSE, right = FALSE)
  near <- tax_audit[tax_audit$shown != "yes" &
                      tax_audit$n_clusters >= max(1, HOSTTAX_MIN_N - 2), ]
  if (nrow(near))
    cat("\nJust below the threshold (lower HOSTTAX_MIN_N to show these): ",
        paste(sprintf("%s(%d)", near$taxon, near$n_clusters), collapse = ", "),
        "\n", sep = "")

  pool_ph  <- function(v) ifelse(v %in% keep_ph, v, "Other phylum")
  ## Pooled orders keep their parent phylum, so a rare chordate order and a
  ## rare arthropod order do not collapse into one anonymous bucket. With a
  ## single "other order" the reader cannot tell which phylum it belongs to,
  ## and only one of the two ever appears if just one phylum has a tail.
  pool_ord <- function(phy, ord) ifelse(!nzchar(ord), "",
                              ifelse(ord %in% keep_ord, ord, other_lab(phy)))

  ## --- phyla -> shapes ---
  ph_levels <- c(sort(keep_ph), if (length(pooled_ph)) "Other phylum")
  shape_map <- setNames(rep(NA_integer_, length(ph_levels)), ph_levels)
  res <- PHYLUM_SHAPES[intersect(names(PHYLUM_SHAPES), ph_levels)]
  shape_map[names(res)] <- as.integer(res)
  if ("Other phylum" %in% ph_levels && is.na(shape_map[["Other phylum"]]))
    shape_map[["Other phylum"]] <- 1L
  need <- ph_levels[is.na(shape_map)]
  if (length(need)) {
    free <- setdiff(PHYLUM_POOL, as.integer(shape_map[!is.na(shape_map)]))
    if (length(free) < length(need)) {
      cat("*** only ", length(free), " distinct symbols for ", length(need),
          " phyla; some reused ***\n", sep = "")
      free <- rep(if (length(free)) free else PHYLUM_POOL, length.out = length(need))
    }
    shape_map[need] <- as.integer(free[seq_along(need)])
  }
  stopifnot(!any(is.na(shape_map)))

  ## --- colours: order where it exists, else phylum ---
  ## Group the orders BY PARENT PHYLUM, then alphabetically within each, and
  ## put each phylum's "other <Phylum> order" bucket at the end of its own
  ## group. A flat alphabetical list interleaves chordate and arthropod orders,
  ## which makes the legend read as one undifferentiated block.
  ord_of_phylum <- function(phy) {
    o <- keep_ord[vapply(keep_ord, function(x) {
      pp <- ord_parent[match(x, names(ord_parent))]
      !is.na(pp) && pp == phy
    }, logical(1))]
    sort(o)
  }
  split_present <- SPLIT_TO_ORDER[SPLIT_TO_ORDER %in% ph_levels]
  ord_levels <- character(0)
  ord_group  <- character(0)   # order level -> phylum, for legend sub-headers
  for (phy in split_present) {
    grp <- ord_of_phylum(phy)
    if (phy %in% pooled_parent) grp <- c(grp, other_lab(phy))
    if (!length(grp)) next
    ord_levels <- c(ord_levels, grp)
    ord_group  <- c(ord_group, setNames(rep(phy, length(grp)), grp))
  }
  ## any order whose parent is not a split phylum (should not happen, but a
  ## surprise in the taxonomy table must not silently vanish from the legend)
  leftover <- setdiff(c(keep_ord,
                        if (length(pooled_parent)) other_lab(pooled_parent)),
                      ord_levels)
  if (length(leftover)) {
    cat("NOTE: order(s) with no split-phylum parent, appended last: ",
        paste(leftover, collapse = ", "), "\n", sep = "")
    ord_levels <- c(ord_levels, sort(leftover))
    ord_group <- c(ord_group, setNames(rep("", length(leftover)), leftover))
  }

  ph_colour_levels <- sort(setdiff(ph_levels, SPLIT_TO_ORDER))
  other_levels <- character(0)   # already folded into ord_levels, per phylum
  colour_levels <- unique(c(ord_levels, ph_colour_levels))
  if (!length(colour_levels)) colour_levels <- "unassigned"

  ## Orders are assigned from the front of the palette and phyla from the back,
  ## so the two kinds of entry stay far apart in colour space while every
  ## individual colour remains distinct.
  n_o <- length(ord_levels)
  n_p <- length(colour_levels) - n_o
  pal <- rep(DISTINCT_PAL, length.out = max(length(colour_levels),
                                            length(DISTINCT_PAL)))
  if (length(colour_levels) > length(DISTINCT_PAL))
    cat("*** ", length(colour_levels), " colour levels but only ",
        length(DISTINCT_PAL), " distinct colours; some repeat. Raise ",
        "HOSTTAX_MIN_N to pool more. ***\n", sep = "")
  cols <- c(if (n_o > 0) pal[seq_len(n_o)],
            if (n_p > 0) rev(pal)[seq_len(n_p)])
  colour_map <- setNames(cols[seq_along(colour_levels)], colour_levels)
  if (!("unassigned" %in% names(colour_map)))
    colour_map <- c(colour_map, unassigned = NOORDER_COL)
  colour_levels <- names(colour_map)
  stopifnot(!any(is.na(colour_map)))

  colour_of <- function(phy, ord) {
    out <- ifelse(nzchar(ord), ord, phy)
    ifelse(out %in% colour_levels, out, "unassigned")
  }

  cat("Symbol colours: ", length(ord_levels), " order(s)",
      if (length(pooled_ord)) " + other order" else "", ", ",
      length(ph_colour_levels), " phylum-level colour(s)\n", sep = "")

  sym <- do.call(rbind, lapply(seq_len(n_tip), function(i) {
    d <- tip_taxa[[i]]
    if (!nrow(d)) return(NULL)
    d$phylum <- pool_ph(d$phylum)
    d$order  <- pool_ord(d$phylum, d$order)
    d <- d[!duplicated(paste(d$phylum, d$order)), , drop = FALSE]
    d <- head(d[order(d$phylum, d$order), , drop = FALSE], HOSTSYM_MAX)
    data.frame(y = tip_y[i],
               x = xmax * (r + HOSTSYM_START + (seq_len(nrow(d)) - 1) * HOSTSYM_GAP),
               phylum = factor(d$phylum, levels = ph_levels),
               order  = factor(colour_of(d$phylum, d$order), levels = colour_levels),
               stringsAsFactors = FALSE)
  }))

  if (!is.null(sym) && nrow(sym)) {
    cat("\n-- symbol colour assignment (this is what is drawn) --\n")
    print(table(colour = as.character(sym$order), useNA = "ifany"))
    if (any(is.na(sym$order)) || any(is.na(sym$phylum)))
      cat("*** NA shape/colour present -- those render grey50 ***\n")
    p <- p +
      geom_point(data = sym, aes(x = x, y = y, shape = phylum, colour = order),
                 size = HOSTSYM_SIZE, inherit.aes = FALSE) +
      scale_shape_manual(values = shape_map, drop = FALSE, guide = "none") +
      scale_colour_manual(values = colour_map, drop = FALSE, guide = "none")
    cat("Host symbols: ", nrow(sym), " drawn across ", length(unique(sym$y)),
        " tips; ", length(ph_levels), " phyla, ", length(ord_levels),
        " orders\n", sep = "")
  }
} else {
  cat("Host symbols: none (no host taxonomy resolved)\n")
}

## ---- legend ------------------------------------------------------------------
## built from the clade list, so nothing can silently go missing
add_legend_box <- function() {
  ## no border box: at 7 inches the frame only steals width from the tree

  ## Two-column flow layout. Entries are queued and then laid out, so adding a
  ## category never requires recomputing y positions by hand -- the previous
  ## hardcoded version silently pushed later blocks off the bottom of the box.
  ## Columns are packed tight: the second starts just past the longest label
  ## rather than at the halfway mark, which previously left a wide empty gutter.
  ## The legend now sits in a WIDE, SHORT band under the tree, so the column
  ## geometry is computed from LEGEND_COLS rather than hardcoded for two.
  ## Row height is in npc of a band LEGEND_IN inches tall: 6 pt with 1.55x
  ## leading is 0.129 in, i.e. 0.129/LEGEND_IN of the band.
  col_w   <- 1 / LEGEND_COLS
  COL_X   <- (seq_len(LEGEND_COLS) - 1) * col_w + 0.006
  SW_DX   <- col_w * 0.055          # swatch offset from the column edge
  TXT_DX  <- col_w * 0.150          # text offset
  ROW_H   <- (LEGEND_FS / 72 * 1.55) / LEGEND_IN
  TOP     <- 1 - ROW_H * 0.55
  BOTTOM  <- ROW_H * 0.35

  st <- new.env()
  st$col <- 1; st$y <- TOP

  st$broke <- FALSE
  nl <- function(n = 1) {
    st$y <- st$y - n * ROW_H
    if (st$y < BOTTOM && st$col < LEGEND_COLS) {
      st$col <- st$col + 1; st$y <- TOP; st$broke <- TRUE
    }
  }
  ## `rows` is how many lines the whole block needs. Breaking to the next
  ## column only when the HEADER does not fit leaves the header stranded at the
  ## bottom of one column with its entries at the top of the next.
  hdr <- function(txt, rows = 2) {
    if (st$y - (rows + 1) * ROW_H < BOTTOM && st$col < LEGEND_COLS) {
      st$col <- st$col + 1; st$y <- TOP
    }
    grid.text(txt, x = COL_X[st$col], y = st$y, just = "left",
              gp = gpar(fontface = "bold", fontsize = LEGEND_FS + 1.5))
    nl()
  }
  note <- function(txt) {
    grid.text(txt, x = COL_X[st$col], y = st$y, just = "left",
              gp = gpar(fontsize = LEGEND_FS - 1, col = "grey35"))
    nl(0.75)
  }
  swatch <- function(txt, fill) {
    grid.rect(x = COL_X[st$col] + SW_DX, y = st$y, width = col_w * 0.075,
              height = ROW_H * 0.62, gp = gpar(fill = fill, col = NA))
    grid.text(txt, x = COL_X[st$col] + TXT_DX, y = st$y, just = "left",
              gp = gpar(fontsize = LEGEND_FS))
    nl()
  }
  glyph <- function(txt, pch, col = "grey15") {
    grid.points(x = unit(COL_X[st$col] + SW_DX, "npc"), y = unit(st$y, "npc"),
                pch = pch, size = unit(0.42, "char"), gp = gpar(col = col))
    grid.text(txt, x = COL_X[st$col] + TXT_DX, y = st$y, just = "left",
              gp = gpar(fontsize = LEGEND_FS))
    nl()
  }

  ## --- species status ---
  hdr("Species status",
      length(ENTITY_LABELS) + 1 + as.integer(any(is.na(tip_entity))))
  for (lv in ENTITY_LABELS) swatch(lv, ENTITY_COLORS[[lv]])
  ## only show "unclassified" if a tip actually is -- an entry for a category
  ## with no members just adds noise
  if (any(is.na(tip_entity))) swatch("unclassified", ENTITY_NA_COL)
  nl(0.4)

  ## --- clades ---
  hdr("Clade", length(clade_pieces))
  for (nm in names(clade_pieces)) {
    ## computed for this tree above: contiguous arcs the reader can actually
    ## count, not the (possibly larger) number of pure subclades
    np <- if (nm %in% names(clade_arcs)) clade_arcs[[nm]] else 0L
    swatch(if (np > 1) sprintf("%s (%d groups)", nm, np) else nm,
           clade_colors[[nm]])
  }
  nl(0.4)

  ## --- environments ---
  hdr("Environment rings", length(ENV_ORDER) + 2)
  for (e in ENV_ORDER) swatch(cap1(e), ENV_COLORS[[e]])
  swatch("Not in this habitat", ABSENT_COL)
  swatch("No habitat record", NODATA_COL)
  nl(0.4)

  ## --- host phylum (shape) ---
  if (length(ph_levels)) {
    hdr("Host phylum (symbol)", length(ph_levels) + 2)
    note(paste0(">= ", HOSTTAX_MIN_N, " clusters; rest pooled"))
    for (ph in ph_levels) {
      pchv <- shape_map[match(ph, names(shape_map))]
      ## a phylum that is NOT split to order carries its own colour, so draw it
      ## here in that colour and leave it out of the colour block below --
      ## otherwise Nematoda, Mollusca etc. appear twice in the legend
      cv <- colour_map[match(ph, names(colour_map))]
      glyph(ph, if (is.na(pchv)) 1L else as.integer(pchv),
            if (is.na(cv)) "grey15" else cv)
    }
    nl(0.4)
  }

  ## --- host order (colour) ---
  if (length(setdiff(colour_levels, c(ph_levels, "unassigned")))) {
    ord_show <- setdiff(colour_levels, c(ph_levels, "unassigned"))
    hdr("Host order (colour)", length(ord_show) + 3)
    note(paste0(paste(SPLIT_TO_ORDER, collapse = " & "), " only;"))
    note(paste0("other phyla take the colour"))
    note("shown above")
    ## only true ORDERS here: phylum-level colours are already shown in the
    ## phylum block, drawn in their own colour
    last_grp <- NA_character_
    for (o in ord_show) {
      ## unname(): ord_group is a NAMED vector, and identical() compares names
      ## too, so without this the sub-header repeats on every row
      grp <- unname(ord_group[match(o, names(ord_group))])
      if (!is.na(grp) && nzchar(grp) && !identical(grp, last_grp)) {
        ## start each group in a fresh column so it is never split
        if (ORDER_GROUP_NEW_COL && !identical(last_grp, NA_character_) &&
            st$col < LEGEND_COLS) {
          st$col <- st$col + 1; st$y <- TOP
        }
        note(paste0(grp, ":"))
        last_grp <- grp
      } else if (!is.na(grp) && nzchar(grp) && st$broke) {
        note(paste0(grp, " (cont.):"))
      }
      st$broke <- FALSE
      cv <- colour_map[match(o, names(colour_map))]
      ## draw each colour entry in the shape it actually appears as on the tree:
      ## an order takes its parent phylum's shape, a phylum-level colour takes
      ## its own shape
      par_ph <- ord_parent[match(o, names(ord_parent))]
      if (is.na(par_ph) && o %in% names(shape_map)) par_ph <- o
      ## "other <Phylum> order" carries its phylum in the label itself
      if (is.na(par_ph)) {
        m <- regmatches(o, regexec("^other (.+) order$", o))[[1]]
        if (length(m) == 2 && m[2] %in% names(shape_map)) par_ph <- m[2]
      }
      pchv <- if (!is.na(par_ph)) shape_map[match(par_ph, names(shape_map))] else NA
      if (o == "unassigned") pchv <- 16L
      glyph(o, if (is.na(pchv)) 16L else as.integer(pchv),
            if (is.na(cv)) "grey50" else cv)
    }
    nl(0.4)
  }

  ## --- footer ---
  ## The "Reading the figure" block was removed: the ring order and symbol
  ## stacking belong in the figure caption, not in the plot itself.
}

## ---- per-tip table: clade, status, environments, and WHY anything is blank ---
## Covers EVERY tip, not just those inside a clade, so a tip missing from the
## figure can be traced. status_reason explains any tip that is neither named
## nor provisional; the three causes are quite different:
##   "not in cluster table"  the accession did not match ATTR_FILE at all
##   "insufficient"          cluster_to_database.py found neither a binomial
##                           nor a reported host, so it is neither named nor
##                           provisional by that script's definition
##   "<other value>"         an entity_type this script does not recognise
tip_of_clade <- rep(NA_character_, n_tip)
for (nm in names(clade_pieces)) {
  ti <- unlist(desc[clade_pieces[[nm]]], use.names = FALSE)
  ti <- ti[ti <= n_tip]
  tip_of_clade[ti] <- nm
}

raw_entity <- if (!is.null(attrs) && "entity_type" %in% names(attrs))
  attrs$entity_type[attr_idx] else rep(NA_character_, n_tip)

status_reason <- ifelse(
  !is.na(tip_entity), "",
  ifelse(is.na(attr_idx), "not in cluster table",
         ifelse(is.na(raw_entity), "entity_type blank in cluster table",
                paste0("entity_type = '", raw_entity, "'"))))

assignments <- data.frame(
  tip_label        = tree$tip.label,
  clade            = tip_of_clade,
  status           = ifelse(is.na(tip_entity), "UNCLASSIFIED", tip_entity),
  entity_type_raw  = raw_entity,
  status_reason    = status_reason,
  matched_cluster_table = !is.na(attr_idx),
  cluster_species  = tip_label_sp,
  hosts            = tip_hosts,
  stringsAsFactors = FALSE
)
for (j in seq_along(ENV_ORDER))
  assignments[[ENV_ORDER[j]]] <- ifelse(is.na(env_mat[, j]), "no data",
                                        ifelse(env_mat[, j], "yes", "no"))
assignments$host_taxa   <- vapply(tip_taxa, paste, character(1), collapse = "; ")
assignments$n_host_taxa <- lengths(tip_taxa)

write.csv(assignments, "clade_assignments.csv", row.names = FALSE)
cat("\n", sum(!is.na(tip_of_clade)), "of", n_tip, "tips assigned to",
    length(clade_pieces), "clades\n")

## ---- why are any tips unclassified? -----------------------------------------
uncl <- assignments[assignments$status == "UNCLASSIFIED", ]
if (nrow(uncl)) {
  cat("\n*** ", nrow(uncl), " tip(s) are neither named nor provisional ***\n",
      sep = "")
  print(as.data.frame(table(reason = uncl$status_reason)), row.names = FALSE)
  cat("\nfirst few:\n")
  print(head(uncl[, c("tip_label", "status_reason", "cluster_species")], 8),
        row.names = FALSE)
  write.csv(uncl, "unclassified_tips.csv", row.names = FALSE)
  cat("full list -> unclassified_tips.csv\n")
  if (any(uncl$status_reason == "not in cluster table")) {
    cat("\n'not in cluster table' means the tip label did not match",
        "centroid_accession.\n  tree e.g.: ",
        paste(head(uncl$tip_label, 3), collapse = ", "), "\n", sep = "")
    if (!is.null(attrs) && !is.na(acc_col))
      cat("  table e.g.: ",
          paste(head(attrs[[acc_col]], 3), collapse = ", "), "\n", sep = "")
  }
} else {
  cat("all tips classified as named or provisional\n")
}

## warn about tips claimed by more than one clade
in_two <- vapply(seq_len(n_tip), function(i)
  sum(vapply(clade_pieces, function(nds)
    i %in% unlist(desc[nds], use.names = FALSE), logical(1))), integer(1))
if (any(in_two > 1))
  cat("*** WARNING:", sum(in_two > 1),
      "tip(s) fall inside more than one clade definition ***\n")


## ---- write figures -----------------------------------------------------------
## Each run gets its own timestamped filenames. Overwriting one fixed name makes
## it impossible to tell whether an open PDF is the current run or the previous
## one, which is a real hazard when iterating on a figure.
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
stamp    <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_pdf  <- file.path(OUT_DIR, sprintf("%s_%s.pdf", OUT_STEM, stamp))
out_png  <- file.path(OUT_DIR, sprintf("%s_%s.png", OUT_STEM, stamp))
out_lab  <- file.path(OUT_DIR, sprintf("%s_%s_TIPLABELS.pdf", OUT_STEM, stamp))

draw_figure <- function() {
  grid.newpage()
  ## tree on top (full width, square), legend band beneath
  pushViewport(viewport(layout = grid.layout(2, 1,
                heights = unit(c(FIG_H - LEGEND_IN, LEGEND_IN), "in"))))
  print(p, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
  add_legend_box()
  popViewport()
  popViewport()
}

## PDF -- cairo_pdf where available, so the point symbols and any non-ASCII
## taxon names render identically to the PNG
if (capabilities("cairo")) {
  cairo_pdf(out_pdf, width = FIG_W, height = FIG_H)
} else {
  pdf(out_pdf, width = FIG_W, height = FIG_H)
}
tryCatch(draw_figure(), finally = dev.off())

if (OUT_PNG) {
  ok <- TRUE
  tryCatch({
    if (capabilities("cairo")) {
      png(out_png, width = FIG_W, height = FIG_H, units = "in", res = PNG_DPI,
          type = "cairo")
    } else {
      png(out_png, width = FIG_W, height = FIG_H, units = "in", res = PNG_DPI)
    }
    draw_figure()
    dev.off()
  }, error = function(e) {
    ok <<- FALSE
    if (length(dev.list())) dev.off()
    cat("PNG not written:", conditionMessage(e), "\n")
  })
  if (ok && !file.exists(out_png)) cat("PNG not written (device produced no file)\n")
}

## ---- second figure: tip names, clade ring only -------------------------------
## Rebuilt from scratch rather than by adding labels to `p`, because `p` already
## carries the fill scale for six rings and the shape/colour scales for the host
## symbols. Starting clean keeps the clade ring the only annotation.
if (OUT_LABELS) {
  p_lab <- ggtree(tree, layout = "circular",
                  branch.length = if (CLADOGRAM) "none" else "branch.length",
                  size = 0.2, colour = BRANCH_COL) +
    theme_tree()

  pdl   <- p_lab$data
  xmaxl <- max(pdl$x, na.rm = TRUE)
  tip_yl <- pdl$y[match(seq_len(n_tip), pdl$node)]

  ## clade ring, immediately outside the tips and inside the labels
  ring_lab <- data.frame(y = tip_yl, x = xmaxl * (1 + RING_GAP),
                         clade = tip_clade, stringsAsFactors = FALSE)
  ring_lab <- ring_lab[!is.na(ring_lab$clade), ]
  if (nrow(ring_lab)) {
    ring_lab$clade <- factor(ring_lab$clade, levels = names(clade_pieces))
    p_lab <- p_lab +
      geom_tile(data = ring_lab, aes(x = x, y = y, fill = clade),
                width = xmaxl * LABEL_RING_W, height = 1, inherit.aes = FALSE) +
      scale_fill_manual(values = clade_colors[names(clade_pieces)],
                        guide = "none")
  }

  ## Tip labels: "Species name (accession)" rather than the bare accession.
  ## The species name is tip_label_sp (cluster_label in
  ## 03_CLUSTER_HOST_ATTRIBUTES.tsv). A tip with no matching row -- an outgroup,
  ## say -- keeps just its accession rather than reading "NA (ACC123.1)".
  lab_df <- data.frame(
    label = tree$tip.label,
    tiplab = ifelse(is.na(tip_label_sp) | !nzchar(tip_label_sp),
                    tree$tip.label,
                    paste0(tip_label_sp, " (", tree$tip.label, ")")),
    stringsAsFactors = FALSE)
  cat("Tip labels: ", sum(lab_df$tiplab != lab_df$label), " of ", n_tip,
      " carry a species name; the rest show the accession alone\n", sep = "")

  p_lab <- p_lab %<+% lab_df +
    geom_tiplab(aes(label = tiplab), size = LABEL_SIZE,
                offset = xmaxl * (RING_GAP + LABEL_RING_W + 0.01)) +
    # Legend drawn manually in a band below (see draw_lab), exactly like the
    # main figure: CoordPolar caps the circle at ~80% of the panel, so a
    # ggplot bottom legend can never hug the tips. LABEL_TREE_MARGIN plays
    # the role TREE_MARGIN plays for the main figure.
    theme(legend.position = "none",
          plot.margin = margin(LABEL_TREE_MARGIN, LABEL_TREE_MARGIN,
                               LABEL_TREE_MARGIN, LABEL_TREE_MARGIN,
                               unit = "in")) +
    xlim(NA, xmaxl * LABEL_XLIM)

  ## Tree in the top band, legend drawn with grid in a short band below --
  ## the same architecture as the main figure, so the legend's distance from
  ## the circle is set here directly instead of fought out of CoordPolar.
  draw_lab <- function() {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(
      2, 1, heights = unit(c(LABEL_FIG_H - LABEL_LEGEND_IN, LABEL_LEGEND_IN),
                           "in"))))
    print(p_lab, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
    pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
    nms  <- names(clade_pieces)
    ncol <- ceiling(length(nms) / 2)
    colw <- unit(LABEL_LEGEND_COLW, "in")
    x0   <- unit(0.5, "npc") - 0.5 * ncol * colw
    for (i in seq_along(nms)) {
      row <- ceiling(i / ncol)                # row-major: first half on top
      col <- i - (row - 1) * ncol
      x <- x0 + (col - 0.5) * colw
      y <- unit(1, "npc") - unit(row * 0.34 - 0.12, "in")
      grid.rect(x = x - unit(0.45, "in"), y = y,
                width = unit(0.2, "in"), height = unit(0.2, "in"),
                gp = gpar(fill = clade_colors[[nms[i]]], col = NA))
      grid.text(nms[i], x = x - unit(0.30, "in"), y = y, just = "left",
                gp = gpar(fontsize = LABEL_LEGEND_FS))
    }
    popViewport(2)
  }
  if (capabilities("cairo")) {
    cairo_pdf(out_lab, width = LABEL_FIG_W, height = LABEL_FIG_H)
  } else {
    pdf(out_lab, width = LABEL_FIG_W, height = LABEL_FIG_H)
  }
  tryCatch(draw_lab(), finally = dev.off())
  cat("Tip-label figure: ", n_tip, " labels at size ", LABEL_SIZE,
      "; clade ring only\n", sep = "")
}

## also draw to the interactive device, if there is one
if (interactive()) draw_figure()

cat("\nWrote:\n")
cat("  ", out_pdf, "\n", sep = "")
if (OUT_PNG && file.exists(out_png))
  cat("  ", out_png, " (", PNG_DPI, " dpi)\n", sep = "")
if (OUT_LABELS && file.exists(out_lab))
  cat("  ", out_lab, " (tip names, clade ring only)\n", sep = "")
cat("   clade_assignments.csv\n")
if (file.exists("host_taxon_counts.csv"))
  cat("   host_taxon_counts.csv (why each phylum/order is shown or pooled)\n")
cat("   clade_fragmentation_report.csv\n")
if (exists("uncl") && nrow(uncl)) cat("   unclassified_tips.csv\n")

## keep a record of which run produced which file
log_line <- data.frame(
  timestamp = stamp, pdf = basename(out_pdf),
  png = if (OUT_PNG && file.exists(out_png)) basename(out_png) else "",
  labels_pdf = if (OUT_LABELS && file.exists(out_lab)) basename(out_lab) else "",
  tree = TREE_FILE, outgroup = OUTGROUP, n_tips = Ntip(tree),
  n_clades = length(clade_pieces),
  stringsAsFactors = FALSE)
log_path <- file.path(OUT_DIR, "figure_run_log.tsv")
write.table(log_line, log_path, sep = "\t", row.names = FALSE,
            col.names = !file.exists(log_path), append = file.exists(log_path),
            quote = FALSE)
cat("   figure_run_log.tsv (appended)\n")
