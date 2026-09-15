#!/usr/bin/env Rscript
###############################################################################
# Figure -- geographic distribution of microsporidia records
#
# Choropleth of the number of species recorded per country, parsed out of the
# free-text Locality column, split by named vs provisional species.
#
# MAP DATA LICENSING
#   Basemap: Natural Earth (naturalearthdata.com), 1:110m Admin 0 countries,
#   supplied by the rnaturalearthdata package. Natural Earth is in the PUBLIC
#   DOMAIN (CC0): "no permission is needed to use Natural Earth. Crediting the
#   authors is unnecessary." Nothing here needs a third-party attribution or a
#   figure permission for publication. A credit line is optional courtesy.
#   The fallback basemap (maps::map_data) is also Natural Earth derived.
#
# WHY A DICTIONARY, NOT A SPLIT
#   Locality strings are free text with the country in any position:
#     "U.S.A. (California)", "(Vancouver Island) Canada", "Russia: Lake Baikal",
#     "Lake Biwa, Shiga Prefecture, Japan", "(TN) USA"
#   Splitting on commas therefore fails constantly. Instead every known country
#   name, historical name, misspelling and subnational unit is searched for
#   anywhere in the string, longest alias first, and each match is blanked out
#   before the next is tried. That is what stops "Nigeria" matching "Niger",
#   "Indiana" matching "India", and "Indian Ocean" matching either.
#
# HISTORICAL NAMES are mapped to a modern country and every such row is listed
# in a QC file, because some of them are genuinely lossy (see AMBIGUOUS below).
#
# Usage:  Rscript fig_locality_map.R [path/to/database.xlsx]
###############################################################################

## ---------------------------------------------------------------- config ----
CFG <- list(
  db_path   = "Microsporidia_Characteristics_Database_merge_pro4.xlsx",
  sheet     = "Actively Updated Masterlist",
  name_col  = "Species Name",
  loc_col   = "Locality",

  # "species"  count each species once per country (a species in 3 countries
  #            adds 1 to each) -- the usual reading of "species reported here"
  # "records"  count database rows instead (identical unless rows duplicate)
  count_unit = "species",

  # "all" | "named" | "provisional"
  subset     = "all",

  projection = "robinson",     # "robinson" | "equirectangular"
  ne_scale   = 50,             # Natural Earth scale: 50 (1:50m) or 110 (1:110m)
                               # 110 drops small states (Malta, Grenada, Singapore)
  bin_breaks = c(0, 1, 2, 5, 10, 25, 50, 100, Inf),
  palette    = "viridis",      # "viridis" | "magma" | "blues"
  label_top  = 0,              # label the N highest countries on the map
  # open-water records, drawn as labelled circles on the map
  #   "basin" roll marginal seas up into Atlantic/Pacific/Indian/Arctic/Southern
  #   "sea"   one circle per named water body (Mediterranean, Red Sea, ...)
  #   "none"  countries only; open-water records still reported in the console
  sea_display = "basin",
  # which records the circles count
  #   "open_water_only" only entries whose Locality resolved to NO country, so
  #                     the circles and the choropleth never count the same
  #                     record twice (default)
  #   "all_mentions"    every entry naming that water body, including coastal
  #                     records that also resolved to a country
  sea_scope   = "open_water_only",
  sea_radius  = 6.5,           # circle size (mm); grows for 3-digit counts
  show_antarctica = FALSE,     # grey Antarctica, for context under the
                               # Southern Ocean circle

  outdir     = "map_output",
  fig_width  = 10,
  fig_height = 5.4
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) CFG$db_path <- args[1]

## -------------------------------------------------------------- packages ----
for (p in c("readxl", "ggplot2")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({ library(readxl); library(ggplot2) })
has_sf   <- requireNamespace("sf", quietly = TRUE) &&
            requireNamespace("rnaturalearthdata", quietly = TRUE)
has_stringi <- requireNamespace("stringi", quietly = TRUE)
if (has_sf) suppressPackageStartupMessages(library(sf))

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

show_tbl <- function(df, n = 25, width = 60, file = NULL) {
  if (nrow(df) == 0) return(invisible(NULL))
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
  if (nrow(df) == 0) return(invisible(NULL))
  write.table(df, file.path(CFG$outdir, file), sep = "\t",
              row.names = FALSE, quote = TRUE, fileEncoding = "UTF-8")
}

## ------------------------------------------------------- species classifier --
# classify_species() comes from classify_species.R -- the ONE canonical
# named/provisional implementation, shared by every figure script, the
# assembler and the flowchart. Do not define a local copy here.
source("classify_species.R")

## ============================================================================
## GAZETTEER
## ============================================================================
# Values are ISO 3166-1 alpha-3, which is what the basemap is joined on.

# --- historical and defunct states -------------------------------------------
# AMBIGUOUS: these lose information and every affected row is written to
# QC_historical_names.tsv so the assignment can be checked case by case.
#   USSR/Soviet Union -> RUS  (could be any of 15 successor states)
#   Czechoslovakia    -> CZE  (could be Slovakia)
#   Yugoslavia        -> SRB  (could be any of 7 successor states)
#   Belgian Congo     -> COD
# WEAK: applied only when nothing else in the string matched, so
# "U.S.S.R. (Moldavia)" resolves to Moldova rather than to Russia AND Moldova.
WEAK <- c("ussr", "soviet union", "sovjet union", "czechoslovakia",
          "czecho slovakia", "yugoslavia")

HISTORICAL <- c(
  "ussr" = "RUS", "soviet union" = "RUS", "sovjet union" = "RUS",
  "czechoslovakia" = "CZE", "czecho slovakia" = "CZE",
  "yugoslavia" = "SRB",
  "east germany" = "DEU", "west germany" = "DEU",
  "german democratic republic" = "DEU",
  "federal republic of germany" = "DEU",
  "belgian congo" = "COD", "belgain congo" = "COD", "zaire" = "COD",
  "burma" = "MMR", "siam" = "THA", "ceylon" = "LKA", "persia" = "IRN",
  "formosa" = "TWN", "rhodesia" = "ZWE", "upper volta" = "BFA",
  "dahomey" = "BEN", "tanganyika" = "TZA", "zanzibar" = "TZA",
  "roumania" = "ROU", "rumania" = "ROU",
  "new hebrides" = "VUT", "dutch east indies" = "IDN"
)

# --- spelling and format variants --------------------------------------------
VARIANTS <- c(
  "usa" = "USA", "united states of america" = "USA", "united states" = "USA",
  "uk" = "GBR", "united kingdom" = "GBR", "great britain" = "GBR",
  "britain" = "GBR",
  "england" = "GBR", "scotland" = "GBR", "wales" = "GBR",
  "northern ireland" = "GBR", "isle of man" = "GBR", "hebrides" = "GBR",
  "thames" = "GBR", "london" = "GBR", "plymouth" = "GBR", "devon" = "GBR",
  "brittas" = "IRL", "eire" = "IRL",
  "russian federation" = "RUS", "siberia" = "RUS", "lake baikal" = "RUS",
  "lake baical" = "RUS", "kamchatka" = "RUS",
  # Soviet-era regional names, matched BEFORE the weak "ussr" alias so the
  # record lands in the right successor state
  "karelia" = "RUS", "kola" = "RUS", "amur" = "RUS", "volga" = "RUS",
  "povolzhye" = "RUS", "ural" = "RUS", "sverdlovsk" = "RUS",
  "krasnodar" = "RUS", "chukotka" = "RUS", "peterhoff" = "RUS",
  "kirov" = "RUS", "novosibirsk" = "RUS", "tomsk" = "RUS",
  "kemerovo" = "RUS", "leningrad" = "RUS", "moscow" = "RUS",
  "st petersberg" = "RUS", "st petersburg" = "RUS", "mary a s s r" = "RUS",
  "moldavia" = "MDA", "vilnius" = "LTU",
  "lvov" = "UKR", "cernovcy" = "UKR", "sevastopol" = "UKR",
  "turkestan" = "UZB", "suchomi" = "GEO", "sukhumi" = "GEO",
  "frankreich" = "FRA",
  "turkiye" = "TUR", "trabzon" = "TUR", "ortahisar" = "TUR", "ordu" = "TUR",
  "venezeula" = "VEN",
  "malaya" = "MYS", "west malaysia" = "MYS", "negri sembilan" = "MYS",
  "rab island" = "HRV",
  "amazon river" = "BRA", "florianopolis" = "BRA",
  "visakhapatnam" = "IND",
  "floridia" = "USA", "ascension" = "SHN",
  "strait of georgia" = "CAN",      # BC, not the country -- must beat "georgia"
  "banyuls" = "FRA", "orsay" = "FRA", "essonne" = "FRA", "concarneau" = "FRA",
  "trebon" = "CZE", "chotebor" = "CZE", "lednice" = "CZE",
  "okhotsk" = "RUS", "vrevo" = "RUS",
  "republic of korea" = "KOR", "south korea" = "KOR",
  "north korea" = "PRK", "korea" = "KOR",
  "ivory coast" = "CIV", "cote divoire" = "CIV",
  "czech republic" = "CZE", "czechia" = "CZE", "bohemia" = "CZE",
  "moravia" = "CZE", "slovakia" = "SVK",
  "holland" = "NLD", "the netherlands" = "NLD",
  "vietnam" = "VNM", "viet nam" = "VNM",
  "hong kong" = "HKG", "macau" = "HKG",
  "reunion" = "FRA", "corsica" = "FRA", "guadeloupe" = "FRA",
  "french polynesia" = "PYF", "new caledonia" = "NCL",
  "sardinia" = "ITA", "sicily" = "ITA",
  "madeira" = "PRT", "azores" = "PRT",
  "canary islands" = "ESP", "balearic" = "ESP",
  "greenland" = "GRL", "faroe" = "DNK",
  "tasmania" = "AUS", "papua new guinea" = "PNG",
  "democratic republic of the congo" = "COD", "dr congo" = "COD",
  "republic of the congo" = "COG",
  "myanmar" = "MMR", "eswatini" = "SWZ", "swaziland" = "SWZ",
  "north macedonia" = "MKD", "macedonia" = "MKD",
  "bosnia" = "BIH", "herzegovina" = "BIH",
  "cape verde" = "CPV", "cabo verde" = "CPV", "ngazidja" = "COM",
  "puerto rico" = "PRI"
)

# Aliases whose meaning flips in a US context. value = ISO used when the string
# does NOT also name the USA.
AMBIGUOUS_US <- c("georgia" = "GEO", "cornwall" = "GBR")

# --- subnational units that appear INSTEAD of the country --------------------
US_STATES <- c(
  "alabama", "alaska", "arizona", "arkansas", "california", "colorado",
  "connecticut", "delaware", "florida", "hawaii", "idaho", "illinois",
  "indiana", "iowa", "kansas", "kentucky", "louisiana", "maine", "maryland",
  "massachusetts", "michigan", "minnesota", "mississippi", "missouri",
  "montana", "nebraska", "nevada", "new hampshire", "new jersey",
  "new mexico", "new york", "north carolina", "north dakota", "ohio",
  "oklahoma", "oregon", "pennsylvania", "rhode island", "south carolina",
  "south dakota", "tennessee", "texas", "utah", "vermont", "virginia",
  "washington", "west virginia", "wisconsin", "wyoming",
  # towns and features that appear without a state
  "gainesville", "biscayne bay", "puget sound", "chesapeake bay",
  "cornwall", "wisconisin",   # Cornwall CT; misspelling of Wisconsin in one row
  "great salt lake", "cape cod", "long island sound"
)
CA_PROVINCES <- c(
  "nova scotia", "ontario", "quebec", "british columbia", "alberta",
  "manitoba", "saskatchewan", "newfoundland", "labrador", "new brunswick",
  "prince edward island", "yukon", "northwest territories", "nunavut",
  "vancouver island", "toronto", "montreal"
)
# Provinces/regions used alone for other countries
OTHER_SUBNATIONAL <- c(
  "jiangsu" = "CHN", "jiangxi" = "CHN", "zhejiang" = "CHN", "guangdong" = "CHN",
  "shandong" = "CHN", "hubei" = "CHN", "henan" = "CHN", "hunan" = "CHN",
  "sichuan" = "CHN", "yunnan" = "CHN", "chongqing" = "CHN", "shenzhen" = "CHN",
  "nanjing" = "CHN", "shanghai" = "CHN", "beijing" = "CHN", "tibet" = "CHN",
  "kerala" = "IND", "assam" = "IND", "andhra pradesh" = "IND",
  "tamil nadu" = "IND", "karnataka" = "IND", "west bengal" = "IND",
  "maharashtra" = "IND", "waltair" = "IND",
  "queensland" = "AUS", "new south wales" = "AUS", "victoria, australia" = "AUS",
  "lower saxony" = "DEU", "bavaria" = "DEU", "saxony" = "DEU",
  "brittany" = "FRA", "normandy" = "FRA", "provence" = "FRA",
  "alsace" = "FRA", "roscoff" = "FRA", "montpellier" = "FRA",
  "sao paulo" = "BRA", "minas gerais" = "BRA", "amazonia" = "BRA",
  "rio de janeiro" = "BRA", "para" = "BRA",
  "buenos aires" = "ARG", "neuquen" = "ARG", "patagonia" = "ARG",
  "hokkaido" = "JPN", "honshu" = "JPN", "kyushu" = "JPN",
  "lake biwa" = "JPN", "ibaraki" = "JPN", "iwate" = "JPN", "shiga" = "JPN",
  "crimea" = "UKR", "kiev" = "UKR", "kyiv" = "UKR",
  "tian shan" = "UZB", "bukharskiy" = "UZB",
  "panama canal" = "PAN"
)

# --- open water: recorded, but not attributable to a country -----------------
SEAS <- c(
  "atlantic" = "Atlantic Ocean", "pacific" = "Pacific Ocean",
  "indian ocean" = "Indian Ocean", "southern ocean" = "Southern Ocean",
  "antarctic ocean" = "Southern Ocean", "arctic ocean" = "Arctic Ocean",
  "mediterranean" = "Mediterranean Sea", "baltic" = "Baltic Sea",
  "north sea" = "North Sea", "black sea" = "Black Sea",
  "red sea" = "Red Sea", "caspian" = "Caspian Sea",
  "bay of bengal" = "Bay of Bengal", "gulf of mexico" = "Gulf of Mexico",
  "sea of japan" = "Sea of Japan", "adriatic" = "Adriatic Sea",
  "arabian gulf" = "Persian Gulf", "persian gulf" = "Persian Gulf",
  "barents sea" = "Barents Sea", "bering sea" = "Bering Sea",
  "weddell sea" = "Southern Ocean", "white sea" = "White Sea",
  "gulf of finland" = "Baltic Sea", "finnish bay" = "Baltic Sea",
  "finish bay" = "Baltic Sea", "sea of azov" = "Black Sea"
)

# Where each water body is drawn, and which ocean basin it rolls up into.
# lon/lat are plotting positions in open water, not centroids.
WATERBODY <- data.frame(
  sea = c("Atlantic Ocean", "Pacific Ocean", "Indian Ocean", "Southern Ocean",
          "Arctic Ocean", "Mediterranean Sea", "Baltic Sea", "North Sea",
          "Black Sea", "White Sea", "Barents Sea", "Red Sea", "Persian Gulf",
          "Bay of Bengal", "Sea of Japan", "Bering Sea", "Gulf of Mexico",
          "Adriatic Sea", "Caspian Sea", "Sea of Okhotsk"),
  basin = c("Atlantic Ocean", "Pacific Ocean", "Indian Ocean", "Southern Ocean",
            "Arctic Ocean", "Atlantic Ocean", "Atlantic Ocean", "Atlantic Ocean",
            "Atlantic Ocean", "Arctic Ocean", "Arctic Ocean", "Indian Ocean",
            "Indian Ocean", "Indian Ocean", "Pacific Ocean", "Pacific Ocean",
            "Atlantic Ocean", "Atlantic Ocean", "Caspian Sea", "Pacific Ocean"),
  lon = c(-32, -150, 76, 5, 5, 17, 19, 3, 34, 39, 42, 38, 51, 89, 135, -177,
          -91, 16, 51, 150),
  lat = c(22, 0, -18, -62, 84, 35, 58, 56, 43, 66, 74, 20, 27, 15, 40, 58,
          25, 43, 42, 55),
  stringsAsFactors = FALSE
)
# basin plotting positions (used when sea_display = "basin")
BASIN_POS <- data.frame(
  basin = c("Atlantic Ocean", "Pacific Ocean", "Indian Ocean",
            "Southern Ocean", "Arctic Ocean", "Caspian Sea"),
  lon = c(-32, -150, 76, 5, 5, 51),
  lat = c(22, 0, -18, -58, 79, 42),
  stringsAsFactors = FALSE
)

# --- strings that carry no geography -----------------------------------------
NOISE <- c("europe", "asia", "africa", "america", "north america",
           "south america", "central america", "west africa", "east africa",
           "caribbean", "caribbeans", "west indies", "transcaucasia",
           "worldwide", "cosmopolitan", "widespread", "wide spread",
           "tropical regions", "north temperate zone", "temperate",
           "laboratory", "lab culture", "insectary", "biological control",
           "unknown", "not stated", "not specified", "not reported",
           "no data", "aquarium", "commercial", "imported", "obtained",
           "various", "several countries", "multiple", "many localities")

# filler that appears in almost every locality string; removed before deciding
# whether leftover text is worth reporting as a possible missed country
GEO_FILLER <- c("lake", "lakes", "river", "rivers", "pool", "pond", "stream",
                "bay", "gulf", "sea", "coast", "coastal", "island", "islands",
                "province", "region", "district", "county", "city", "town",
                "village", "near", "vicinity", "valley", "mountain", "mnt",
                "range", "north", "south", "east", "west", "northern",
                "southern", "eastern", "western", "middle", "central",
                "upper", "lower", "state", "prefecture", "territory",
                "strait", "straits", "sound", "inlet", "estuary", "canal",
                "zone", "area", "basin", "reservoir", "farm", "forest",
                "marsh", "swamp", "beach", "harbour", "harbor", "port",
                "and", "the", "of", "in", "at", "from", "etc", "presumably",
                "study", "isolate", "locality", "stated", "abstract",
                "accessible", "precise", "obtained", "collected", "samples")

## ---------------------------------------------------- build the gazetteer ---
build_gazetteer <- function(world_names) {
  # every country name the basemap knows, plus the tables above
  aliases <- world_names                       # named vector: alias -> ISO3
  add <- function(a, iso) {
    a <- tolower(deaccent(a))
    keep <- !(a %in% names(aliases))
    aliases[a[keep]] <<- iso[keep]
  }
  add(names(HISTORICAL), unname(HISTORICAL))
  add(names(VARIANTS), unname(VARIANTS))
  add(US_STATES, rep("USA", length(US_STATES)))
  add(CA_PROVINCES, rep("CAN", length(CA_PROVINCES)))
  add(names(OTHER_SUBNATIONAL), unname(OTHER_SUBNATIONAL))
  # Alias keys go through the same normalisation as the locality text, or
  # hyphenated and abbreviated names ("Guinea-Bissau", "St. Kitts and Nevis")
  # can never match once punctuation has been stripped from the text.
  names(aliases) <- trimws(gsub("\\s+", " ",
                                gsub("[^a-z0-9 ]+", " ", names(aliases))))
  aliases <- aliases[nzchar(names(aliases)) & !duplicated(names(aliases))]
  # longest first, so "Papua New Guinea" wins over "Guinea" and
  # "Indian Ocean" is consumed before "India" can match
  aliases[order(-nchar(names(aliases)))]
}

# Split one locality string into ISO3 codes, sea names and leftovers.
# Dotted abbreviations must be collapsed BEFORE punctuation is stripped: once
# "U.S.A." becomes "u s a" no alias can match it, and gluing single letters
# afterwards mangles "U.S.A. (A permanent pond...)" into "usaa". This was the
# single most common way a locality silently failed to resolve.
ABBREV <- c(
  "u\\.\\s?s\\.\\s?s\\.\\s?r\\.?" = " ussr ",
  "u\\.\\s?s\\.\\s?a\\.?"            = " usa ",
  "u\\.\\s?s\\.(?![a-z])"             = " usa ",
  "u\\.\\s?k\\.?"                     = " uk ",
  "a\\.\\s?s\\.\\s?s\\.\\s?r\\.?" = " russia "
)
expand_abbrev <- function(s) {
  for (i in seq_along(ABBREV))
    s <- gsub(names(ABBREV)[i], ABBREV[i], s, perl = TRUE)
  s
}

# One segment of a locality string. The weak-alias rule is deliberately scoped
# to the SEGMENT, not the whole cell: see parse_locality below.
parse_segment <- function(txt, gaz) {
  s <- paste0(" ", tolower(deaccent(txt)), " ")
  s <- expand_abbrev(s)                      # before punctuation is stripped
  s <- gsub("[^a-z0-9 ]+", " ", s)           # punctuation -> space, keeps words
  s <- gsub("\\s+", " ", s)

  seas <- character(0)
  for (i in order(-nchar(names(SEAS)))) {
    rx <- paste0("\\b", names(SEAS)[i], "\\b")
    if (grepl(rx, s, perl = TRUE)) {
      seas <- c(seas, unname(SEAS)[i])
      s <- gsub(rx, " ", s, perl = TRUE)     # consume, so "indian ocean" != India
    }
  }

  in_us <- grepl("\\b(usa|united states)\\b", s, perl = TRUE)

  iso <- character(0)
  weak_hits <- character(0)
  for (i in seq_along(gaz)) {
    alias <- names(gaz)[i]
    rx <- paste0("\\b", alias, "\\b")
    if (!grepl(rx, s, perl = TRUE)) next
    if (alias %in% WEAK) {                   # defer: USSR, Czechoslovakia, ...
      weak_hits <- c(weak_hits, unname(gaz)[i])
      next
    }
    hit <- unname(gaz)[i]
    if (alias %in% names(AMBIGUOUS_US))
      # a US state matched earlier in the string also counts as US context:
      # "Puget Sound; Strait of Georgia" is not the country Georgia
      hit <- if (in_us || "USA" %in% iso) "USA" else AMBIGUOUS_US[[alias]]
    iso <- c(iso, hit)
    s <- gsub(rx, " ", s, perl = TRUE)
  }
  # a weak alias only counts if the string named no other country
  if (length(iso) == 0 && length(weak_hits)) {
    iso <- weak_hits
    for (a in intersect(names(gaz), WEAK))
      s <- gsub(paste0("\\b", a, "\\b"), " ", s, perl = TRUE)
  }

  for (w in c(NOISE, WEAK)) s <- gsub(paste0("\\b", w, "\\b"), " ", s, perl = TRUE)
  s <- gsub("\\b[a-z]{1,2}\\b", " ", s, perl = TRUE)   # stray initials
  informative <- s
  for (w in GEO_FILLER)
    informative <- gsub(paste0("\\b", w, "\\b"), " ", informative, perl = TRUE)
  informative <- gsub("\\b[0-9]+\\b", " ", informative, perl = TRUE)

  list(iso = unique(iso), seas = unique(seas),
       leftover = trimws(gsub("\\s+", " ", s)),
       informative = trimws(gsub("\\s+", " ", informative)))
}

# A locality cell can list several places separated by ";". Each is parsed on
# its own and the results are unioned, because a species genuinely occurs in
# every country listed and must be counted in all of them.
#
# Scoping the weak-alias rule to the segment is what makes this work. Applied
# to the whole cell, "...; Czechoslovakia (near Chotebor) [Weiser]; Canada..."
# loses Czechia entirely: Canada and France matched, so the weak "czechoslovakia"
# alias was suppressed. Per segment, the Czechoslovakia segment matches nothing
# else, so the weak alias fires and the record is counted in Czechia too --
# while "U.S.S.R. (Moldavia)", a single segment, still resolves to Moldova
# alone rather than to both Moldova and Russia.
parse_locality <- function(txt, gaz) {
  segs <- trimws(strsplit(as.character(txt), ";")[[1]])
  segs <- segs[nzchar(segs)]
  if (length(segs) == 0) segs <- as.character(txt)
  parts <- lapply(segs, parse_segment, gaz = gaz)
  unresolved_seg <- segs[vapply(parts, function(p)
    length(p$iso) == 0 && length(p$seas) == 0 && nzchar(p$informative), logical(1))]
  list(
    iso  = unique(unlist(lapply(parts, `[[`, "iso"))),
    seas = unique(unlist(lapply(parts, `[[`, "seas"))),
    leftover = trimws(paste(vapply(parts, `[[`, character(1), "leftover"),
                            collapse = " ")),
    informative = trimws(paste(vapply(parts, `[[`, character(1), "informative"),
                               collapse = " ")),
    unresolved_segments = unresolved_seg
  )
}

## ------------------------------------------------------------------ load ----
if (!file.exists(CFG$db_path)) stop("database not found: ", CFG$db_path, call. = FALSE)
raw <- read_excel(CFG$db_path, sheet = CFG$sheet, col_types = "text",
                  .name_repair = "minimal")
message(sprintf("read %d rows x %d columns from sheet '%s'",
                nrow(raw), ncol(raw), CFG$sheet))

dat <- data.frame(
  row_xlsx     = seq_len(nrow(raw)) + 1L,
  species_name = trimws(get_col(raw, CFG$name_col)),
  locality     = get_col(raw, CFG$loc_col),
  stringsAsFactors = FALSE
)
dat$locality[is.na(dat$locality)] <- ""

cat("\n================ LOCALITY PARSING REPORT ================\n")
hdr <- which(norm_key(dat$species_name) == norm_key(CFG$name_col))
if (length(hdr)) {
  qc("!! %d repeated HEADER row(s) inside the data (xlsx row %s) -- dropped\n",
     length(hdr), paste(dat$row_xlsx[hdr], collapse = ", "))
  dat <- dat[-hdr, ]
}
dat <- dat[!is.na(dat$species_name) & nzchar(dat$species_name), ]
dat$category <- classify_species(dat$species_name)

if (CFG$subset != "all") {
  n0 <- nrow(dat)
  dat <- dat[dat$category == CFG$subset, ]
  qc("subset = %s: %d of %d entries kept\n", CFG$subset, nrow(dat), n0)
}

## ------------------------------------------------------------- basemap -----
if (has_sf) {
  world <- sf::st_as_sf(if (CFG$ne_scale == 50) rnaturalearthdata::countries50
                        else rnaturalearthdata::countries110)
  # Natural Earth leaves iso_a3 as "-99" for a few states (France, Norway,
  # Kosovo, N. Cyprus, Somaliland); adm0_a3 is populated for all of them.
  world$iso3 <- ifelse(is.na(world$iso_a3) | world$iso_a3 == "-99",
                       world$adm0_a3, world$iso_a3)
  if (!CFG$show_antarctica)
    world <- world[world$iso3 != "ATA", ]     # no records, and it distorts badly
  world_names <- c(
    setNames(world$iso3, tolower(deaccent(world$name))),
    setNames(world$iso3, tolower(deaccent(world$name_long))),
    setNames(world$iso3, tolower(deaccent(world$admin))),
    setNames(world$iso3, tolower(deaccent(world$sovereignt))),
    setNames(world$iso3, tolower(deaccent(world$formal_en)))
  )
  world_names <- world_names[!is.na(names(world_names)) & nzchar(names(world_names))]
  world_names <- world_names[!duplicated(names(world_names))]
} else {
  stop("sf + rnaturalearthdata are required for the public-domain basemap.\n",
       "  install.packages(c('sf','rnaturalearthdata'))\n",
       "  (Ubuntu: apt install r-cran-sf r-cran-rnaturalearthdata)", call. = FALSE)
}

gaz <- build_gazetteer(world_names)
qc("gazetteer: %d aliases -> %d countries\n", length(gaz), length(unique(gaz)))

## ------------------------------------------------------------ parse all -----
parsed <- lapply(dat$locality, parse_locality, gaz = gaz)
dat$n_iso    <- vapply(parsed, function(p) length(p$iso), integer(1))
dat$iso_list <- vapply(parsed, function(p) paste(p$iso, collapse = ";"), character(1))
dat$seas     <- vapply(parsed, function(p) paste(p$seas, collapse = ";"), character(1))
dat$leftover <- vapply(parsed, function(p) p$leftover, character(1))
dat$informative <- vapply(parsed, function(p) p$informative, character(1))

n_blank <- sum(!nzchar(trimws(dat$locality)))
qc("entries                     : %d\n", nrow(dat))
qc("  no Locality recorded      : %d (%.1f%%)\n", n_blank, 100 * n_blank / nrow(dat))
has_loc <- nzchar(trimws(dat$locality))
qc("  Locality recorded         : %d\n", sum(has_loc))
qc("    resolved to >=1 country : %d\n", sum(dat$n_iso > 0))
qc("    open water only         : %d\n", sum(dat$n_iso == 0 & nzchar(dat$seas)))
nogeo <- dat[has_loc & dat$n_iso == 0 & !nzchar(dat$seas) & !nzchar(dat$informative), ]
unresolved <- dat[has_loc & dat$n_iso == 0 & !nzchar(dat$seas) & nzchar(dat$informative), ]
qc("    no usable geography     : %d  (\"Cosmopolitan\", \"Europe\", lab stock)\n",
   nrow(nogeo))
qc("    UNRESOLVED              : %d\n", nrow(unresolved))
if (nrow(nogeo))
  write_qc(nogeo[, c("row_xlsx", "species_name", "locality")],
           "QC_no_usable_geography.tsv")
if (nrow(unresolved)) {
  qc("\nunresolved localities (add these to the gazetteer):\n")
  u <- as.data.frame(table(locality = trimws(unresolved$locality)),
                     stringsAsFactors = FALSE)
  u <- u[order(-u$Freq), c("locality", "Freq")]
  names(u)[2] <- "n_entries"
  show_tbl(u, n = 25, file = "QC_unresolved_localities.tsv")
  write_qc(u, "QC_unresolved_localities.tsv")
}

# A locality can list several places separated by ";". Reporting every row with
# leftover text just lists village names. What actually matters is a SEGMENT
# that resolved to nothing while the row as a whole resolved -- that is where a
# second country is silently lost.
bad_seg <- vapply(parsed, function(p)
  paste(p$unresolved_segments, collapse = " | "), character(1))
sel <- which(dat$n_iso > 0 & nzchar(bad_seg))
seg_report <- if (length(sel))
  data.frame(row_xlsx = dat$row_xlsx[sel], species_name = dat$species_name[sel],
             locality = dat$locality[sel], matched = dat$iso_list[sel],
             unresolved_segment = bad_seg[sel], stringsAsFactors = FALSE) else NULL
if (!is.null(seg_report) && nrow(seg_report)) {
  qc("\n%d multi-part localities have a segment that resolved to nothing\n",
     nrow(seg_report))
  qc("   (a second country may be missing from these rows):\n")
  show_tbl(seg_report[, c("locality", "matched", "unresolved_segment")],
           n = 12, file = "QC_unresolved_segments.tsv")
  write_qc(seg_report, "QC_unresolved_segments.tsv")
}

## ---------------------------------------------------- historical warnings ---
hist_rx <- paste0("(", paste(names(HISTORICAL), collapse = "|"), ")")
hist_rows <- dat[grepl(hist_rx, tolower(deaccent(dat$locality)), perl = TRUE), ]
if (nrow(hist_rows)) {
  qc("\n%d entries use a historical or defunct country name:\n", nrow(hist_rows))
  h <- as.data.frame(table(locality = trimws(hist_rows$locality)),
                     stringsAsFactors = FALSE)
  h <- h[order(-h$Freq), ]; names(h)[2] <- "n_entries"
  show_tbl(h, n = 12, file = "QC_historical_names.tsv")
  qc("   USSR->Russia, Czechoslovakia->Czechia, Yugoslavia->Serbia are LOSSY:\n")
  qc("   the record may belong to any successor state. Check before publishing.\n")
  write_qc(hist_rows[, c("row_xlsx", "species_name", "locality", "iso_list")],
           "QC_historical_names.tsv")
}

## -------------------------------------------------------------- counting ----
long <- do.call(rbind, lapply(which(dat$n_iso > 0), function(i) {
  data.frame(species_name = dat$species_name[i],
             category = dat$category[i],
             iso3 = strsplit(dat$iso_list[i], ";")[[1]],
             stringsAsFactors = FALSE)
}))

if (CFG$count_unit == "species") {
  long <- long[!duplicated(long[, c("species_name", "iso3")]), ]
}
counts <- as.data.frame(table(iso3 = long$iso3), stringsAsFactors = FALSE)
names(counts)[2] <- "n"
by_cat <- as.data.frame(table(iso3 = long$iso3, category = long$category),
                        stringsAsFactors = FALSE)
by_cat <- reshape(by_cat, idvar = "iso3", timevar = "category", direction = "wide")
names(by_cat) <- sub("^Freq\\.", "", names(by_cat))
counts <- merge(counts, by_cat, by = "iso3", all.x = TRUE)
counts$country <- world$name[match(counts$iso3, world$iso3)]
counts <- counts[order(-counts$n), ]

cat("\n================ SPECIES PER COUNTRY ================\n")
qc("%d countries, %d species-country pairs\n\n", nrow(counts), sum(counts$n))
show_tbl(counts[, c("country", "iso3", "n",
                    intersect(c("named", "provisional"), names(counts)))],
         n = 25, file = "counts_by_country.tsv")
write_qc(counts, "counts_by_country.tsv")

# Open water is counted the same way countries are: one species counted once
# per water body. Under the default scope only entries that resolved to NO
# country are counted, so a circle and a filled country never represent the
# same record -- "South Africa; Gulf of Mexico" is a South African record, not
# an Atlantic one. Set sea_scope = "all_mentions" to count coastal records in
# both places instead.
sea_sel <- if (identical(CFG$sea_scope, "all_mentions"))
  which(nzchar(dat$seas)) else which(nzchar(dat$seas) & dat$n_iso == 0)
n_coastal <- sum(nzchar(dat$seas) & dat$n_iso > 0)
sea_long <- do.call(rbind, lapply(sea_sel, function(i) {
  data.frame(species_name = dat$species_name[i], category = dat$category[i],
             sea = strsplit(dat$seas[i], ";")[[1]], stringsAsFactors = FALSE)
}))
sea_counts <- data.frame()
if (!is.null(sea_long) && nrow(sea_long)) {
  if (CFG$count_unit == "species")
    sea_long <- sea_long[!duplicated(sea_long[, c("species_name", "sea")]), ]
  sea_long$basin <- WATERBODY$basin[match(sea_long$sea, WATERBODY$sea)]
  unknown_sea <- unique(sea_long$sea[is.na(sea_long$basin)])
  if (length(unknown_sea)) {
    warning("water body missing from WATERBODY table, not drawn: ",
            paste(unknown_sea, collapse = ", "), call. = FALSE)
    sea_long$basin[is.na(sea_long$basin)] <- sea_long$sea[is.na(sea_long$basin)]
  }

  unit <- if (CFG$sea_display == "basin") "basin" else "sea"
  if (unit == "basin")
    sea_long <- sea_long[!duplicated(sea_long[, c("species_name", "basin")]), ]
  sea_counts <- as.data.frame(table(unit = sea_long[[unit]]),
                              stringsAsFactors = FALSE)
  names(sea_counts) <- c("water_body", "n")
  pos <- if (unit == "basin") BASIN_POS else WATERBODY[, c("sea", "lon", "lat")]
  names(pos)[1] <- "water_body"
  sea_counts <- merge(sea_counts, pos[, c("water_body", "lon", "lat")],
                      by = "water_body", all.x = TRUE)
  sea_counts <- sea_counts[order(-sea_counts$n), ]

  cat("\n---- open-water records ----\n")
  qc("scope: %s\n", CFG$sea_scope)
  if (identical(CFG$sea_scope, "open_water_only") && n_coastal)
    qc("%d further entries name a water body but also resolved to a country;\n  they are counted in the choropleth only.\n", n_coastal)
  qc("%d species-water body pairs across %d water bodies%s\n\n",
     sum(sea_counts$n), nrow(sea_counts),
     if (CFG$sea_display == "none") " (not drawn: sea_display = \"none\")"
     else " (drawn as circles)")
  print(sea_counts[, c("water_body", "n")], row.names = FALSE, right = FALSE)
  detail <- as.data.frame(table(sea = sea_long$sea, basin = sea_long$basin),
                          stringsAsFactors = FALSE)
  detail <- detail[detail$Freq > 0, ]; names(detail)[3] <- "n"
  write_qc(detail[order(-detail$n), ], "counts_by_sea.tsv")
  if (any(is.na(sea_counts$lon)))
    warning("no plotting position for: ",
            paste(sea_counts$water_body[is.na(sea_counts$lon)], collapse = ", "),
            call. = FALSE)
}

## ------------------------------------------------------------------ plot ----
brk <- CFG$bin_breaks
lab <- c()
for (i in seq_len(length(brk) - 1)) {
  lo <- brk[i] + 1; hi <- brk[i + 1]
  lab <- c(lab, if (is.infinite(hi)) sprintf("%d+", lo)
                else if (lo == hi) sprintf("%d", lo)
                else sprintf("%d\u2013%d", lo, hi))
}
counts$bin <- cut(counts$n, breaks = brk, labels = lab, right = TRUE)

world$n <- counts$n[match(world$iso3, counts$iso3)]
world$bin <- cut(world$n, breaks = brk, labels = lab, right = TRUE)

# A binned scale, not a continuous gradient: the counts span 1 to >100 and are
# heavily right-skewed, so on a linear ramp every country except Russia looks
# identical. Bins make "1 record" and "50 records" separable by eye.
pal_fun <- switch(CFG$palette,
  viridis = function(n) grDevices::hcl.colors(n, "viridis", rev = TRUE),
  magma   = function(n) grDevices::hcl.colors(n, "Inferno", rev = TRUE),
  blues   = function(n) grDevices::hcl.colors(n, "Blues 3", rev = TRUE),
  function(n) grDevices::hcl.colors(n, "viridis", rev = TRUE))
fill_cols <- setNames(pal_fun(length(lab)), lab)

crs_out <- if (CFG$projection == "robinson") "+proj=robin +lon_0=0" else "EPSG:4326"
# Polygons that straddle the antimeridian (Russia, Fiji, NZ) smear right across
# the map once reprojected, unless they are cut at 180 degrees first.
world_p <- sf::st_make_valid(world)
world_p <- suppressWarnings(sf::st_wrap_dateline(
  world_p, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180")))
world_p <- sf::st_transform(world_p, crs_out)

p <- ggplot(world_p) +
  geom_sf(aes(fill = bin), colour = "grey35", linewidth = 0.12) +
  scale_fill_manual(values = fill_cols, na.value = "grey93", drop = FALSE,
                    name = sprintf("%s reported",
                                   if (CFG$count_unit == "species") "Species"
                                   else "Records"),
                    guide = guide_legend(nrow = 1, label.position = "bottom",
                                         keywidth = unit(1.5, "lines"),
                                         keyheight = unit(0.55, "lines"))) +
  theme_void(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 9, vjust = 0.9),
    legend.text = element_text(size = 8),
    legend.margin = margin(t = 2),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(4, 6, 4, 6)
  )

## ---- open-water circles ----------------------------------------------------
# Same binned fill as the countries, so one legend reads for both: a circle in
# the same colour as a country means the same number of species.
if (CFG$sea_display != "none" && nrow(sea_counts) &&
    any(!is.na(sea_counts$lon))) {
  sc <- sea_counts[!is.na(sea_counts$lon), ]
  sc$bin <- cut(sc$n, breaks = brk, labels = lab, right = TRUE)
  pts <- sf::st_as_sf(sc, coords = c("lon", "lat"), crs = 4326)
  pts <- sf::st_transform(pts, crs_out)
  xy <- sf::st_coordinates(pts)
  sc$x <- xy[, 1]; sc$y <- xy[, 2]
  # dark fills need light text
  sc$txt <- ifelse(as.integer(sc$bin) > ceiling(length(lab) / 2), "white", "grey10")

  # the circle has to hold its own label, so it grows with the digit count
  sc$radius <- CFG$sea_radius + 1.6 * (nchar(as.character(sc$n)) - 1)
  pts$radius <- sc$radius

  p <- p +
    geom_sf(data = pts, aes(fill = bin), shape = 21, colour = "grey20",
            size = pts$radius, stroke = 0.5, show.legend = FALSE) +
    geom_text(data = sc, aes(x = x, y = y, label = n), inherit.aes = FALSE,
              colour = sc$txt, size = 2.9, fontface = "bold")
  qc("\ncircles drawn for %d water bodies (%s level)\n",
     nrow(sc), CFG$sea_display)
}

if (CFG$label_top > 0) {
  top <- head(counts[order(-counts$n), ], CFG$label_top)
  cen <- suppressWarnings(sf::st_centroid(world_p[world_p$iso3 %in% top$iso3, ]))
  xy <- sf::st_coordinates(cen)
  lbl <- data.frame(x = xy[, 1], y = xy[, 2], n = cen$n)
  p <- p + geom_text(data = lbl, aes(x = x, y = y, label = n),
                     inherit.aes = FALSE, size = 2.6, fontface = "bold",
                     colour = "white")
}

ggsave(file.path(CFG$outdir, "fig_locality_map.pdf"), p,
       width = CFG$fig_width, height = CFG$fig_height,
       device = if (capabilities("cairo")) cairo_pdf else pdf)
ggsave(file.path(CFG$outdir, "fig_locality_map.png"), p,
       width = CFG$fig_width, height = CFG$fig_height, dpi = 600)

write_qc(dat[, c("row_xlsx", "species_name", "category", "locality",
                 "iso_list", "seas", "leftover")], "per_entry_localities.tsv")

qc("\nwritten to %s/\n  fig_locality_map.pdf / .png\n  counts_by_country.tsv\n  counts_by_sea.tsv\n  per_entry_localities.tsv\n  QC_*.tsv\n",
   CFG$outdir)
qc("\nBasemap: Natural Earth 1:110m, public domain (CC0). No attribution or\n")
qc("permission is required to publish this figure.\n")
