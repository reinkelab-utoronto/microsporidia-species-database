#!/usr/bin/env python3
"""
Host taxonomy and environment lookup for the microsporidia attribute database.

Reads the database .xlsx directly, pulls every host name out of the free-text
host columns, resolves each one to a full taxonomic lineage and to a habitat
(marine / brackish / freshwater / terrestrial), and writes tidy TSVs ready for
plotting in R.

SOURCES
  taxonomy     GBIF -> Catalogue of Life -> WoRMS -> NCBI Taxonomy
  environment  WoRMS environment flags -> GBIF habitats -> COL -> optional
               Claude API. Derived from the HOST ORGANISM ONLY: the database's
               own "Host Environment" column is recorded alongside for
               comparison but does not feed the answer, because it is attached
               to a microsporidia row rather than to a host, so the same host
               can carry different values in different rows. Pass
               --use-db-environment to let it override.
  NCBI is new: it is the only source with good coverage of insect and nematode
  names that GBIF backbone matching misses, so it runs last and rescues the
  long tail. It needs an email address (NCBI policy) and honours an API key.

WHAT CHANGED FROM THE PREVIOUS SCRIPT
  * input is the spreadsheet, not a hand-made species list
  * host cells are parsed properly: "Aedes cantans (type host); Aedes flavescens"
    is two hosts, "(= Ochlerotatus crinifer)" is a synonym to retry with, and
    "D. pulex" is expanded from the genus named earlier in the same cell
  * results are cached on disk, so a re-run costs nothing and an interrupted run
    resumes where it stopped
  * output is a tidy long table (one row per microsporidia species x host) plus
    a host lookup table and a QC file of everything that failed to resolve
  * every host gets a coarse host_group (Insecta, Crustacea, Fish, ...) derived
    from the lineage, which is what the figures actually plot

USAGE
  python host_taxonomy_environment.py --db database.xlsx --email you@utoronto.ca
  python host_taxonomy_environment.py --db database.xlsx --email ... --limit 50
  python host_taxonomy_environment.py --db database.xlsx --offline   # cache only
"""

import argparse
import csv
import json
import os
import re
import sys
import time
import unicodedata
import xml.etree.ElementTree as ET
from collections import Counter, OrderedDict

# THE canonical named/provisional rule -- one implementation, shared with
# every figure script via classify_species.R and with the other Python
# scripts via this module. Do not define a local copy; edit
# classify_species.py instead.
from classify_species import classify_microsporidia

try:
    import requests
except ImportError:
    sys.exit("requests required: pip install requests")

# --------------------------------------------------------------------------- #
GBIF_BASE = "https://api.gbif.org/v1"
WORMS_API = "https://www.marinespecies.org/rest"
COL_API = "https://api.catalogueoflife.org"
COL_DATASET = "3LR"
NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
ENVIRONMENTS = ["marine", "brackish", "freshwater", "terrestrial"]

# --------------------------------------------------------------------------- #
# host name parsing
# --------------------------------------------------------------------------- #
# Parenthetical content that is a role or a note, not a name.
ROLE_WORDS = re.compile(
    r"^(?:type host|intermediate host|definitive host|paratype|natural host"
    r"|experimental|primary host|secondary host|alternate host|vector"
    r"|\d+|[a-z]\d?|type \d|host \d|\?|unknown|likely|possibly|probably"
    r"|should be .*|according to .*|misidentified.*|as .*)$", re.I)

# "(= Ochlerotatus crinifer)" / "(also Paranucleospora theridion)"
SYNONYM_RE = re.compile(r"^(?:=|also|syn\.?|synonym of)\s*(.+)$", re.I)

# not a taxon name at all
NON_TAXON = re.compile(
    r"^(?:unknown|not (?:stated|reported|known|specified)|n/?a|none|various"
    r"|several|many|unidentified|undetermined|environmental|host|hosts"
    r"|\?+|-+)$", re.I)

# --- hyperparasite context ---------------------------------------------------
# Microsporidia frequently infect other parasites, and the database records the
# whole chain in one cell:
#     "Sycia inopinata (in the polychaete annelid Audouinia tentaculata)"
#     "Heterophyes heterophyes (parasite in Liza ramada)"
# The MICROSPORIDIAN's host is the OUTER name; the organism inside the
# parenthetical is that host's own host and must not be counted as a
# microsporidian host. It is captured separately so the chain is not lost.
# Subgenus in parentheses: "Lecudina (Ophioidina) sp." is ONE name, not a name
# plus a note. The subgenus is recorded and removed so the lookup sees
# "Lecudina sp." rather than a string no database can match.
SUBGENUS = re.compile(r"^([A-Z][a-z-]+)\s*\(([A-Z][a-z-]+)\)\s*(.*)$")

# Hyperparasite context that is NOT in parentheses:
#   'Lecudina (Ophioidina) sp." in Lumbriconereis latreilli and L. zonata'
# Everything from " in <Genus species>" onwards says where the host lives.
BARE_CONTEXT = re.compile(
    r"[,;\"]?\s*\b(?:hyperparasitic\s+)?(?:in|of|from|within)\b\s+"
    r"(?:the\s+)?(?:[a-z-]+\s+){0,4}([A-Z][a-z]{2,}\s+[a-z][a-z-]{2,}).*$")

# A fragment that does not START with a capitalised genus is a DESCRIPTION, not
# a name: 'A gregarine, "form en comete," in the polychaete annelid Capitella
# capitata'. Mining the only binomial out of it yields the HOST'S host, which is
# how that entry became an annelid. These resolve to nothing by design.
# NB: the lowercase-start test is a separate isupper() check, NOT "^[a-z]" with
# re.I -- under IGNORECASE that class also matches capitals, which classified
# every name as a description.
DESCRIPTIVE = re.compile(r"^(?:a|an|the|unnamed|unidentified|undescribed|"
                         r"gregarine|form|parasitic)\b", re.I)


def tsv_safe(value):
    """
    Strip the three characters that break a TSV: the tab itself, newlines, and
    the double quote. A single stray quote in a host cell makes R's read.delim
    swallow every following line until it finds a partner, which silently drops
    most of the table. Quoting the field instead only helps if every reader
    agrees on the quoting convention, and they do not -- so the characters are
    removed at the source.
    """
    if value is None:
        return ""
    s = str(value).replace("\t", " ").replace("\r", " ").replace("\n", " ")
    return s.replace('"', "'")


def is_description(text):
    return bool(text) and (not text[0].isupper() or bool(DESCRIPTIVE.match(text)))

HYPER_CONTEXT = re.compile(
    r"^(?:hyperparasit\w*|parasit\w*|found|living|encysted|larva\w*|cyst\w*)?"
    r"\s*(?:in|of|within|inside|from)\b", re.I)
BINOMIAL_IN = re.compile(r"\b([A-Z][a-z]{2,})\s+([a-z][a-z-]{2,})\b")
HYPER_WORD = re.compile(r"hyperparasit|parasit(?:e|ic)\b|gregarin|in gut of", re.I)

ABBREV_GENUS = re.compile(r"^([A-Z])\.\s*([a-z][a-z-]+)$")


def deaccent(text):
    s = unicodedata.normalize("NFKD", str(text))
    return "".join(c for c in s if not unicodedata.combining(c))


# --- name repairs before lookup ---------------------------------------------
# Each of these was derived from the unresolved list, and each is reported so a
# repair can be audited rather than trusted.
CONFUSABLE = [
    # a capital I typed for a lower-case l, and a capital O for a zero-like o:
    # "Daphnia Iongispina", "EuboreIlia plebeja"
    (re.compile(r"(?<=[a-z])I(?=[a-z])"), "l"),
    (re.compile(r"(?<=[a-z])O(?=[a-z])"), "o"),
]
# "Ageneotettix d. deorum" is a trinomial with the subspecies rank abbreviated;
# the binomial is Genus + the final epithet.
ABBREV_MIDDLE = re.compile(r"^([A-Z][a-z-]{2,})\s+[a-z]\.\s+([a-z-]{3,})$")
# "S.mauritia acronyctoides" -- the space after an abbreviated genus is missing
MISSING_SPACE = re.compile(r"^([A-Z])\.([a-z-]{3,})")
# "Daphnia longispina/galeata complex" -- an aggregate; the first name is the
# one that can be looked up
SLASH_AGG = re.compile(r"^([A-Z][a-z-]+\s+[a-z-]+)/[a-z-]+(\s+complex)?$", re.I)


def repair_name(name):
    """Return (repaired_name, what_was_done). Empty note means untouched."""
    n = " ".join(str(name).split())
    notes = []
    # Confusables are only repaired in the EPITHET, never the genus: Ixodes and
    # Oscheius legitimately start with those letters, but "Iongispina" is a
    # capital I standing in for an l.
    toks = n.split()
    fixed = []
    for i, t in enumerate(toks):
        if i > 0 and re.fullmatch(r"[IO][a-z-]{2,}", t):
            fixed.append({"I": "l", "O": "o"}[t[0]] + t[1:])
        elif i > 0:
            fixed.append(re.sub(r"(?<=[a-z])I(?=[a-z])", "l", t))
        else:
            fixed.append(t)
    if fixed != toks:
        notes.append("confusable letter")
        n = " ".join(fixed)
    m = MISSING_SPACE.match(n)
    if m:
        n = f"{m.group(1)}. {m.group(2)}" + n[m.end():]
        notes.append("missing space after initial")
    m = SLASH_AGG.match(n)
    if m:
        n = m.group(1)
        notes.append("aggregate reduced to first name")
    m = ABBREV_MIDDLE.match(n)
    if m:
        n = f"{m.group(1)} {m.group(2)}"
        notes.append("abbreviated subspecies dropped")
    # "Dissosteira Carolina" / "Halyomorpha Halys" -- a capitalised epithet
    m = re.match(r"^([A-Z][a-z-]{2,})\s+([A-Z][a-z-]{2,})$", n)
    if m and m.group(2).lower() not in ("sp", "spp"):
        n = f"{m.group(1)} {m.group(2).lower()}"
        notes.append("epithet decapitalised")
    return n, "; ".join(notes)


def clean_epithet(name):
    """Trim trailing junk that stops a name matching in any database."""
    n = re.sub(r"\s+", " ", deaccent(name)).strip(" .,;:")
    n = re.sub(r"\b(?:sp|spp)\.?\s*(?:nov\.?)?\s*\d*$", "sp.", n, flags=re.I)
    n = re.sub(r"\b(?:cf|aff|nr)\.?\s+", "", n, flags=re.I)
    n = re.sub(r"\s+(?:sensu lato|s\.l\.|sensu stricto|s\.s\.)$", "", n, flags=re.I)
    # strip a taxonomic authority: "Nosema bombycis Nageli, 1857"
    n = re.sub(r",?\s*\(?[A-Z][a-z]+(?:\s*&\s*[A-Z][a-z]+)?\)?,?\s*\d{4}$", "", n)
    return n.strip()



# A cell may separate hosts with commas instead of semicolons:
#   "Bombus honshuensis, Bombus diversus diversus, Bombus ussurensis"
# Splitting on every comma would break "Nosema bombycis Nageli, 1857" and
# "Aedes sp., type host", so a comma split is only taken when at least two of
# the resulting fragments look like binomials and none is a bare year.
BINOMIAL_FRAG = re.compile(r"^[A-Z][a-z-]{2,}\s+[a-z][a-z-]{2,}")
YEAR_FRAG = re.compile(r"^\(?\d{4}\)?$")


def split_hosts(text):
    parts = [p for p in re.split(r";", str(text)) if p.strip()]
    out = []
    for p in parts:
        if "," not in p:
            out.append(p)
            continue
        frags = [f.strip() for f in p.split(",") if f.strip()]
        looks = sum(bool(BINOMIAL_FRAG.match(f)) for f in frags)
        has_year = any(YEAR_FRAG.match(f) for f in frags)
        out.extend(frags if (looks >= 2 and not has_year) else [p])
    return out


def parse_host_cell(cell):
    """
    One host cell -> list of dicts:
      raw       the fragment as written
      name      the name to look up
      rank      'species' | 'genus' | 'other'
      role      'type host', 'intermediate host', ... or ''
      synonyms  alternate names to retry with if `name` fails
    """
    if cell is None:
        return []
    text = str(cell).strip()
    if not text:
        return []

    out = []
    seen_genera = {}          # first letter -> genus named earlier in this cell
    for frag in split_hosts(text):
        frag = frag.strip()
        if not frag:
            continue

        # pull the parentheticals out and sort them into roles vs synonyms
        role, synonyms = "", []
        secondary, context, repair_note = "", "", ""
        for inner in re.findall(r"\(([^)]*)\)", frag):
            inner = inner.strip()
            m = SYNONYM_RE.match(inner)
            if m:
                syn = clean_epithet(m.group(1))
                if syn:
                    synonyms.append(syn)
            elif HYPER_CONTEXT.match(inner) or HYPER_WORD.search(inner):
                # "(in the polychaete annelid Audouinia tentaculata)" names the
                # HOST'S host, not a second microsporidian host
                b = BINOMIAL_IN.search(inner)
                if b:
                    secondary = f"{b.group(1)} {b.group(2)}"
                context = inner.lower()
            elif ROLE_WORDS.match(inner):
                role = role or inner.lower()
                # "(should be Simulium aureum, according to ...)" is a role note
                # that also names a better binomial: keep it as a fallback
                m2 = re.search(r"\b([A-Z][a-z]+ [a-z-]{3,})", inner)
                if m2:
                    synonyms.append(clean_epithet(m2.group(1)))
            elif re.match(r"^[A-Z][a-z]+ [a-z-]+", inner):
                synonyms.append(clean_epithet(inner))       # common name or alt
            # anything else (e.g. "Swamp Guppy") is dropped
        # strip a subgenus first, keeping it on record
        subgenus = ""
        sg = SUBGENUS.match(frag.strip())
        if sg:
            subgenus = sg.group(2)
            frag_use = (sg.group(1) + " " + sg.group(3)).strip()
        else:
            frag_use = frag
        # unparenthesised "... in Genus species ..." names the host's host
        bc = BARE_CONTEXT.search(frag_use)
        if bc:
            if not secondary:
                secondary = bc.group(1)
            context = context or bc.group(0).strip().lower()
            frag_use = frag_use[:bc.start()]
        base = re.sub(r"\([^)]*\)", " ", frag_use)
        base = clean_epithet(base).strip(' "\'')
        if is_description(base):
            # a description, not a name: do not mine a binomial out of it
            out.append({"raw": frag, "name": "", "rank": "unidentified",
                        "role": role, "synonyms": [], "subgenus": subgenus,
                        "secondary_organism": secondary, "context": context,
                    "subgenus": subgenus, "repair": repair_note})
            continue
        if not base or NON_TAXON.match(base):
            continue

        # "D. pulex" after "Daphnia magna" in the same cell
        # "D. pulex" resolves against any genus named earlier in the cell,
        # not only the immediately preceding one
        m = ABBREV_GENUS.match(base)
        if m and m.group(1) in seen_genera:
            base = f"{seen_genera[m.group(1)]} {m.group(2)}"

        base, repair_note = repair_name(base)
        toks = base.split()
        if len(toks) >= 2 and toks[1].lower() in ("sp.", "sp", "spp.", "spp"):
            rank, lookup = "genus", toks[0]
        elif len(toks) >= 2 and re.match(r"^[a-z-]+$", toks[1]):
            rank, lookup = "species", " ".join(toks[:2])
        elif len(toks) == 1 and re.match(r"^[A-Z][a-z]+$", toks[0]):
            rank, lookup = "genus", toks[0]
        else:
            rank, lookup = "other", base
        if rank in ("species", "genus"):
            g = lookup.split()[0]
            seen_genera.setdefault(g[0], g)

        out.append({"raw": frag, "name": lookup, "rank": rank,
                    "role": role, "synonyms": synonyms,
                    "secondary_organism": secondary, "context": context,
                    "subgenus": subgenus, "repair": repair_note})
    return out


# --------------------------------------------------------------------------- #
# cache
# --------------------------------------------------------------------------- #
class Cache:
    """Disk-backed lookup cache. Nothing is fetched twice, ever."""

    def __init__(self, path):
        self.path = path
        self.data = {}
        self.dirty = 0
        if path and os.path.exists(path):
            try:
                with open(path, encoding="utf-8") as fh:
                    self.data = json.load(fh)
                sys.stderr.write(f"cache: {len(self.data)} name(s) loaded from {path}\n")
            except Exception as exc:                        # noqa: BLE001
                sys.stderr.write(f"cache unreadable ({exc}), starting fresh\n")

    def get(self, name):
        return self.data.get(name)

    def put(self, name, record):
        self.data[name] = record
        self.dirty += 1
        if self.dirty >= 25:
            self.flush()

    def flush(self):
        if not self.path or not self.dirty:
            return
        tmp = self.path + ".tmp"
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(self.data, fh, ensure_ascii=False)
        os.replace(tmp, self.path)
        self.dirty = 0


# --------------------------------------------------------------------------- #
# HTTP
# --------------------------------------------------------------------------- #
class Http:
    def __init__(self, delay=0.34, retries=3, timeout=15):
        self.s = requests.Session()
        self.s.headers.update({"User-Agent": "microsporidia-host-lookup/2.0",
                               "Accept": "application/json"})
        self.delay = delay
        self.retries = retries
        self.timeout = timeout
        self.calls = Counter()
        self.errors = Counter()

    def get(self, url, params=None, want="json", source=""):
        self.calls[source] += 1
        for attempt in range(self.retries):
            try:
                r = self.s.get(url, params=params, timeout=self.timeout)
                if r.status_code == 404:
                    return None
                r.raise_for_status()
                time.sleep(self.delay)
                return r.json() if want == "json" else r.text
            except Exception as exc:                        # noqa: BLE001
                self.errors[f"{source}:{type(exc).__name__}"] += 1
                if attempt == self.retries - 1:
                    return None
                time.sleep(2 ** attempt)
        return None


# --------------------------------------------------------------------------- #
# taxonomy sources
# --------------------------------------------------------------------------- #

# --- phylum synonymy across sources ------------------------------------------
# The four sources do not use the same phylum names, so the same organism is
# filed under different phyla depending on which source answered -- which
# splits a real group across two bars and undercounts both.
#   GBIF puts gregarines in Chromista / MYZOZOA; NCBI and COL call the same
#   organisms APICOMPLEXA. Myzozoa also contains dinoflagellates, so the
#   rename is applied only when the class is a genuinely apicomplexan one.
# phylum_normalised is written as an extra column; the raw value is kept.
PHYLUM_ALIASES = {
    # NCBI reports these amoebozoan classes in its phylum slot
    "discosea": "Amoebozoa", "tubulinea": "Amoebozoa", "evosea": "Amoebozoa",
    "conosa": "Amoebozoa", "lobosa": "Amoebozoa", "longamoebia": "Amoebozoa",
    "ectoprocta": "Bryozoa",
    "craniata": "Chordata",
    "vertebrata": "Chordata",
    "urochordata": "Chordata",
    "aschelminthes": "Nematoda",
    "priapulida": "Cephalorhyncha",
    "kinorhyncha": "Cephalorhyncha",
    "loricifera": "Cephalorhyncha",
    "rhombozoa": "Dicyemida",
    "pentastomida": "Arthropoda",
    "myzostomida": "Annelida",
    "echiura": "Annelida",
    "sipuncula": "Annelida",
    "apicomplexa": "Apicomplexa",
}
APICOMPLEXAN_CLASSES = {"conoidasida", "gregarinomorphea", "aconoidasida",
                        "gregarinasina", "coccidia", "hematozoa",
                        "sporozoasida", "gregarinia"}
# COL frequently returns no class at all for gregarines, only an order, so the
# order is tested too -- otherwise "Miozoa / <blank> / Eugregarinida" fails the
# class test and stays under a name no other source uses.
APICOMPLEXAN_ORDERS = {"eugregarinida", "eugregarinorida", "archigregarinida",
                       "neogregarinorida", "blastogregarinorida",
                       "agamococcidiorida", "eucoccidiorida", "adeleorina",
                       "haemosporida", "piroplasmida", "cryptosporidiida",
                       "gregarinida"}


def normalise_phylum(phylum, klass="", order=""):
    p = (phylum or "").strip()
    if not p:
        return ""
    key = p.lower()
    # GBIF says Myzozoa, COL says Miozoa, NCBI and others say Apicomplexa --
    # three names for one group. Both aliases resolve, but only when the class
    # OR the order is genuinely apicomplexan, since Myzozoa/Miozoa also hold
    # the dinoflagellates.
    if key in ("myzozoa", "miozoa"):
        if (klass or "").strip().lower() in APICOMPLEXAN_CLASSES or \
           (order or "").strip().lower() in APICOMPLEXAN_ORDERS:
            return "Apicomplexa"
        return p
    return PHYLUM_ALIASES.get(key, p)


def norm_tax(d):
    """Keep only the seven ranks, blank-filled."""
    return {r: (d.get(r) or "") for r in RANKS}


def tax_completeness(tax):
    return sum(1 for r in RANKS if tax and tax.get(r))


def gbif_lookup(http, name, rank):
    """
    NOTE on fuzzy matching. strict=false lets GBIF correct a misspelled name,
    which is usually what we want -- "Euprocitis" finds Euproctis -- but it
    silently answers a different question from the one asked. matchType and the
    name GBIF actually matched are therefore captured and reported, so a fuzzy
    hit is visible rather than invisible.
    """
    m = http.get(f"{GBIF_BASE}/species/match", {"name": name, "strict": "false"},
                 source="gbif")
    if not m or not m.get("usageKey"):
        return None, None, {}
    if m.get("matchType") == "NONE":
        return None, None, {}
    meta = {"match_type": m.get("matchType", ""),
            "match_confidence": m.get("confidence", ""),
            "matched_scientific_name": m.get("scientificName") or m.get("canonicalName") or "",
            "taxonomic_status": m.get("status", "")}
    key = m["usageKey"]
    data = http.get(f"{GBIF_BASE}/species/{key}", source="gbif") or m
    tax = norm_tax(data)

    env = set()
    hab = http.get(f"{GBIF_BASE}/species/{key}/habitats", source="gbif")
    for h in (hab or {}).get("results", []) or []:
        env |= env_words(h.get("habitat", ""))
    if not env:
        prof = http.get(f"{GBIF_BASE}/species/{key}/speciesProfiles", source="gbif")
        for p in (prof or {}).get("results", []) or []:
            if p.get("habitat"):
                env |= env_words(p["habitat"])
            if p.get("marine"):
                env.add("marine")
            if p.get("freshwater"):
                env.add("freshwater")
            if p.get("terrestrial"):
                env.add("terrestrial")
    return tax, sorted(env), meta


def worms_lookup(http, name, rank):
    rec = http.get(f"{WORMS_API}/AphiaRecordsByMatchNames",
                   {"scientificnames[]": name, "marine_only": "false"},
                   source="worms")
    if not rec:
        return None, None, {}
    first = rec[0]
    if isinstance(first, list):
        first = first[0] if first else None
    if not first:
        return None, None, {}
    # WoRMS reports exact / like / phonetic / near_1 / near_2
    meta = {"match_type": first.get("match_type", ""),
            "matched_scientific_name": first.get("scientificname", ""),
            "taxonomic_status": first.get("status", "")}

    tax = norm_tax({
        "kingdom": first.get("kingdom"), "phylum": first.get("phylum"),
        "class": first.get("class"), "order": first.get("order"),
        "family": first.get("family"), "genus": first.get("genus"),
        "species": first.get("scientificname") if first.get("rank") == "Species" else "",
    })
    # WoRMS environment flags are authoritative and need no scraping
    env = set()
    for flag, label in (("isMarine", "marine"), ("isBrackish", "brackish"),
                        ("isFreshwater", "freshwater"),
                        ("isTerrestrial", "terrestrial")):
        if first.get(flag) == 1:
            env.add(label)
    return tax, sorted(env), meta


def col_lookup(http, name, rank):
    res = http.get(f"{COL_API}/dataset/{COL_DATASET}/nameusage/search",
                   {"q": name, "limit": 5, "content": "SCIENTIFIC_NAME"},
                   source="col")
    hits = (res or {}).get("result") or []
    if not hits:
        return None, None, {}
    usage = hits[0]
    meta = {"matched_scientific_name":
            (usage.get("usage", {}).get("name", {}) or {}).get("scientificName", ""),
            "taxonomic_status": (usage.get("usage", {}) or {}).get("status", "")}
    tax = {}
    for c in usage.get("classification") or []:
        r = (c.get("rank") or "").lower()
        if r in RANKS:
            tax[r] = c.get("name")
    if usage.get("usage", {}).get("label"):
        tax.setdefault("species", usage["usage"].get("name", {}).get("scientificName"))
    env = set()
    for key in ("environments", "environment"):
        for e in usage.get(key) or []:
            env |= env_words(str(e))
    return (norm_tax(tax) if tax else None), sorted(env), meta


def ncbi_lookup(http, name, rank, email, api_key):
    """
    NCBI Taxonomy. The fallback that rescues insect and nematode names the
    other three miss. esearch for a TaxID, then efetch the full lineage.
    """
    common = {"db": "taxonomy", "email": email, "tool": "microsporidia-host-lookup"}
    if api_key:
        common["api_key"] = api_key
    q = dict(common, term=f'"{name}"[Scientific Name]', retmode="json")
    res = http.get(f"{NCBI_BASE}/esearch.fcgi", q, source="ncbi")
    ids = ((res or {}).get("esearchresult") or {}).get("idlist") or []
    match_type = "exact"
    if not ids:      # retry as a free-text search
        q = dict(common, term=name, retmode="json")
        res = http.get(f"{NCBI_BASE}/esearch.fcgi", q, source="ncbi")
        ids = ((res or {}).get("esearchresult") or {}).get("idlist") or []
        match_type = "free_text"
    if not ids:
        return None, None, {}

    xml = http.get(f"{NCBI_BASE}/efetch.fcgi", dict(common, id=ids[0], retmode="xml"),
                   want="text", source="ncbi")
    if not xml:
        return None, None, {}
    try:
        root = ET.fromstring(xml)
    except ET.ParseError:
        return None, None, {}
    node = root.find("Taxon")
    if node is None:
        return None, None, {}
    meta = {"match_type": match_type,
            "matched_scientific_name": node.findtext("ScientificName") or ""}

    tax = {}
    for t in node.findall("./LineageEx/Taxon"):
        r = (t.findtext("Rank") or "").lower()
        if r in RANKS:
            tax[r] = t.findtext("ScientificName")
    own_rank = (node.findtext("Rank") or "").lower()
    if own_rank in RANKS:
        tax[own_rank] = node.findtext("ScientificName")
    # NCBI uses superkingdom/clade rather than kingdom for some lineages
    if not tax.get("kingdom"):
        for t in node.findall("./LineageEx/Taxon"):
            if (t.findtext("Rank") or "").lower() in ("superkingdom", "clade") \
                    and t.findtext("ScientificName") in ("Metazoa", "Viridiplantae",
                                                         "Fungi", "Eukaryota"):
                tax["kingdom"] = t.findtext("ScientificName")
                break
    return (norm_tax(tax) if tax else None), None, meta   # NCBI has no habitats


# --------------------------------------------------------------------------- #
# environment
# --------------------------------------------------------------------------- #
def env_words(text):
    """Free text -> canonical environment terms."""
    t = str(text).lower()
    out = set()
    if re.search(r"marine|salt\s*water|saltwater|sea\b|ocean|deep-sea|pelagic"
                 r"|benthic|reef|estuar", t):
        out.add("marine")
    if "brackish" in t or "estuar" in t:
        out.add("brackish")
    if re.search(r"fresh\s*water|freshwater|limnetic|lacustrine|riverine|pond"
                 r"|lake|stream", t):
        out.add("freshwater")
    if re.search(r"terrestrial|soil|land\b|forest|compost|tree-bark|domestic", t):
        out.add("terrestrial")
    return out


# Coarse group for plotting, from the lineage. First match wins.
HOST_GROUP_RULES = [
    ("Insecta",          {"class": ["Insecta", "Hexapoda"]}),
    ("Crustacea",        {"class": ["Malacostraca", "Branchiopoda", "Copepoda",
                                    "Maxillopoda", "Ostracoda", "Hexanauplia",
                                    "Thecostraca", "Cephalocarida", "Ichthyostraca"],
                          "phylum": ["Crustacea"], "subphylum": ["Crustacea"]}),
    ("Arachnida & myriapods", {"class": ["Arachnida", "Chilopoda", "Diplopoda",
                                         "Merostomata", "Pycnogonida"]}),
    ("Fish",             {"class": ["Actinopterygii", "Chondrichthyes",
                                    "Sarcopterygii", "Actinopteri", "Teleostei",
                                    "Myxini", "Petromyzonti", "Cephalaspidomorphi",
                                    "Hyperoartia", "Elasmobranchii", "Holocephali"]}),
    ("Mammals",          {"class": ["Mammalia"]}),
    ("Birds",            {"class": ["Aves"]}),
    ("Reptiles & amphibians", {"class": ["Reptilia", "Amphibia", "Squamata",
                                         "Testudines", "Lepidosauria"]}),
    ("Nematoda",         {"phylum": ["Nematoda", "Nematomorpha"]}),
    ("Annelida",         {"phylum": ["Annelida"]}),
    ("Mollusca",         {"phylum": ["Mollusca"]}),
    ("Platyhelminthes",  {"phylum": ["Platyhelminthes"]}),
    ("Rotifera & bryozoa", {"phylum": ["Rotifera", "Bryozoa", "Gastrotricha",
                                       "Entoprocta", "Tardigrada"]}),
    ("Cnidaria & sponges", {"phylum": ["Cnidaria", "Porifera", "Ctenophora"]}),
    ("Echinodermata",    {"phylum": ["Echinodermata", "Chaetognatha"]}),
    ("Other invertebrates", {"kingdom": ["Metazoa", "Animalia"]}),
    ("Protists",         {"kingdom": ["Chromista", "Protozoa", "Protista",
                                      "Eukaryota", "Harosa", "Discoba"]}),
    ("Fungi",            {"kingdom": ["Fungi"]}),
    ("Plants & algae",   {"kingdom": ["Plantae", "Viridiplantae"]}),
]


def host_group(tax):
    if not tax:
        return "Unresolved"
    for label, rules in HOST_GROUP_RULES:
        for rank, names in rules.items():
            val = (tax.get(rank) or "").strip()
            if val and any(val.lower() == n.lower() for n in names):
                return label
    if any(tax.get(r) for r in RANKS):
        return "Other / unplaced"
    return "Unresolved"


# --------------------------------------------------------------------------- #
class ClaudeClient:
    """Optional last-resort environment call."""

    def __init__(self, api_key, model="claude-sonnet-4-5"):
        self.key = api_key
        self.model = model

    def classify(self, name, tax):
        lineage = ", ".join(f"{r}: {tax[r]}" for r in RANKS if tax and tax.get(r))
        prompt = (
            f"Organism: {name}\n"
            f"{('Lineage: ' + lineage) if lineage else ''}\n\n"
            "Which habitats does this organism live in? Choose from: marine, "
            "brackish, freshwater, terrestrial. Reply with ONLY a JSON object: "
            '{\"environments\": [...], \"confidence\": 0-100}. '
            "If you are not confident, return an empty list."
        )
        try:
            r = requests.post(
                "https://api.anthropic.com/v1/messages",
                headers={"x-api-key": self.key, "anthropic-version": "2023-06-01",
                         "content-type": "application/json"},
                json={"model": self.model, "max_tokens": 200,
                      "messages": [{"role": "user", "content": prompt}]},
                timeout=30)
            r.raise_for_status()
            txt = "".join(b.get("text", "") for b in r.json().get("content", []))
            txt = re.sub(r"```json|```", "", txt).strip()
            obj = json.loads(txt)
            if int(obj.get("confidence", 0)) < 70:
                return [], obj.get("confidence", 0)
            envs = [e for e in obj.get("environments", []) if e in ENVIRONMENTS]
            return envs, obj.get("confidence", 0)
        except Exception:                                   # noqa: BLE001
            return [], 0


# --------------------------------------------------------------------------- #
def _edits(a, b):
    """Levenshtein distance, short strings only."""
    if abs(len(a) - len(b)) > 2:
        return 99
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


def looks_like_a_taxon(name):
    """
    "dog", "laboratory", "hamsters" are not scientific names. Sending them to a
    fuzzy search invites a nonsense hit -- COL answers "dog" with the fossil
    brachiopod genus Dogdospirifer, and it is a COMPLETE seven-rank lineage, so
    a completeness-based chooser prefers it over every correct blank.
    """
    n = name.strip()
    if not n or not n[0].isupper():
        return False
    if len(n) < 4:
        return False
    return bool(re.match(r"^[A-Z][a-z-]+(?:\s+[a-z-]+)?$", n))


# Microsporidia infect animals and protists, occasionally fungi. A lineage that
# lands in Plantae is a HOMONYM: the same genus name is legal in both botany and
# zoology, so COL answers "Limnetis" (a clam shrimp) with a grass and
# "Hymenolepis" (a tapeworm) with a daisy. Both were silently plausible-looking
# rows. Any source returning one of these kingdoms is rejected and the next
# source is used instead.
IMPLAUSIBLE_KINGDOMS = {"plantae", "viridiplantae", "bacteria", "archaea",
                        "viruses", "monera"}


def kingdom_plausible(tax):
    k = (tax or {}).get("kingdom", "").strip().lower()
    if k and k in IMPLAUSIBLE_KINGDOMS:
        return False, f"kingdom {tax['kingdom']} is not a microsporidian host"
    return True, ""


# Authorities and subgenera are part of a returned scientific name but not part
# of the NAME: GBIF answers a genus-level query with "Simulium Latreille, 1802"
# and COL with "Aedes (Isoaedes)". Comparing token 2 against the epithet then
# reads "Latreille" as the epithet and rejects a perfectly good match. This was
# by far the largest single cause of unresolved hosts.
AUTHORITY = re.compile(r"\s*\(?\b[A-Z][A-Za-z.'-]*\b\.?(?:\s*&\s*[A-Z][A-Za-z.'-]*)*"
                       r"\s*,?\s*\d{4}\)?\s*$")


def strip_authority(name):
    """'Simulium Latreille, 1802' -> 'Simulium'; 'Aedes (Isoaedes)' -> 'Aedes'."""
    n = re.sub(r"\([^)]*\)", " ", str(name))          # subgenus / (Author, year)
    n = AUTHORITY.sub("", n)
    toks = n.split()
    if not toks:
        return ""
    # keep the genus plus only genuinely lower-case epithets; an author surname
    # is capitalised, an epithet never is
    out = [toks[0]]
    for t in toks[1:]:
        if re.fullmatch(r"[a-z][a-z-]+", t):
            out.append(t)
        else:
            break
    return " ".join(out[:2])


# Latin gender endings change with the genus: complana/complanum, aceri/aceris,
# unipuncta/unipunctata, thumi/thummi. The stem is what matters.
GENDER_END = re.compile(r"(us|um|a|is|e|i|ii|ae|orum|arum)$")


def _stem(epithet):
    e = GENDER_END.sub("", epithet)
    return e if len(e) >= 3 else epithet


def name_agrees(query, returned, rank):
    """
    Does the lineage the source returned belong to the name asked for?

    Returns True (accept), False (reject) or None (provisional: a two-edit
    genus difference, accepted later only if a second source corroborates it).
    """
    if not returned:
        return True, ""
    q = re.sub(r"[^A-Za-z ]", " ", deaccent(query)).split()
    r = strip_authority(deaccent(returned)).split()
    if not q or not r:
        return True, ""
    gq, gr = q[0].lower(), r[0].lower()
    tol = 1 if min(len(gq), len(gr)) >= 6 else 0
    d = _edits(gq, gr)
    if d > tol:
        if d <= 2 and min(len(gq), len(gr)) >= 6:
            return None, f"genus differs by {d}: asked {q[0]}, got {r[0]}"
        return False, f"genus differs: asked {q[0]}, got {r[0]}"
    # A genus-only answer to a species query is a HIGHER-RANK match, not a
    # wrong one: the genus lineage is still correct down to genus.
    if rank == "species" and len(q) > 1 and len(r) < 2:
        return True, ""
    if rank == "species" and len(q) > 1 and len(r) > 1:
        eq, er = q[1].lower(), r[1].lower()
        if eq != er and _stem(eq) != _stem(er):
            de = _edits(eq, er)
            if de > (2 if min(len(eq), len(er)) >= 6 else 1):
                return False, f"epithet differs: asked {q[1]}, got {r[1]}"
    return True, ""


def consensus_rank(per_source, rank):
    """
    Majority value for one rank across the accepted sources.

    Also reports whether the vote was TIED. A 1-1 split has no majority, and
    resolving it by whichever source happened to be queried first is a coin
    toss dressed up as a decision -- see --on-tie.
    """
    # Vote on the NORMALISED phylum. Myzozoa (GBIF), Miozoa (COL) and
    # Apicomplexa (NCBI) are one group under three names, and Discosea or
    # Tubulinea (NCBI) sit inside Amoebozoa -- counting the raw strings turned
    # agreement into a tie and threw the lineage away.
    def _val(d):
        v = d["tax"].get(rank, "")
        if rank == "phylum":
            v = normalise_phylum(v, d["tax"].get("class", ""),
                                 d["tax"].get("order", ""))
        return v
    vals = [(src, _val(d)) for src, d in per_source.items()
            if d["accepted"] and d["tax"] and d["tax"].get(rank)]
    if not vals:
        return "", [], [], False
    tally = Counter(v for _, v in vals)
    ranked = tally.most_common()
    top, n_top = ranked[0]
    tied = len(ranked) > 1 and ranked[1][1] == n_top
    agree = [s for s, v in vals if v == top]
    disagree = [f"{s}:{v}" for s, v in vals if v != top]
    return top, agree, disagree, tied


def resolve_host(http, name, rank, synonyms, args, claude=None):
    """
    Ask EVERY source, keep every answer, then decide by agreement.

    The previous version kept whichever lineage had the most ranks filled. That
    is a completeness contest, not a correctness one, and a confident wrong
    answer beats three correct blanks every time -- which is how "dog" became a
    brachiopod. Now:
      1. sources that answered about a different organism are rejected outright
      2. the phylum is decided by majority vote among the survivors
      3. the lineage is taken from a source that voted with the majority,
         completeness breaking ties
      4. every disagreement at phylum or order is recorded for review
    """
    rec = {"lookup_name": name, "rank": rank, "matched_name": "",
           "taxonomy_source": "", "env_sources": {}, "n_ranks": 0,
           "match_type": "", "match_confidence": "",
           "matched_scientific_name": "", "taxonomic_status": "",
           "phylum_tied": "",
           "n_sources_queried": 0, "n_sources_found": 0, "sources_found": "",
           "sources_rejected": "", "phylum_agreement": "", "phylum_dispute": "",
           "order_dispute": "", "flag": ""}
    for r in RANKS:
        rec[r] = ""

    if not looks_like_a_taxon(name):
        rec["flag"] = "not a scientific name; not looked up"
        rec["host_group"] = "Unresolved"
        return rec

    sources = []
    if not args.no_gbif:
        sources.append(("GBIF", gbif_lookup))
    if not args.no_col:
        sources.append(("COL", col_lookup))
    if not args.no_worms:
        sources.append(("WoRMS", worms_lookup))
    if not args.no_ncbi:
        sources.append(("NCBI", lambda h, n, rk: ncbi_lookup(
            h, n, rk, args.email, args.ncbi_api_key)))

    candidates = [name] + [s for s in synonyms if s and s != name]
    per_source, env_by_source = OrderedDict(), OrderedDict()
    used_name = name

    for cand in candidates:
        for src_name, fn in sources:
            if src_name in per_source and per_source[src_name]["accepted"]:
                continue                      # already have a good answer
            tax, env, meta = fn(http, cand, rank)
            meta = meta or {}
            rec["n_sources_queried"] += 1
            got = (meta.get("matched_scientific_name") or "").strip()
            ok, why = name_agrees(cand, got, rank)
            # NCBI's free-text search resolves COMMON names -- "Atlantic cod"
            # comes back as Gadus morhua, which the name check then rejects
            # because the genus does not match. Those are accepted, but always
            # flagged, because a free-text hit can also be loose ("Bamboo rat"
            # -> Rattus).
            if (ok is not True and src_name == "NCBI"
                    and meta.get("match_type") == "free_text"
                    and not args.no_vernacular and got):
                ok, why = True, ""
                meta["vernacular"] = True
            if ok and not args.allow_any_kingdom:
                ok, why = kingdom_plausible(tax)
            entry = {"tax": tax or {}, "meta": meta,
                     "accepted": bool(tax) and ok is True,
                     "provisional": bool(tax) and ok is None,
                     "reject_reason": "" if ok is True else why,
                     "query": cand, "returned": got}
            prev = per_source.get(src_name)
            if prev is None or (entry["accepted"] and not prev["accepted"]):
                per_source[src_name] = entry
            if env and ok:
                env_by_source.setdefault(src_name, sorted(set(env)))
        if any(d["accepted"] for d in per_source.values()):
            used_name = cand
            break

    # Promote provisional (two-edit) matches that a second source corroborates:
    # the same corrected scientific name coming back independently is far better
    # evidence than the edit distance.
    prov = {k: v for k, v in per_source.items() if v.get("provisional")}
    if prov:
        seen = Counter()
        for k, v in per_source.items():
            nm = re.sub(r"[^a-z ]", "", (v.get("returned") or "").lower()).strip()
            if nm and (v["accepted"] or v.get("provisional")):
                seen[" ".join(nm.split()[:2])] += 1
        for k, v in prov.items():
            nm = re.sub(r"[^a-z ]", "", (v.get("returned") or "").lower()).strip()
            if seen[" ".join(nm.split()[:2])] >= 2:
                v["accepted"] = True
                v["reject_reason"] = ""
                v["corroborated"] = True

    accepted = {k: v for k, v in per_source.items() if v["accepted"]}
    rejected = {k: v for k, v in per_source.items()
                if not v["accepted"] and v["tax"]}
    rec["n_sources_found"] = len(accepted)
    rec["sources_found"] = ";".join(accepted)
    rec["sources_rejected"] = ";".join(
        f"{k}({v['returned']}: {v['reject_reason']})" for k, v in rejected.items())

    if accepted:
        phy, agree, dis, tied = consensus_rank(per_source, "phylum")
        rec["phylum_agreement"] = (f"{len(agree)}/{len(agree) + len(dis)}"
                                   if (agree or dis) else "")
        rec["phylum_dispute"] = ";".join(dis)
        rec["phylum_tied"] = "yes" if tied else ""
        _, o_agree, o_dis, o_tied = consensus_rank(per_source, "order")
        rec["order_dispute"] = ";".join(o_dis)

        # prefer a source that voted with the majority phylum; completeness
        # only breaks ties among those
        pool = [k for k in accepted if not phy or
                accepted[k]["tax"].get("phylum", "") in ("", phy)]
        if not pool:
            pool = list(accepted)
        best = max(pool, key=lambda k: (tax_completeness(accepted[k]["tax"]),
                                        -list(accepted).index(k)))
        btax, bmeta = accepted[best]["tax"], accepted[best]["meta"]
        rec.update({r: btax.get(r, "") for r in RANKS})
        rec["taxonomy_source"] = best
        rec["n_ranks"] = tax_completeness(btax)
        rec["matched_name"] = accepted[best]["query"]
        rec["match_type"] = bmeta.get("match_type", "")
        rec["match_confidence"] = bmeta.get("match_confidence", "")
        rec["matched_scientific_name"] = bmeta.get("matched_scientific_name", "")
        rec["taxonomic_status"] = bmeta.get("taxonomic_status", "")

        # Should a disputed host be USED, or only reported? Neither answer is
        # right for every purpose, so it is a switch:
        #   include  take the majority and flag it (default; a 2-1 majority is
        #            usually the correct call and dropping it loses real data)
        #   exclude  withhold the lineage entirely until it has been checked, so
        #            a disputed host cannot reach a figure. It is still listed.
        # A TIE has no majority at all and is withheld unless --on-tie
        # source-order is given, because otherwise the winner is decided by
        # which source was queried first.
        withhold = ""
        if tied and args.on_tie == "blank":
            withhold = "TIE: no majority phylum, lineage withheld"
        elif dis and args.disputed == "exclude":
            withhold = "DISPUTED: lineage withheld (--disputed exclude)"
        if withhold:
            for r_ in RANKS:
                rec[r_] = ""
            rec["taxonomy_source"] = ""
            rec["n_ranks"] = 0
            rec["flag"] = withhold + "; " + ";".join(dis)
            rec["env_sources"] = env_by_source
            rec["host_group"] = "Unresolved"
            rec["name_changed"] = ""
            return rec

        flags = []
        if any(v["meta"].get("vernacular") for v in accepted.values()):
            flags.append("resolved from a common name -- verify")
        if any(v.get("corroborated") for v in accepted.values()):
            flags.append("spelling corrected, corroborated by 2+ sources")
        if tied:
            flags.append("TIE broken by source order")
        if dis:
            flags.append("PHYLUM DISAGREEMENT")
        if o_dis:
            flags.append("order disagreement")
        if len(accepted) == 1:
            flags.append("single source only")
        if rejected:
            flags.append("a source answered about a different organism")
        rec["flag"] = "; ".join(flags)
    else:
        rec["flag"] = ("all sources rejected: " + rec["sources_rejected"]
                       if rejected else "no source had this name")

    rec["env_sources"] = env_by_source
    if claude and not env_by_source and accepted:
        envs, conf = claude.classify(name, {r: rec[r] for r in RANKS})
        if envs:
            rec["env_sources"]["Claude"] = envs
            rec["claude_confidence"] = conf
    rec["host_group"] = host_group({r: rec[r] for r in RANKS})

    got = (rec.get("matched_scientific_name") or "").strip()
    rec["name_changed"] = ""
    if got:
        a = re.sub(r"[^a-z ]", "", (rec["matched_name"] or name).lower()).split()
        b = re.sub(r"[^a-z ]", "", got.lower()).split()
        if a[:2] != b[:2]:
            rec["name_changed"] = got
    return rec


def unify_env(db_env, rec_env_sources, use_db=False):
    """
    Environment is a property of the HOST ORGANISM, so it comes from the host
    lookups only: the union of what the sources report, WoRMS flags first.

    The database's "Host Environment" column is deliberately NOT used by
    default. It is recorded per microsporidia row, not per host, so the same
    host resolves to different environments depending on which parasite row it
    appears in -- and it is populated for only a small minority of rows, which
    would make the figure inconsistent between hosts that happen to have a
    curated value and hosts that do not. It is written to a separate column so
    the two can be compared. --use-db-environment restores the old precedence.
    """
    if use_db and db_env:
        return sorted(db_env), "database"
    order = ["WoRMS", "GBIF", "COL", "Claude"]
    out, used = set(), []
    for src in order:
        vals = rec_env_sources.get(src)
        if vals:
            out |= set(vals)
            used.append(src)
    for src, vals in rec_env_sources.items():
        if src not in order and vals:
            out |= set(vals)
            used.append(src)
    return sorted(out), "+".join(used)


# --------------------------------------------------------------------------- #
def nzchar_py(x):
    return bool(x) and bool(str(x).strip())


def norm_key(text):
    return re.sub(r"[^a-z0-9]", "", deaccent(str(text)).lower())


def find_col(columns, wanted):
    """Header lookup on a normalised key (micro sign vs mu, spacing, case)."""
    for c in columns:
        if norm_key(c) == norm_key(wanted):
            return c
    return None


def main():
    p = argparse.ArgumentParser(
        description="Resolve microsporidia host names to taxonomy and habitat.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--db", help="attribute database .xlsx")
    p.add_argument("--file", help="alternative input: one host name per line")
    p.add_argument("--sheet", default="Actively Updated Masterlist")
    p.add_argument("--species-col", default="Species Name")
    p.add_argument("--host-cols", nargs="+",
                   default=["Natural Host(s)", "Experimental Host(s)"])
    p.add_argument("--env-col", default="Host Environment")
    p.add_argument("--outdir", default="host_taxonomy_output")
    p.add_argument("--cache", default="host_lookup_cache.json")
    p.add_argument("--email", default=os.environ.get("NCBI_EMAIL", ""),
                   help="required by NCBI; also used as a contact header")
    p.add_argument("--ncbi-api-key", default=os.environ.get("NCBI_API_KEY", ""))
    p.add_argument("--claude-api-key", default=os.environ.get("ANTHROPIC_API_KEY", ""))
    p.add_argument("--use-db-environment", action="store_true",
                   help="let the spreadsheet's Host Environment column override "
                        "the host lookups (off by default: environment is taken "
                        "from the host organism only)")
    p.add_argument("--no-vernacular", action="store_true",
                   help="do not accept NCBI free-text hits on common names "
                        "such as 'Atlantic cod'; they are flagged when used")
    p.add_argument("--allow-any-kingdom", action="store_true",
                   help="keep lineages that place a host in Plantae, Bacteria "
                        "etc.; by default these are rejected as genus homonyms")
    p.add_argument("--disputed", choices=["include", "exclude"],
                   default="include",
                   help="what to do when the sources disagree on phylum or "
                        "order: 'include' uses the majority and flags it, "
                        "'exclude' withholds the lineage until you have checked "
                        "it (the host is reported either way)")
    p.add_argument("--on-tie", choices=["blank", "source-order"],
                   default="blank",
                   help="a tied vote has no majority; 'blank' withholds the "
                        "lineage, 'source-order' takes whichever source was "
                        "queried first")
    p.add_argument("--refresh", default="",
                   help="comma-separated cache entries to discard and re-fetch "
                        "before running, e.g. 'unresolved,rejected'. Choose "
                        "from: all, unresolved, single-source, rejected, "
                        "disputed, tied, changed-name, implausible, flagged. "
                        "Names the parser now spells differently are re-fetched "
                        "automatically and need no selector.")
    p.add_argument("--refresh-names", action="store_true",
                   help="refetch cache entries written before the name-match "
                        "check was added")
    p.add_argument("--limit", type=int, help="only look up the first N new hosts")
    p.add_argument("--offline", action="store_true",
                   help="use the cache only, make no network calls")
    p.add_argument("--delay", type=float, default=0.34, help="seconds between calls")
    for src in ("gbif", "col", "worms", "ncbi"):
        p.add_argument(f"--no-{src}", action="store_true", help=f"skip {src.upper()}")
    args = p.parse_args()

    # Every input path is checked before any work starts. A missing file is an
    # error, never a silent fall-back to a partial run.
    for label, path in [("database", getattr(args, "db", None)),
                        ("host list", getattr(args, "file", None)),
                        ("compare table", getattr(args, "compare", None))]:
        if path and not os.path.exists(path):
            sys.exit(f"{label} not found: {path}")
    for m in (getattr(args, "meta", None) or []):
        if not os.path.exists(m):
            sys.exit(f"--meta file not found: {m}")

    if not args.db and not args.file:
        p.error("give --db (spreadsheet) or --file (name list)")
    if not args.offline and not args.no_ncbi and not args.email:
        p.error("NCBI requires --email (or set NCBI_EMAIL); or pass --no-ncbi")

    os.makedirs(args.outdir, exist_ok=True)
    out = lambda f: os.path.join(args.outdir, f)

    # ---- read input -------------------------------------------------------
    pairs = []          # (microsporidia species, category, host dict, db_env)
    if args.db:
        try:
            import pandas as pd
        except ImportError:
            sys.exit("pandas required: pip install pandas openpyxl")
        df = pd.read_excel(args.db, sheet_name=args.sheet, dtype=str).fillna("")
        sp_col = find_col(df.columns, args.species_col)
        if not sp_col:
            sys.exit(f"column not found: {args.species_col}\n  have: {list(df.columns)}")
        host_cols = [c for c in (find_col(df.columns, h) for h in args.host_cols) if c]
        if not host_cols:
            sys.exit(f"none of the host columns found: {args.host_cols}")
        env_col = find_col(df.columns, args.env_col)
        sys.stderr.write(f"host columns: {host_cols}\n")

        n_hdr = 0
        for _, row in df.iterrows():
            sp = str(row[sp_col]).strip()
            if not sp:
                continue
            if norm_key(sp) == norm_key(args.species_col):
                n_hdr += 1
                continue                     # header row pasted into the data
            db_env = env_words(row[env_col]) if env_col else set()
            for hc in host_cols:
                for h in parse_host_cell(row[hc]):
                    h["host_column"] = hc
                    h["cell"] = str(row[hc]).strip()
                    pairs.append((sp, classify_microsporidia(sp), h, db_env))
        if n_hdr:
            sys.stderr.write(f"note: skipped {n_hdr} repeated header row(s)\n")
    else:
        with open(args.file, encoding="utf-8") as fh:
            for line in fh:
                for h in parse_host_cell(line.strip()):
                    h["host_column"] = "input file"
                    pairs.append(("", "", h, set()))

    # ---- expand abbreviated genera across the whole database ---------------
    # "C. tarsalis" is unresolvable on its own, but the database elsewhere names
    # Culex, Culiseta and Cotylurus. Expansion is only accepted when exactly ONE
    # full genus in the entire host list starts with that letter and is already
    # paired with that epithet somewhere, or when only one candidate genus
    # exists at all; anything ambiguous is left alone and reported.
    full_genera = Counter()
    genus_species = set()
    for _, _, h, _ in pairs:
        t = h["name"].split()
        if len(t) >= 1 and re.fullmatch(r"[A-Z][a-z-]{2,}", t[0]):
            full_genera[t[0]] += 1
            if len(t) > 1:
                genus_species.add((t[0], t[1].lower()))
    expanded = []
    for _, _, h, _ in pairs:
        m = re.fullmatch(r"([A-Z])\.?\s+([a-z-]{3,})", h["name"])
        if not m:
            continue
        initial, epi = m.group(1), m.group(2).lower()
        cands = [g for g in full_genera if g.startswith(initial)]
        exact = [g for g in cands if (g, epi) in genus_species]
        pick = exact[0] if len(exact) == 1 else (cands[0] if len(cands) == 1 else None)
        if pick:
            expanded.append((h["name"], f"{pick} {epi}",
                             "unique match" if len(exact) == 1 else "only candidate"))
            h["abbrev_expanded_from"] = h["name"]
            h["name"] = f"{pick} {epi}"
        else:
            h["abbrev_ambiguous"] = ";".join(sorted(cands)[:6])
    if expanded:
        sys.stderr.write(f"expanded {len(expanded)} abbreviated genus name(s) "
                         f"using genera named elsewhere in the database\n")
        for a, b, why in expanded[:8]:
            sys.stderr.write(f"    {a:<22} -> {b}  ({why})\n")
    amb = [h for _, _, h, _ in pairs if h.get("abbrev_ambiguous")]
    if amb:
        sys.stderr.write(f"{len(amb)} abbreviated name(s) were ambiguous and "
                         f"left as they are\n")

    unique_names = OrderedDict()
    for _, _, h, _ in pairs:
        unique_names.setdefault(h["name"], h)
    print(f"{len(pairs)} species-host pairs, {len(unique_names)} unique host names")

    rank_counts = Counter(h["rank"] for _, _, h, _ in pairs)
    print("  host name ranks: " + ", ".join(f"{k} {v}" for k, v in rank_counts.items()))

    # ---- look up ----------------------------------------------------------
    cache = Cache(args.cache)
    http = Http(delay=args.delay)
    claude = ClaudeClient(args.claude_api_key) if args.claude_api_key else None
    # ---- selective cache invalidation --------------------------------------
    # The cache is keyed by HOST NAME and stores a finished record, so:
    #   * a name the parser now spells differently is a NEW key and is fetched
    #     automatically -- nothing to do
    #   * a name that is unchanged keeps whatever the OLD resolution rules
    #     decided, which is wrong for anything the new rules would now accept
    #     or reject
    #   * phylum_normalised is computed when the tables are written, not when
    #     the lookup happens, so taxonomy renames need no refetch at all
    # Only the affected entries need discarding, which is far cheaper than
    # deleting the cache.
    SELECTORS = {
        "unresolved":   lambda r: not r.get("n_ranks"),
        "single-source": lambda r: r.get("n_sources_found") == 1,
        "rejected":     lambda r: bool(r.get("sources_rejected")),
        "disputed":     lambda r: bool(r.get("phylum_dispute")),
        "tied":         lambda r: r.get("phylum_tied") == "yes",
        "changed-name": lambda r: bool(r.get("name_changed")),
        "implausible":  lambda r: (r.get("kingdom", "") or "").lower()
                                  in IMPLAUSIBLE_KINGDOMS,
        "flagged":      lambda r: bool(r.get("flag")),
    }
    if args.refresh:
        wanted = [w.strip() for w in args.refresh.split(",") if w.strip()]
        bad = [w for w in wanted if w not in SELECTORS and w != "all"]
        if bad:
            sys.exit(f"--refresh: unknown selector(s) {bad}. "
                     f"Choose from: all, {', '.join(sorted(SELECTORS))}")
        if "all" in wanted:
            drop = set(cache.data)
        else:
            drop = {n for n, r in cache.data.items() if isinstance(r, dict)
                    and any(SELECTORS[w](r) for w in wanted)}
        print(f"\n--refresh {args.refresh}: discarding {len(drop)} of "
              f"{len(cache.data)} cached entries; the rest are reused")
        for w in wanted:
            if w != "all":
                n = sum(1 for n_, r in cache.data.items()
                        if isinstance(r, dict) and SELECTORS[w](r))
                print(f"  {n:>6}  {w}")
        for n in drop:
            cache.data.pop(n, None)
        cache.dirty += 1

    # entries written by an older version lack the name-match fields; they are
    # still usable, but the fuzzy-match report will be blank for them unless
    # they are refetched
    stale = [n for n, r in cache.data.items()
             if isinstance(r, dict) and "n_sources_found" not in r]
    if stale:
        sys.stderr.write(
            f"note: {len(stale)} cached entr{'y' if len(stale)==1 else 'ies'} "
            f"predate the cross-source consensus check and cannot be "
            f"re-examined; pass --refresh-names to refetch them\n")
        if args.refresh_names:
            for n in stale:
                cache.data.pop(n, None)
            cache.dirty += 1
    todo = [n for n in unique_names if cache.get(n) is None]
    print(f"  cached {len(unique_names) - len(todo)}, to look up {len(todo)}")
    if args.offline:
        print("  --offline: no network calls; uncached hosts stay unresolved")
        todo = []
    if args.limit:
        todo = todo[:args.limit]
        print(f"  --limit {args.limit}: stopping after that many")

    for i, name in enumerate(todo, 1):
        h = unique_names[name]
        sys.stderr.write(f"[{i}/{len(todo)}] {name} ... ")
        rec = resolve_host(http, name, h["rank"], h.get("synonyms", []), args, claude)
        cache.put(name, rec)
        sys.stderr.write(f"{rec['taxonomy_source'] or 'unresolved'} "
                         f"({rec['n_ranks']}/7 ranks), "
                         f"{rec['host_group']}\n")
    cache.flush()

    # ---- assemble ---------------------------------------------------------
    host_rows, long_rows = [], []
    for name, h in unique_names.items():
        rec = cache.get(name) or {"lookup_name": name, "rank": h["rank"],
                                  "taxonomy_source": "", "n_ranks": 0,
                                  "env_sources": {}, "host_group": "Unresolved"}
        envs, env_src = unify_env(set(), rec.get("env_sources", {}))
        row = {"host_name": name, "rank": rec.get("rank", ""),
               "n_sources_found": rec.get("n_sources_found", 0),
               "sources_found": rec.get("sources_found", ""),
               "sources_rejected": rec.get("sources_rejected", ""),
               "phylum_agreement": rec.get("phylum_agreement", ""),
               "phylum_dispute": rec.get("phylum_dispute", ""),
               "phylum_tied": rec.get("phylum_tied", ""),
               "order_dispute": rec.get("order_dispute", ""),
               "flag": rec.get("flag", ""),
               "matched_name": rec.get("matched_name", ""),
               "matched_scientific_name": rec.get("matched_scientific_name", ""),
               "match_type": rec.get("match_type", ""),
               "match_confidence": rec.get("match_confidence", ""),
               "taxonomic_status": rec.get("taxonomic_status", ""),
               "name_changed": rec.get("name_changed", ""),
               "taxonomy_source": rec.get("taxonomy_source", ""),
               "ranks_resolved": rec.get("n_ranks", 0),
               "host_group": rec.get("host_group", "Unresolved"),
               "environment": "; ".join(envs), "environment_source": env_src}
        row.update({r: rec.get(r, "") for r in RANKS})
        row["phylum_normalised"] = normalise_phylum(rec.get("phylum", ""),
                                                    rec.get("class", ""),
                                                    rec.get("order", ""))
        for e in ENVIRONMENTS:
            row[f"is_{e}"] = int(e in envs)
        host_rows.append(row)

    host_by_name = {r["host_name"]: r for r in host_rows}
    n_db_agree = n_db_conflict = n_db_only = 0
    for sp, cat, h, db_env in pairs:
        base = host_by_name[h["name"]]
        rec = cache.get(h["name"]) or {}
        host_envs, host_src = unify_env(set(), rec.get("env_sources", {}))
        envs, env_src = ((sorted(db_env), "database")
                         if args.use_db_environment and db_env
                         else (host_envs, host_src))
        # the curated column is compared against the HOST-DERIVED value, not
        # against whatever ended up in `environment`, so the check still means
        # something when --use-db-environment is on
        if db_env:
            if not host_envs:
                n_db_only += 1
            elif set(db_env) & set(host_envs):
                n_db_agree += 1
            else:
                n_db_conflict += 1
        r = {"microsporidia_species": sp, "microsporidia_category": cat,
             "host_raw": h["raw"], "host_name": h["name"], "host_rank": h["rank"],
             "host_role": h["role"], "host_column": h["host_column"],
             "host_group": base["host_group"],
             "secondary_organism": h.get("secondary_organism", ""),
             "hyperparasite_context": h.get("context", ""),
             "cell": h.get("cell", ""),
             "taxonomy_source": base["taxonomy_source"],
             "environment": "; ".join(envs), "environment_source": env_src,
             "host_derived_environment": "; ".join(host_envs),
             "db_host_environment": "; ".join(sorted(db_env)),
             "db_env_agrees": ("" if not db_env else
                               "no_host_env" if not host_envs else
                               "yes" if set(db_env) & set(host_envs) else "CONFLICT")}
        r.update({r_: base[r_] for r_ in RANKS})
        r["phylum_normalised"] = base.get("phylum_normalised", "")
        # The QC verdict travels with EVERY row of the main table, not only in
        # the separate disagreement file. needs_review is the one column to
        # filter or facet on; the rest say why.
        r["n_sources_found"] = base.get("n_sources_found", 0)
        r["phylum_agreement"] = base.get("phylum_agreement", "")
        r["phylum_dispute"] = base.get("phylum_dispute", "")
        r["taxonomy_flag"] = base.get("flag", "")
        reasons = []
        if base.get("phylum_dispute"):
            reasons.append("phylum disagreement")
        if base.get("phylum_tied") == "yes":
            reasons.append("tied vote")
        if base.get("sources_rejected"):
            reasons.append("a source answered about another organism")
        if base.get("name_changed"):
            reasons.append("source matched a different name")
        if "spelling corrected" in str(base.get("flag", "")):
            reasons.append("spelling corrected")
        if "common name" in str(base.get("flag", "")):
            reasons.append("resolved from a common name")
        if base.get("n_sources_found") == 1:
            reasons.append("single source")
        if base.get("n_sources_found") == 0 and nzchar_py(base.get("host_name")):
            reasons.append("unresolved")
        r["needs_review"] = "; ".join(reasons)
        for e in ENVIRONMENTS:
            r[f"is_{e}"] = int(e in envs)
        long_rows.append(r)

    def write_tsv(path, rows, cols=None):
        if not rows:
            return
        cols = cols or list(rows[0].keys())
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                               extrasaction="ignore", lineterminator="\n")
            w.writeheader()
            w.writerows([{k: tsv_safe(v) for k, v in r.items()}
                         for r in rows])

    write_tsv(out("01_species_host_long.tsv"), long_rows)
    write_tsv(out("02_host_lookup.tsv"), host_rows)

    # a fuzzy or near match means the database answered about a different
    # spelling from the one asked for. Usually a corrected typo, occasionally a
    # different organism entirely, so all of them are written out for review.
    fuzzy = [r for r in host_rows
             if r.get("name_changed")
             or str(r.get("match_type", "")).upper() in
                ("FUZZY", "LIKE", "PHONETIC", "NEAR_1", "NEAR_2", "FREE_TEXT",
                 "HIGHERRANK")]
    write_qc_rows = fuzzy
    unresolved = [r for r in host_rows if r["ranks_resolved"] == 0]
    write_tsv(out("03_UNRESOLVED_hosts.tsv"), unresolved)
    no_env = [r for r in host_rows if not r["environment"]]
    write_tsv(out("04_NO_environment.tsv"), no_env)
    write_tsv(out("07_fuzzy_name_matches.tsv"), write_qc_rows)

    disagree = [r for r in host_rows
                if r.get("phylum_dispute") or r.get("order_dispute")]
    write_tsv(out("08_TAXONOMY_DISAGREEMENTS.tsv"), disagree)
    single = [r for r in host_rows if r.get("n_sources_found") == 1]
    write_tsv(out("09_single_source_hosts.tsv"), single)
    wrong_org = [r for r in host_rows if r.get("sources_rejected")]
    write_tsv(out("10_rejected_wrong_organism.tsv"), wrong_org)

    # ---- cells naming more than one organism -------------------------------
    # Two situations, and they need different treatment:
    #   A. the second organism sits INSIDE a parenthetical after "in"/"parasite
    #      in" -- it is the host's host, already excluded from the counts, and
    #      recorded in secondary_organism.
    #   B. two organisms separated by ";" where one is a parasite of the other
    #      ("Enterocystis rhithrogenae (hyperparasitic); Rhithrogena
    #      semicolorata"). Both are counted as hosts per the database's own
    #      convention, but which one the microsporidian actually infects is not
    #      recoverable from the text, so every such cell is flagged.
    # Obligate-parasite phyla only. Nematoda and Ciliophora are deliberately
    # absent: most nematode and ciliate hosts here are free-living (C. elegans,
    # Paramecium), so including them would flag hundreds of ordinary entries and
    # bury the real cases.
    PARASITE_PHYLA = {"Platyhelminthes", "Apicomplexa", "Acanthocephala",
                      "Myxozoa", "Nematomorpha"}
    by_cell = OrderedDict()
    for r in long_rows:
        by_cell.setdefault((r["microsporidia_species"], r.get("host_column", ""),
                            r.get("cell", "")), []).append(r)
    multi = []
    for (sp, hc, cell), rows_ in by_cell.items():
        names = [x["host_name"] for x in rows_]
        phyla = {x.get("phylum", "") for x in rows_ if x.get("phylum")}
        secondaries = [x["secondary_organism"] for x in rows_
                       if x.get("secondary_organism")]
        ctx = [x["hyperparasite_context"] for x in rows_
               if x.get("hyperparasite_context")]
        reasons = []
        if secondaries:
            reasons.append("host's host named in parentheses (excluded)")
        if len(names) > 1 and (phyla & PARASITE_PHYLA) and len(phyla) > 1:
            reasons.append("a parasite phylum AND another phylum listed "
                           "together: which is the microsporidian's host?")
        if len(names) > 1 and any(HYPER_WORD.search(c) for c in ctx):
            reasons.append("hyperparasite wording alongside a second organism")
        if not reasons:
            continue
        multi.append({
            "microsporidia_species": sp, "host_column": hc, "cell": cell,
            "hosts_counted": "; ".join(names),
            "phyla": "; ".join(sorted(phyla)),
            "secondary_organisms_excluded": "; ".join(sorted(set(secondaries))),
            "why_flagged": " | ".join(reasons)})
    write_tsv(out("11_MULTIPLE_ORGANISMS.tsv"), multi)

    # per-microsporidia-species summary, one row per species
    per_sp = OrderedDict()
    for r in long_rows:
        k = r["microsporidia_species"]
        d = per_sp.setdefault(k, {"microsporidia_species": k,
                                  "microsporidia_category": r["microsporidia_category"],
                                  "n_hosts": 0, "host_groups": set(),
                                  "phyla": set(), "environments": set()})
        d["n_hosts"] += 1
        d["host_groups"].add(r["host_group"])
        if r["phylum"]:
            d["phyla"].add(r["phylum"])
        d["environments"] |= {e for e in r["environment"].split("; ") if e}
    summary = []
    for d in per_sp.values():
        summary.append({
            "microsporidia_species": d["microsporidia_species"],
            "microsporidia_category": d["microsporidia_category"],
            "n_hosts": d["n_hosts"],
            "n_host_groups": len(d["host_groups"] - {"Unresolved"}),
            "host_groups": "; ".join(sorted(d["host_groups"])),
            "host_phyla": "; ".join(sorted(d["phyla"])),
            "environments": "; ".join(sorted(d["environments"])) or "unknown",
            "n_environments": len(d["environments"]),
        })
    write_tsv(out("05_per_microsporidia_species.tsv"), summary)
    conflicts = [r for r in long_rows if r["db_env_agrees"] == "CONFLICT"]
    write_tsv(out("06_env_conflicts.tsv"), conflicts)

    # ---- report -----------------------------------------------------------
    print("\n================ RESOLUTION SUMMARY ================")
    print(f"unique host names           : {len(host_rows)}")
    res = sum(1 for r in host_rows if r["ranks_resolved"] > 0)
    print(f"  resolved to a lineage     : {res} ({100*res/max(1,len(host_rows)):.1f}%)")
    print(f"  UNRESOLVED                : {len(unresolved)}")
    print(f"  with an environment       : {len(host_rows)-len(no_env)}")
    print(f"  no environment            : {len(no_env)}")

    if fuzzy:
        changed = [r for r in fuzzy if r.get("name_changed")]
        print(f"\nnames the sources did NOT match exactly : {len(fuzzy)}")
        print(f"  of which a DIFFERENT binomial came back : {len(changed)}")
        print("  (see 07_fuzzy_name_matches.tsv; most are corrected "
              "misspellings, but check)")
        for r in changed[:10]:
            print(f"    {r['host_name']:<34} -> {r['name_changed']} "
                  f"[{r['taxonomy_source']} {r.get('match_type','')}]")

    print("\nhow many of the four sources recognised each host name:")
    nf = Counter(r.get("n_sources_found", 0) for r in host_rows)
    for k in sorted(nf):
        print(f"  {nf[k]:>6}  {k} source(s)"
              + ("   <- highest risk: nothing to cross-check" if k == 1 else "")
              + ("   <- unresolved" if k == 0 else ""))
    print("\nper source, how often it supplied a usable lineage:")
    sc = Counter()
    for r in host_rows:
        for srcname in str(r.get("sources_found", "")).split(";"):
            if srcname:
                sc[srcname] += 1
    for k, v in sc.most_common():
        print(f"  {v:>6}  {k}")

    tied_rows = [r for r in host_rows if r.get("phylum_tied") == "yes"]
    withheld = [r for r in host_rows if str(r.get("flag", "")).startswith(
        ("TIE:", "DISPUTED:"))]
    if tied_rows or withheld:
        print(f"\ntied votes (no majority phylum) : {len(tied_rows)}")
        print(f"lineages WITHHELD from the output: {len(withheld)}"
              f"   [--disputed {args.disputed}, --on-tie {args.on_tie}]")
        if withheld:
            print("  these hosts have no phylum in the tables, so they cannot "
                  "reach a figure until reviewed:")
            for r in withheld[:10]:
                print(f"     {r['host_name']:<32} {r['flag'][:80]}")

    nr = sum(1 for r in long_rows if r.get("needs_review"))
    print(f"\nrows of 01_species_host_long.tsv marked needs_review: {nr} of "
          f"{len(long_rows)} ({100 * nr / max(1, len(long_rows)):.1f}%)")
    rc = Counter()
    for r in long_rows:
        for reason in str(r.get("needs_review", "")).split("; "):
            if reason:
                rc[reason] += 1
    for k, v in rc.most_common():
        print(f"  {v:>6}  {k}")

    if disagree:
        print(f"\n!! {len(disagree)} host name(s) where the sources DISAGREE on "
              f"phylum or order")
        pd_ = [r for r in disagree if r.get("phylum_dispute")]
        print(f"   {len(pd_)} of them disagree at PHYLUM -- check these first "
              f"(08_TAXONOMY_DISAGREEMENTS.tsv)")
        for r in pd_[:12]:
            print(f"     {r['host_name']:<32} chose {r.get('phylum','?'):<16} "
                  f"({r.get('taxonomy_source','')}), others said "
                  f"{r.get('phylum_dispute','')}")
    if wrong_org:
        print(f"\n{len(wrong_org)} host(s) had at least one source answer about a "
              f"DIFFERENT organism (rejected, 10_rejected_wrong_organism.tsv):")
        for r in wrong_org[:8]:
            print(f"     {r['host_name']:<32} {r['sources_rejected'][:90]}")
    if single:
        print(f"\n{len(single)} host(s) were placed by ONE source only, so nothing "
              f"cross-checks them (09_single_source_hosts.tsv)")

    n_sec = sum(1 for r in long_rows if r.get("secondary_organism"))
    if n_sec or multi:
        print(f"\nHYPERPARASITE / MULTI-ORGANISM CELLS")
        print(f"  {n_sec:>6}  host rows where a SECOND organism was named inside "
              f"parentheses")
        print( "          (the host's own host: excluded from the host counts, "
               "kept in secondary_organism)")
        amb = [r for r in multi if "which is the microsporidian" in r["why_flagged"]]
        print(f"  {len(multi):>6}  cells flagged for review "
              f"(11_MULTIPLE_ORGANISMS.tsv)")
        print(f"  {len(amb):>6}  of those list a parasite phylum AND another "
              f"phylum -- genuinely ambiguous")
        for r in multi[:8]:
            print(f"     {r['microsporidia_species']:<34} {r['cell'][:64]}")

    print("\ntaxonomy source:")
    for k, v in Counter(r["taxonomy_source"] or "(none)" for r in host_rows).most_common():
        print(f"  {v:>6}  {k}")
    renamed = [r for r in host_rows
               if r.get("phylum_normalised") and r.get("phylum")
               and r["phylum_normalised"] != r["phylum"]]
    if renamed:
        print(f"\nphylum names normalised across sources: {len(renamed)} host(s)")
        for k, v in Counter(f"{r['phylum']} -> {r['phylum_normalised']}"
                            for r in renamed).most_common():
            print(f"  {v:>6}  {k}")
        print("  use phylum_normalised for figures; the raw value is kept in phylum")

    print("\nhost group (unique host names):")
    for k, v in Counter(r["host_group"] for r in host_rows).most_common():
        print(f"  {v:>6}  {k}")
    print("\nenvironment (unique host names, a host may have several):")
    ec = Counter()
    for r in host_rows:
        for e in r["environment"].split("; "):
            if e:
                ec[e] += 1
    for k, v in ec.most_common():
        print(f"  {v:>6}  {k}")
    if http.calls:
        print("\nAPI calls: " + ", ".join(f"{k} {v}" for k, v in http.calls.items()))
    if http.errors:
        print("API errors: " + ", ".join(f"{k} {v}" for k, v in http.errors.most_common(6)))
    if n_db_agree or n_db_conflict or n_db_only:
        print("\ncuration check -- spreadsheet 'Host Environment' vs host lookup:")
        print(f"  {n_db_agree:>6}  agree")
        print(f"  {n_db_conflict:>6}  CONFLICT (see db_env_agrees in 01_..., "
              f"and 06_env_conflicts.tsv)")
        print(f"  {n_db_only:>6}  spreadsheet has a value, host lookup does not")
        print("  the spreadsheet column " +
              ("OVERRIDES `environment`" if args.use_db_environment else
               "does NOT feed `environment` (--use-db-environment to change)"))

    if unresolved:
        print(f"\nfirst unresolved host names (see 03_UNRESOLVED_hosts.tsv):")
        for r in unresolved[:15]:
            print(f"  {r['host_name']}  [{r['rank']}]")

    print(f"\nwritten to {args.outdir}/")
    print("  01_species_host_long.tsv        one row per microsporidia species x host")
    print("      ^ carries needs_review / taxonomy_flag / phylum_agreement /")
    print("        phylum_dispute / n_sources_found, so a suspect row can be")
    print("        spotted or filtered without opening the QC files")
    print("  02_host_lookup.tsv              one row per unique host name")
    print("  03_UNRESOLVED_hosts.tsv         no lineage from any source")
    print("  04_NO_environment.tsv           lineage but no habitat")
    print("  05_per_microsporidia_species.tsv one row per microsporidia species")
    print("  06_env_conflicts.tsv            spreadsheet env disagrees with the host")
    print("  07_fuzzy_name_matches.tsv       source matched a different spelling")
    print("  08_TAXONOMY_DISAGREEMENTS.tsv   sources disagree on phylum or order")
    print("  09_single_source_hosts.tsv      only one source placed the host")
    print("  10_rejected_wrong_organism.tsv  a source answered about another taxon")
    print("  11_MULTIPLE_ORGANISMS.tsv       cell names >1 organism; host may be ambiguous")


if __name__ == "__main__":
    main()
