#!/usr/bin/env python3
"""
Predict the habitat of every host in the microsporidia database from taxonomic
knowledge, with no API calls.

WHAT THIS IS
  A curated genus -> habitat table plus family/order-level rules, written out
  by hand. It is offline knowledge, not a database lookup: nothing here was
  fetched from GBIF, WoRMS, COL or NCBI. Treat it as a strong prior to be
  checked, not as a source of record. Every assignment carries a confidence
  and a basis so the weak ones can be found and verified.

  It is complementary to host_taxonomy_environment.py, not a replacement:
    * this script covers the long tail of insect and nematode names that the
      habitat databases simply do not hold
    * the API script gives lineages, authority-backed habitats, and evidence
  Run both and compare. The comparison is the useful artefact.

THE LIFE-STAGE CONVENTION (matters for this dataset)
  Mosquitoes, blackflies, chironomids, mayflies and caddisflies are aquatic as
  larvae and terrestrial as adults. Microsporidian infection in these groups is
  overwhelmingly acquired and expressed in the aquatic larval stage, so:
    environment_primary = freshwater
    environment_all     = freshwater; terrestrial
    life_stage_note     = "aquatic larva, terrestrial adult"
  Use environment_primary for a single-habitat figure; use environment_all if
  the figure should show every habitat the host occupies. Anadromous and
  catadromous fish (Salmo, Oncorhynchus, Anguilla, Alosa, Osmerus) are handled
  the same way, primary = the habitat where infection is usually reported.

CONFIDENCE
  high    the genus has one unambiguous habitat (Daphnia, Bombyx, Penaeus)
  medium  the genus is mostly one habitat but has exceptions, or the call rests
          on a family-level rule rather than the genus itself
  low     a guess from the name alone; verify before using

USAGE
  python predict_host_environment.py --db database.xlsx
  python predict_host_environment.py --db database.xlsx --compare api_output.tsv
"""

import argparse
import csv
import os
import re
import sys
import unicodedata
from collections import Counter, OrderedDict

# THE canonical named/provisional rule -- one shared implementation.
# Do not define a local copy; edit classify_species.py instead.
from classify_species import classify_microsporidia

# --------------------------------------------------------------------------- #
# host-name parsing, identical rules to host_taxonomy_environment.py
# --------------------------------------------------------------------------- #
ROLE_WORDS = re.compile(
    r"^(?:type host|intermediate host|definitive host|paratype|natural host"
    r"|experimental|primary host|secondary host|alternate host|vector"
    r"|\d+|[a-z]\d?|type \d|host \d|\?|unknown|likely|possibly|probably"
    r"|should be .*|according to .*|misidentified.*|as .*)$", re.I)
SYNONYM_RE = re.compile(r"^(?:=|also|syn\.?|synonym of)\s*(.+)$", re.I)
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
    n = re.sub(r"\s+", " ", deaccent(name)).strip(" .,;:")
    n = re.sub(r"\b(?:sp|spp)\.?\s*(?:nov\.?)?\s*\d*$", "sp.", n, flags=re.I)
    n = re.sub(r"\b(?:cf|aff|nr)\.?\s+", "", n, flags=re.I)
    n = re.sub(r"\s+(?:sensu lato|s\.l\.|sensu stricto|s\.s\.)$", "", n, flags=re.I)
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
    if cell is None:
        return []
    text = str(cell).strip()
    if not text:
        return []
    out, seen_genera = [], {}
    for frag in split_hosts(text):
        frag = frag.strip()
        if not frag:
            continue
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
        subgenus = ""
        sg = SUBGENUS.match(frag.strip())
        if sg:
            subgenus = sg.group(2)
            frag_use = (sg.group(1) + " " + sg.group(3)).strip()
        else:
            frag_use = frag
        bc = BARE_CONTEXT.search(frag_use)
        if bc:
            if not secondary:
                secondary = bc.group(1)
            context = context or bc.group(0).strip().lower()
            frag_use = frag_use[:bc.start()]
        base = clean_epithet(re.sub(r"\([^)]*\)", " ", frag_use)).strip(' "\'')
        if is_description(base):
            out.append({"raw": frag, "name": "", "rank": "unidentified",
                        "role": role, "synonyms": [], "subgenus": subgenus,
                        "secondary_organism": secondary, "context": context,
                        "repair": ""})
            continue
        if not base or NON_TAXON.match(base):
            continue
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
        out.append({"raw": frag, "name": lookup, "rank": rank, "role": role,
                    "synonyms": synonyms})
    return out


# --------------------------------------------------------------------------- #
# GROUPS. (primary environment, all environments, coarse host group, note)
# --------------------------------------------------------------------------- #
FW = "freshwater"
MA = "marine"
BR = "brackish"
TE = "terrestrial"
AQ_LARVA = "aquatic larva, terrestrial adult"
DIADROM = "migrates between fresh and salt water"

GROUPS = {
    # ---- insects with aquatic larvae ------------------------------------- #
    "mosquito": (FW, [FW, TE], "Insecta", AQ_LARVA, "high"),
    "blackfly": (FW, [FW, TE], "Insecta", AQ_LARVA, "high"),
    "chironomid": (FW, [FW, TE], "Insecta", AQ_LARVA, "high"),
    "mayfly": (FW, [FW, TE], "Insecta", AQ_LARVA, "high"),
    "caddisfly": (FW, [FW, TE], "Insecta", AQ_LARVA, "high"),
    "odonate": (FW, [FW, TE], "Insecta", AQ_LARVA, "high"),
    "aquatic_bug": (FW, [FW], "Insecta", "aquatic throughout", "high"),
    "biting_midge": (FW, [FW, TE], "Insecta", AQ_LARVA, "medium"),
    # ---- fully terrestrial insects --------------------------------------- #
    "lepidoptera": (TE, [TE], "Insecta", "", "high"),
    "hymenoptera": (TE, [TE], "Insecta", "", "high"),
    "coleoptera": (TE, [TE], "Insecta", "", "high"),
    "orthoptera": (TE, [TE], "Insecta", "", "high"),
    "diptera_terr": (TE, [TE], "Insecta", "", "high"),
    "hemiptera": (TE, [TE], "Insecta", "", "high"),
    "other_insect": (TE, [TE], "Insecta", "", "medium"),
    # ---- crustaceans ------------------------------------------------------ #
    "cladoceran": (FW, [FW], "Crustacea", "", "high"),
    "copepod_fw": (FW, [FW], "Crustacea", "", "high"),
    "copepod_marine": (MA, [MA], "Crustacea", "", "high"),
    "amphipod_fw": (FW, [FW], "Crustacea", "", "high"),
    "amphipod_baikal": (FW, [FW], "Crustacea", "Lake Baikal endemic", "high"),
    "amphipod_marine": (MA, [MA, BR], "Crustacea", "", "high"),
    "amphipod_brackish": (BR, [BR, FW], "Crustacea", "", "medium"),
    "decapod_marine": (MA, [MA, BR], "Crustacea", "", "high"),
    "decapod_fw": (FW, [FW], "Crustacea", "", "high"),
    "decapod_brackish": (BR, [BR, MA, FW], "Crustacea", "", "medium"),
    "isopod_fw": (FW, [FW], "Crustacea", "", "high"),
    "isopod_terr": (TE, [TE], "Crustacea", "", "high"),
    "brine_shrimp": (BR, [BR, MA], "Crustacea", "hypersaline", "high"),
    "other_crustacean_fw": (FW, [FW], "Crustacea", "", "medium"),
    "other_crustacean_marine": (MA, [MA], "Crustacea", "", "medium"),
    # ---- fish -------------------------------------------------------------#
    "fish_marine": (MA, [MA], "Fish", "", "high"),
    "fish_fw": (FW, [FW], "Fish", "", "high"),
    "fish_diadromous": (FW, [FW, MA], "Fish", DIADROM, "high"),
    "fish_brackish": (BR, [BR, MA, FW], "Fish", "", "medium"),
    # ---- other animals ---------------------------------------------------- #
    "mammal": (TE, [TE], "Mammals", "", "high"),
    "bird": (TE, [TE], "Birds", "", "high"),
    "reptile": (TE, [TE], "Reptiles & amphibians", "", "high"),
    "amphibian": (FW, [FW, TE], "Reptiles & amphibians", "aquatic larva", "high"),
    "nematode_terr": (TE, [TE], "Nematoda", "soil", "high"),
    "nematode_marine": (MA, [MA], "Nematoda", "", "medium"),
    "annelid_fw": (FW, [FW], "Annelida", "", "high"),
    "annelid_terr": (TE, [TE], "Annelida", "", "high"),
    "annelid_marine": (MA, [MA], "Annelida", "", "high"),
    "mollusc_marine": (MA, [MA], "Mollusca", "", "high"),
    "mollusc_fw": (FW, [FW], "Mollusca", "", "high"),
    "mollusc_terr": (TE, [TE], "Mollusca", "", "high"),
    "bryozoan_fw": (FW, [FW], "Rotifera & bryozoa", "", "high"),
    "rotifer_fw": (FW, [FW], "Rotifera & bryozoa", "", "high"),
    "mite_terr": (TE, [TE], "Arachnida & myriapods", "", "high"),
    "collembola": (TE, [TE], "Insecta", "soil", "high"),
    "gregarine_marine": (MA, [MA], "Protists", "parasite of marine invertebrates",
                         "medium"),
    "protist_fw": (FW, [FW], "Protists", "", "medium"),
    "protist_marine": (MA, [MA], "Protists", "", "medium"),
    "helminth_fw": (FW, [FW], "Platyhelminthes", "freshwater life cycle", "medium"),
    "helminth_terr": (TE, [TE], "Platyhelminthes", "terrestrial life cycle",
                      "medium"),
    "nemertean_marine": (MA, [MA], "Other invertebrates", "", "medium"),
    "cnidarian_marine": (MA, [MA], "Cnidaria & sponges", "", "high"),
    "echinoderm_marine": (MA, [MA], "Echinodermata", "", "high"),
    "tunicate_marine": (MA, [MA], "Other invertebrates", "", "high"),
}

# --------------------------------------------------------------------------- #
# GENUS -> group. Hand-written. Ordered by how often each appears in the sheet.
# --------------------------------------------------------------------------- #
GENUS = {}


def add(group, *genera):
    for g in genera:
        GENUS[g.lower()] = group


add("mosquito", "Aedes", "Culex", "Anopheles", "Ochlerotatus", "Culiseta",
    "Mansonia", "Psorophora", "Aedeomyia", "Orthopodomyia", "Uranotaenia",
    "Coquillettidia", "Toxorhynchites", "Wyeomyia", "Armigeres", "Tripteroides",
    "Haemagogus", "Sabethes", "Deinocerites", "Culicella", "Theobaldia",
    "Ochleratatus", "Finlaya", "Verrallina", "Aedimorphus")
add("blackfly", "Simulium", "Odagmia", "Wilhelmia", "Wilhemia", "Prosimulium",
    "Eusimulium", "Gigantodax", "Tetisimulium", "Cnephia", "Byssodon",
    "Metacnephia", "Boophthora", "Schoenbaueria", "Stegopterna", "Twinnia",
    "Greniera", "Helodon", "Austrosimulium", "Nevermannia", "Hellichiella")
add("chironomid", "Chironomus", "Camptochironomus", "Endochironomus",
    "Kiefferulus", "Tanypus", "Orthocladius", "Cricotopus", "Corynoneura",
    "Microtendipes", "Polypedilum", "Glyptotendipes", "Procladius",
    "Cryptochironomus", "Einfeldia", "Limnochironomus", "Paratanytarsus",
    "Tanytarsus", "Chaoborus", "Culicoides", "Ceratopogon", "Dixa", "Psychoda",
    "Sciomyza")
add("mayfly", "Ephemera", "Cloeon", "Baetis", "Rhithrogena", "Ephemerella",
    "Centroptilum", "Caenis", "Leptophlebia", "Habrophlebia", "Siphlonurus",
    "Ecdyonurus", "Heptagenia", "Paraleptophlebia")
add("caddisfly", "Polycentropus", "Hydropsyche", "Limnephilus", "Rhyacophila",
    "Phryganea", "Trichostegia", "Agapetus", "Anabolia", "Halesus",
    "Potamophylax", "Chaetopteryx", "Glyphotaelius", "Molanna", "Micrasema")
add("odonate", "Aeschna", "Aeshna", "Libellula", "Coenagrion", "Ischnura",
    "Sympetrum", "Anax")
add("aquatic_bug", "Nepa", "Notonecta", "Corixa", "Gerris", "Ranatra", "Naucoris")

add("lepidoptera", "Bombyx", "Pieris", "Malacosoma", "Antheraea", "Heliothis",
    "Spodoptera", "Lymantria", "Galleria", "Plodia", "Ephestia", "Agrotis",
    "Choristoneura", "Euproctis", "Hyphantria", "Colias", "Ostrinia", "Chilo",
    "Crambus", "Trichoplusia", "Carpocapsa", "Cydia", "Paramyelois",
    "Pseudaletia", "Vanessa", "Operophtera", "Tortrix", "Cacoecia",
    "Thaumetopoea", "Eriogaster", "Aporia", "Sesia", "Nonagria", "Gnorimoschema",
    "Estigmene", "Phryganidia", "Stilpnotia", "Leucoma", "Eilema", "Porthetria",
    "Nygmia", "Helicoverpa", "Mamestra", "Plutella", "Laphygma", "Loxostege",
    "Argyresthia", "Cheimatobia", "Lithosa", "Lithosia", "Callimorpha", "Sphinx",
    "Hyalophora", "Dione", "Argyrogramma", "Manduca", "Danaus", "Papilio",
    "Maruca", "Diatraea", "Anagasta", "Cadra", "Corcyra", "Sitotroga",
    "Phthorimaea", "Tuta", "Earias", "Achroia", "Hypera", "Pectinophora",
    "Prodenia", "Barathra", "Panolis", "Bupalus", "Hyponomeuta", "Yponomeuta",
    "Zeiraphera", "Dendrolimus", "Orgyia", "Pandemis", "Adoxophyes",
    "Homona", "Ectomyelois", "Amyelois", "Pyrausta", "Crocidolomia",
    "Lobesia", "Grapholita", "Epiphyas", "Chrysodeixis", "Autographa",
    "Agraulis", "Junonia", "Heliconius", "Aglais", "Melitaea", "Euphydryas",
    "Boloria", "Actias", "Samia", "Philosamia", "Attacus")
add("hymenoptera", "Apis", "Bombus", "Solenopsis", "Neodiprion", "Pristiphora",
    "Macrocentrus", "Vespula", "Formica", "Nosema", "Melipona", "Osmia",
    "Megachile", "Lasius", "Myrmica", "Camponotus", "Monomorium", "Linepithema",
    "Atta", "Acromyrmex", "Cotesia", "Aphidius", "Trichogramma", "Diprion",
    "Athalia", "Cephus", "Sirex", "Polistes", "Vespa", "Dolichovespula",
    "Anoplius", "Wasp")
add("coleoptera", "Scolytus", "Ips", "Leptinotarsa", "Tribolium", "Oryzaephilus",
    "Dinoderus", "Anthonomus", "Otiorrhynchus", "Otiorhynchus", "Listronotus",
    "Phyllotreta", "Nisotra", "Geotrupes", "Hylobius", "Blaps", "Dermestes",
    "Sitophilus", "Tenebrio", "Trogoderma", "Callosobruchus", "Popillia",
    "Melolontha", "Phyllophaga", "Diabrotica", "Meligethes", "Dendroctonus",
    "Tomicus", "Pityogenes", "Hypothenemus", "Xyleborus", "Cryptolestes",
    "Rhyzopertha", "Lasioderma", "Stegobium", "Alphitobius", "Coccinella",
    "Harmonia", "Adalia", "Carabus", "Pterostichus", "Curculio", "Sitona",
    "Ceutorhynchus", "Bruchus", "Acanthoscelides", "Anomala", "Costelytra",
    "Heteronychus", "Oryctes", "Agriotes", "Limonius", "Cylas")
add("orthoptera", "Melanoplus", "Locusta", "Schistocerca", "Anabrus",
    "Chorthippus", "Acrida", "Gryllus", "Locustana", "Metator", "Phoetatiotes",
    "Camnula", "Aulocara", "Dissosteira", "Oedaleus", "Calliptamus",
    "Dociostaurus", "Nomadacris", "Zonocerus", "Gryllidae", "Acheta",
    "Gryllotalpa", "Conocephalus", "Tettigonia", "Chortoicetes", "Austroicetes")
add("diptera_terr", "Drosophila", "Calliphora", "Musca", "Mucsa", "Muscina",
    "Delia", "Sciara", "Rhynchosciara", "Bibio", "Tabanus", "Hybomitra",
    "Glossina", "Lucilia", "Sarcophaga", "Stomoxys", "Ceratitis", "Bactrocera",
    "Rhagoletis", "Phormia", "Chrysomya", "Hylemya", "Pegomya", "Liriomyza",
    "Contarinia", "Mayetiola", "Lycoriella", "Bradysia", "Megaselia",
    "Phlebotomus", "Lutzomyia", "Fannia", "Hydrotaea", "Haematobia",
    "Dasyneura", "Sitodiplosis", "Chrysops", "Atylotus")
add("hemiptera", "Euschistus", "Halyomorpha", "Nezara", "Lygus", "Aphis",
    "Myzus", "Bemisia", "Trialeurodes", "Nilaparvata", "Sogatella",
    "Dysdercus", "Oncopeltus", "Chinavia", "Podisus", "Perillus", "Adelphocoris",
    "Rhodnius", "Triatoma", "Cimex", "Pediculus", "Pthirus", "Blattella",
    "Periplaneta", "Blatta", "Supella", "Xanthocaecilius", "Polypsocus",
    "Psocus", "Liposcelis", "Thrips", "Frankliniella")
add("other_insect", "Reticulitermes", "Coptotermes", "Termes", "Kalotermes",
    "Zootermopsis", "Nasutitermes", "Odontotermes", "Microtermes",
    "Forficula", "Ctenocephalides", "Xenopsylla", "Nosopsyllus", "Pulex",
    "Chrysoperla", "Chrysopa", "Myrmeleon", "Panorpa", "Sialis")
add("collembola", "Lepidocyrtus", "Tomocerus", "Folsomia", "Orchesella",
    "Isotoma", "Sinella", "Hypogastrura", "Onychiurus", "Sminthurus",
    "Entomobrya", "Protaphorura")

add("cladoceran", "Daphnia", "Ceriodaphnia", "Simocephalus", "Moina", "Bosmina",
    "Chydorus", "Diaphanosoma", "Scapholeberis", "Alona", "Sida", "Polyphemus",
    "Leptodora", "Limnetis", "Lepidurus", "Branchinecta", "Streptocephalus",
    "Eubranchipus", "Triops", "Macrothrix", "Ilyocryptus", "Acroperus",
    "Pleuroxus", "Graptoleberis", "Eurycercus", "Bythotrephes", "Cercopagis")
add("copepod_fw", "Cyclops", "Macrocyclops", "Acanthocyclops", "Megacyclops",
    "Eucyclops", "Mesocyclops", "Thermocyclops", "Diaptomus", "Apocyclops",
    "Paracyclops", "Metacyclops", "Microcyclops", "Tropocyclops",
    "Ectocyclops", "Diacyclops", "Cyclopina", "Eudiaptomus", "Acanthodiaptomus",
    "Arctodiaptomus", "Heliodiaptomus", "Neodiaptomus", "Attheyella",
    "Canthocamptus", "Bryocamptus", "Nitokra", "Candona", "Cypris",
    "Cyclocypris", "Notodromas", "Ilyocypris")
add("copepod_marine", "Lepeophtheirus", "Caligus", "Calanus", "Acartia",
    "Temora", "Centropages", "Paracalanus", "Oithona", "Pseudocalanus",
    "Euterpina", "Tisbe", "Tigriopus", "Harpacticus", "Microsetella",
    "Corycaeus", "Oncaea", "Metridia", "Euchaeta", "Labidocera")
add("brine_shrimp", "Artemia")

add("amphipod_baikal", "Brandtia", "Eulimnogammarus", "Gmelinoides",
    "Micruropus", "Acanthogammarus", "Crypturopus", "Pallasea", "Pallaseopsis",
    "Linevichella", "Dorogostaiskia", "Asprogammarus", "Odontogammarus",
    "Macrohectopus", "Poekilogammarus", "Hyalellopsis", "Ommatogammarus",
    "Garjajewia", "Abyssogammarus", "Parapallasea", "Boeckaxelia",
    "Carinogammarus", "Eucarinogammarus", "Spinacanthus", "Brachyuropus",
    "Plesiogammarus", "Homocerisca", "Cryptoropus", "Axelboeckia")
add("amphipod_fw", "Gammarus", "Crangonyx", "Niphargus", "Synurella",
    "Hyalella", "Paracalliope", "Diporeia", "Pontoporeia", "Rivulogammarus",
    "Gammarellus", "Stygobromus", "Bactrurus", "Chaetogammarus")
add("amphipod_marine", "Orchestia", "Melita", "Corophium", "Talitrus",
    "Ampelisca", "Caprella", "Jassa", "Gammaropsis", "Hyale", "Amphithoe",
    "Talorchestia", "Platorchestia", "Marinogammarus", "Elasmopus",
    "Leptocheirus", "Monoporeia", "Themisto", "Parathemisto", "Hyperia")
add("amphipod_brackish", "Dikerogammarus", "Pontogammarus", "Echinogammarus",
    "Obesogammarus", "Chelicorophium", "Gmelina", "Niphargoides")

add("decapod_marine", "Penaeus", "Litopenaeus", "Metapenaeus",
    "Farfantepenaeus", "Fenneropenaeus", "Marsupenaeus", "Panulirus",
    "Homarus", "Cancer", "Carcinus", "Callinectes", "Crangon", "Pandalus",
    "Palaemon", "Nephrops", "Maja", "Hyas", "Chionoecetes", "Paralithodes",
    "Pagurus", "Eriphia", "Necora", "Liocarcinus", "Portunus", "Scylla",
    "Sicyonia", "Trachypenaeus", "Xiphopenaeus", "Pleoticus", "Solenocera",
    "Nematocarcinus", "Munida", "Galathea", "Upogebia", "Callianassa",
    "Emerita", "Uca", "Sesarma", "Panopeus", "Rhithropanopeus", "Hemigrapsus",
    "Eriocheir", "Charybdis", "Chaceon", "Lithodes", "Sclerocrangon")
add("decapod_fw", "Cherax", "Procambarus", "Faxonius", "Orconectes", "Astacus",
    "Pacifastacus", "Austropotamobius", "Caridina", "Neocaridina",
    "Atyephira", "Atya", "Paratya", "Cambarus", "Cambaroides", "Parastacus",
    "Euastacus", "Aegla", "Potamon", "Dilocarcinus")
add("decapod_brackish", "Palaemonetes", "Macrobrachium", "Exopalaemon",
    "Leander", "Nika", "Athanas", "Alpheus")
add("isopod_fw", "Asellus", "Proasellus", "Caecidotea", "Lirceus")
add("isopod_terr", "Armadillidium", "Porcellio", "Oniscus", "Trichoniscus",
    "Philoscia", "Porcellionides", "Cylisticus")
add("other_crustacean_marine", "Balanus", "Semibalanus", "Chthamalus",
    "Sacculina", "Idotea", "Sphaeroma", "Limnoria", "Gnathia", "Cyathura",
    "Nebalia", "Euphausia", "Meganyctiphanes", "Mysis", "Neomysis",
    "Praunus", "Leptomysis", "Siriella")
add("other_crustacean_fw", "Branchiura", "Argulus", "Lernaea", "Ergasilus")

add("fish_marine", "Gadus", "Merluccius", "Lophius", "Sparus", "Dentex",
    "Sardinella", "Theragra", "Melanogrammus", "Pleuronectes", "Glyptocephalus",
    "Limanda", "Solea", "Synaptura", "Lutjanus", "Nemipterus", "Pagrus",
    "Seriola", "Thunnus", "Epinephelus", "Saurida", "Callionymus",
    "Ctenolabrus", "Cyclopterus", "Hippocampus", "Clupea", "Sardina",
    "Engraulis", "Scomber", "Trachurus", "Micromesistius", "Pollachius",
    "Molva", "Trisopterus", "Hippoglossus", "Reinhardtius", "Platichthys",
    "Scophthalmus", "Psetta", "Dicentrarchus", "Diplodus", "Boops",
    "Mullus", "Trigla", "Chelidonichthys", "Zeus", "Conger", "Raja",
    "Squalus", "Scyliorhinus", "Sebastes", "Anarhichas", "Cyclothone",
    "Myctophum", "Benthosema", "Maurolicus", "Argentina", "Mallotus",
    "Ammodytes", "Zoarces", "Lycodes", "Gymnocanthus", "Myoxocephalus",
    "Hexagrammos", "Pleurogrammus", "Paralichthys", "Verasper", "Kareius",
    "Acanthopagrus", "Rachycentron", "Trachinotus", "Caranx", "Decapterus",
    "Sillago", "Gymnothorax", "Muraena", "Mugil", "Liza", "Chelon")
add("fish_fw", "Hyphessobrycon", "Hemigrammus", "Brachydanio", "Danio",
    "Xiphophorus", "Lepomis", "Esox", "Silurus", "Abramis", "Leuciscus",
    "Lucioperca", "Sander", "Acerina", "Gymnocephalus", "Cottus", "Lota",
    "Culter", "Oreochromis", "Tilapia", "Carassius", "Cyprinus", "Rutilus",
    "Perca", "Barbus", "Tinca", "Scardinius", "Blicca", "Alburnus",
    "Chondrostoma", "Phoxinus", "Gobio", "Aspius", "Vimba", "Pelecus",
    "Ictalurus", "Ameiurus", "Micropterus", "Pomoxis", "Ambloplites",
    "Notropis", "Pimephales", "Catostomus", "Moxostoma", "Umbra",
    "Thymallus", "Coregonus", "Stizostedion", "Channa", "Clarias",
    "Heteropneustes", "Labeo", "Catla", "Cirrhinus", "Poecilia",
    "Micropoecilia", "Gambusia", "Xenotoca", "Astyanax", "Colossoma",
    "Piaractus", "Pygocentrus", "Serrasalmus", "Corydoras", "Pterophyllum",
    "Symphysodon", "Betta", "Trichogaster", "Macropodus", "Paracheirodon",
    "Nothobranchius", "Aphyosemion", "Brachyhypopomus", "Electrophorus",
    "Gymnotus", "Misgurnus", "Cobitis", "Nemacheilus", "Rhodeus",
    "Pseudorasbora", "Hypophthalmichthys", "Ctenopharyngodon", "Megalobrama",
    "Parabramis", "Elopichthys", "Erythroculter", "Hemiculter", "Xenocypris")
add("fish_diadromous", "Oncorhynchus", "Salmo", "Salvelinus", "Plecoglossus",
    "Osmerus", "Anguilla", "Alosa", "Acipenser", "Huso", "Petromyzon",
    "Lampetra", "Coregonus", "Salvethymus", "Hucho", "Parahucho")
add("fish_brackish", "Gasterosteus", "Pungitius", "Neogobius", "Fundulus",
    "Trypauchen", "Boleophthalmus", "Colomesus", "Syngnathus", "Pomatoschistus",
    "Knipowitschia", "Proterorhinus", "Benthophilus", "Clupeonella",
    "Atherina", "Menidia", "Cyprinodon", "Aphanius", "Periophthalmus")

add("mammal", "Homo", "Oryctolagus", "Canis", "Vulpes", "Saimiri", "Ondatra",
    "Mus", "Rattus", "Macaca", "Sus", "Bos", "Ovis", "Capra", "Equus",
    "Felis", "Mesocricetus", "Cricetulus", "Cavia", "Chinchilla", "Meriones",
    "Apodemus", "Microtus", "Myodes", "Clethrionomys", "Sorex", "Talpa",
    "Erinaceus", "Nyctereutes", "Mustela", "Neovison", "Martes", "Meles",
    "Procyon", "Cervus", "Capreolus", "Rangifer", "Alces", "Odocoileus",
    "Lepus", "Sciurus", "Marmota", "Spermophilus", "Callithrix", "Papio",
    "Pan", "Gorilla", "Pongo", "Lemur", "Eptesicus", "Myotis", "Pipistrellus")
add("bird", "Gallus", "Melopsittacus", "Eclectus", "Chalcopsitta", "Agapornis",
    "Anas", "Anser", "Meleagris", "Columba", "Passer", "Sturnus", "Corvus",
    "Struthio", "Coturnix", "Psittacus", "Amazona", "Ara", "Nymphicus",
    "Carduelis", "Serinus", "Taeniopygia", "Falco", "Buteo", "Accipiter",
    "Larus", "Sterna", "Phalacrocorax", "Ciconia", "Grus", "Fulica")
add("reptile", "Pogona", "Lacerta", "Tropidonotus", "Natrix", "Podarcis",
    "Anolis", "Chamaeleo", "Python", "Boa", "Testudo", "Trachemys",
    "Chelonia", "Crocodylus", "Alligator", "Gekko", "Eublepharis")
add("amphibian", "Bufo", "Rana", "Xenopus", "Ambystoma", "Triturus", "Lissotriton",
    "Pelophylax", "Hyla", "Salamandra", "Notophthalmus", "Lithobates")

add("nematode_terr", "Caenorhabditis", "Oscheius", "Panagrellus", "Neoaplectana",
    "Steinernema", "Heterorhabditis", "Pristionchus", "Rhabditis",
    "Acrobeloides", "Plectus", "Aphelenchus", "Ditylenchus", "Meloidogyne",
    "Heterodera", "Globodera", "Pratylenchus", "Bursaphelenchus",
    "Panagrolaimus", "Diploscapter", "Mesorhabditis", "Poikilolaimus",
    "Turbatrix", "Aphelenchoides", "Trichinella", "Ascaris", "Toxocara",
    "Strongyloides", "Haemonchus", "Trichuris", "Enterobius", "Necator")
add("nematode_marine", "Sabatieria", "Monhystera", "Enoplus", "Desmodora",
    "Chromadora", "Anisakis", "Pseudoterranova", "Contracaecum")

add("annelid_fw", "Limnodrilus", "Tubifex", "Lumbriculus", "Chaetogaster",
    "Stylaria", "Nais", "Dero", "Aulodrilus", "Branchiura", "Erpobdella",
    "Helobdella", "Piscicola", "Hirudo", "Glossiphonia")
add("annelid_terr", "Pheretima", "Lumbricus", "Eisenia", "Aporrectodea",
    "Allolobophora", "Dendrobaena", "Octolasion", "Fridericia", "Enchytraeus",
    "Amynthas", "Metaphire")
add("annelid_marine", "Pygospio", "Polydora", "Capitella", "Capitellides",
    "Spio", "Nereis", "Hediste", "Arenicola", "Lanice", "Sabella",
    "Serpula", "Pomatoceros", "Potamoceros", "Owenia", "Scoloplos",
    "Travisia", "Ophelia", "Glycera", "Nephtys", "Harmothoe", "Lepidonotus",
    "Terebella", "Amphitrite", "Chaetopterus", "Myxicola", "Protula")

add("mollusc_marine", "Mytilus", "Tellina", "Crassostrea", "Ostrea", "Ruditapes",
    "Venerupis", "Cerastoderma", "Macoma", "Mya", "Ensis", "Pecten",
    "Placopecten", "Argopecten", "Haliotis", "Littorina", "Nucella",
    "Buccinum", "Nassarius", "Hydrobia", "Patella", "Mercenaria",
    "Aplysia", "Loligo", "Sepia", "Octopus", "Illex", "Todarodes")
add("mollusc_fw", "Lymnaea", "Physa", "Biomphalaria", "Bulinus", "Planorbis",
    "Radix", "Stagnicola", "Anodonta", "Unio", "Dreissena", "Corbicula",
    "Pisidium", "Sphaerium", "Viviparus", "Bithynia", "Valvata", "Ancylus",
    "Potamopyrgus", "Melanoides")
add("mollusc_terr", "Deroceras", "Arion", "Helix", "Cornu", "Cepaea",
    "Succinea", "Limax", "Achatina", "Theba", "Oxychilus", "Discus")

add("bryozoan_fw", "Plumatella", "Fredericella", "Cristatella", "Pectinatella",
    "Lophopus", "Hyalinella")
add("rotifer_fw", "Brachionus", "Asplanchna", "Keratella", "Philodina",
    "Rotaria", "Euchlanis", "Lecane", "Epiphanes")
add("mite_terr", "Phytoseiulus", "Tyrophagus", "Damaeus", "Phtiracarus",
    "Varroa", "Falculifer", "Acarus", "Dermatophagoides", "Tetranychus",
    "Panonychus", "Oribatula", "Steganacarus", "Hypoaspis", "Sarcoptes",
    "Psoroptes", "Ixodes", "Dermacentor", "Rhipicephalus", "Amblyomma",
    "Araneus", "Pardosa", "Lycosa", "Tegenaria", "Linyphia", "Erigone",
    "Lithobius", "Julus", "Polydesmus", "Glomeris")

add("gregarine_marine", "Lecudina", "Selenidium", "Polyrhabdina", "Ancora",
    "Gonospora", "Urospora", "Cephaloidophora", "Nematopsis", "Porospora",
    "Pterospora", "Difficilina", "Lankesteria", "Ganymedes", "Ophioidina")
add("protist_fw", "Thecamoeba", "Amoeba", "Acanthamoeba", "Naegleria",
    "Paramecium", "Tetrahymena", "Euglena", "Chlamydomonas", "Stylonychia",
    "Euplotes", "Vorticella", "Ciliatosporidium", "Platyophrya", "Colpoda",
    "Spirostomum", "Blepharisma", "Actinophrys", "Difflugia", "Arcella")
add("protist_marine", "Vannella", "Marteilia", "Bodo", "Cafeteria",
    "Oxyrrhis", "Noctiluca", "Gymnodinium", "Alexandrium", "Thraustochytrium",
    "Labyrinthula", "Haplosporidium", "Bonamia", "Perkinsus", "Minchinia")

add("helminth_fw", "Echinostoma", "Echinoparyphium", "Cotylurus", "Cercaria",
    "Telorchis", "Diplostomum", "Schistosoma", "Fasciola", "Opisthorchis",
    "Clonorchis", "Ligula", "Diphyllobothrium", "Proteocephalus",
    "Caryophyllaeus", "Bothriocephalus", "Acanthocephalus",
    "Pomphorhynchus", "Neoechinorhynchus", "Acanthocephaloides",
    "Diplacanthus", "Echinorhynchus")
add("helminth_terr", "Hymenolepis", "Taenia", "Echinococcus", "Moniezia",
    "Dicrocoelium", "Dugesia", "Schmidtea", "Girardia", "Polycelis",
    "Planaria", "Macrostomum")
add("nemertean_marine", "Maculaura", "Lineus", "Cerebratulus", "Tetrastemma",
    "Amphiporus", "Prostoma")
add("cnidarian_marine", "Aurelia", "Hydra", "Obelia", "Metridium", "Actinia",
    "Anemonia", "Halichondria", "Suberites", "Haliclona", "Ephydatia",
    "Spongilla", "Sycon", "Leucosolenia")
add("echinoderm_marine", "Asterias", "Echinus", "Paracentrotus",
    "Strongylocentrotus", "Ophiura", "Amphiura", "Holothuria", "Cucumaria",
    "Psolus", "Antedon")
add("tunicate_marine", "Ciona", "Botryllus", "Molgula", "Styela", "Oikopleura",
    "Salpa")


# --- second pass, added after reviewing the unmatched list ----------------- #
add("lepidoptera", "Phalera", "Cactoblastis", "Archips", "Danais", "Tyria",
    "Orthosia", "Ceratomia", "Plusia", "Diacrisia", "Cnaphalocrocis",
    "Laspeyresia", "Chloridea", "Oncopera", "Anisota", "Alsophila",
    "Symmerista", "Notodonta", "Dioryctria", "Trachea", "Clysia", "Petrova",
    "Hyphantia", "Loxagrotis", "Peridroma", "Spilosoma", "Arctia", "Zeuzera",
    "Cossus", "Acronicta", "Euxoa", "Scotia", "Heliothrips")
add("coleoptera", "Hippodamia", "Epilachna", "Pissodes", "Chaetocnema",
    "Gastrophysa", "Polygramma", "Trox", "Xyloterus", "Pityokteines",
    "Agrilus", "Crepidodera", "Pyrrhalta", "Xanthogaleruca", "Anisoplia",
    "Gracilia", "Lasioderma", "Coccinula", "Phaedon", "Psylliodes",
    "Aphthona", "Longitarsus", "Altica", "Oulema", "Lema")
add("hymenoptera", "Trichiocampus", "Campoletis", "Cremastus", "Perisierola",
    "Bracon", "Andrena", "Hemichroa", "Dahlbominus", "Pikonema", "Arge",
    "Nylanderia", "Habrobracon", "Apanteles", "Exeristes", "Itoplectis",
    "Pteromalus", "Nasonia", "Melittobia")
add("orthoptera", "Pyrgomorpha", "Romalea", "Poecilimon", "Acrididae",
    "Tetrigidae", "Opeia", "Amphitornus", "Ageneotettix", "Phlibostroma",
    "Arphia", "Encoptolophus", "Hesperotettix", "Hypochlora", "Cordillacris",
    "Drepanopterna", "Spharagemon", "Trachyrhachys", "Covasacris", "Tristira",
    "Obuchivia", "Aeropedellus", "Boopedon", "Eritettix", "Mermiria")
add("diptera_terr", "Dacus", "Homalomyia", "Anastrepha", "Scatophaga",
    "Protophormia", "Cochliomyia")
add("other_insect", "Leucophaea", "Blatella", "Macrotermes", "Uncitermes",
    "Hylobittacus", "Chloroperla", "Nemoura", "Perla", "Isoperla", "Leuctra")
add("mosquito", "Manosonia", "Psorophopa", "Mochlonyx", "Corethrella")
add("blackfly", "Chelocnetha", "Cnetha", "Odagnia", "Sulcicnephia")
add("chironomid", "Ablabesmyia", "Arctopelopia", "Psectrotanipus",
    "Psectrotanypus", "Ptychoptera", "Ptyahoptera", "Diapteris", "Helodes",
    "Tipula", "Dicranota", "Pericoma")
add("caddisfly", "Plectrocnemia", "Holocentropus", "Sericostoma", "Micropterna",
    "Anabolin", "Athripsodes", "Mystacides", "Oecetis", "Brachycentrus")
add("mayfly", "Ameletus", "Emphemera", "Cinygma", "Isonychia")
add("odonate", "Tramea", "Calopteryx", "Lestes", "Enallagma")
add("aquatic_bug", "Sigara", "Velia", "Aphelocheirus", "Ilyocoris")

add("cladoceran", "Ceriodaphniae", "Simocephaus", "Holopedium", "Latona",
    "Drepanothrix", "Streblocerus")
add("copepod_fw", "Cyclopoida", "Neocyclops", "Metacyclopina")
add("decapod_marine", "Eupagurus", "Metacarcinus", "Xiphopeneus", "Acetes",
    "Metanephrops", "Sicyonella", "Parapenaeus", "Aristeus", "Plesionika")
add("decapod_fw", "Cambarellus", "Paranephrops", "Cherax")
add("amphipod_fw", "Pontogamtnarus", "Obesogammarus", "Anisogammarus")
add("copepod_marine", "Lepeoptheirus")

add("fish_marine", "Opsanus", "Anarrhichas", "Dorosoma", "Macrourus",
    "Julis", "Carangoides", "Cephalopholis", "Serranus", "Pseudopleuronectes",
    "Parophrys", "Microgadus", "Ophiodon", "Anoplopoma", "Trichiurus",
    "Hyperoplus", "Cepola", "Motella", "Symphodus", "Taurulus", "Crenilabrus",
    "Crenilabpus", "Leiostomus", "Chloroscombrus", "Cypselurus", "Entelurus",
    "Myrophis", "Cymatogaster", "Gobius", "Clevelandia", "Vincentia",
    "Vanstraelenia", "Dicologolossa", "Trisopterus", "Zeugopterus")
add("fish_fw", "Hasemania", "Puntius", "Apistogramma", "Oryzias", "Siluris",
    "Mylopharyngodon", "Lebistes", "Platypoecilus", "Molliensia", "Colisa",
    "Hemichromis", "Pseudocrenilabrus", "Pangasius", "Gymnorhamphichthys",
    "Girardinus", "Pimphales", "Savelinus", "Prosopium", "Apeltes",
    "Eucyclogobius", "Planiliza", "Hypseleotris", "Melanotaenia")
add("fish_diadromous", "Coregonus", "Stenodus")

add("mite_terr", "Rhysotritia", "Carabodes", "Microtritia", "Xenillus",
    "Anoetus", "Opilio", "Xysticus", "Anocentor", "Damaeus", "Oppia",
    "Scheloribates", "Tectocepheus", "Nothrus", "Hermannia")
add("other_crustacean_fw", "Limnochares", "Limmochares", "Piona", "Hydrachna",
    "Unionicola")
add("annelid_fw", "Rhyacodrilus", "Rhynchelmis", "Limmodrilus", "Limnodriulus",
    "Ilodrilus", "Lumbriaulus", "Slavina", "Opistocysta", "Criodrilus",
    "Stylodrilus", "Aeolosoma")
add("annelid_marine", "Scolelepis", "Spirobutschliella", "Sycia",
    "Amphicteis", "Melinna", "Pectinaria")
add("mollusc_fw", "Actinonaias", "Ptychobranchus", "Tropidiscus", "Planorbid",
    "Gyraulus", "Segmentina", "Aplexa", "Lampsilis", "Elliptio", "Villosa")
add("mollusc_terr", "Schizophyllum", "Zonitoides", "Vitrina")
add("helminth_fw", "Posthodiplostomum", "Uvulifer", "Echinostama", "Fasaiola",
    "Allocreadium", "Encyclometra", "Distomum", "Spelotrema", "Plerocercoides",
    "Xiphidiocercaria", "Xiphidiocercariae", "Yamagutisentis", "Dorchis",
    "Apatemon", "Tylodelphys", "Ichthyocotylurus", "Sanguinicola")
add("protist_fw", "Saccamoeba", "Monocystis", "Gregarina", "Gamocystis",
    "Dolyocystis", "Cirrophrya", "Vahlkampfia", "Mayorella", "Cochliopodium")
add("protist_marine", "Ceratomyxa", "Myxobolus", "Myxidium", "Kudoa",
    "Sphaerospora", "Thelohanellus", "Henneguya", "Zschokkella")
add("nematode_terr", "Procephalobus", "Rhabditella", "Procamallanus",
    "Cephalobus", "Zeldia", "Chiloplacus")
add("nematode_marine", "Theristus", "Turbanella", "Kinorhynchus", "Phoronis")
add("reptile", "Mabuya", "Podracis", "Podarcis", "Chalcides", "Emys", "Caretta")
add("mammal", "Callicebus", "Chinchilla", "Nyctalus")
add("lepidoptera", "Abraxas", "Acronycta", "Apopia", "Cheimatobia")
add("other_insect", "Ctenocephalus", "Ctenocephalides")
add("myriapod", "Pachyiulus", "Xenobolus")
GROUPS["myriapod"] = (TE, [TE], "Arachnida & myriapods", "", "high")

# --------------------------------------------------------------------------- #
# Name-shape rules, applied only when the genus is not in the table above.
# These are weaker and are always reported as medium or low confidence.
# --------------------------------------------------------------------------- #
SUFFIX_RULES = [
    (re.compile(r"gammarus$|gammarellus$", re.I), "amphipod_fw", "low"),
    (re.compile(r"cyclops$|cyclopina$", re.I), "copepod_fw", "medium"),
    (re.compile(r"daphnia$|daphnid", re.I), "cladoceran", "medium"),
    (re.compile(r"simulium$|simuliid", re.I), "blackfly", "medium"),
    (re.compile(r"chironomus$|chironom", re.I), "chironomid", "medium"),
    (re.compile(r"penaeus$", re.I), "decapod_marine", "medium"),
    (re.compile(r"cambarus$|astacus$", re.I), "decapod_fw", "medium"),
    (re.compile(r"^culex|^aedes|^anopheles", re.I), "mosquito", "medium"),
]

def _edits(a, b):
    """Levenshtein distance, small strings only."""
    if abs(len(a) - len(b)) > 2:
        return 99
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1,
                           prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


def fuzzy_genus(genus):
    """
    Host genera in the sheet are frequently misspelled -- Chroistoneura,
    Euprocitis, Antherea, Lepeoptheirus, Limmodrilus. Rather than lose those
    rows, a genus within two edits of a known one (same first letter, at least
    six characters) is treated as that genus and reported as a suspected typo.
    """
    g = genus.lower()
    if len(g) < 6:
        return None, None
    best, best_d = None, 3
    for known in GENUS:
        if known[0] != g[0] or abs(len(known) - len(g)) > 2:
            continue
        d = _edits(g, known)
        if d < best_d:
            best, best_d = known, d
    return (best, best_d) if best else (None, None)


NOT_A_TAXON = re.compile(
    r"^(?:amphipod|gammaridean|laboratory|labatory|hamsters|unidentified"
    r"|and|wasp|host|mice|rats|fish|insect|isopod|shrimp|crab|snail"
    r"|[A-Z]\.?)$", re.I)


def predict(host_name):
    """Return (primary, all, group, note, confidence, basis)."""
    if not host_name or not host_name.strip():
        # a fragment that was a description rather than a name
        return ("", [], "Unresolved", "", "none", "no host name parsed")
    if NOT_A_TAXON.match(host_name.strip()):
        return ("", [], "Unresolved", "", "none", "not a taxon name")
    genus = host_name.split()[0]
    key = genus.lower()
    if key in GENUS:
        grp = GENUS[key]
        prim, alls, hg, note, conf = GROUPS[grp]
        return (prim, alls, hg, note, conf, f"genus table ({grp})")
    for rx, grp, conf in SUFFIX_RULES:
        if rx.search(genus):
            prim, alls, hg, note, _ = GROUPS[grp]
            return (prim, alls, hg, note, conf, f"name pattern ({grp})")
    near, dist = fuzzy_genus(genus)
    if near and dist == 1:
        # one edit is safe enough to act on: Euprocitis -> Euproctis
        grp = GENUS[near]
        prim, alls, hg, note, _ = GROUPS[grp]
        return (prim, alls, hg, note, "medium",
                f"fuzzy match to {near.capitalize()} (1 edit) - LIKELY MISSPELLING")
    if near:
        # two edits is NOT acted on. "Abraxas" (a moth) is two edits from
        # "Abramis" (a fish), which would have put a terrestrial lepidopteran
        # in freshwater. The suggestion is reported for review, no habitat is
        # assigned, and the entry stays in the unknown pile.
        return ("", [], "Unresolved", "", "none",
                f"possible misspelling of {near.capitalize()} "
                f"({dist} edits) - NOT APPLIED, check manually")
    return ("", [], "Unresolved", "", "none", "genus not in table")


# --------------------------------------------------------------------------- #
def norm_key(text):
    return re.sub(r"[^a-z0-9]", "", deaccent(str(text)).lower())


def find_col(columns, wanted):
    for c in columns:
        if norm_key(c) == norm_key(wanted):
            return c
    return None


def main():
    p = argparse.ArgumentParser(
        description="Predict host habitat from taxonomic knowledge, offline.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--db", required=True, help="attribute database .xlsx")
    p.add_argument("--sheet", default="Actively Updated Masterlist")
    p.add_argument("--species-col", default="Species Name")
    p.add_argument("--host-cols", nargs="+",
                   default=["Natural Host(s)", "Experimental Host(s)"])
    p.add_argument("--env-col", default="Host Environment")
    p.add_argument("--outdir", default="predicted_environment")
    p.add_argument("--compare", help="01_species_host_long.tsv from the API "
                                     "script, to check the two against each other")
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

    try:
        import pandas as pd
    except ImportError:
        sys.exit("pandas required: pip install pandas openpyxl")
    os.makedirs(args.outdir, exist_ok=True)
    out = lambda f: os.path.join(args.outdir, f)

    df = pd.read_excel(args.db, sheet_name=args.sheet, dtype=str).fillna("")
    sp_col = find_col(df.columns, args.species_col)
    host_cols = [c for c in (find_col(df.columns, h) for h in args.host_cols) if c]
    env_col = find_col(df.columns, args.env_col)
    if not sp_col or not host_cols:
        sys.exit("could not find the species or host columns")

    pairs = []
    for _, row in df.iterrows():
        sp = str(row[sp_col]).strip()
        if not sp or norm_key(sp) == norm_key(args.species_col):
            continue
        db_env = str(row[env_col]).strip().lower() if env_col else ""
        for hc in host_cols:
            for h in parse_host_cell(row[hc]):
                pairs.append((sp, classify_microsporidia(sp), h, hc, db_env))

    # fragments that were descriptions rather than names parse to "" -- they are
    # kept as rows (so the pair is not lost) but never looked up as a host
    # Same database-wide expansion of abbreviated genera as
    # host_taxonomy_environment.py, so the two scripts keep identical host
    # lists and their outputs still join on host_name.
    full_genera, genus_species = Counter(), set()
    for _, _, h, _, _ in pairs:
        t = h["name"].split()
        if t and re.fullmatch(r"[A-Z][a-z-]{2,}", t[0]):
            full_genera[t[0]] += 1
            if len(t) > 1:
                genus_species.add((t[0], t[1].lower()))
    n_exp = 0
    for _, _, h, _, _ in pairs:
        m2 = re.fullmatch(r"([A-Z])\.?\s+([a-z-]{3,})", h["name"])
        if not m2:
            continue
        cands = [g for g in full_genera if g.startswith(m2.group(1))]
        exact = [g for g in cands if (g, m2.group(2).lower()) in genus_species]
        pick = exact[0] if len(exact) == 1 else (cands[0] if len(cands) == 1 else None)
        if pick:
            h["name"] = f"{pick} {m2.group(2).lower()}"
            n_exp += 1
    if n_exp:
        print(f"expanded {n_exp} abbreviated genus name(s)")

    hosts = OrderedDict()
    for _, _, h, _, _ in pairs:
        if h["name"]:
            hosts.setdefault(h["name"], h["rank"])

    host_rows = []
    for name, rank in hosts.items():
        prim, alls, hg, note, conf, basis = predict(name)
        r = {"host_name": name, "host_rank": rank,
             "genus": name.split()[0] if name.split() else "",
             "environment_primary": prim,
             "environment_all": "; ".join(alls),
             "host_group": hg, "life_stage_note": note,
             "confidence": conf, "basis": basis}
        for e in (FW, BR, MA, TE):
            r["is_" + e] = int(e in alls)
        host_rows.append(r)
    by_name = {r["host_name"]: r for r in host_rows}
    blank = {"host_group": "Unresolved", "environment_primary": "",
             "environment_all": "", "life_stage_note": "", "confidence": "none",
             "basis": "no host name parsed",
             **{"is_" + e: 0 for e in (FW, BR, MA, TE)}}

    long_rows = []
    for sp, cat, h, hc, db_env in pairs:
        b = by_name.get(h["name"], blank)
        r = {"microsporidia_species": sp, "microsporidia_category": cat,
             "host_raw": h["raw"], "host_name": h["name"],
             "host_rank": h["rank"], "host_role": h["role"], "host_column": hc,
             "host_group": b["host_group"],
             "environment_primary": b["environment_primary"],
             "environment_all": b["environment_all"],
             "life_stage_note": b["life_stage_note"],
             "confidence": b["confidence"], "basis": b["basis"],
             "db_host_environment": db_env}
        for e in (FW, BR, MA, TE):
            r["is_" + e] = b["is_" + e]
        long_rows.append(r)

    def write_tsv(path, rows):
        if not rows:
            return
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t",
                               extrasaction="ignore", lineterminator="\n")
            w.writeheader()
            w.writerows([{k: tsv_safe(v) for k, v in r.items()}
                         for r in rows])

    write_tsv(out("01_species_host_predicted.tsv"), long_rows)
    write_tsv(out("02_host_predicted.tsv"), host_rows)

    # The genus -> habitat table itself, so it can be deposited and audited
    # rather than existing only inside this script. One row per genus, whether
    # or not that genus appears in the current database.
    genus_rows = []
    for g in sorted(GENUS):
        grp = GENUS[g]
        prim, alls, hg, note, conf = GROUPS[grp]
        genus_rows.append({
            "genus": g.capitalize(), "group": grp,
            "environment_primary": prim, "environment_all": "; ".join(alls),
            "host_group": hg, "life_stage_note": note, "confidence": conf})
    write_tsv(out("00_genus_habitat_table.tsv"), genus_rows)
    print(f"\ngenus habitat table: {len(genus_rows)} genera, "
          f"{len(set(GROUPS))} habitat groups -> 00_genus_habitat_table.tsv")
    used = {r["genus"].lower() for r in host_rows if r.get("genus")}
    print(f"  of these, {len(used & set(GENUS))} were matched by a host in "
          f"this database")
    unknown = [r for r in host_rows if r["confidence"] == "none"]
    write_tsv(out("03_UNKNOWN_hosts.tsv"), unknown)
    lowconf = [r for r in host_rows if r["confidence"] == "low"]
    write_tsv(out("04_low_confidence.tsv"), lowconf)
    typos = [r for r in host_rows if "MISSPELLING" in r["basis"]]
    maybe_typos = [r for r in host_rows if "NOT APPLIED" in r["basis"]]
    write_tsv(out("06_possible_misspellings_unapplied.tsv"), maybe_typos)
    write_tsv(out("05_suspected_misspellings.tsv"), typos)

    # ---- report -----------------------------------------------------------
    print(f"{len(pairs)} species-host pairs, {len(host_rows)} unique host names\n")
    placed = [r for r in host_rows if r["environment_primary"]]
    print(f"predicted            : {len(placed)}/{len(host_rows)} host names "
          f"({100*len(placed)/len(host_rows):.1f}%)")
    pl_pairs = sum(1 for r in long_rows if r["environment_primary"])
    print(f"                       {pl_pairs}/{len(long_rows)} species-host pairs "
          f"({100*pl_pairs/len(long_rows):.1f}%)")
    print("\nconfidence (unique host names):")
    for k, v in Counter(r["confidence"] for r in host_rows).most_common():
        print(f"  {v:>6}  {k}")
    print("\nprimary environment (unique host names):")
    for k, v in Counter(r["environment_primary"] or "(unknown)"
                        for r in host_rows).most_common():
        print(f"  {v:>6}  {k}")
    print("\nprimary environment (species-host pairs):")
    for k, v in Counter(r["environment_primary"] or "(unknown)"
                        for r in long_rows).most_common():
        print(f"  {v:>6}  {k}")
    print("\nhost group (species-host pairs):")
    for k, v in Counter(r["host_group"] for r in long_rows).most_common():
        print(f"  {v:>6}  {k}")

    if typos:
        print(f"\n{len(typos)} host genera look like MISSPELLINGS of a known "
              f"genus (see 05_suspected_misspellings.tsv):")
        for r in typos[:15]:
            print(f"  {r['host_name']:<38} {r['basis']}")
        if len(typos) > 15:
            print(f"  ... {len(typos) - 15} more")
        print("  worth correcting in the spreadsheet: they also break any "
              "host-name join")
    if maybe_typos:
        print(f"\n{len(maybe_typos)} further genera are 2 edits from a known one "
              f"but were NOT applied\n  (too risky - see "
              f"06_possible_misspellings_unapplied.tsv):")
        for r in maybe_typos[:8]:
            print(f"  {r['host_name']:<38} {r['basis']}")

    if unknown:
        print(f"\n{len(unknown)} host names have no prediction "
              f"(see 03_UNKNOWN_hosts.tsv); most frequent genera:")
        gc = Counter(r["genus"] for r in unknown)
        for g, n in gc.most_common(20):
            print(f"  {n:>4}  {g}")

    # ---- optional comparison with the API run ------------------------------
    if args.compare and not os.path.exists(args.compare):
        sys.exit(f"--compare file not found: {args.compare}\n"
                 "  it is written by host_taxonomy_environment.py as\n"
                 "  <outdir>/01_species_host_long.tsv or 02_host_lookup.tsv")
    if args.compare:
        api = {}
        with open(args.compare, encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                envs = {e for e in row.get("environment", "").split("; ") if e}
                if envs:
                    api[row["host_name"]] = envs
        agree = conflict = api_only = pred_only = 0
        rows = []
        for r in host_rows:
            a = api.get(r["host_name"])
            pset = {e for e in r["environment_all"].split("; ") if e}
            if not a and not pset:
                continue
            if not a:
                pred_only += 1
                verdict = "prediction only"
            elif not pset:
                api_only += 1
                verdict = "API only"
            elif a & pset:
                agree += 1
                verdict = "agree"
            else:
                conflict += 1
                verdict = "CONFLICT"
            rows.append({"host_name": r["host_name"],
                         "predicted": r["environment_all"],
                         "api": "; ".join(sorted(a)) if a else "",
                         "verdict": verdict, "confidence": r["confidence"],
                         "basis": r["basis"], "host_group": r["host_group"]})
        write_tsv(out("05_prediction_vs_api.tsv"), rows)
        print(f"\nprediction vs API run ({args.compare}):")
        print(f"  {agree:>6}  agree (overlap in at least one habitat)")
        print(f"  {conflict:>6}  CONFLICT -- check these first")
        print(f"  {api_only:>6}  API only")
        print(f"  {pred_only:>6}  prediction only (the gap this script fills)")

        # ---- how well did each confidence level actually do? --------------
        # Scored only on host names where BOTH have a habitat, since a name the
        # API could not place says nothing about whether the prediction was
        # right. "Correct" = the two overlap in at least one habitat.
        conf_stats = {}
        for r in host_rows:
            a = api.get(r["host_name"])
            pset = {e for e in r["environment_all"].split("; ") if e}
            if not a or not pset:
                continue
            d = conf_stats.setdefault(r["confidence"],
                                      {"n": 0, "ok": 0, "exact": 0})
            d["n"] += 1
            if a & pset:
                d["ok"] += 1
            if a == pset:
                d["exact"] += 1
        if conf_stats:
            print("\naccuracy by confidence level (host names scored by both):")
            print(f"  {'confidence':<12}{'n':>6}{'overlap':>10}{'exact':>10}")
            order = ["high", "medium", "low", "none"]
            tot = {"n": 0, "ok": 0, "exact": 0}
            for c in order + [k for k in conf_stats if k not in order]:
                d = conf_stats.get(c)
                if not d:
                    continue
                for k in tot:
                    tot[k] += d[k]
                print(f"  {c:<12}{d['n']:>6}{100*d['ok']/d['n']:>9.1f}%"
                      f"{100*d['exact']/d['n']:>9.1f}%")
            if tot["n"]:
                print(f"  {'ALL':<12}{tot['n']:>6}{100*tot['ok']/tot['n']:>9.1f}%"
                      f"{100*tot['exact']/tot['n']:>9.1f}%")
            print("  overlap = share at least one habitat; exact = identical set")
            print("  NB this scores the prediction against the APIs, which are")
            print("  themselves incomplete -- a 'conflict' can be the API's error.")

        # accuracy broken down by basis, to show which rules are weak
        basis_stats = {}
        for r in host_rows:
            a = api.get(r["host_name"])
            pset = {e for e in r["environment_all"].split("; ") if e}
            if not a or not pset:
                continue
            key = re.sub(r"\(.*", "", r["basis"]).strip() or r["basis"]
            if "MISSPELLING" in r["basis"]:
                key = "fuzzy match (misspelling)"
            d = basis_stats.setdefault(key, {"n": 0, "ok": 0})
            d["n"] += 1
            if a & pset:
                d["ok"] += 1
        if basis_stats:
            print("\naccuracy by rule type:")
            for k, d in sorted(basis_stats.items(), key=lambda kv: -kv[1]["n"]):
                print(f"  {k:<28}{d['n']:>6}{100*d['ok']/d['n']:>9.1f}%")

        # which host groups the prediction gets wrong most often
        grp_stats = {}
        for r in host_rows:
            a = api.get(r["host_name"])
            pset = {e for e in r["environment_all"].split("; ") if e}
            if not a or not pset:
                continue
            d = grp_stats.setdefault(r["host_group"], {"n": 0, "bad": 0})
            d["n"] += 1
            if not (a & pset):
                d["bad"] += 1
        worst = [(k, v) for k, v in grp_stats.items() if v["bad"]]
        if worst:
            print("\nhost groups with conflicts:")
            for k, v in sorted(worst, key=lambda kv: -kv[1]["bad"]):
                print(f"  {k:<28}{v['bad']:>4} of {v['n']:<5} "
                      f"({100*v['bad']/v['n']:.0f}% conflict)")

    print(f"\nwritten to {args.outdir}/")
    print("  01_species_host_predicted.tsv  one row per microsporidia species x host")
    print("  00_genus_habitat_table.tsv     the genus -> habitat table itself")
    print("  02_host_predicted.tsv          one row per unique host name")
    print("  03_UNKNOWN_hosts.tsv           genus not in the table, no prediction")
    print("  04_low_confidence.tsv          predicted from the name alone, verify")
    print("  05_suspected_misspellings.tsv  1 edit from a known genus, applied")
    print("  06_possible_misspellings_unapplied.tsv  2 edits, flagged only")
    print("\nNOTE: these habitats come from taxonomic knowledge written into this")
    print("script, not from any database. Verify the low-confidence rows, and")
    print("prefer the API result where the two disagree.")


if __name__ == "__main__":
    main()
