#!/usr/bin/env python3
"""
Turn 18S sequence clusters into three clearly-labelled outputs.

OUTPUT 1  01_ALREADY_IN_DATABASE.tsv
          Clusters that match a species already in your attribute database.
          Nothing to add taxonomically, BUT each row lists the hosts,
          localities and accessions the cluster carries that your database
          row does NOT -- i.e. new host records for known species.

OUTPUT 2  02_NEW_ENTRIES_TO_ADD.xlsx / .tsv
          Clusters NOT in your attribute database, in your exact column
          order, ready to paste in. One row per cluster. Includes both:
            - named species     (cluster carries a real binomial)
            - provisional species (no binomial, but has an 18S sequence and
                                   a reported host) -- given a stable ID
                                   like MSP-C0042 so you can refer to them.

OUTPUT 3  03_CLUSTER_HOST_ATTRIBUTES.tsv  <- read by every downstream step
          Every cluster, named and provisional, with all hosts, localities
          and accessions. This is the tree-annotation table: join it to your
          tree on cluster_label.

          all_hosts here is the UNION of the GenBank metadata and, for a
          cluster matched to a database entry, that entry's curated
          Natural Host(s). Without the merge a host recorded only in the
          curated database never reaches the tree, the environment rings or
          the host-taxonomy symbols, all of which read this file. Provenance
          is preserved in hosts_from_genbank / hosts_from_database /
          host_source, and --no-merge-database-hosts restores the
          GenBank-only view.

          Curated attributes go to ONE cluster per database row. Where a
          species name spans several clusters, only the cluster containing
          that row's voucher accessions inherits its curated hosts; the
          others keep the GenBank records of their own members. Without this
          a handful of literature records would be replicated across dozens
          of tips, inflating apparent host breadth. attribute_basis and
          owns_database_record record which applies to each cluster.

Plus two supporting files:
  04_TREE_TIP_LABELS.tsv     centroid -> display label mapping for the tree
  07_CENTROIDS_WITH_HOSTS.fasta  the centroid sequence of every cluster that has
                             at least one reported host -- the tree-building
                             input, filtered so no host-less cluster appears
  05_FLAGGED_FOR_REVIEW.tsv  clusters needing a human decision

Usage
-----
  python cluster_to_database.py \
      --uc clusters99.uc \
      --meta trimmed.annot.fasta.meta.tsv \
      --db Microsporidia_Characteristics_Database.xlsx \
      --outdir cluster_output
"""

import argparse
import csv
import os
import re
import sys
import unicodedata
import difflib
from collections import Counter, defaultdict
from datetime import datetime

# ---------------------------------------------------------------------------
# database layout
# ---------------------------------------------------------------------------
DB_COLUMNS = [
    "Timestamp", "Species Name", "Date Identified (year)", "Natural Host(s)",
    "Experimental Host(s)", "Host Environment",
    "Host Life Stage during Infection", "Site of Infection", "Transmission",
    "Spore Length Average (\u00b5m)", "Spore Width Average (\u00b5m)",
    "Calculated Volume (\u00b5m\u00b3) (see methods)",
    "Spore Shape (Class; Condition)", "Locality", "Nucleus",
    "Measured Polar Tubule Length Max (\u03bcm)",
    "Measured Polar Tubule Length Min (\u03bcm)",
    "Measured Polar Tubule (\u03bcm)",
    "Calculated Polar Tubule (\u03bcm) (see methods)",
    "Polar Tubule Coils Range", "Polar Tubule Coils Average",
    "18S Accession #", "Has the genome been sequenced?", "Important Remarks",
    "References",
]

ACC_RE = re.compile(r"^[A-Z]{1,6}_?\d{5,9}(\.\d+)?$")

# names that carry no taxonomic information. The prefix list is IMPORTED
# from the shared canonical classifier (classify_species.py), so it can
# never drift from the rule the figure scripts and the flowchart use.
from classify_species import UNINFORMATIVE_RX, classify_microsporidia
UNINFORMATIVE = re.compile(UNINFORMATIVE_RX, re.I)


# ---------------------------------------------------------------------------
# normalisation
# ---------------------------------------------------------------------------
def year_of(text):
    """Earliest plausible 4-digit year in a date-ish string, or ''."""
    if not text:
        return ""
    yrs = [int(y) for y in re.findall(r"(?<!\d)(1[6-9]\d{2}|20\d{2})(?!\d)", str(text))]
    yrs = [y for y in yrs if 1600 <= y <= datetime.now().year]
    return str(min(yrs)) if yrs else ""


def acc_key(a):
    """
    Canonical form for comparing accessions.

    RefSeq accessions carry an underscore (XR_552278) that the curated database
    sometimes omits (XR552278). Comparing the raw strings makes those two
    different accessions, so a database record can silently match no cluster at
    all. Dropping the underscore makes the two forms equal. Ordinary GenBank
    accessions contain no underscore, so nothing else is affected.
    """
    return str(a).upper().replace("_", "").split(".")[0]


def norm_acc(text):
    """Normalise one accession token; return '' if it isn't accession-like."""
    t = str(text).strip().upper().split(".")[0]
    return t if ACC_RE.match(t) or ACC_RE.match(t + ".1") else ""


def split_accessions(cell):
    """Pull every accession out of a free-text database cell."""
    if not cell:
        return []
    s = str(cell).upper()
    # repair accessions written with an internal space, e.g. "MG 708238"
    s = re.sub(r"\b([A-Z]{1,6})\s+(\d{5,9})\b", r"\1\2", s)
    out = []
    for tok in re.split(r"[,;/|\s]+", s):
        tok = tok.strip("()[].")
        a = norm_acc(tok)
        if a:
            out.append(a)
    return out


# nomenclatural acts and other trailing apparatus that is not part of the name
NOM_ACTS = re.compile(
    r"\b(?:"
    r"gen\.?\s*(?:et|and|&)?\s*(?:sp|spp)\.?\s*nov\.?"
    r"|(?:sp|spp|ssp|subsp|var|comb|stat|fam|ord|nom)\.?\s*nov\.?"
    r"|nov\.?\s*(?:gen|sp|comb)\.?"
    r"|n\.\s*(?:gen|sp)\.?"
    r"|emend\.?|sensu\s+(?:lato|stricto)|s\.\s*[ls]\."
    r"|incertae\s+sedis"
    r")", re.I)

QUALIFIERS = re.compile(r"\b(?:cf|aff|nr|ex)\.?\s+", re.I)

# infraspecific / culture apparatus: drop the marker and whatever follows it
INFRA = re.compile(
    r"\b(?:var|subsp|ssp|forma|morphotype|genotype|serotype|strain|isolate"
    r"|clone|haplotype|lineage|type)\.?\s*\S*", re.I)

# trailing authority, e.g. "Nageli, 1857" or "(Balbiani, 1882)" or "Sprague & Vavra 1976"
AUTHORITY = re.compile(
    r"\s*\(?\b[A-Z][A-Za-z\u00c0-\u024f'\u2019-]{2,}"
    r"(?:\s*(?:&|and|et)\s*[A-Z][A-Za-z\u00c0-\u024f'\u2019-]{2,})*"
    r"(?:\s+et\s+al\.?)?[,\s]+(?:1[5-9]\d{2}|20\d{2})\b\)?")

SYN_PREFIX = re.compile(
    r"^(?:=+|syn\.?|synonym|synonymous\s+with|formerly|previously|was|now|also"
    r"|aka|a\.k\.a\.?|orig\.?|basionym|renamed|reassigned)"
    r"\s*:?\s*", re.I)


def strip_accents(text):
    d = unicodedata.normalize("NFKD", str(text))
    return "".join(c for c in d if not unicodedata.combining(c))


def alt_names(raw):
    """
    Pull alternative names out of a database cell. A parenthetical counts as a
    synonym only if it actually looks like a binomial -- "Ameson nelsoni
    (Perezia nelsoni)" yields Perezia nelsoni, while "(siberian isolate)",
    "(18S partial)" and "(?)" yield nothing.
    """
    out = []
    text = str(raw)
    for inner in re.findall(r"\(([^)]*)\)", text):
        cand = SYN_PREFIX.sub("", inner.strip())
        toks = re.sub(r"[^A-Za-z\s.]", " ", strip_accents(cand)).split()
        if (len(toks) >= 2 and len(toks[0]) > 2
                and toks[0][:1].isupper() and toks[1][:1].islower()
                and toks[1].rstrip(".") not in ("sp", "spp")):
            out.append(" ".join(toks[:2]))
    # names joined outside parentheses: "Nosema apis = Vairimorpha apis"
    for part in re.split(r"\s*=+\s*|\s+syn\.?\s+", re.sub(r"\([^)]*\)", " ", text))[1:]:
        toks = re.sub(r"[^A-Za-z\s.]", " ", strip_accents(part)).split()
        if (len(toks) >= 2 and len(toks[0]) > 2
                and toks[0][:1].isupper() and toks[1][:1].islower()):
            out.append(" ".join(toks[:2]))
    return out


def norm_name(name, keep_all_tokens=False):
    """
    Normalise a species name to a comparable key.

    Handles: accents, trailing authority and year, nomenclatural acts
    ("sp. nov.", "gen. et sp. nov.", "n. sp."), cf./aff. qualifiers,
    infraspecific and strain apparatus, and parenthetical asides. The key is
    genus + specific epithet only, so "Metchnikovella dobrovolskiji sp. nov."
    and "Metchnikovella dobrovolskiji" collapse to the same string.
    """
    if not name:
        return ""
    s = AUTHORITY.sub(" ", str(name))
    s = strip_accents(s).lower()
    s = re.sub(r"\(.*?\)", " ", s)
    s = NOM_ACTS.sub(" ", s)
    s = QUALIFIERS.sub(" ", s)
    s = INFRA.sub(" ", s)
    s = re.sub(r"\bspp?\.?(?=\s|$)", " sp ", s)
    s = re.sub(r"[^a-z0-9\s]", " ", s)
    tokens = [t for t in s.split() if t not in ("nov", "gen", "et", "and")]
    if keep_all_tokens:
        return " ".join(tokens)
    return " ".join(tokens[:2])


def name_keys(raw):
    """Every key a database entry should be findable under."""
    keys = []
    primary = norm_name(raw)
    if primary:
        keys.append(primary)
    for a in alt_names(raw):
        k = norm_name(a)
        if k and k not in keys:
            keys.append(k)
    return keys


def abbrev_key(key):
    """'nematocida parisii' -> 'n parisii', to match 'N. parisii'."""
    parts = key.split()
    if len(parts) == 2 and len(parts[0]) > 1 and parts[1] != "sp":
        return f"{parts[0][0]} {parts[1]}"
    return ""


def name_rank(name):
    """0 = usable binomial, 1 = 'Genus sp.' or placeholder, 2 = uninformative.

    Aligned with the canonical R classifier used by the figure scripts and
    the flowchart: a leading bare "sp" is never a genus, and a purely numeric
    second token ("Microsporidium 1") is a placeholder, not an epithet -- both
    count as rank 1, so such clusters come out PROVISIONAL, exactly as the
    same names do in every panel.
    """
    n = norm_name(name)
    if not n or UNINFORMATIVE.match(n):
        return 2
    parts = n.split()
    if parts[0] == "sp":
        return 1
    if len(parts) >= 2 and parts[1] != "sp" and not parts[1].isdigit():
        return 0
    return 1


def split_multi(cell, extra_sep=False):
    """Split a database cell holding several values."""
    if not cell:
        return []
    parts = re.split(r"\s*;\s*|\n", str(cell))
    if extra_sep:
        parts = [q for p in parts for q in re.split(r"\s*,\s*", p)]
    return [p.strip() for p in parts if p.strip()]


def dedupe(values):
    """Case-insensitive de-duplication, preserving first-seen spelling."""
    seen, out = set(), []
    for v in values:
        v = " ".join(str(v).split())
        k = v.lower()
        if v and k not in seen:
            seen.add(k)
            out.append(v)
    return out


# ---------------------------------------------------------------------------
# inputs
# ---------------------------------------------------------------------------
def clean_label(label):
    """Strip vsearch size annotations and accession versions from a .uc label."""
    lab = label.split(";")[0].strip()
    lab = re.sub(r";?size=\d+;?", "", lab)
    return lab


def read_uc(path):
    """Return {cluster_id: [member_label, ...]} from a vsearch .uc file."""
    clusters = defaultdict(list)
    centroid_of = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10 or f[0] == "C":
                continue
            cid, query, target = f[1], clean_label(f[8]), clean_label(f[9])
            if f[0] == "S":
                clusters[cid].append(query)
                centroid_of[cid] = query
            elif f[0] == "H":
                clusters[cid].append(query)
                centroid_of.setdefault(cid, target)
    return dict(clusters), centroid_of


def read_meta(path):
    """Return {accession_no_version: {field: value}} from the metadata TSV."""
    meta = {}
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            acc = row.get("accession", "")
            key = norm_acc(acc) or acc.split(".")[0].upper()
            if key:
                meta[key] = row
    return meta


def read_db(path, sheets):
    """Return (rows, accession_set, {normalised_name: [row_index, ...]})."""
    try:
        import pandas as pd
    except ImportError:
        sys.exit("pandas is required to read the database (pip install pandas openpyxl)")

    frames = []
    xl = pd.ExcelFile(path)
    wanted = sheets or xl.sheet_names
    for sh in wanted:
        if sh not in xl.sheet_names:
            sys.stderr.write(f"  note: sheet {sh!r} not in workbook, skipping\n")
            continue
        df = xl.parse(sh).fillna("")
        df["__sheet__"] = sh
        frames.append(df)
    if not frames:
        sys.exit("No usable sheets found in the database workbook.")

    rows, acc_set = [], set()
    by_name, by_abbrev, by_epithet = defaultdict(list), defaultdict(list), defaultdict(list)
    for df in frames:
        for _, r in df.iterrows():
            rec = {k: ("" if v == "" else str(v)) for k, v in r.items()}
            idx = len(rows)
            rows.append(rec)
            for a in split_accessions(rec.get("18S Accession #", "")):
                acc_set.add(acc_key(a))
            for k in name_keys(rec.get("Species Name", "")):
                by_name[k].append(idx)
                ab = abbrev_key(k)
                if ab:
                    by_abbrev[ab].append(idx)
                parts = k.split()
                if len(parts) == 2 and parts[1] != "sp":
                    by_epithet[parts[1]].append(idx)
    return (rows, acc_set, dict(by_name), dict(by_abbrev), dict(by_epithet))


# ---------------------------------------------------------------------------
# cluster aggregation
# ---------------------------------------------------------------------------
def pick_species_name(names):
    """Choose one display name for the cluster; return (name, all_distinct)."""
    distinct = dedupe(n for n in names if n)
    if not distinct:
        return "", []
    best = min(name_rank(n) for n in distinct)
    candidates = [n for n in distinct if name_rank(n) == best]
    # Among names for the same species, prefer the cleanest rendering: the one
    # that is closest to a bare binomial, so "Nematocida parisii" wins over
    # "Nematocida parisii ERTm1". Frequency breaks remaining ties.
    counts = Counter(n for n in names if n in candidates)
    chosen = sorted(candidates,
                    key=lambda n: (len(n.split()) > 2, len(n), -counts[n]))[0]
    return chosen, distinct


def aggregate_cluster(cid, members, centroid, meta):
    infos = []
    for m in members:
        key = norm_acc(m) or m.split(".")[0].upper()
        infos.append(meta.get(key, {"accession": m}))

    names = [i.get("organism", "") for i in infos]
    chosen, distinct_names = pick_species_name(names)
    # Count distinct SPECIES, not distinct name strings. "Nematocida parisii
    # ERTm1" and "Nematocida parisii" are one species written two ways, so
    # collapse to the normalised genus+epithet key before counting.
    binomial_keys = dedupe(norm_name(n) for n in distinct_names
                           if name_rank(n) == 0)
    binomials = binomial_keys

    hosts = dedupe(h for i in infos for h in split_multi(i.get("host", "")))
    lab_hosts = dedupe(h for i in infos for h in split_multi(i.get("lab_host", "")))
    localities = dedupe(i.get("geo_loc_name", "") for i in infos)
    countries = dedupe(i.get("country", "") for i in infos)
    coll_years = dedupe(i.get("collection_year", "") for i in infos)
    pub_years = dedupe(i.get("pub_year", "") for i in infos)
    # submitted_year is the date the sequence was uploaded to GenBank. Every
    # record has one; pub_year exists only for published records, which is why
    # relying on publication alone leaves clusters dateless.
    sub_years = dedupe(i.get("submitted_year", "") for i in infos)
    # last resort: the LOCUS line date (last modified) -- an upper bound on the
    # submission date, never earlier than it.
    rec_years = dedupe(year_of(i.get("record_date", "")) for i in infos)
    sources = dedupe(i.get("isolation_source", "") for i in infos)
    strains = dedupe(v for i in infos for v in
                     (i.get("isolate", ""), i.get("strain", "")))
    pmids = dedupe(p for i in infos for p in split_multi(i.get("pubmed", ""), True))
    refs = dedupe(i.get("reference", "") for i in infos)
    accs = [i.get("accession", "") or m for i, m in zip(infos, members)]

    # per-accession provenance: which record reported which host/locality.
    # This is what makes an "ambiguous" multi-species cluster resolvable --
    # a host belongs to the species named on the record that reported it.
    per_rec = []
    for i, mem in zip(infos, members):
        hs = dedupe(h for h in split_multi(i.get("host", "")))
        if not hs and not i.get("geo_loc_name"):
            continue
        per_rec.append({
            "accession": i.get("accession", "") or mem,
            "genbank_organism": i.get("organism", ""),
            "hosts": "; ".join(hs),
            "locality": i.get("geo_loc_name", ""),
        })

    return {
        "cluster_id": cid,
        "centroid": centroid,
        "per_record": per_rec,
        "n_members": len(members),
        "n_with_host": sum(1 for i in infos if i.get("host")),
        "n_with_locality": sum(1 for i in infos if i.get("geo_loc_name")),
        "species_name": chosen,
        "n_distinct_names": len(distinct_names),
        "n_binomials": len(binomials),
        "distinct_species": " | ".join(binomial_keys),
        "all_names": " | ".join(distinct_names),
        "hosts": "; ".join(hosts),
        "n_hosts": len(hosts),
        "experimental_hosts": "; ".join(lab_hosts),
        "localities": "; ".join(localities),
        "countries": "; ".join(countries),
        "n_countries": len(countries),
        "collection_years": "; ".join(sorted(coll_years)),
        "publication_years": "; ".join(sorted(pub_years)),
        "submission_years": "; ".join(sorted(sub_years)),
        "record_years": "; ".join(sorted(rec_years)),
        "isolation_sources": "; ".join(sources),
        "isolates_strains": "; ".join(strains),
        "accessions": "; ".join(accs),
        "pubmed": ", ".join(pmids),
        "references": " | ".join(refs),
    }


# ---------------------------------------------------------------------------
# database matching
# ---------------------------------------------------------------------------
def match_cluster(cl, db_acc_set, db_by_name, db_by_abbrev, db_by_epithet,
                  db_rows, fuzzy_cutoff=0.9):
    """
    Match a cluster to existing database rows in tiers, most reliable first.
    Returns (status, matched_by, matched_row_indices).

    status is "existing" for confident matches and "probable" for the two
    heuristic tiers, which are held back for review rather than being treated
    as either new or already-present.
    """
    # 1. shared accession -- the only unambiguous signal
    hits = set()
    for a in split_accessions(cl["accessions"]):
        if acc_key(a) in db_acc_set:
            for i, row in enumerate(db_rows):
                if acc_key(a) in {acc_key(x) for x in
                                  split_accessions(row.get("18S Accession #", ""))}:
                    hits.add(i)
    if hits:
        return "existing", "accession", sorted(hits)

    name = cl["species_name"]
    if not name or name_rank(name) != 0:
        return "new", "", []
    key = norm_name(name)
    genus, epithet = (key.split() + [""])[:2]

    # 2. exact key, including synonyms indexed from parentheses
    if key in db_by_name:
        return "existing", "species name", db_by_name[key]

    # 3. abbreviated genus in either direction
    ab = abbrev_key(key)
    if ab and ab in db_by_abbrev:
        return "existing", "abbreviated genus", db_by_abbrev[ab]
    if len(genus) == 1 and f"{genus} {epithet}" in db_by_abbrev:
        return "existing", "abbreviated genus", db_by_abbrev[f"{genus} {epithet}"]

    # 4. same epithet under a different genus -- the reassignment signature
    if epithet and epithet in db_by_epithet:
        return "probable", "same epithet, different genus", db_by_epithet[epithet]

    # 5. near-identical spelling within the same genus (nelsoni / nelsonii)
    same_genus = [k for k in db_by_name if k.split()[:1] == [genus]]
    close = difflib.get_close_matches(key, same_genus, n=3, cutoff=fuzzy_cutoff)
    if close:
        idxs = sorted({i for k in close for i in db_by_name[k]})
        return "probable", f"close spelling to {close[0]}", idxs

    return "new", "", []


def novel_values(cluster_values, db_rows, idxs, column, extra_sep=False):
    """Values the cluster has that none of the matched database rows record."""
    have = set()
    for i in idxs:
        for v in split_multi(db_rows[i].get(column, ""), extra_sep):
            have.add(v.lower())
    out = []
    for v in split_multi(cluster_values):
        low = v.lower()
        if low and not any(low in h or h in low for h in have):
            out.append(v)
    return out


# ---------------------------------------------------------------------------
# output
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# cluster classification
# ---------------------------------------------------------------------------
def classify_cluster(cl, min_hosts_for_provisional=1):
    """
    Decide what kind of entity a cluster represents.

      "named"        the cluster carries a real binomial
      "provisional"  no binomial, but it has an 18S sequence and a reported
                     host -- a real biological entity awaiting a name
      "insufficient" no binomial and no host: an 18S sequence with no
                     associated biology, nothing to enter

    Returns (kind, reason).
    """
    if cl["n_binomials"] >= 1:
        return "named", ""
    # use the merged host list where it exists, so a cluster known only from
    # the curated database still counts as provisional
    hosts = cl.get("hosts_merged", cl.get("hosts", ""))
    n_hosts = len([h for h in hosts.split("; ") if h])
    if n_hosts >= min_hosts_for_provisional:
        return "provisional", ""
    return "insufficient", "no binomial and no reported host"


def provisional_label(cl, prefix):
    """A stable, human-readable name for an unnamed cluster."""
    genera = [n.split()[0].capitalize() for n in
              (cl["all_names"].split(" | ") if cl["all_names"] else [])
              if n and not UNINFORMATIVE.match(n.lower())]
    genus = Counter(genera).most_common(1)[0][0] if genera else ""
    tag = f"{prefix}-C{int(cl['cluster_id']):04d}" if cl["cluster_id"].isdigit() \
        else f"{prefix}-C{cl['cluster_id']}"
    return f"{genus} sp. {tag}" if genus else f"Microsporidia sp. {tag}"


# ---------------------------------------------------------------------------
# database row construction
# ---------------------------------------------------------------------------
# Candidate dates for 'Date Identified (year)', in the order --year-from auto
# tries them. GenBank carries no species-description date, so all four are
# proxies: publication is the closest to a description year, submission is the
# only one that is always present, collection and record date are weaker.
YEAR_SOURCES = [
    ("publication", "publication_years", "publication year"),
    ("submission", "submission_years", "GenBank submission year"),
    ("collection", "collection_years", "sample collection year"),
    ("record", "record_years", "GenBank record date"),
]


def resolve_cluster_year(cl, year_from):
    """Return (year_string, basis_label, source_key). '' if nothing available."""
    if year_from == "none":
        return "", "", ""
    order = YEAR_SOURCES if year_from == "auto" else \
        ([t for t in YEAR_SOURCES if t[0] == year_from] +
         [t for t in YEAR_SOURCES if t[0] != year_from])
    for key, field, label in order:
        val = str(cl.get(field, "") or "").strip()
        if val:
            return val, label, key
    return "", "", ""


def build_db_row(cl, timestamp, year_from, kind, display_name):
    row = {c: "" for c in DB_COLUMNS}
    row["Timestamp"] = timestamp
    row["Species Name"] = display_name
    row["Natural Host(s)"] = cl["hosts"]
    row["Experimental Host(s)"] = cl["experimental_hosts"]
    row["Locality"] = cl["localities"]
    row["18S Accession #"] = cl["accessions"]

    year, basis, basis_key = resolve_cluster_year(cl, year_from)
    row["Date Identified (year)"] = year
    cl["_date_basis"] = basis_key or "none"

    remarks = []
    if kind == "provisional":
        remarks.append(
            "PROVISIONAL SPECIES: defined by 18S sequence cluster, no formal "
            "species name in GenBank")
    remarks.append(f"cluster {cl['cluster_id']}: {cl['n_members']} sequences, "
                   f"centroid {cl['centroid']}")
    for label, key in [("isolates/strains", "isolates_strains"),
                       ("isolation source", "isolation_sources"),
                       ("collection year", "collection_years"),
                       ("GenBank submission year", "submission_years")]:
        if cl.get(key):
            remarks.append(f"{label}: {cl[key]}")
    # provenance of the date cell: none of these is a description year, so say
    # which proxy was used. Downstream figure scripts key off this tag.
    if year and basis:
        # the cell may hold several years ("2019; 2021"); the tag carries the
        # earliest, and no semicolon, so it survives the "; ".join below
        remarks.append(f"date basis: {basis} {year_of(year) or year}")
    if cl["all_names"]:
        remarks.append(f"GenBank organism names in cluster: {cl['all_names']}")
    row["Important Remarks"] = "; ".join(remarks)

    refs = []
    if cl["references"]:
        refs.append(cl["references"])
    if cl["pubmed"]:
        refs.append("PMID: " + cl["pubmed"])
    refs.append("GenBank " + cl["accessions"])
    row["References"] = "; ".join(refs)
    return row


def write_tsv(path, fieldnames, rows):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t",
                           extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def write_xlsx(path, rows, note):
    try:
        from openpyxl import Workbook
        from openpyxl.styles import Alignment, Font, PatternFill
        from openpyxl.utils import get_column_letter
    except ImportError:
        sys.stderr.write("openpyxl not installed -- .xlsx skipped, .tsv written\n")
        return False
    wb = Workbook()
    ws = wb.active
    ws.title = "New entries"
    ws.append(DB_COLUMNS)
    for c in ws[1]:
        c.font = Font(name="Arial", size=10, bold=True)
        c.alignment = Alignment(wrap_text=True, vertical="top")
        c.fill = PatternFill("solid", fgColor="E8E8E8")
    for r in rows:
        ws.append([r.get(c, "") for c in DB_COLUMNS])
    for r in ws.iter_rows(min_row=2):
        for c in r:
            c.font = Font(name="Arial", size=10)
            c.alignment = Alignment(vertical="top")
            if isinstance(c.value, str):
                c.number_format = "@"
    widths = {"Species Name": 32, "Natural Host(s)": 34, "Locality": 26,
              "Important Remarks": 60, "References": 45, "Timestamp": 20,
              "18S Accession #": 26}
    for i, col in enumerate(DB_COLUMNS, start=1):
        ws.column_dimensions[get_column_letter(i)].width = widths.get(col, 14)
    ws.freeze_panes = "A2"

    info = wb.create_sheet("README")
    for i, line in enumerate(note, start=1):
        info.cell(row=i, column=1, value=line).font = Font(name="Arial", size=10)
    info.column_dimensions["A"].width = 100
    wb.save(path)
    return True


# ---------------------------------------------------------------------------
# merging curated database values into a cluster's aggregate
# ---------------------------------------------------------------------------
def db_values(db_rows, idxs, column, extra_sep=False):
    """Every distinct value of `column` across the matched database rows."""
    out = []
    for i in idxs:
        out.extend(split_multi(db_rows[i].get(column, ""), extra_sep))
    return dedupe(out)


def merge_sources(genbank, database):
    """
    Union two "; "-joined lists, case-insensitively, keeping GenBank's spelling
    where both have the same value. Returns (merged, n_only_db).
    """
    gb = split_multi(genbank)
    db = split_multi(database)
    seen = {v.lower() for v in gb}
    extra = [v for v in db if v.lower() not in seen]
    return dedupe(gb + extra), len(extra)


# ---------------------------------------------------------------------------
# filtered centroid FASTA
# ---------------------------------------------------------------------------
def read_fasta(path):
    """Yield (header_without_gt, sequence)."""
    head, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if head is not None:
                    yield head, "".join(seq)
                head, seq = line[1:], []
            else:
                seq.append(line.strip())
    if head is not None:
        yield head, "".join(seq)


def fasta_id(header):
    """
    First whitespace-delimited token, with vsearch's decorations removed:
      "240|OK155998.1 desc"  -> "OK155998.1"
      "OK155998.1;size=12;"  -> "OK155998.1"
    """
    tok = header.split()[0] if header.split() else header
    tok = re.sub(r"^[0-9]+\|", "", tok)          # cluster-number prefix
    tok = re.sub(r";size=\d+;?", "", tok)        # abundance annotation
    return tok.strip(";")


def write_filtered_centroids(src, dest, keep, label_mode, line_width=60):
    """
    Copy the centroids of clusters in `keep` into a new FASTA.

    Matching is tried on the exact id, then with the version suffix dropped,
    so a mismatch in how the accession was written cannot silently empty the
    output.
    """
    keep_base = {k.split(".")[0]: v for k, v in keep.items()}
    written, seen = 0, set()
    with open(dest, "w") as out:
        for header, seq in read_fasta(src):
            fid = fasta_id(header)
            name = keep.get(fid)
            if name is None:
                name = keep_base.get(fid.split(".")[0])
            if name is None:
                continue
            seen.add(fid)
            if label_mode == "display":
                hdr = re.sub(r"[\s,:;()\[\]'\"|/\\]+", "_", name)
                hdr = f"{hdr} [centroid={fid}]"
            else:
                hdr = f"{fid} [cluster_label={name}]"
            out.write(f">{hdr}\n")
            for i in range(0, len(seq), line_width):
                out.write(seq[i:i + line_width] + "\n")
            written += 1
    return written, seen


def main():
    p = argparse.ArgumentParser(
        description="Turn 18S clusters into database entries and a tree "
                    "attribute table.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--uc", required=True, help="vsearch .uc cluster file")
    p.add_argument("--meta", required=True,
                   help="metadata TSV from annotate_fasta_ncbi.py")
    p.add_argument("--db", required=True, help="attribute database .xlsx")
    p.add_argument("--db-sheets", nargs="*", default=None,
                   help="sheets to check (default: every sheet)")
    p.add_argument("--outdir", default="cluster_output")
    p.add_argument("--timestamp", default=None)
    p.add_argument("--year-from",
                   choices=["none", "auto", "publication", "submission",
                            "collection", "record"],
                   default="auto",
                   help="what to put in 'Date Identified (year)'. 'auto' tries "
                        "publication year, then GenBank submission year, then "
                        "collection year, then record date -- so a cluster with "
                        "no publication is dated rather than left blank. Naming "
                        "one source tries it first and falls back to the rest; "
                        "'none' leaves the column empty. Whichever was used is "
                        "written to Important Remarks as 'date basis: ...'")
    p.add_argument("--provisional-prefix", default="MSP",
                   help="prefix for provisional species IDs, e.g. MSP-C0042")
    p.add_argument("--no-provisional", action="store_true",
                   help="only propose clusters that carry a real binomial")
    p.add_argument("--require-host-for-provisional", type=int, default=1,
                   help="minimum reported hosts for an unnamed cluster to "
                        "count as a provisional species (0 to include all)")
    p.add_argument("--fuzzy-cutoff", type=float, default=0.9)
    p.add_argument("--probable-as-new", action="store_true",
                   help="treat heuristic name matches as new rather than "
                        "holding them for review")
    p.add_argument("--min-members", type=int, default=1)
    p.add_argument("--no-merge-database-hosts", dest="merge_db_hosts",
                   action="store_false", default=True,
                   help="do NOT merge curated hosts/localities from the "
                        "database into 03_CLUSTER_HOST_ATTRIBUTES.tsv; keep "
                        "that file as a GenBank-only view")
    p.add_argument("--centroids",
                   help="FASTA of cluster centroids (vsearch --centroids "
                        "output). When given, a filtered copy containing only "
                        "clusters with at least --min-hosts host(s) is written "
                        "to 07_CENTROIDS_WITH_HOSTS.fasta")
    p.add_argument("--min-hosts", type=int, default=1,
                   help="minimum distinct reported hosts for a cluster's "
                        "centroid to be kept in the filtered FASTA")
    p.add_argument("--fasta-label", choices=["accession", "display"],
                   default="accession",
                   help="header for the filtered FASTA: the original centroid "
                        "accession, or the cluster's display name")
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    out = lambda f: os.path.join(args.outdir, f)
    timestamp = args.timestamp or datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    clusters, centroids = read_uc(args.uc)
    meta = read_meta(args.meta)
    (db_rows, db_acc_set, db_by_name, db_by_abbrev,
     db_by_epithet) = read_db(args.db, args.db_sheets)
    sys.stderr.write(
        f"{len(clusters)} clusters / {sum(len(v) for v in clusters.values())} sequences\n"
        f"{len(meta)} accessions with metadata\n"
        f"database: {len(db_rows)} rows, {len(db_acc_set)} accessions, "
        f"{len(db_by_name)} name keys\n\n")

    attributes, tip_labels, attrib_rows = [], [], []
    n_augmented = n_hosts_added = 0
    # accession -> the database rows that list it
    acc_to_rows = defaultdict(list)
    for i, row in enumerate(db_rows):
        for a in split_accessions(row.get("18S Accession #", "")):
            # key on acc_key() so a RefSeq accession written XR_552278 here
            # matches the XR552278 form used when this dict is queried below;
            # without this every RefSeq lookup missed and the per-record host
            # attribution fell back to a weaker rule
            acc_to_rows[acc_key(a)].append(i)
    keep_centroids = {}      # centroid accession -> display name, for the FASTA
    already, to_add, flagged, undated = [], [], [], []
    counts = Counter()
    date_basis = Counter()

    ## ---- PASS 1: assign each database row's curated data to ONE cluster -----
    ## A species name can be spread over many clusters (Vittaforma corneae spans
    ## 32 here). Letting every one of them inherit the same curated hosts would
    ## replicate a handful of literature records across dozens of tips, which
    ## inflates apparent host breadth and pseudo-replicates the host-specificity
    ## analysis. Only the cluster that actually CONTAINS a row's voucher
    ## accessions has a demonstrated link to that record, so it alone inherits
    ## the curated attributes.
    ##
    ## Winner per database row: the cluster holding the most of that row's
    ## accessions; ties go to the cluster holding the accession listed FIRST in
    ## the row's "18S Accession #" cell. Naming is unaffected -- a cluster still
    ## takes the species name by accession OR by name, exactly as before.
    ordered_cluster_ids = [c for c in sorted(clusters,
                           key=lambda c: int(c) if c.isdigit() else c)
                           if len(clusters[c]) >= args.min_members]
    cl_cache, row_owner = {}, {}
    row_hits = defaultdict(dict)          # db row -> {cid: n accessions shared}
    for cid in ordered_cluster_ids:
        cl = aggregate_cluster(cid, clusters[cid], centroids.get(cid, ""), meta)
        cl_cache[cid] = cl
        cl_accs = {acc_key(a) for a in split_accessions(cl["accessions"])}
        for i, row in enumerate(db_rows):
            shared = cl_accs & {acc_key(a) for a in
                                split_accessions(row.get("18S Accession #", ""))}
            if shared:
                row_hits[i][cid] = len(shared)

    for i, hits in row_hits.items():
        best = max(hits.values())
        tied = [c for c, n in hits.items() if n == best]
        if len(tied) == 1:
            row_owner[i] = tied[0]
        else:
            # tie: the cluster holding the accession listed first in the cell
            order = split_accessions(db_rows[i].get("18S Accession #", ""))
            pick = None
            for a in order:
                for c in tied:
                    if acc_key(a) in {acc_key(x) for x in
                                      split_accessions(cl_cache[c]["accessions"])}:
                        pick = c; break
                if pick: break
            row_owner[i] = pick or tied[0]

    n_tied = sum(1 for i, h in row_hits.items() if len(h) > 1)
    if n_tied:
        sys.stderr.write(
            f"{n_tied} database row(s) have accessions in more than one cluster; "
            f"curated attributes were given to the cluster holding the most of "
            f"them\n")

    ## ---- PASS 2: build the outputs -----------------------------------------
    for cid in ordered_cluster_ids:
        members = clusters[cid]
        cl = cl_cache[cid]
        status, matched_by, idxs = match_cluster(
            cl, db_acc_set, db_by_name, db_by_abbrev, db_by_epithet, db_rows,
            args.fuzzy_cutoff)
        if status == "probable" and args.probable_as_new:
            status, matched_by, idxs = "new", "", []

        cl["matched_species"] = "; ".join(
            dedupe(db_rows[i].get("Species Name", "") for i in idxs))

        # ---- merge curated database values into the cluster aggregate -----
        # 03_CLUSTER_HOST_ATTRIBUTES.tsv is what every downstream step reads
        # (tree symbols, environment rings, host taxonomy). Built from GenBank
        # alone it silently drops hosts that are only in the curated database,
        # so for a matched cluster those are merged in here. Provenance is kept
        # in separate columns so the two sources stay distinguishable.
        ## only the rows this cluster OWNS contribute curated attributes
        owned = [i for i in idxs if row_owner.get(i) == cid]
        db_hosts = db_locs = []
        if owned and args.merge_db_hosts:
            db_hosts = db_values(db_rows, owned, "Natural Host(s)")
            db_locs = db_values(db_rows, owned, "Locality")
        cl["owns_db_rows"] = "; ".join(
            dedupe(db_rows[i].get("Species Name", "") for i in owned))
        if owned:
            basis = "genbank + curated (owns the database record)"
        elif idxs:
            claimed = [i for i in idxs if i in row_owner]
            basis = ("genbank only (another cluster holds the database "
                     "accessions)" if claimed else
                     "genbank only (matched by name; the database record's "
                     "accessions are in no cluster)")
        else:
            basis = "genbank only (no database match)"
        cl["attribute_basis"] = basis
        gb_hosts, gb_locs = cl["hosts"], cl["localities"]
        merged_hosts, n_host_only_db = merge_sources(gb_hosts, "; ".join(db_hosts))
        merged_locs, n_loc_only_db = merge_sources(gb_locs, "; ".join(db_locs))
        cl["hosts_merged"] = "; ".join(merged_hosts)
        cl["localities_merged"] = "; ".join(merged_locs)
        cl["hosts_from_genbank"] = gb_hosts
        cl["hosts_from_database"] = "; ".join(db_hosts)
        cl["n_hosts_merged"] = len(merged_hosts)
        cl["n_hosts_only_in_database"] = n_host_only_db
        if n_host_only_db or n_loc_only_db:
            n_augmented += 1
            n_hosts_added += n_host_only_db
        # Classification must run AFTER the merge and on the MERGED host list.
        # Running it first, on GenBank hosts only, made a cluster whose hosts
        # are all curated come out "insufficient" -- yet the same cluster
        # passed the --min-hosts FASTA filter, which does use the merged count.
        # The tree then contained tips the status ring could not classify.
        kind, why = classify_cluster(cl, args.require_host_for_provisional)

        cl["host_source"] = (
            "genbank+database" if (gb_hosts and db_hosts) else
            "database only" if db_hosts else
            "genbank only" if gb_hosts else "none")

        display = (cl["species_name"] if kind == "named"
                   else provisional_label(cl, args.provisional_prefix))

        # entity_type is the tree's status ring. It must agree with the
        # canonical named/provisional rule applied to the display name, or the
        # tree would disagree with every other figure. name_rank() answers a
        # different question (is this GenBank string a usable binomial?), so
        # the two are checked against each other rather than merged.
        if kind in ("named", "provisional") and \
                classify_microsporidia(display) != kind:
            sys.stderr.write(
                f"WARNING cluster {cid}: entity_type={kind} but the canonical "
                f"classifier says {classify_microsporidia(display)} for "
                f"{display!r}\n")

        # ---- decide the disposition -------------------------------------
        if status == "existing":
            disposition = "already in database"
        elif status == "probable":
            disposition = "probable match - review"
        elif kind == "insufficient":
            disposition = "insufficient data"
        elif kind == "provisional" and args.no_provisional:
            disposition = "provisional - excluded by flag"
        elif cl["n_binomials"] > 1:
            disposition = "conflicting names - review"
        else:
            disposition = "ADD TO DATABASE"
        counts[disposition] += 1

        # ---- OUTPUT 3: attributes for every cluster ----------------------
        attributes.append({
            "cluster_label": display,
            "cluster_id": cid,
            "entity_type": kind,
            "in_attribute_database": "yes" if status == "existing" else "no",
            "disposition": disposition,
            "centroid_accession": cl["centroid"],
            "n_sequences": cl["n_members"],
            "n_sequences_with_host": cl["n_with_host"],
            "n_distinct_hosts": cl["n_hosts_merged"],
            "all_hosts": cl["hosts_merged"],
            "host_source": cl["host_source"],
            "attribute_basis": cl["attribute_basis"],
            "owns_database_record": cl["owns_db_rows"],
            "hosts_from_genbank": cl["hosts_from_genbank"],
            "hosts_from_database": cl["hosts_from_database"],
            "n_hosts_only_in_database": cl["n_hosts_only_in_database"],
            "n_distinct_countries": cl["n_countries"],
            "all_countries": cl["countries"],
            "all_localities": cl["localities_merged"],
            "all_accessions": cl["accessions"],
            "genbank_organism_names": cl["all_names"],
            "experimental_hosts": cl["experimental_hosts"],
            "isolation_sources": cl["isolation_sources"],
            "collection_years": cl["collection_years"],
            "matched_database_species": cl.get("matched_species", ""),
        })
        # Keep the FASTA in step with the classifier: a cluster that cannot be
        # classified has no place in the tree, or it appears as an
        # unclassifiable tip.
        if (cl["n_hosts_merged"] >= args.min_hosts and cl["centroid"]
                and kind != "insufficient"):
            keep_centroids[cl["centroid"]] = display

        tip_labels.append({
            "centroid_accession": cl["centroid"],
            "tree_tip_label": re.sub(r"[\s,:;()\[\]'\"|/\\]+", "_", display),
            "display_name": display,
            "entity_type": kind,
            "n_hosts": cl["n_hosts_merged"],
            "primary_host": (cl["hosts_merged"].split("; ")[0]
                             if cl["hosts_merged"] else ""),
            "all_hosts": cl["hosts_merged"],
        })

        # ---- OUTPUT 1 / 2 / review --------------------------------------
        if status in ("existing", "probable"):
            # ---- attribute each NEW host/locality to one database row -------
            db_h = {v.lower() for v in db_values(db_rows, idxs, "Natural Host(s)")}
            db_l = {v.lower() for v in db_values(db_rows, idxs, "Locality")}
            for rec in cl.get("per_record", []):
                new_h = [h for h in split_multi(rec["hosts"])
                         if h.lower() not in db_h]
                new_l = ([rec["locality"]]
                         if rec["locality"] and rec["locality"].lower() not in db_l
                         else [])
                if not new_h and not new_l:
                    continue
                # the accession that reported this host: which database row is it?
                own = acc_to_rows.get(
                    acc_key((split_accessions(rec["accession"]) or [""])[0]), [])
                if len(own) == 1:
                    conf, target = "high (accession is in this row)", own
                elif len(own) > 1:
                    conf, target = f"ambiguous ({len(own)} rows list this accession)", own
                elif len(idxs) == 1:
                    conf, target = "high (cluster matches one row)", idxs
                else:
                    # fall back to the database row whose species name matches
                    # the organism on THIS record
                    k = norm_name(rec["genbank_organism"])
                    byname = [i for i in idxs
                              if norm_name(db_rows[i].get("Species Name", "")) == k]
                    if len(byname) == 1:
                        conf, target = "medium (organism name matches row)", byname
                    else:
                        conf, target = "UNRESOLVED - assign by hand", idxs
                attrib_rows.append({
                    "cluster_id": cid,
                    "cluster_label": display,
                    "accession": rec["accession"],
                    "genbank_organism": rec["genbank_organism"],
                    "new_hosts": "; ".join(new_h),
                    "new_locality": "; ".join(new_l),
                    "confidence": conf,
                    "assign_to_species": "; ".join(dedupe(
                        db_rows[i].get("Species Name", "") for i in target)),
                    "n_candidate_rows": len(target),
                })

            already.append({
                "cluster_label": display,
                "cluster_id": cid,
                "match_confidence": ("confirmed" if status == "existing"
                                     else "PROBABLE - CHECK THIS"),
                "matched_how": matched_by,
                "database_species_matched": cl["matched_species"],
                "n_sequences": cl["n_members"],
                "NEW_hosts_not_in_database":
                    "; ".join(novel_values(cl["hosts"], db_rows, idxs,
                                           "Natural Host(s)")),
                "NEW_localities_not_in_database":
                    "; ".join(novel_values(cl["localities"], db_rows, idxs,
                                           "Locality")),
                "NEW_accessions_not_in_database": "; ".join(
                    a for a in split_accessions(cl["accessions"])
                    # compare on acc_key(): db_acc_set holds underscore-stripped
                    # keys, so testing the raw accession reported every RefSeq
                    # accession (XR_...) as new even when the database had it
                    if acc_key(a) not in db_acc_set),
                "all_hosts_in_cluster": cl["hosts"],
                "all_accessions_in_cluster": cl["accessions"],
            })
        elif disposition == "ADD TO DATABASE":
            new_row = build_db_row(cl, timestamp, args.year_from, kind, display)
            to_add.append(new_row)
            date_basis[cl.get("_date_basis", "none")] += 1
            if not new_row["Date Identified (year)"] and args.year_from != "none":
                undated.append({
                    "cluster_label": display,
                    "cluster_id": cid,
                    "n_sequences": cl["n_members"],
                    "accessions": cl["accessions"],
                    "why": "no publication, submission, collection or record "
                           "date on any accession in the cluster",
                })

        needs_review = (disposition.endswith("review")
                        or kind == "insufficient"
                        or cl["n_binomials"] > 1)
        if needs_review:
            if cl["n_binomials"] > 1 and not why and \
                    disposition == "already in database":
                why = (f"matched database as {cl['matched_species']!r}, but the "
                       f"cluster holds {cl['n_binomials']} distinct species "
                       f"({cl['distinct_species']}) -- possible lumping, "
                       f"chimera, or a misannotated accession")
            flagged.append({
                "cluster_label": display,
                "cluster_id": cid,
                "why_flagged": why or {
                    "probable match - review":
                        f"heuristic name match ({matched_by}) to "
                        f"{cl['matched_species']} -- confirm before merging",
                    "conflicting names - review":
                        f"{cl['n_binomials']} distinct species in one "
                        f"cluster: {cl['distinct_species']}",
                }.get(disposition, disposition),
                "n_sequences": cl["n_members"],
                "distinct_species_in_cluster": cl["distinct_species"],
                "genbank_organism_names_verbatim": cl["all_names"],
                "hosts": cl["hosts"],
                "localities": cl["localities"],
                "accessions": cl["accessions"],
            })

    # ---- write ----------------------------------------------------------
    if already:
        write_tsv(out("01_ALREADY_IN_DATABASE.tsv"), list(already[0].keys()), already)
    write_tsv(out("02_NEW_ENTRIES_TO_ADD.tsv"), DB_COLUMNS, to_add)
    n_prov = sum(1 for r in to_add if "PROVISIONAL" in r["Important Remarks"])
    write_xlsx(out("02_NEW_ENTRIES_TO_ADD.xlsx"), to_add, [
        "OUTPUT 2: new entries to add to the attribute database.",
        "",
        "Columns match your masterlist exactly -- copy rows straight in.",
        "One row per SEQUENCE CLUSTER, never one per strain or accession.",
        "",
        f"{len(to_add) - n_prov} named species (cluster carries a binomial).",
        f"{n_prov} provisional species (no binomial, but has 18S + a host);",
        "  these are named like 'Nosema sp. MSP-C0042' and are marked",
        "  PROVISIONAL SPECIES in the Important Remarks column.",
        "",
        "Blank columns are spore attributes that must come from the",
        "literature -- GenBank does not carry them.",
    ])
    if attributes:
        write_tsv(out("03_CLUSTER_HOST_ATTRIBUTES.tsv"),
                  list(attributes[0].keys()), attributes)
        write_tsv(out("04_TREE_TIP_LABELS.tsv"),
                  list(tip_labels[0].keys()), tip_labels)
    if flagged:
        write_tsv(out("05_FLAGGED_FOR_REVIEW.tsv"), list(flagged[0].keys()), flagged)

    # ---- OUTPUT 6: new host records, attributed to a single database row ----
    # 01_ALREADY_IN_DATABASE.tsv reports new hosts per CLUSTER, which is
    # ambiguous when a cluster spans several database species. Here each new
    # host is traced back to the accession that reported it, and that
    # accession is matched to one database row -- so the record can be filed
    # against the right species instead of being duplicated across all of them.
    if attrib_rows:
        write_tsv(out("06_NEW_RECORDS_BY_ACCESSION.tsv"),
                  list(attrib_rows[0].keys()), attrib_rows)

    # ---- OUTPUT 7: centroids of clusters that have a host --------------------
    n_fasta = 0
    if args.centroids:
        if not os.path.isfile(args.centroids):
            sys.stderr.write(
                f"\nWARNING: --centroids file not found: {args.centroids}\n"
                "         07_CENTROIDS_WITH_HOSTS.fasta not written\n")
        else:
            dest = out("07_CENTROIDS_WITH_HOSTS.fasta")
            n_fasta, matched = write_filtered_centroids(
                args.centroids, dest, keep_centroids, args.fasta_label)
            missing = set(keep_centroids) - matched
            sys.stderr.write(
                f"\nFiltered centroid FASTA: {n_fasta} of "
                f"{len(keep_centroids)} host-bearing clusters written\n")
            if missing:
                sys.stderr.write(
                    f"  WARNING: {len(missing)} centroid(s) not found in "
                    f"{args.centroids}\n"
                    f"           e.g. {', '.join(sorted(missing)[:5])}\n"
                    "           (check for vsearch 'NNN|' prefixes or version "
                    "suffixes)\n")
            write_tsv(out("07_CENTROIDS_WITH_HOSTS.tsv"),
                      ["centroid_accession", "cluster_label", "in_fasta"],
                      [{"centroid_accession": a, "cluster_label": n,
                        "in_fasta": "yes" if a in matched else "no"}
                       for a, n in sorted(keep_centroids.items())])

    # ---- report ---------------------------------------------------------
    if undated:
        write_tsv(out("06_NO_DATE.tsv"), list(undated[0].keys()), undated)

    n_named = sum(1 for a in attributes if a["entity_type"] == "named")
    n_prov_all = sum(1 for a in attributes if a["entity_type"] == "provisional")
    n_aug = sum(1 for r in already if r["NEW_hosts_not_in_database"]
                or r["NEW_localities_not_in_database"])
    if args.merge_db_hosts:
        sys.stderr.write(
            f"Merged curated database hosts into {n_augmented} cluster(s); "
            f"{n_hosts_added} host record(s) that GenBank lacks are now in "
            f"{out('03_CLUSTER_HOST_ATTRIBUTES.tsv')}\n")
    else:
        sys.stderr.write(
            "NOTE: --no-merge-database-hosts set; output 3 is GenBank-only\n")
    sys.stderr.write(
        f"{len(attributes)} clusters: {n_named} named, {n_prov_all} provisional, "
        f"{len(attributes) - n_named - n_prov_all} insufficient data\n\n"
        f"OUTPUT 1  {out('01_ALREADY_IN_DATABASE.tsv')}\n"
        f"          {len(already)} clusters already in your database; "
        f"{n_aug} carry hosts or localities it lacks\n\n"
        f"OUTPUT 2  {out('02_NEW_ENTRIES_TO_ADD.xlsx')}\n"
        f"          {len(to_add)} rows to add "
        f"({len(to_add) - n_prov} named, {n_prov} provisional)\n\n"
        f"OUTPUT 3  {out('03_CLUSTER_HOST_ATTRIBUTES.tsv')}\n"
        f"          {len(attributes)} clusters with all hosts, for tree mapping\n"
        f"          (tip labels: {out('04_TREE_TIP_LABELS.tsv')})\n\n")
    if args.centroids:
        n_drop = len(attributes) - len(keep_centroids)
        sys.stderr.write(
            f"OUTPUT 7  {out('07_CENTROIDS_WITH_HOSTS.fasta')}\n"
            f"          {n_fasta} centroid sequences "
            f"(clusters with >= {args.min_hosts} host); "
            f"{n_drop} host-less cluster(s) excluded\n\n")
    if flagged:
        sys.stderr.write(f"REVIEW    {out('05_FLAGGED_FOR_REVIEW.tsv')} "
                         f"({len(flagged)} clusters)\n\n")
    if args.year_from != "none":
        labels = {k: lab for k, _f, lab in YEAR_SOURCES}
        labels["none"] = "NO DATE AVAILABLE"
        sys.stderr.write("Date Identified (year) filled from:\n")
        for k, v in date_basis.most_common():
            sys.stderr.write(f"  {v:>6}  {labels.get(k, k)}\n")
        if date_basis.get("none"):
            sys.stderr.write(f"          -> listed in {out('06_NO_DATE.tsv')}\n")
        sys.stderr.write(
            "  (none of these is a species-description year; the proxy used is\n"
            "   recorded per row as 'date basis: ...' in Important Remarks)\n\n")

    sys.stderr.write("Breakdown:\n")
    for k, v in counts.most_common():
        sys.stderr.write(f"  {v:>6}  {k}\n")


if __name__ == "__main__":
    main()
