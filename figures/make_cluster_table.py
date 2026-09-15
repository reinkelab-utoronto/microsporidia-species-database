#!/usr/bin/env python3
"""
Build a readable supplemental table describing the 18S sequence clusters.

03_CLUSTER_HOST_ATTRIBUTES.tsv holds everything the pipeline needs, including
diagnostic columns that only matter while debugging. This selects and renames
the columns a reader wants, drops clusters that never entered the analysis, and
sorts by clade then name.

Usage:
  python make_cluster_table.py --attrs cluster_output/03_CLUSTER_HOST_ATTRIBUTES.tsv \\
      --assignments clade_assignments.csv --out TableS2_clusters.tsv
"""
import argparse, csv, os, sys

COLS = [
    ("cluster_label",        "Cluster / species"),
    ("clade",                "Clade"),
    ("entity_type",          "Status"),
    ("in_attribute_database","In attribute database"),
    ("n_sequences",          "Sequences in cluster"),
    ("centroid_accession",   "Representative accession"),
    ("n_distinct_hosts",     "Hosts (n)"),
    ("all_hosts",            "Hosts"),
    ("host_source",          "Host record source"),
    ("attribute_basis",      "Attribute basis"),
    ("owns_database_record", "Curated record inherited from"),
    ("all_countries",        "Countries"),
    ("collection_years",     "Collection year(s)"),
    ("all_accessions",       "All accessions"),
]

p = argparse.ArgumentParser()
p.add_argument("--attrs", default="cluster_output/03_CLUSTER_HOST_ATTRIBUTES.tsv")
p.add_argument("--assignments", help="clade_assignments.csv from the tree script, "
                                     "to add the clade column")
p.add_argument("--out", default="TableS2_clusters.tsv")
p.add_argument("--only-in-tree", action="store_true",
               help="keep only clusters that have at least one host, i.e. those "
                    "that entered the tree")
a = p.parse_args()

rows = list(csv.DictReader(open(a.attrs, newline="", encoding="utf-8"), delimiter="\t"))
print(f"{len(rows)} clusters in {a.attrs}")

clade = {}
if a.assignments and os.path.exists(a.assignments):
    for r in csv.DictReader(open(a.assignments, newline="", encoding="utf-8")):
        k = r.get("tip_label") or r.get("centroid_accession")
        if k and r.get("clade"): clade[k] = r["clade"]
    print(f"clade assignments for {len(clade)} tips")

out = []
for r in rows:
    if a.only_in_tree and r.get("n_distinct_hosts", "0") in ("", "0"):
        continue
    r["clade"] = clade.get(r.get("centroid_accession", ""), "")
    out.append({new: r.get(old, "") for old, new in COLS})

out.sort(key=lambda x: (x["Clade"] or "zzz", x["Cluster / species"]))
with open(a.out, "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=[n for _, n in COLS], delimiter="\t",
                       lineterminator="\n")
    w.writeheader(); w.writerows(out)

print(f"wrote {len(out)} rows -> {a.out}")
dups = {}
for x in out: dups.setdefault(x["Cluster / species"], []).append(x)
rep = {k: v for k, v in dups.items() if len(v) > 1}
if rep:
    print(f"\n*** {len(rep)} name(s) appear on more than one cluster -- these are "
          f"species split by the clustering,\n    and if any were also added as new "
          f"database entries they will be duplicated there:")
    for k, v in list(rep.items())[:8]:
        print(f"    {k}  ({len(v)} clusters)")
