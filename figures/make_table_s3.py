#!/usr/bin/env python3
"""Build Table S3: per-host taxonomy and environment, COMBINING the taxonomic
databases (host_taxonomy_environment.py) with the genus-level habitat
predictions (predict_host_environment.py) -- the table the S3 legend promises.

Neither pipeline writes this combined view on its own: 02_host_lookup.tsv is
databases-only and 02_host_predicted.tsv is genus-table-only. The combination
rule here is the SAME one the figures use (api_first, as in
fig_environment_venn's source_rule and the tree's HOST_ENV_RULE): the
database-derived environment where one exists, the genus-table prediction
otherwise.

Usage:
  python make_table_s3.py \
      --lookup host_taxonomy_output/02_host_lookup.tsv \
      --pred   predicted_environment/02_host_predicted.tsv \
      --out    TableS3_host_taxonomy_environment.tsv
"""
import argparse
import csv
import sys

# (source, source column, published column label)
COLS = [
    ("L", "host_name",            "Host"),
    ("L", "rank",                 "Name rank"),
    ("L", "kingdom",              "Kingdom"),
    ("L", "phylum_normalised",    "Phylum"),
    ("L", "phylum",               "Phylum (as returned)"),
    ("L", "class",                "Class"),
    ("L", "order",                "Order"),
    ("L", "family",               "Family"),
    ("L", "genus",                "Genus"),
    ("L", "species",              "Species"),
    ("L", "taxonomy_source",      "Taxonomy source"),
    ("L", "n_sources_found",      "Sources agreeing (n)"),
    ("L", "phylum_agreement",     "Phylum agreement"),
    ("L", "phylum_dispute",       "Phylum dispute"),
    ("L", "flag",                 "Taxonomy flag"),
    ("L", "environment",          "Environment (taxonomic databases)"),
    ("L", "environment_source",   "Environment source (databases)"),
    ("P", "environment_all",      "Environment (genus table)"),
    ("P", "environment_primary",  "Primary environment (genus table)"),
    ("P", "confidence",           "Genus-table confidence"),
    ("P", "basis",                "Genus-table basis"),
    ("P", "life_stage_note",      "Life-stage note"),
    ("C", "env_combined",         "Environment used in analyses"),
    ("C", "env_used_source",      "Environment used from"),
]


def read_tsv(path, key="host_name"):
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    out = {}
    for r in rows:
        out.setdefault(r[key].strip(), r)
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--lookup", default="host_taxonomy_output/02_host_lookup.tsv")
    p.add_argument("--pred", default="predicted_environment/02_host_predicted.tsv")
    p.add_argument("--out", default="TableS3_host_taxonomy_environment.tsv")
    a = p.parse_args()

    look = read_tsv(a.lookup)
    pred = read_tsv(a.pred)
    # description-only fragments parse to an empty host name; they are kept in
    # the per-pair tables but have no place in a per-host supplemental table
    look.pop("", None); pred.pop("", None)
    hosts = sorted(set(look) | set(pred), key=str.lower)

    only_l = [h for h in hosts if h not in pred]
    only_p = [h for h in hosts if h not in look]
    if only_l or only_p:
        sys.stderr.write(
            f"NOTE: host lists differ between the two pipelines: "
            f"{len(only_l)} only in the database lookup, {len(only_p)} only in "
            f"the prediction. The parsers are meant to be identical -- "
            f"inspect these:\n")
        for h in (only_l + only_p)[:10]:
            sys.stderr.write(f"    {h}\n")

    out_rows, n_db, n_pred, n_none = [], 0, 0, 0
    for h in hosts:
        L = look.get(h, {})
        P = pred.get(h, {})
        db_env = (L.get("environment") or "").strip()
        pr_env = (P.get("environment_all") or "").strip()
        if db_env:
            comb, src = db_env, "taxonomic databases"
            n_db += 1
        elif pr_env:
            comb, src = pr_env, "genus table"
            n_pred += 1
        else:
            comb, src = "", ""
            n_none += 1
        C = {"env_combined": comb, "env_used_source": src}
        row = {}
        for where, col, label in COLS:
            src_d = {"L": L, "P": P, "C": C}[where]
            row[label] = src_d.get(col, "")
        row["Host"] = h
        out_rows.append(row)

    labels = [label for _, _, label in COLS]
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=labels, delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(out_rows)
    sys.stderr.write(
        f"wrote {a.out}: {len(out_rows)} hosts "
        f"({n_db} environment from the databases, {n_pred} from the genus "
        f"table, {n_none} unplaced)\n")


if __name__ == "__main__":
    main()
