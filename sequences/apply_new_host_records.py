#!/usr/bin/env python3
"""
Apply the NEW_* columns from 01_ALREADY_IN_DATABASE.tsv back into the
attribute database.

cluster_to_database reports, for every cluster that matched an existing
database entry, the hosts / localities / accessions the cluster carries that
the database row does not. Those are new records for described species, and
until they are written back they exist only in the report.

This script does that write-back:

  * matches each report row to database rows BY ACCESSION (the same key the
    cluster script matched on), falling back to normalised species name
  * appends the missing values to the existing cell, de-duplicated, keeping
    the database's own spelling where both have a value
  * writes a NEW .xlsx -- the input file is never modified
  * writes a change log with one line per edited cell, so every change can be
    checked or reversed

Rows that match more than one database entry are NOT edited. That case is
genuinely ambiguous (which of the entries gains the host?) and guessing would
silently corrupt curated data, so they are written to a review file instead.

Usage
-----
  # see what would change, without writing
  python apply_new_host_records.py --db database.xlsx \\
      --report cluster_output/01_ALREADY_IN_DATABASE.tsv --dry-run

  # write the merged file
  python apply_new_host_records.py --db database.xlsx \\
      --report cluster_output/01_ALREADY_IN_DATABASE.tsv \\
      --out database_merged.xlsx
"""

import argparse
import csv
import os
import re
import sys
import unicodedata
from collections import defaultdict

ACC_RE = re.compile(r"^[A-Z]{1,6}_?\d{5,9}(\.\d+)?$")

# report column -> database column
FIELD_MAP = {
    "NEW_hosts_not_in_database": "Natural Host(s)",
    "NEW_localities_not_in_database": "Locality",
    "NEW_accessions_not_in_database": "18S Accession #",
}


def split_accessions(cell):
    if not cell:
        return []
    s = str(cell).upper()
    s = re.sub(r"\b([A-Z]{1,6})\s+(\d{5,9})\b", r"\1\2", s)
    out = []
    for tok in re.split(r"[,;/|\s]+", s):
        tok = tok.strip("()[].")
        base = tok.split(".")[0]
        if ACC_RE.match(tok) or ACC_RE.match(base + ".1"):
            out.append(base)
    return out


def norm_name(name):
    if not name:
        return ""
    s = unicodedata.normalize("NFKD", str(name))
    s = "".join(c for c in s if not unicodedata.combining(c)).lower()
    s = re.sub(r"\(.*?\)", " ", s)
    s = re.sub(r"\b(n|sp|spp|gen|comb|stat|nov|et|cf|aff|var|subsp)\.?\b", " ", s)
    s = re.sub(r"[^a-z0-9\s]", " ", s)
    return " ".join(s.split()[:2])


def split_multi(cell):
    if not cell:
        return []
    parts = re.split(r"\s*;\s*|\n", str(cell))
    return [p.strip() for p in parts if p.strip()]


def append_values(existing, new_values):
    """Append new_values to existing, case-insensitively de-duplicated."""
    cur = split_multi(existing)
    seen = {v.lower() for v in cur}
    added = [v for v in new_values if v.lower() not in seen]
    if not added:
        return existing, []
    merged = cur + added
    return "; ".join(merged), added


def main():
    p = argparse.ArgumentParser(
        description="Write new host/locality/accession records back into the "
                    "attribute database.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--db", required=True, help="attribute database .xlsx")
    p.add_argument("--report", required=True,
                   help="01_ALREADY_IN_DATABASE.tsv from cluster_to_database")
    p.add_argument("--sheet", default="Actively Updated Masterlist")
    p.add_argument("--out", help="output .xlsx (default: <db>_merged.xlsx)")
    p.add_argument("--log", help="change log TSV (default: <out>.changes.tsv)")
    p.add_argument("--fields", nargs="+", default=list(FIELD_MAP),
                   help="which NEW_* columns to apply")
    p.add_argument("--dry-run", action="store_true",
                   help="report what would change; write nothing")
    p.add_argument("--allow-multi-match", action="store_true",
                   help="apply to EVERY matching row when a report row matches "
                        "more than one (default: skip and flag for review)")
    args = p.parse_args()

    try:
        from openpyxl import load_workbook
    except ImportError:
        sys.exit("openpyxl required: pip install openpyxl")

    ## default output goes to the CURRENT directory, not beside the input --
    ## the input may sit on a read-only mount
    out_path = args.out or (
        os.path.basename(re.sub(r"\.xlsx$", "", args.db)) + "_merged.xlsx")
    log_path = args.log or re.sub(r"\.xlsx$", "", out_path) + ".changes.tsv"

    ## load the workbook with openpyxl (not pandas) so every other sheet,
    ## column, formula and bit of formatting survives untouched
    wb = load_workbook(args.db)
    if args.sheet not in wb.sheetnames:
        sys.exit(f"sheet {args.sheet!r} not in {args.db} "
                 f"(have: {', '.join(wb.sheetnames)})")
    ws = wb[args.sheet]

    header = [c.value for c in ws[1]]
    col_of = {h: i + 1 for i, h in enumerate(header) if h}
    for f in args.fields:
        dbcol = FIELD_MAP[f]
        if dbcol not in col_of:
            sys.exit(f"database has no column {dbcol!r} (needed for {f})")
    acc_col = col_of.get("18S Accession #")
    sp_col = col_of.get("Species Name")
    if not acc_col or not sp_col:
        sys.exit("database needs both 'Species Name' and '18S Accession #'")

    ## index the database rows
    by_acc, by_name = defaultdict(list), defaultdict(list)
    for r in range(2, ws.max_row + 1):
        for a in split_accessions(ws.cell(r, acc_col).value):
            by_acc[a].append(r)
        n = norm_name(ws.cell(r, sp_col).value)
        if n:
            by_name[n].append(r)
    sys.stderr.write(f"database: {ws.max_row - 1} rows, {len(by_acc)} accessions\n")

    changes, review = [], []
    n_rows_edited = set()

    with open(args.report, newline="", encoding="utf-8") as fh:
        for rep in csv.DictReader(fh, delimiter="\t"):
            pending = {f: split_multi(rep.get(f, "")) for f in args.fields}
            if not any(pending.values()):
                continue

            ## match on accession first -- the same key the cluster script used
            rows, how = set(), ""
            for a in split_accessions(rep.get("all_accessions_in_cluster", "")):
                rows.update(by_acc.get(a, []))
            if rows:
                how = "accession"
            else:
                for nm in split_multi(rep.get("database_species_matched", "")):
                    rows.update(by_name.get(norm_name(nm), []))
                if rows:
                    how = "species name"

            if not rows:
                review.append({**rep, "problem": "no matching database row"})
                continue
            if len(rows) > 1 and not args.allow_multi_match:
                review.append({
                    **rep, "problem": f"matches {len(rows)} database rows "
                    f"(rows {', '.join(map(str, sorted(rows)))}) -- "
                    "ambiguous which entry gains the record"})
                continue

            for row in sorted(rows):
                for f, vals in pending.items():
                    if not vals:
                        continue
                    c = col_of[FIELD_MAP[f]]
                    before = ws.cell(row, c).value or ""
                    after, added = append_values(before, vals)
                    if not added:
                        continue
                    changes.append({
                        "excel_row": row,
                        "species_name": ws.cell(row, sp_col).value,
                        "column": FIELD_MAP[f],
                        "matched_by": how,
                        "added": "; ".join(added),
                        "before": str(before)[:200],
                        "after": str(after)[:200],
                    })
                    n_rows_edited.add(row)
                    if not args.dry_run:
                        ws.cell(row, c).value = after

    ## ---- report -----------------------------------------------------------
    sys.stderr.write(
        f"\n{len(changes)} cell edit(s) across {len(n_rows_edited)} row(s)\n")
    by_col = defaultdict(int)
    for c in changes:
        by_col[c["column"]] += 1
    for k, v in sorted(by_col.items()):
        sys.stderr.write(f"  {v:>5}  {k}\n")
    if review:
        sys.stderr.write(f"{len(review)} report row(s) NOT applied "
                         "(see the review file)\n")

    if changes:
        with open(log_path, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(changes[0].keys()),
                               delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(changes)
    if review:
        rp = re.sub(r"\.xlsx$", "", out_path) + ".needs_review.tsv"
        with open(rp, "w", newline="", encoding="utf-8") as fh:
            cols = list(review[0].keys())
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                               extrasaction="ignore", lineterminator="\n")
            w.writeheader()
            w.writerows(review)
        sys.stderr.write(f"  review -> {rp}\n")

    if args.dry_run:
        sys.stderr.write("\nDRY RUN -- nothing written. "
                         f"Change log: {log_path if changes else '(none)'}\n")
        return

    wb.save(out_path)
    sys.stderr.write(f"\nWrote {out_path}\n")
    if changes:
        sys.stderr.write(f"Change log: {log_path}\n")
    sys.stderr.write("The input file was not modified.\n")


if __name__ == "__main__":
    main()
