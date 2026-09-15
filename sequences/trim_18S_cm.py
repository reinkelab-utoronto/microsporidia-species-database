#!/usr/bin/env python3
"""
Trim 18S sequences to the region matching a covariance model (default: the
Rfam microsporidia SSU model RF02542, 1312 columns).

Differences from the simple version:
  * safe with duplicate / awkward sequence IDs (sequences are renamed to
    serial IDs for the cmsearch run and mapped back afterwards)
  * preserves the original FASTA description, appending trim annotation
  * filters on MODEL COVERAGE as well as trimmed sequence length, and writes
    the model coordinates into the header so downstream clustering can
    require overlap in a common coordinate system
  * applies a significance threshold (inclusion flag or E-value)
  * chains colinear hits instead of keeping only the single best one
  * writes a TSV log of every input sequence and why it was kept or dropped

Usage
-----
  python trim_18S_cm.py in.fasta out.fasta --cm_file RF02542.cm \
      --min_model_cov 300 --min_length 250 --cpu 8

  # then, e.g., keep only sequences overlapping the V4 region
  # (model coords roughly 550-950 for RF02542 -- check against your own data)
"""

import argparse
import csv
import os
import shutil
import subprocess
import sys
import tempfile

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# tblout column indices (Infernal 1.1.x cmsearch --tblout), verified against 1.1.5:
# 0 target  1 acc  2 query  3 acc  4 mdl  5 mdl_from  6 mdl_to  7 seq_from
# 8 seq_to  9 strand  10 trunc  11 pass  12 gc  13 bias  14 score  15 E-value
# 16 inc  17+ description
C_TARGET, C_MDL_FROM, C_MDL_TO = 0, 5, 6
C_SEQ_FROM, C_SEQ_TO, C_STRAND = 7, 8, 9
C_TRUNC, C_SCORE, C_EVALUE, C_INC = 10, 14, 15, 16


def run_cmsearch(cm_file, fasta, tblout, cpu, extra_args):
    cmd = ["cmsearch", "--tblout", tblout, "--noali", "-o", os.devnull]
    if cpu:
        cmd += ["--cpu", str(cpu)]
    cmd += list(extra_args) + [cm_file, fasta]
    sys.stderr.write("Running: " + " ".join(cmd) + "\n")
    proc = subprocess.run(cmd, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        sys.exit(f"cmsearch failed (exit {proc.returncode}):\n{proc.stderr}")


def parse_tblout(path, max_evalue, require_inc):
    """Return {target: [hit, ...]} for hits passing the significance filter."""
    hits = {}
    n_seen = n_kept = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) <= C_INC:
                continue
            n_seen += 1

            try:
                evalue = float(parts[C_EVALUE])
                score = float(parts[C_SCORE])
            except ValueError:
                continue

            if require_inc and parts[C_INC] != "!":
                continue
            if max_evalue is not None and evalue > max_evalue:
                continue
            n_kept += 1

            seq_from, seq_to = int(parts[C_SEQ_FROM]), int(parts[C_SEQ_TO])
            hits.setdefault(parts[C_TARGET], []).append({
                "mdl_from": int(parts[C_MDL_FROM]),
                "mdl_to": int(parts[C_MDL_TO]),
                "start": min(seq_from, seq_to),   # 1-based inclusive
                "end": max(seq_from, seq_to),
                "strand": parts[C_STRAND],
                "trunc": parts[C_TRUNC],
                "score": score,
                "evalue": evalue,
            })
    sys.stderr.write(
        f"cmsearch reported {n_seen} hits; {n_kept} passed the significance filter\n"
    )
    return hits


def chain_hits(hit_list, max_gap, overlap_tol=50):
    """
    Group hits into colinear chains and return the highest-scoring chain.

    Two hits chain if they are on the same strand, adjacent in the sequence
    (gap <= max_gap), and genuinely ADVANCE through the model -- the next hit
    must start at or beyond the end of the previous one (allowing overlap_tol
    positions of overlap). Requiring advancement rather than mere ordering is
    what stops two tandem rRNA copies on one contig from being merged into a
    single span.
    """
    best = None
    for strand in ("+", "-"):
        same = sorted([h for h in hit_list if h["strand"] == strand],
                      key=lambda h: h["start"])
        chain = []
        for h in same:
            if not chain:
                chain = [h]
                continue
            prev = chain[-1]
            gap = h["start"] - prev["end"] - 1
            if strand == "+":
                advances = h["mdl_from"] >= prev["mdl_to"] - overlap_tol
            else:
                advances = h["mdl_to"] <= prev["mdl_from"] + overlap_tol
            if gap <= max_gap and advances:
                chain.append(h)
            else:
                best = better(best, chain)
                chain = [h]
        best = better(best, chain)
    return best


def better(best, chain):
    if not chain:
        return best
    total = sum(h["score"] for h in chain)
    if best is None or total > sum(h["score"] for h in best):
        return chain
    return best


def model_coverage(chain):
    """Number of distinct model positions covered by a chain of hits."""
    spans = sorted((min(h["mdl_from"], h["mdl_to"]),
                    max(h["mdl_from"], h["mdl_to"])) for h in chain)
    covered, cur_s, cur_e = 0, None, None
    for s, e in spans:
        if cur_s is None:
            cur_s, cur_e = s, e
        elif s <= cur_e + 1:
            cur_e = max(cur_e, e)
        else:
            covered += cur_e - cur_s + 1
            cur_s, cur_e = s, e
    if cur_s is not None:
        covered += cur_e - cur_s + 1
    return covered


def main():
    p = argparse.ArgumentParser(
        description="Trim 18S sequences to a covariance model region.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("input_fasta")
    p.add_argument("output_fasta")
    p.add_argument("--cm_file", default="RF02542.cm", help="path to the CM file")
    p.add_argument("-l", "--min_length", type=int, default=0,
                   help="minimum length of the trimmed sequence (bp)")
    p.add_argument("-c", "--min_model_cov", type=int, default=0,
                   help="minimum number of model positions covered "
                        "(RF02542 is 1312 columns long)")
    p.add_argument("--model_length", type=int, default=1312,
                   help="model length, used only for the coverage fraction "
                        "reported in the header and log")
    p.add_argument("--max_evalue", type=float, default=1e-5,
                   help="drop hits above this E-value")
    p.add_argument("--require_inc", action="store_true",
                   help="additionally require the inclusion flag '!'")
    p.add_argument("--max_gap", type=int, default=300,
                   help="max sequence gap (bp) between hits that still chain")
    p.add_argument("--cpu", type=int, default=None, help="threads for cmsearch")
    p.add_argument("--cmsearch_arg", action="append", default=[],
                   help="extra argument passed through to cmsearch "
                        "(repeatable, e.g. --cmsearch_arg=--cut_ga)")
    p.add_argument("--log", help="TSV log of every input sequence "
                                 "(default: <output>.trimlog.tsv)")
    p.add_argument("--no_description", action="store_true",
                   help="drop the original description instead of keeping it")
    p.add_argument("--line_width", type=int, default=60)
    args = p.parse_args()

    if not os.path.isfile(args.cm_file):
        sys.exit(f"Error: CM file {args.cm_file} not found.")
    if shutil.which("cmsearch") is None:
        sys.exit("Error: 'cmsearch' not found. Install Infernal and put it on your PATH.")

    log_path = args.log or args.output_fasta + ".trimlog.tsv"

    records = list(SeqIO.parse(args.input_fasta, "fasta"))
    if not records:
        sys.exit(f"No sequences found in {args.input_fasta}")

    dups = len(records) - len({r.id for r in records})
    if dups:
        sys.stderr.write(
            f"NOTE: {dups} duplicate sequence ID(s) present; handled safely "
            "via internal serial IDs, but you may want to deduplicate.\n"
        )

    with tempfile.TemporaryDirectory() as tmp:
        # rename to serial IDs so duplicate or awkward names can't collide
        tmp_fa = os.path.join(tmp, "input.fasta")
        with open(tmp_fa, "w") as fh:
            for i, rec in enumerate(records):
                fh.write(f">s{i}\n{str(rec.seq)}\n")

        tblout = os.path.join(tmp, "hits.tbl")
        run_cmsearch(args.cm_file, tmp_fa, tblout, args.cpu, args.cmsearch_arg)
        hits = parse_tblout(tblout, args.max_evalue, args.require_inc)

    kept, log_rows = [], []
    for i, rec in enumerate(records):
        key = f"s{i}"
        row = {"input_id": rec.id, "input_length": len(rec.seq), "status": "",
               "reason": "", "n_hits": len(hits.get(key, [])), "strand": "",
               "trunc": "", "seq_from": "", "seq_to": "", "trimmed_length": "",
               "mdl_from": "", "mdl_to": "", "model_cov": "", "score": "",
               "evalue": ""}

        if key not in hits:
            row["status"], row["reason"] = "dropped", "no significant CM hit"
            log_rows.append(row)
            continue

        chain = chain_hits(hits[key], args.max_gap)
        start = min(h["start"] for h in chain) - 1          # to 0-based
        end = max(h["end"] for h in chain)
        strand = chain[0]["strand"]
        mdl_lo = min(min(h["mdl_from"], h["mdl_to"]) for h in chain)
        mdl_hi = max(max(h["mdl_from"], h["mdl_to"]) for h in chain)
        cov = model_coverage(chain)

        subseq = rec.seq[start:end]
        if strand == "-":
            subseq = subseq.reverse_complement()

        row.update({
            "strand": strand, "trunc": chain[0]["trunc"],
            "seq_from": start + 1, "seq_to": end, "trimmed_length": len(subseq),
            "mdl_from": mdl_lo, "mdl_to": mdl_hi, "model_cov": cov,
            "score": f"{sum(h['score'] for h in chain):.1f}",
            "evalue": f"{min(h['evalue'] for h in chain):.2g}",
        })

        if len(subseq) < args.min_length:
            row["status"] = "dropped"
            row["reason"] = f"trimmed length {len(subseq)} < {args.min_length}"
            log_rows.append(row)
            continue
        if cov < args.min_model_cov:
            row["status"] = "dropped"
            row["reason"] = f"model coverage {cov} < {args.min_model_cov}"
            log_rows.append(row)
            continue

        tags = (f"[cm_from={mdl_lo}] [cm_to={mdl_hi}] "
                f"[cm_cov={cov}/{args.model_length}] [strand={strand}] "
                f"[trunc={chain[0]['trunc']}]")
        original = "" if args.no_description else \
            rec.description[len(rec.id):].strip()
        desc = (original + " " + tags).strip()

        kept.append(SeqRecord(subseq, id=rec.id, description=desc))
        row["status"], row["reason"] = "kept", ""
        log_rows.append(row)

    with open(args.output_fasta, "w") as out:
        for rec in kept:
            out.write(f">{rec.id} {rec.description}\n".replace(" \n", "\n"))
            s = str(rec.seq)
            if args.line_width > 0:
                for j in range(0, len(s), args.line_width):
                    out.write(s[j:j + args.line_width] + "\n")
            else:
                out.write(s + "\n")

    with open(log_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(log_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(log_rows)

    n_drop = len(records) - len(kept)
    sys.stderr.write(
        f"\nKept {len(kept)}/{len(records)} sequences -> {args.output_fasta}\n"
        f"Dropped {n_drop} (see {log_path} for per-sequence reasons)\n"
    )
    if kept:
        covs = sorted(int(r["model_cov"]) for r in log_rows if r["status"] == "kept")
        sys.stderr.write(
            f"Model coverage of kept sequences: min {covs[0]}, "
            f"median {covs[len(covs) // 2]}, max {covs[-1]} "
            f"of {args.model_length}\n"
        )


if __name__ == "__main__":
    main()
