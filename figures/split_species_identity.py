#!/usr/bin/env python3
"""
Species that fall into more than one 98 % cluster: how far apart are their
clusters?  (calculations only; fig_split_species_identity.R draws the figure)

For every species label that owns two or more clusters in Table S2, every
pair of its clusters is compared: each sequence of one cluster is aligned to
each sequence of the other and the identities are averaged. That gives one
between-cluster identity per cluster pair; the mean of those is the species
value. Within-cluster identity (all pairs of sequences sharing a cluster) and
each minor cluster's identity to the species' largest cluster come from the
same alignments.

WHICH SEQUENCES
  Sequences are assigned to clusters by accession, using the "All accessions"
  column of Table S2. The FASTA can be the full set of trimmed sequences that
  went into the clustering, or only the centroids. Every centroid is always
  used; additional members enter the averages only if they are at least
  --min-len bp. The floor matters because identity ignores terminal gaps: a
  700 bp amplicon scores ~99.7 % against two centroids that differ from each
  other only in the variable ends it never covers, so without a floor partial
  sequences vote "identical" on regions they have not seen. Cluster sizes
  (for labels and point sizes) still count every sequence in Table S2.

IDENTITY
  vsearch --allpairs_global with --iddef (default 2: matching columns /
  alignment columns, terminal gaps excluded), the same definition the
  clustering used, so the values sit on the same scale as the threshold.

OUTPUT (tab-separated, in --outdir, prefix --stem)
  <stem>_cluster_pairs.tsv   one row per cluster pair: mean/min/max identity,
                             sequences used on each side, sizes from Table S2
  <stem>_within_cluster.tsv  one row per cluster with >= 2 sequences used
  <stem>_vs_largest.tsv      each minor cluster against the species' largest
  <stem>_species.tsv         one row per species: the "average of averages"
  <stem>_clusters.tsv        one row per cluster: size, centroid, centroid
                             length, sequences used, short-centroid flag
  <stem>_matrix_<species>.tsv  cluster x cluster identity matrix for each
                             --detail species (diagonal = 100), for the
                             heatmap panels

Usage
  python3 split_species_identity.py TableS2_clusters_2.tsv trimmed.fasta
      [--min-len 1000] [--threshold 98] [--iddef 2] [--threads 2]
      [--detail "Vittaforma corneae" "Nosema bombycis" | --detail-n 2]
      [--outdir figure_output] [--stem fig_split_species_identity]
"""

import argparse
import csv
import os
import re
import statistics
import subprocess
import sys
import tempfile
from collections import defaultdict
from itertools import combinations


# ----------------------------------------------------------------- inputs --
def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("table", help="Table S2 (tab-separated, one row per cluster)")
    ap.add_argument("fasta", help="trimmed sequences (all members) or centroids")
    ap.add_argument("--min-len", type=int, default=1000,
                    help="members shorter than this are not used in identity "
                         "averages; centroids are always used (default 1000)")
    ap.add_argument("--threshold", type=float, default=98.0,
                    help="clustering identity, for the counts reported (default 98)")
    ap.add_argument("--iddef", type=int, default=2,
                    help="vsearch --iddef used for clustering (default 2)")
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--min-clusters", type=int, default=2,
                    help="a species needs at least this many clusters (default 2)")
    ap.add_argument("--detail", nargs="*", default=None,
                    help="species to write cluster x cluster matrices for")
    ap.add_argument("--detail-n", type=int, default=2,
                    help="if --detail is not given, the N species with the most "
                         "clusters (default 2)")
    ap.add_argument("--vsearch", default="vsearch")
    ap.add_argument("--outdir", default="figure_output")
    ap.add_argument("--stem", default="fig_split_species_identity")
    return ap.parse_args()


def strip_ver(acc):
    return re.sub(r"\.\d+$", "", acc)


def read_table(path):
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    need = ["Cluster / species", "Status", "Sequences in cluster",
            "Representative accession", "All accessions"]
    missing = [c for c in need if c not in rows[0]]
    if missing:
        sys.exit(f"Table S2 is missing column(s): {', '.join(missing)}")
    clusters = {}
    acc2cl = {}
    for r in rows:
        cl = r["Representative accession"].strip()
        if cl in clusters:
            sys.exit(f"representative accession {cl} occurs twice; cannot be a cluster id")
        clusters[cl] = dict(cluster=cl, species=r["Cluster / species"].strip(),
                            status=r["Status"].strip(),
                            n_seq=int(r["Sequences in cluster"]))
        for a in r["All accessions"].split(";"):
            a = strip_ver(a.strip())
            if a and a not in acc2cl:          # first listing wins
                acc2cl[a] = cl
    return clusters, acc2cl


def read_fasta(path):
    seqs = {}
    name = None
    buf = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(re.sub(r"[^A-Za-z]", "", line).upper())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


# --------------------------------------------------------------- vsearch --
def allpairs(vsearch, seqs, iddef, threads):
    """all unordered pairs -> dict[(a, b)] = identity (percent)"""
    with tempfile.TemporaryDirectory(prefix="allpairs_") as td:
        fa = os.path.join(td, "in.fasta")
        out = os.path.join(td, "pairs.tsv")
        with open(fa, "w") as fh:
            for k, s in seqs.items():
                fh.write(f">{k}\n{s}\n")
        cmd = [vsearch, "--allpairs_global", fa, "--acceptall",
               "--iddef", str(iddef), "--threads", str(threads),
               "--userout", out, "--userfields", "query+target+id", "--quiet"]
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL)
        except FileNotFoundError:
            sys.exit(f"vsearch not found at '{vsearch}' (use --vsearch)")
        res = {}
        with open(out) as fh:
            for line in fh:
                q, t, idp = line.rstrip("\n").split("\t")
                res[(q, t)] = float(idp)
        return res


def mean(x):
    return sum(x) / len(x)


# ------------------------------------------------------------------ main --
def main():
    a = parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    clusters, acc2cl = read_table(a.table)
    seqs = read_fasta(a.fasta)

    # sequence -> cluster, matching on the versionless accession
    seq_cluster = {}
    for k in seqs:
        cl = acc2cl.get(strip_ver(k))
        if cl:
            seq_cluster[k] = cl
    print(f"{len(clusters)} clusters; {len(seqs)} sequences in FASTA, "
          f"{len(seq_cluster)} match Table S2")
    if not seq_cluster:
        sys.exit("no FASTA id matches an accession in Table S2")

    # centroid sequence of each cluster, if present
    by_key = {strip_ver(k): k for k in seq_cluster}
    centroid_seq = {cl: by_key.get(strip_ver(cl)) for cl in clusters}
    for cl, c in clusters.items():
        c["centroid_len"] = len(seqs[centroid_seq[cl]]) if centroid_seq[cl] else None

    per_cluster = defaultdict(int)
    for cl in seq_cluster.values():
        per_cluster[cl] += 1
    centroid_only = all(v == 1 for v in per_cluster.values()) and \
        all(strip_ver(k) in {strip_ver(c) for c in clusters} for k in seq_cluster)
    print("sequence file provides " +
          ("ONE sequence per cluster (centroids only)" if centroid_only else
           f"cluster members; members < {a.min_len} bp are excluded from the "
           "averages, centroids always kept"))

    species_clusters = defaultdict(list)
    for cl, c in clusters.items():
        species_clusters[c["species"]].append(cl)
    multi = sorted(sp for sp, cls in species_clusters.items()
                   if len(cls) >= a.min_clusters)
    print(f"{len(multi)} species in >= {a.min_clusters} clusters")

    pairs_rows, within_rows, largest_rows, species_rows = [], [], [], []
    used_per_cluster = defaultdict(int)
    matrices = {}
    n_dropped_total = 0

    for sp in multi:
        cl_sp = species_clusters[sp]
        members = {k: cl for k, cl in seq_cluster.items() if cl in cl_sp}
        have = sorted({cl for cl in members.values()})
        if len(have) < 2:
            print(f"  {sp}: only {len(have)} of {len(cl_sp)} clusters have a "
                  "sequence in the FASTA; skipped")
            continue
        # length floor
        use = {}
        n_short = 0
        for k, cl in members.items():
            if k == centroid_seq.get(cl) or len(seqs[k]) >= a.min_len:
                use[k] = cl
            else:
                n_short += 1
        n_dropped_total += n_short
        for k, cl in use.items():
            used_per_cluster[cl] += 1
        ident = allpairs(a.vsearch, {k: seqs[k] for k in use}, a.iddef, a.threads)

        between = defaultdict(list)
        within = defaultdict(list)
        low = 0
        for (q, t), idp in ident.items():
            if idp < 50:
                low += 1
            cq, ct = use[q], use[t]
            if cq == ct:
                within[cq].append(idp)
            else:
                between[tuple(sorted((cq, ct)))].append(idp)
        if low:
            print(f"  {sp}: {low} sequence pair(s) below 50 % identity -- "
                  "reverse strand or wrong region?")

        n_used = defaultdict(int)
        for cl in use.values():
            n_used[cl] += 1
        largest = max(cl_sp, key=lambda c: clusters[c]["n_seq"])

        pair_means = []
        for (c1, c2), v in sorted(between.items()):
            n1, n2 = clusters[c1]["n_seq"], clusters[c2]["n_seq"]
            m = mean(v)
            pair_means.append(m)
            pairs_rows.append(dict(
                species=sp, cluster_1=c1, cluster_2=c2, n_seq_1=n1, n_seq_2=n2,
                smaller=min(n1, n2), n_used_1=n_used[c1], n_used_2=n_used[c2],
                seq_pairs=len(v), identity_mean=round(m, 3),
                identity_min=min(v), identity_max=max(v)))
            if largest in (c1, c2):
                minor = c2 if c1 == largest else c1
                largest_rows.append(dict(
                    species=sp, largest=largest,
                    n_seq_largest=clusters[largest]["n_seq"], cluster=minor,
                    n_seq=clusters[minor]["n_seq"], n_used=n_used[minor],
                    identity_mean=round(m, 3)))
        for cl, v in sorted(within.items()):
            within_rows.append(dict(species=sp, cluster=cl, seq_pairs=len(v),
                                    identity_mean=round(mean(v), 3)))

        species_rows.append(dict(
            species=sp, status=clusters[cl_sp[0]]["status"],
            n_clusters=len(cl_sp), n_clusters_compared=len(have),
            n_sequences_table=sum(clusters[c]["n_seq"] for c in cl_sp),
            n_sequences_used=len(use), n_sequences_below_min_len=n_short,
            largest_cluster=largest, cluster_pairs=len(pair_means),
            identity_mean=round(mean(pair_means), 3),
            identity_sd=round(statistics.stdev(pair_means), 3)
            if len(pair_means) > 1 else "",
            identity_min=min(pair_means), identity_max=max(pair_means),
            pairs_above_threshold=sum(m >= a.threshold for m in pair_means),
            pairs_below_95=sum(m < 95 for m in pair_means),
            pairs_below_90=sum(m < 90 for m in pair_means)))
        print(f"  {sp:38s} {len(have):2d} clusters  {len(use):4d} seqs used "
              f"({n_short:3d} < {a.min_len} bp dropped)  {len(pair_means):4d} "
              f"cluster pairs  mean {mean(pair_means):.1f}%  "
              f"[{min(pair_means):.1f}-{max(pair_means):.1f}]")

        # matrix for the detail species (diagonal = a cluster against itself)
        matrices[sp] = (have, {k: mean(v) for k, v in between.items()})

    if not pairs_rows:
        sys.exit("no species with >= 2 clusters represented in the FASTA")

    # ------------------------------------------------------------- writing --
    def write(name, rows, fields):
        path = os.path.join(a.outdir, f"{a.stem}_{name}.tsv")
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                               lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        return path

    written = [
        write("cluster_pairs", pairs_rows, list(pairs_rows[0])),
        write("within_cluster", within_rows,
              ["species", "cluster", "seq_pairs", "identity_mean"]),
        write("vs_largest", largest_rows, list(largest_rows[0])),
        write("species", species_rows, list(species_rows[0]))]

    cluster_rows = []
    for cl, c in clusters.items():
        if c["species"] not in multi:
            continue
        cluster_rows.append(dict(
            species=c["species"], cluster=cl, n_seq=c["n_seq"],
            centroid=cl, centroid_len=c["centroid_len"] or "",
            centroid_short=(int(c["centroid_len"] < a.min_len)
                            if c["centroid_len"] else ""),
            n_used=used_per_cluster.get(cl, 0)))
    written.append(write("clusters", cluster_rows, list(cluster_rows[0])))

    detail = a.detail if a.detail else \
        [r["species"] for r in sorted(species_rows, key=lambda r: -r["n_clusters"])][:a.detail_n]
    for sp in detail:
        if sp not in matrices:
            print(f"  no matrix for '{sp}' (not a multi-cluster species in the data)")
            continue
        have, bt = matrices[sp]
        slug = re.sub(r"[^A-Za-z0-9]+", "_", sp).strip("_")
        path = os.path.join(a.outdir, f"{a.stem}_matrix_{slug}.tsv")
        with open(path, "w") as fh:
            fh.write("cluster\t" + "\t".join(have) + "\n")
            for c1 in have:
                vals = []
                for c2 in have:
                    if c1 == c2:
                        vals.append("100")
                    else:
                        v = bt.get(tuple(sorted((c1, c2))))
                        vals.append("" if v is None else f"{v:.3f}")
                fh.write(c1 + "\t" + "\t".join(vals) + "\n")
        written.append(path)

    # ------------------------------------------------- numbers behind it --
    thr = a.threshold
    print("\n================ NUMBERS BEHIND THE FIGURE ================")
    print(f"  {'sequence file':46s} {os.path.basename(a.fasta)}")
    print(f"  {'comparison unit':46s} " +
          ("centroid vs centroid" if centroid_only else
           f"members >= {a.min_len} bp (centroids always included)"))
    print(f"  {'identity':46s} vsearch --allpairs_global --iddef {a.iddef}")
    print(f"  {'species plotted':46s} {len(species_rows)}")
    print(f"  {'clusters those species occupy':46s} "
          f"{sum(r['n_clusters'] for r in species_rows)}")
    print(f"  {'cluster pairs':46s} {len(pairs_rows)}")
    print(f"  {'sequences used / dropped as short':46s} "
          f"{sum(r['n_sequences_used'] for r in species_rows)} / {n_dropped_total}")
    print(f"  {'short centroids (< min-len)':46s} "
          f"{sum(1 for r in cluster_rows if r['centroid_short'] == 1)}")
    sm = [r["identity_mean"] for r in species_rows]
    print(f"  {'median species value':46s} {statistics.median(sm):.1f}%")
    print(f"  {'range of species values':46s} {min(sm):.1f}% - {max(sm):.1f}%")
    print(f"  {'cluster pairs at or above ' + str(thr) + '%':46s} "
          f"{sum(r['identity_mean'] >= thr for r in pairs_rows)} of {len(pairs_rows)}")
    print(f"  {'cluster pairs below ' + str(thr - 3) + '%':46s} "
          f"{sum(r['identity_mean'] < thr - 3 for r in pairs_rows)} of {len(pairs_rows)}")
    sing = [r["identity_mean"] for r in largest_rows if r["n_seq"] == 1]
    more = [r["identity_mean"] for r in largest_rows if r["n_seq"] > 1]
    if sing and more:
        print(f"  {'median identity to largest cluster':46s} "
              f"{statistics.median(sing):.1f}% (singletons) vs "
              f"{statistics.median(more):.1f}% (>= 2 seqs)")
    if within_rows:
        print(f"  {'median within-cluster identity':46s} "
              f"{statistics.median(r['identity_mean'] for r in within_rows):.1f}% "
              f"({len(within_rows)} clusters with >= 2 sequences used)")
    # the fractions quoted in the text for the detail species, with the
    # class counts used in the heatmap legend
    print("\n  detail species: cluster pairs by identity class")
    edges = [90, 95, 98, 99]
    for sp in detail:
        if sp not in matrices:
            continue
        v = [r["identity_mean"] for r in pairs_rows if r["species"] == sp]
        n = len(v)
        cls = [sum(1 for x in v if x < 90),
               sum(1 for x in v if 90 <= x < 95),
               sum(1 for x in v if 95 <= x < 98),
               sum(1 for x in v if 98 <= x < 99),
               sum(1 for x in v if x >= 99)]
        below95 = cls[0] + cls[1]
        print(f"    {sp}: {n} cluster pairs; <90 = {cls[0]}, 90-95 = {cls[1]}, "
              f"95-98 = {cls[2]}, 98-99 = {cls[3]}, >=99 = {cls[4]}")
        print(f"      below 95%: {below95} of {n} = {100 * below95 / n:.0f}%;  "
              f"at or above {thr:g}%: {cls[3] + cls[4]} of {n} = "
              f"{100 * (cls[3] + cls[4]) / n:.0f}%")
    print("\n  per species, most similar first")
    for r in sorted(species_rows, key=lambda r: -r["identity_mean"]):
        print(f"    {r['species']:38s} {r['n_clusters']:2d} clusters "
              f"{r['n_sequences_table']:4d} seqs  mean {r['identity_mean']:5.1f}%  "
              f"min {r['identity_min']:5.1f}%  max {r['identity_max']:5.1f}%  "
              f"largest {r['largest_cluster']}")
    print("\nwritten")
    for p in written:
        print("  " + p)


if __name__ == "__main__":
    main()
