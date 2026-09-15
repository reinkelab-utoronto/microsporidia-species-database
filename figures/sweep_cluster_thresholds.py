#!/usr/bin/env python3
"""
Sweep vsearch identity and length/coverage thresholds, scoring each combination
against species names as ground truth.

Two competing failure modes are counted, exactly as you framed them:

  split_species   binomials that appear in more than one cluster
                  (over-splitting; maximal at 100% identity)
  mixed_clusters  clusters containing more than one binomial
                  (over-lumping; grows as identity falls)

Only full binomials are used for scoring -- "Genus sp.", "uncultured ...",
and unnamed sequences are clustered but excluded from the metrics.

Because raw counts shrink as the length threshold discards sequences, the
fractions (frac_species_split, frac_clusters_mixed) are the comparable
quantities across the grid, and a V-measure is reported as a single objective:

  completeness  formal version of "this species is all in one cluster"
  homogeneity   formal version of "this cluster is all one species"
  v_measure     their harmonic mean -- maximised at the trade-off point

Usage
-----
  python sweep_cluster_thresholds.py \
      --fasta trimmed.fasta --meta trimmed.annot.fasta.meta.tsv \
      --identities 0.97 0.975 0.98 0.985 0.99 0.995 1.0 \
      --min-covs 0 300 500 700 900 \
      --outdir sweep --threads 8
"""

import argparse
import csv
import math
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict

ACC_RE = re.compile(r"^[A-Z]{1,6}_?\d{5,9}(\.\d+)?$")
UNINFORMATIVE = re.compile(
    r"^(uncultured|unidentified|unclassified|environmental|microsporidian?"
    r"|fungal|eukaryot\w*)\b", re.I)


# ---------------------------------------------------------------------------
def norm_acc(text):
    t = str(text).strip().upper().split(".")[0]
    return t if ACC_RE.match(t + ".1") else t


def binomial(name):
    """Return a normalised binomial, or '' if the name is not one."""
    if not name:
        return ""
    s = str(name).lower()
    s = re.sub(r"\(.*?\)", " ", s)
    s = re.sub(r"\b(cf|aff|nr)\.?\s+", " ", s)
    s = re.sub(r"[^a-z0-9\s]", " ", s)
    s = re.sub(r"\s+", " ", s).strip()
    if not s or UNINFORMATIVE.match(s):
        return ""
    parts = s.split()
    if len(parts) < 2 or parts[1] in ("sp", "spp"):
        return ""
    return " ".join(parts[:2])          # genus + specific epithet only


def read_fasta(path):
    """Yield (id, description, sequence)."""
    ident = desc = None
    seq = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if ident is not None:
                    yield ident, desc, "".join(seq)
                head = line[1:]
                ident = head.split()[0] if head.split() else head
                desc = head[len(ident):].strip()
                seq = []
            else:
                seq.append(line.strip())
    if ident is not None:
        yield ident, desc, "".join(seq)


def get_length_value(desc, seq, field):
    """Length metric for filtering: model coverage if available, else bp."""
    if field == "cm_cov":
        m = re.search(r"cm_cov=(\d+)", desc)
        if m:
            return int(m.group(1))
        return len(seq)          # fall back for sequences lacking the tag
    return len(seq)


def read_species(meta_path):
    """Return {accession_no_version: binomial}."""
    out = {}
    with open(meta_path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            acc = norm_acc(row.get("accession", ""))
            b = binomial(row.get("organism", ""))
            if acc and b:
                out[acc] = b
    return out


def read_uc(path):
    """Return ({member_label: cluster_id}, {member_label: pct_id_to_centroid})."""
    assign, pid = {}, {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10 or f[0] == "C":
                continue
            label = f[8].split(";")[0].strip()
            assign[label] = f[1]
            try:
                pid[label] = float(f[3])
            except ValueError:
                pid[label] = 100.0          # seeds report '*'
    return assign, pid


def load_synonyms(path):
    """
    Read a synonym map: two tab- or comma-separated names per line, meaning
    'these are the same species'. Returns {normalised_name: canonical_name}.
    Lines starting with # are ignored.
    """
    if not path:
        return {}
    groups = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = [binomial(p) for p in re.split(r"[\t,;]", line.strip())]
            parts = [p for p in parts if p]
            if len(parts) < 2:
                continue
            canon = sorted(parts)[0]
            root = groups.get(canon, canon)
            for p in parts:
                groups[p] = root
    # resolve chains
    for k in list(groups):
        seen = set()
        v = groups[k]
        while v in groups and groups[v] != v and v not in seen:
            seen.add(v)
            v = groups[v]
        groups[k] = v
    return groups


def classify_mixture(names, synonyms):
    """
    Guess why a cluster holds more than one binomial. The distinctive signature
    of a genus reassignment is the SAME epithet under DIFFERENT genera
    (Nosema necatrix / Vairimorpha necatrix).
    """
    canon = {synonyms.get(n, n) for n in names}
    if len(canon) == 1:
        return "known synonyms"
    genera = {n.split()[0] for n in names}
    epithets = {n.split()[1] for n in names if len(n.split()) > 1}
    if len(genera) > 1 and len(epithets) == 1:
        return "likely genus reassignment (same epithet)"
    if len(genera) > 1 and len(epithets) < len(names):
        return "partial epithet overlap (possible reassignment)"
    if len(genera) == 1:
        return "congeners (same genus, different epithet)"
    return "different genus and epithet"


# ---------------------------------------------------------------------------
# metrics
# ---------------------------------------------------------------------------
def entropy(counts, n):
    return -sum((c / n) * math.log(c / n) for c in counts if c > 0)


def cluster_metrics(pairs):
    """
    pairs: list of (cluster_id, species). Returns the count-based metrics you
    specified plus homogeneity / completeness / V-measure and adjusted Rand.
    """
    n = len(pairs)
    if n == 0:
        return None

    by_species = defaultdict(set)
    by_cluster = defaultdict(set)
    contingency = Counter()
    for cid, sp in pairs:
        by_species[sp].add(cid)
        by_cluster[cid].add(sp)
        contingency[(cid, sp)] += 1

    split_species = sum(1 for s, cs in by_species.items() if len(cs) > 1)
    mixed_clusters = sum(1 for c, ss in by_cluster.items() if len(ss) > 1)

    cluster_sizes = Counter(cid for cid, _ in pairs)
    species_sizes = Counter(sp for _, sp in pairs)

    h_c = entropy(species_sizes.values(), n)     # species entropy
    h_k = entropy(cluster_sizes.values(), n)     # cluster entropy
    h_c_given_k = -sum(
        (cnt / n) * math.log(cnt / cluster_sizes[cid])
        for (cid, sp), cnt in contingency.items() if cnt > 0)
    h_k_given_c = -sum(
        (cnt / n) * math.log(cnt / species_sizes[sp])
        for (cid, sp), cnt in contingency.items() if cnt > 0)

    homogeneity = 1.0 if h_c == 0 else 1 - h_c_given_k / h_c
    completeness = 1.0 if h_k == 0 else 1 - h_k_given_c / h_k
    v = (0.0 if homogeneity + completeness == 0 else
         2 * homogeneity * completeness / (homogeneity + completeness))

    # adjusted Rand index
    comb2 = lambda x: x * (x - 1) / 2
    sum_ij = sum(comb2(c) for c in contingency.values())
    sum_i = sum(comb2(c) for c in cluster_sizes.values())
    sum_j = sum(comb2(c) for c in species_sizes.values())
    total = comb2(n)
    expected = sum_i * sum_j / total if total else 0
    max_index = (sum_i + sum_j) / 2
    ari = 0.0 if max_index == expected else (sum_ij - expected) / (max_index - expected)

    return {
        "n_scored_seqs": n,
        "n_species": len(by_species),
        "n_clusters_with_binomials": len(by_cluster),
        "split_species": split_species,
        "frac_species_split": round(split_species / len(by_species), 4),
        "mixed_clusters": mixed_clusters,
        "frac_clusters_mixed": round(mixed_clusters / len(by_cluster), 4),
        "homogeneity": round(homogeneity, 4),
        "completeness": round(completeness, 4),
        "v_measure": round(v, 4),
        "adjusted_rand": round(ari, 4),
    }


# ---------------------------------------------------------------------------
def run_vsearch(fasta, ident, uc_path, threads, iddef, cmd_name, extra):
    cmd = ["vsearch", f"--{cmd_name}", fasta, "--id", str(ident),
           "--iddef", str(iddef), "--uc", uc_path,
           "--threads", str(threads), "--quiet",
           "--minseqlength", "1"]
    cmd += extra
    proc = subprocess.run(cmd, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        sys.exit(f"vsearch failed:\n{proc.stderr}")


def main():
    p = argparse.ArgumentParser(
        description="Sweep identity and length thresholds against species labels.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--fasta", required=True, help="trimmed FASTA")
    p.add_argument("--meta", required=True, help="metadata TSV with an organism column")
    p.add_argument("--identities", nargs="+", type=float,
                   default=[0.96, 0.97, 0.975, 0.98, 0.985, 0.99, 0.995, 1.0])
    p.add_argument("--min-covs", nargs="+", type=int,
                   default=[0, 300, 500, 700, 900],
                   help="minimum model coverage (or bp) thresholds to test")
    p.add_argument("--length-field", choices=["cm_cov", "length"], default="cm_cov",
                   help="filter on the [cm_cov=...] header tag or on raw bp")
    p.add_argument("--cluster-cmd", default="cluster_fast",
                   choices=["cluster_fast", "cluster_size", "cluster_smallmem"])
    p.add_argument("--iddef", type=int, default=2,
                   help="vsearch identity definition (2 ignores terminal gaps)")
    p.add_argument("--vsearch-arg", action="append", default=[],
                   help="extra argument passed to vsearch (repeatable)")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--outdir", default="sweep")
    p.add_argument("--keep-uc", action="store_true",
                   help="keep every .uc file instead of only the best one")
    p.add_argument("--synonyms",
                   help="TSV/CSV of synonymous name pairs (two names per line) "
                        "to collapse before scoring")
    p.add_argument("--no-plot", action="store_true")
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    species = read_species(args.meta)
    records = list(read_fasta(args.fasta))
    sys.stderr.write(
        f"{len(records)} sequences; {len(species)} have a usable binomial "
        f"({len(set(species.values()))} distinct species)\n"
    )
    if not species:
        sys.exit("No binomials found in the metadata -- nothing to score against.")

    synonyms = load_synonyms(args.synonyms)
    if synonyms:
        sys.stderr.write(f"{len(set(synonyms.values()))} synonym groups loaded "
                         f"covering {len(synonyms)} names\n")
    canon_species = {a: synonyms.get(b, b) for a, b in species.items()}

    results, detail_rows, cause_rows = [], [], []
    for cov in sorted(args.min_covs):
        subset = [(i, d, s) for i, d, s in records
                  if get_length_value(d, s, args.length_field) >= cov]
        if len(subset) < 2:
            sys.stderr.write(f"cov>={cov}: only {len(subset)} sequences, skipping\n")
            continue

        sub_fa = os.path.join(args.outdir, f"subset_cov{cov}.fasta")
        with open(sub_fa, "w") as fh:
            for i, d, s in subset:
                fh.write(f">{i}\n{s}\n")

        n_lab = sum(1 for i, _, _ in subset if norm_acc(i) in canon_species)
        n_sp = len({canon_species[norm_acc(i)] for i, _, _ in subset
                    if norm_acc(i) in canon_species})
        sys.stderr.write(
            f"cov>={cov}: {len(subset)} seqs, {n_lab} labelled, {n_sp} species\n")

        for ident in sorted(args.identities):
            uc = os.path.join(args.outdir, f"uc_cov{cov}_id{ident}.uc")
            run_vsearch(sub_fa, ident, uc, args.threads, args.iddef,
                        args.cluster_cmd, args.vsearch_arg)
            assign, pid = read_uc(uc)

            pairs = [(cid, canon_species[norm_acc(lab)])
                     for lab, cid in assign.items()
                     if norm_acc(lab) in canon_species]
            m = cluster_metrics(pairs)
            if m is None:
                continue
            m.update({
                "min_cov": cov, "identity": ident,
                "n_input_seqs": len(subset),
                "n_clusters_total": len(set(assign.values())),
                "frac_seqs_retained": round(len(subset) / len(records), 4),
                "frac_species_retained": round(n_sp / len(set(canon_species.values())), 4),
            })
            # per-cluster record of every cluster holding >1 binomial
            members_by_cluster = defaultdict(list)
            for lab, cid in assign.items():
                a = norm_acc(lab)
                if a in species:
                    members_by_cluster[cid].append((lab, species[a]))
            cause_counts = Counter()
            for cid, mem in sorted(members_by_cluster.items()):
                names = sorted({n for _, n in mem})
                if len(names) < 2:
                    continue
                cause = classify_mixture(names, synonyms)
                cause_counts[cause] += 1
                per_name = "; ".join(
                    f"{n} (n={sum(1 for _, x in mem if x == n)})" for n in names)
                pids = [pid[l] for l, _ in mem if l in pid]
                detail_rows.append({
                    "min_cov": cov, "identity": ident, "cluster_id": cid,
                    "n_members": len(mem), "n_species": len(names),
                    "cause": cause,
                    "species_in_cluster": per_name,
                    "genera": "; ".join(sorted({n.split()[0] for n in names})),
                    "min_pct_id_to_centroid": round(min(pids), 1) if pids else "",
                    "accessions": "; ".join(l for l, _ in mem),
                })
            for cause, k in cause_counts.items():
                cause_rows.append({"min_cov": cov, "identity": ident,
                                   "cause": cause, "n_clusters": k})
            m["mixed_reassignment_like"] = (
                cause_counts["likely genus reassignment (same epithet)"]
                + cause_counts["partial epithet overlap (possible reassignment)"]
                + cause_counts["known synonyms"])
            m["mixed_genuine"] = m["mixed_clusters"] - m["mixed_reassignment_like"]
            results.append(m)
            sys.stderr.write(
                f"   id={ident:<6} clusters={m['n_clusters_total']:<5} "
                f"split_species={m['split_species']:<4} "
                f"mixed_clusters={m['mixed_clusters']:<4} "
                f"V={m['v_measure']}\n")
            if not args.keep_uc:
                os.remove(uc)

    if not results:
        sys.exit("No results produced.")

    cols = ["min_cov", "identity", "n_input_seqs", "frac_seqs_retained",
            "n_scored_seqs", "n_species", "frac_species_retained",
            "n_clusters_total", "n_clusters_with_binomials",
            "split_species", "frac_species_split",
            "mixed_clusters", "frac_clusters_mixed",
            "mixed_reassignment_like", "mixed_genuine",
            "homogeneity", "completeness", "v_measure", "adjusted_rand"]
    tsv = os.path.join(args.outdir, "threshold_sweep.tsv")
    with open(tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(results)

    if detail_rows:
        dcols = list(detail_rows[0].keys())
        with open(os.path.join(args.outdir, "mixed_clusters_detail.tsv"),
                  "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=dcols, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(detail_rows)
        with open(os.path.join(args.outdir, "mixed_cluster_causes.tsv"),
                  "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["min_cov", "identity", "cause",
                                               "n_clusters"], delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(cause_rows)

        # candidate synonym pairs: same epithet, different genus, seen anywhere
        cand = Counter()
        for d in detail_rows:
            if "reassignment" in d["cause"]:
                names = tuple(sorted(
                    re.sub(r"\s*\(n=\d+\)", "", s).strip()
                    for s in d["species_in_cluster"].split(";")))
                cand[names] += 1
        if cand:
            with open(os.path.join(args.outdir, "synonym_candidates.tsv"),
                      "w", newline="") as fh:
                w = csv.writer(fh, delimiter="\t", lineterminator="\n")
                w.writerow(["name_a", "name_b_and_more", "n_grid_cells_seen"])
                for names, k in cand.most_common():
                    w.writerow([names[0], "; ".join(names[1:]), k])

    best_v = max(results, key=lambda r: r["v_measure"])
    best_ari = max(results, key=lambda r: r["adjusted_rand"])
    best_sum = min(results, key=lambda r: r["frac_species_split"] + r["frac_clusters_mixed"])

    # Pareto front: no other setting is better on both error axes at once
    pareto = [r for r in results if not any(
        o["frac_species_split"] <= r["frac_species_split"]
        and o["frac_clusters_mixed"] <= r["frac_clusters_mixed"]
        and (o["frac_species_split"] < r["frac_species_split"]
             or o["frac_clusters_mixed"] < r["frac_clusters_mixed"])
        for o in results)]
    with open(os.path.join(args.outdir, "pareto_front.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(sorted(pareto, key=lambda r: r["frac_species_split"]))

    fmt = lambda r: (f"min_cov={r['min_cov']} id={r['identity']} -> "
                     f"{r['split_species']} split species "
                     f"({r['frac_species_split']:.1%}), "
                     f"{r['mixed_clusters']} mixed clusters "
                     f"({r['frac_clusters_mixed']:.1%}), "
                     f"V={r['v_measure']}, retains "
                     f"{r['frac_species_retained']:.0%} of species")
    sys.stderr.write(
        f"\nWrote {tsv}\n"
        f"Best V-measure:      {fmt(best_v)}\n"
        f"Best adjusted Rand:  {fmt(best_ari)}\n"
        f"Min summed error:    {fmt(best_sum)}\n"
        f"Pareto front has {len(pareto)} settings "
        f"(see {os.path.join(args.outdir, 'pareto_front.tsv')})\n"
    )
    if cause_rows:
        tot = Counter()
        for r in cause_rows:
            tot[r["cause"]] += r["n_clusters"]
        sys.stderr.write("\nMixed-cluster causes summed over the whole grid:\n")
        for cause, k in tot.most_common():
            sys.stderr.write(f"  {k:>5}  {cause}\n")
        sys.stderr.write(
            f"Per-cluster listing: "
            f"{os.path.join(args.outdir, 'mixed_clusters_detail.tsv')}\n")

    if not args.no_plot:
        make_plots(results, args.outdir)


def make_plots(results, outdir):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        sys.stderr.write("matplotlib not installed -- plots skipped\n")
        return

    covs = sorted({r["min_cov"] for r in results})
    idents = sorted({r["identity"] for r in results})

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))

    for cov in covs:
        rows = sorted([r for r in results if r["min_cov"] == cov],
                      key=lambda r: r["identity"])
        x = [r["identity"] for r in rows]
        axes[0].plot(x, [r["frac_species_split"] for r in rows], "o-",
                     label=f"cov>={cov}")
        axes[0].plot(x, [r["frac_clusters_mixed"] for r in rows], "s--",
                     color=axes[0].lines[-1].get_color(), alpha=0.6)
        axes[1].plot(x, [r["v_measure"] for r in rows], "o-", label=f"cov>={cov}")

    axes[0].set_xlabel("identity threshold")
    axes[0].set_ylabel("fraction")
    axes[0].set_title("solid: species split across clusters\n"
                      "dashed: clusters with >1 species")
    axes[0].legend(fontsize=7)

    axes[1].set_xlabel("identity threshold")
    axes[1].set_ylabel("V-measure")
    axes[1].set_title("V-measure (higher is better)")
    axes[1].legend(fontsize=7)

    sc = axes[2].scatter([r["frac_species_split"] for r in results],
                         [r["frac_clusters_mixed"] for r in results],
                         c=[r["identity"] for r in results], cmap="viridis", s=45)
    for r in results:
        axes[2].annotate(f"{r['identity']}", (r["frac_species_split"],
                                              r["frac_clusters_mixed"]),
                         fontsize=5, alpha=0.6)
    axes[2].set_xlabel("fraction of species split")
    axes[2].set_ylabel("fraction of clusters mixed")
    axes[2].set_title("trade-off (lower-left is better)")
    plt.colorbar(sc, ax=axes[2], label="identity")

    plt.tight_layout()
    path = os.path.join(outdir, "threshold_sweep.png")
    plt.savefig(path, dpi=150)
    sys.stderr.write(f"Plots: {path}\n")


if __name__ == "__main__":
    main()
