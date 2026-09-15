# Microsporidia species attribute database — analysis code

Code for data curation, analysis, and figure generation for:

> \*\*An expanded microsporidia species attribute database documents the extensive ecological and phenotypic diversity of these parasites.\*\*
> Chen C, Alassal A, Reinke AW. \*(citation to be added on publication)\*

The database itself (1,795 microsporidia species and their attributes) is
distributed as Table S1 of the paper, and the most recent version is
maintained at [https://www.reinkelab.org/microsporidiaspecies](https://www.reinkelab.org/microsporidiaspecies), where new
species descriptions can also be submitted.

This repository contains everything needed to go from the curated database
spreadsheet and GenBank 18S rRNA sequences to every figure and supplemental
table in the paper. Scripts are self-contained (no package structure): each
one states its inputs, writes its outputs into its own directory, and prints
a QC report of every parsing decision, exclusion, and repair it makes.

## Overview of the pipeline

```
GenBank 18S sequences                       Curated database (.xlsx)
        |                                            |
  trim\_18S\_cm.py  (Infernal/RF02542)                 |
        |                                            |
  annotate\_fasta\_ncbi.py  (host/locality/date        |
        |                  metadata via E-utilities) |
        |                                            |
  sweep\_cluster\_thresholds.py  -> Fig S2           |
        |                                            |
  VSEARCH clustering (98% id, >=600 bp)              |
        |                                            |
  cluster\_to\_database.py  <------------------------+   matches clusters to
        |            |                               |   database entries;
        |            +--> apply\_new\_host\_records.py+   new entries -> Fig S1
        |                                            |
        |                          host\_taxonomy\_environment.py  (GBIF/COL/
        |                          WoRMS/NCBI lookups)      |
        |                          predict\_host\_environment.py (genus table,
        |                          offline)                 |     -> Table S4
        |                                 |                 |
        |                        make\_table\_s3.py -> Table S3, Fig S3
        |                                 |
        +--> MUSCLE + trimAl + IQ-TREE    +--> R figure scripts -> Fig 1A-F
        |    (external; see Methods)      |
        +--> prune\_tree.R                 +--> assemble\_figure.R -> Fig 1
        +--> plot\_18S\_tree\_fragmented\_clades.R -> Fig 2, Fig S4
        +--> split\_species\_identity.py -> fig\_split\_species\_identity.R -> Fig S5
        |
  make\_cluster\_table.py -> Table S2
  microsporidia\_database\_flowchart.R -> Fig S1 and all database counts
```

## Which script makes which figure/table

|Output|Script(s)|
|-|-|
|Fig 1A (discovery by decade)|`fig1\_discovery\_by\_decade.R`|
|Fig 1B (geography)|`fig\_locality\_map.R`|
|Fig 1C (host phyla)|`fig\_host\_phylum.R`|
|Fig 1D (host environments, Venn)|`fig\_environment\_venn.R`|
|Fig 1E (tissues infected)|`fig\_tissue\_sites.R`|
|Fig 1F (spore size/shape/coils)|`fig\_spore\_size\_shape.R`|
|Fig 1 composite|`assemble\_figure.R` (runs the six panel scripts)|
|Fig 2 / Fig S4 (18S phylogeny)|`prune\_tree.R` → `plot\_18S\_tree\_fragmented\_clades.R`|
|Fig S1 (database assembly flowchart)|`microsporidia\_database\_flowchart.R`|
|Fig S2 (clustering threshold sweep)|`sweep\_cluster\_thresholds.py` → `fig\_cluster\_sweep.R`|
|Fig S3 (habitat concordance)|`fig\_habitat\_concordance.R`|
|Fig S5 (within-species 18S identity)|`split\_species\_identity.py` → `fig\_split\_species\_identity.R`|
|Table S2 (cluster table)|`make\_cluster\_table.py`|
|Table S3 (host taxonomy + environment)|`make\_table\_s3.py`|
|Table S4 (genus-level habitat table)|`predict\_host\_environment.py` (`00\_genus\_habitat\_table.tsv`)|

## The named/provisional classifier

Whether a database entry counts as a **named** species (valid binomial) or a
**provisional** one (sp., isolate code, sequence-defined cluster, …) is
decided by a single rule implemented once per language:

* `classify\_species.R` — sourced by every R figure script, the composite
assembler, and the flowchart
* `classify\_species.py` — imported by every Python script

## Requirements

**Python 3.9+** with `pandas`, `openpyxl`, `requests`, `biopython`
(`pip install pandas openpyxl requests biopython`).

**R 4.3+** with `readxl`, `ggplot2`, `patchwork`, `scales`, `sf`,
`rnaturalearthdata`, `ape`, `ggtree`, `treeio`, `tidytree`, `dplyr`;
optional but recommended: `ggtext` (colored bar-label breakdowns),
`stringi` (accent handling).

**External tools** (versions used in the paper): Infernal 1.1.4
(`cmsearch`), VSEARCH 2.21.1, MUSCLE 5.1, trimAl 1.5, IQ-TREE 3.0.1.

**Inputs not included in this repository:** the curated database
spreadsheet (Table S1 / reinkelab.org), the original 2021 database
(mBio Table S1, `mbio.0149021st001.xlsx`, for the flowchart's provenance
accounting), and the Rfam RF02542 covariance model for 18S trimming.

An NCBI Entrez email (and optionally an API key) is required for
`annotate\_fasta\_ncbi.py` and `host\_taxonomy\_environment.py`.

## Running the pipeline

Every script prints a usage block in its header; the outline below shows the
order and the main handoffs. Exact parameters for the published analysis are
given in the paper's Methods.

**1. Sequence processing and clustering**

```
python trim\_18S\_cm.py all\_18S.fasta trimmed.fasta --cm\_file RF02542.cm --cpu 8
python annotate\_fasta\_ncbi.py trimmed.fasta -o trimmed.annot.fasta --email you@example.org
python sweep\_cluster\_thresholds.py --fasta trimmed.fasta \\
    --meta trimmed.annot.fasta.meta.tsv --outdir sweep        # -> Fig S2 input
vsearch --cluster\_fast trimmed.fasta --id 0.98 --iddef 2 \\
    --minseqlength 600 --uc clusters98.uc --centroids centroids.fasta
```

**2. Cluster ↔ database integration**

```
python cluster\_to\_database.py --uc clusters98.uc \\
    --meta trimmed.annot.fasta.meta.tsv --db database.xlsx \\
    --centroids centroids.fasta --outdir cluster\_output
python apply\_new\_host\_records.py --db database.xlsx \\
    --report cluster\_output/01\_ALREADY\_IN\_DATABASE.tsv --dry-run
```

`cluster\_to\_database.py` writes new database rows
(`02\_NEW\_ENTRIES\_TO\_ADD.xlsx`), the tree-annotation table
(`03\_CLUSTER\_HOST\_ATTRIBUTES.tsv`), and the host-bearing centroid FASTA
(`07\_CENTROIDS\_WITH\_HOSTS.fasta`) that the phylogeny is built from.

**3. Host taxonomy and environment**

```
python host\_taxonomy\_environment.py --db database.xlsx --email you@example.org
python predict\_host\_environment.py  --db database.xlsx \\
    --compare host\_taxonomy\_output/01\_species\_host\_long.tsv
python make\_table\_s3.py
```

Lookups are cached (`host\_lookup\_cache.json`), so re-runs only fetch hosts
not seen before; `--offline` uses the cache alone.

**4. Phylogeny**

Align `07\_CENTROIDS\_WITH\_HOSTS.fasta` with MUSCLE, trim with trimAl
(`-gt 0.1`), infer with IQ-TREE (`-m MFP -B 1000 --alrt 1000 --bnni`), then:

```
Rscript prune\_tree.R  full.treefile  cluster\_output/03\_CLUSTER\_HOST\_ATTRIBUTES.tsv
Rscript plot\_18S\_tree\_fragmented\_clades.R      # Fig 2 + labelled Fig S4
python split\_species\_identity.py                  # Fig S5 inputs
Rscript fig\_split\_species\_identity.R
```

**5. Database figures**

Each Fig 1 panel script runs standalone against the spreadsheet or the
host-taxonomy tables; the composite runs them all and reports every count
behind the figure:

```
Rscript assemble\_figure.R                       # Fig 1 (A-F)
Rscript microsporidia\_database\_flowchart.R      # Fig S1 + counts audit
Rscript fig\_habitat\_concordance.R                 # Fig S3
Rscript fig\_cluster\_sweep.R sweep/threshold\_sweep.tsv   # Fig S2
python make\_cluster\_table.py --assignments clade\_assignments.csv   # Table S2
```

## Design notes

* **Nothing fails silently.** Free-text parsing (localities, tissues, spore
measurements, host names) reports every cell it could not resolve to
`QC\_\*.tsv` files next to each figure's output, and unusual repairs (e.g.
coil ranges that Excel converted to dates, misspelled host genera) are
listed so they can be fixed at source.
* **Counts are printed, not recomputed.** Each figure script prints the
numbers behind the figure (`figure\_numbers.tsv`,
`flowchart\_counts.csv`, …); the manuscript quotes those outputs directly.
* **Provenance lives in the spreadsheet.** Additions and merges are marked
in the `Important Remarks` column (`\[ADDED …]`, `\[TRANSFERRED/ADDED …]`,
`\[MERGED …]`, `PROVISIONAL SPECIES: defined by 18S sequence cluster`,
`date basis: …`), and the flowchart script reconstructs the database's
assembly history from those markers.
* Hyperparasite chains ("X (direct host); Y (hyperhost)") are recorded in
full but only the microsporidian's direct host is counted in analyses.

## AI assistance

As declared in the manuscript, large language models (Anthropic Claude and
OpenAI ChatGPT) were used to help draft and revise these scripts and to
extract species attributes from the literature. All code and extracted data
were reviewed by the authors, and the genus-level habitat table (Table S4)
and every automated parsing decision carry confidence/basis annotations for
auditing.

## Contact

Aaron W. Reinke — aaron.reinke@utoronto.ca — Reinke Lab, Department of
Molecular Genetics, University of Toronto. [https://www.reinkelab.org](https://www.reinkelab.org)

