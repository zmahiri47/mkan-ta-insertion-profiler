# mkan-ta-insertion-profiler

Exploration and description of transposon (Tn) insertion distribution in
*Mycobacterium kansasii*.

This work supports:

> Budell, W. C., Germain, G. A., Janisch, N., McKie-Krisberg, Z.,
> Jayaprakash, A. D., Resnick, A. E., & Quadri, L. E. N. (2020).
> *Transposon mutagenesis in Mycobacterium kansasii links a small RNA gene
> to colony morphology and biofilm formation and identifies 9,885
> intragenic insertions that do not compromise colony outgrowth.*
> MicrobiologyOpen, 9(4), e988.
> https://onlinelibrary.wiley.com/doi/full/10.1002/mbo3.988

## What's in this repo

| File | Description |
|---|---|
| `ta_profiler.py` | CLI tool that classifies every TA transposon insertion site as **intragenic** (inside a gene) or **intergenic** (between two genes), against a GFF gene annotation. |
| `README.md` | This file. |

`ta_profiler.py` is a cleaned-up, runnable rewrite of the project's original
`whole-script.py`, which was a raw export of an interactive Jupyter session
used to develop the analysis. That original script hard-coded local file
names, split GFF attributes on a fixed number of `;` characters, referenced
variables that only existed in the author's live notebook kernel (e.g.
hand-edited intermediate CSVs), and could not be run end-to-end as a
standalone script. This rewrite reproduces the same underlying analysis —
annotate every TA site as intragenic/intergenic and report the relevant
gene context — as a single, documented, parallelized command-line tool.

## How it works

For every TA transposon insertion coordinate, the script asks:

1. **Does this position fall inside a gene's start/end span?**
   If yes → the site is **intragenic**. The tool reports the gene, its
   locus tag, strand, the insertion's position relative to the gene start,
   and what percentage of the way through the gene the insertion falls.

2. **If not, what are the closest flanking genes?**
   The tool finds the nearest upstream gene (the one whose `end` is
   closest to, but not past, the insertion) and the nearest downstream
   gene (the one whose `start` is closest to, but not before, the
   insertion), and reports both, along with the distance from each.

Every insertion site is independent of every other, so this per-site
classification is [embarrassingly parallel](https://en.wikipedia.org/wiki/Embarrassingly_parallel) —
`ta_profiler.py` takes advantage of that by splitting the insertion list
into chunks and distributing them across a pool of worker processes with
Python's `concurrent.futures.ProcessPoolExecutor`.

## Requirements

- Python 3.9+
- [pandas](https://pandas.pydata.org/)

```bash
pip install pandas
```

## Input file formats

### `--gff`: gene annotation (TSV, no header)

Standard 9-column GFF layout:

```
seqname  source  feature  start  end  score  strand  frame  attribute
```

The `attribute` column should contain semicolon-separated `key=value`
pairs, e.g.:

```
ID=gene-MKAN_RS27780;Dbxref=GeneID:29696832;Name=MKAN_RS27780;gbkey=Gene;gene_biotype=protein_coding;locus_tag=MKAN_RS27780;old_locus_tag=MKAN_29095
```

By default, rows with `feature` equal to `region` or `CDS` are excluded
before classification (see `--exclude-features` below), since `region`
describes the whole replicon and `CDS` duplicates its parent `gene` row
for this analysis.

### `--hits`: TA insertion sites (CSV)

A single column of integer genomic coordinates, one per line, with or
without a header row:

```
Plas_Tnhits
720
946
1075
1226
```

## Usage

```bash
python3 ta_profiler.py --gff annotation.tsv --hits ta_hits.csv
```

### Full option list

```
python3 ta_profiler.py --help
```

| Option | Description |
|---|---|
| `--gff PATH` | **(required)** GFF-style annotation TSV file. |
| `--hits PATH` | **(required)** CSV file of TA insertion site coordinates. |
| `--hits-column NAME` | Column name to read positions from, if `--hits` has a header. Defaults to the first column. |
| `--outdir DIR` | Output directory (created if missing). Default: current directory. |
| `--prefix NAME` | Filename prefix for outputs. Default: `ta_profile`. |
| `--exclude-features [F ...]` | GFF `feature` values to drop before classification. Default: `region CDS`. Pass with no values to disable filtering. |
| `--workers N` | Number of parallel worker processes. Default: `1` (serial). Use `--workers 0` to auto-detect and use all CPU cores. |
| `--chunk-size N` | TA sites per worker task when running in parallel. Default: `200`. |
| `-v`, `--verbose` | Enable debug-level logging. |
| `--version` | Print the tool version and exit. |

### Examples

Basic run, serial:

```bash
python3 ta_profiler.py --gff annotation.tsv --hits ta_hits.csv
```

Parallel run using 8 worker processes, custom output directory:

```bash
python3 ta_profiler.py --gff annotation.tsv --hits ta_hits.csv \
    --outdir results/ --workers 8
```

Auto-detect all available CPU cores:

```bash
python3 ta_profiler.py --gff annotation.tsv --hits ta_hits.csv --workers 0
```

TA hit file has a header column named `Tnhits`:

```bash
python3 ta_profiler.py --gff annotation.tsv --hits ta_hits.csv \
    --hits-column Tnhits
```

Keep `CDS` rows in the annotation (only exclude `region`):

```bash
python3 ta_profiler.py --gff annotation.tsv --hits ta_hits.csv \
    --exclude-features region
```

Verbose logging, to watch classification progress on a large genome:

```bash
python3 ta_profiler.py --gff annotation.tsv --hits ta_hits.csv -v
```

## Output

Running the tool produces two files in `--outdir`:

- **`<prefix>.csv`** — one row per TA site / gene relationship, with
  columns:

  | Column | Meaning |
  |---|---|
  | `pos` | TA insertion coordinate |
  | `category` | `intragenic` or `intergenic` |
  | `id1`, `locus_tag1`, `old_locus_tag1`, `start1`, `end1`, `strand1` | The containing gene (intragenic) or nearest upstream gene (intergenic) |
  | `rel_pos_on_gene1`, `gene1_len`, `pct_of_gene1` | Insertion position relative to gene1's start, gene1's length, and the insertion's fractional position through gene1 (intragenic only) |
  | `dist_from_gene1_end` | Distance from gene1's end to the insertion (intergenic only) |
  | `id2`, `locus_tag2`, `old_locus_tag2`, `start2`, `end2`, `strand2` | Nearest downstream gene (intergenic only) |
  | `dist_to_gene2_start` | Distance from the insertion to gene2's start (intergenic only) |

- **`<prefix>_summary.txt`** — a short plain-text summary of how many
  distinct TA sites were classified and how many fell into each category.

## Notes / limitations

- Insertion sites that overlap more than one annotated gene (e.g.
  overlapping genes on opposite strands) produce one output row per
  overlapping gene, all sharing the same `pos`.
- Insertion sites beyond the first or last annotated gene on a replicon
  will have only an upstream or only a downstream gene reported (the
  other side's columns will be empty).
- This tool performs coordinate classification only; it does not itself
  call transposon insertion sites from sequencing reads — `--hits` is
  expected to already be a list of called TA insertion coordinates.
