#!/usr/bin/env python3
"""
ta_profiler.py
==============

TA (transposon) insertion site profiler for *Mycobacterium kansasii* Tn-seq data.

WHAT THIS SCRIPT DOES
----------------------

Given:
  1. A GFF-style gene annotation table (TSV) with the standard 9 GFF columns
     (seqname, source, feature, start, end, score, strand, frame, attribute).
  2. A list of TA-dinucleotide transposon insertion site coordinates.

... the script classifies every insertion site as either:
  * INTRAGENIC -- falling inside the start/end span of one or more annotated
    genes, in which case its position relative to the gene start and its
    percentage-through-the-gene are computed.
  * INTERGENIC -- falling between two genes, in which case the flanking
    (upstream/downstream) genes and the distances to each are reported.

and writes a single tidy CSV describing every insertion.

USAGE
-----
    python3 ta_profiler.py --gff annotation.tsv --hits ta_hits.csv --outdir results/

Run ``python3 ta_profiler.py --help`` for the full option list.
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional

import pandas as pd

__version__ = "1.0.0"

LOG = logging.getLogger("ta_profiler")

# Standard 9-column GFF field names, in order.
GFF_COLUMNS = [
    "seqname", "source", "feature", "start", "end",
    "score", "strand", "frame", "attribute",
]

# Default features to drop from the annotation before insertion-site lookup.
# 'region' rows describe the whole replicon (not a gene) and 'CDS' rows are
# redundant with their parent 'gene'/'pseudogene' row for this analysis.
DEFAULT_EXCLUDED_FEATURES = ("region", "CDS")


# --------------------------------------------------------------------------- #
# Data loading / parsing
# --------------------------------------------------------------------------- #

def parse_gff_attributes(attribute_str: str) -> dict:
    """
    Parse a GFF ``attribute`` column value into a dict.

    GFF attributes are semicolon-separated ``key=value`` pairs, e.g.::

        ID=gene-MKAN_RS27780;Dbxref=GeneID:29696832;Name=MKAN_RS27780;
        gbkey=Gene;gene_biotype=protein_coding;locus_tag=MKAN_RS27780;
        old_locus_tag=MKAN_29095
    """
    fields: dict = {}
    if not isinstance(attribute_str, str):
        return fields
    for chunk in attribute_str.split(";"):
        chunk = chunk.strip()
        if not chunk:
            continue
        if "=" in chunk:
            key, _, value = chunk.partition("=")
            fields[key.strip()] = value.strip()
        else:
            # Malformed / flag-only attribute -- keep it, but namespaced.
            fields.setdefault("_unparsed", []).append(chunk)
    return fields


def load_gff(gff_path: Path, excluded_features: Iterable[str]) -> pd.DataFrame:
    """
    Load a GFF-style annotation TSV and expand its attribute column.

    Returns a DataFrame with the original 9 GFF columns plus one extra
    column per distinct attribute key found in the file (e.g. ``ID``,
    ``Dbxref``, ``Name``, ``gbkey``, ``gene_biotype``, ``locus_tag``,
    ``old_locus_tag``, ``product`` ...), restricted to rows whose
    ``feature`` is not in ``excluded_features``.
    """
    LOG.info("Loading GFF annotation from %s", gff_path)
    genes = pd.read_csv(gff_path, sep="\t", names=GFF_COLUMNS, header=None)
    LOG.info("Loaded %d annotation rows (%d unique feature types)",
              len(genes), genes["feature"].nunique())

    keep_mask = ~genes["feature"].isin(set(excluded_features))
    genes = genes.loc[keep_mask].reset_index(drop=True)
    LOG.info("Kept %d rows after excluding features %s",
              len(genes), sorted(set(excluded_features)))

    # Expand the free-form attribute column into real columns.
    parsed = genes["attribute"].apply(parse_gff_attributes)
    attr_df = pd.json_normalize(parsed)
    genes = pd.concat([genes.reset_index(drop=True), attr_df], axis=1)

    # Normalize numeric columns / common prefixed identifiers.
    genes["start"] = pd.to_numeric(genes["start"], errors="coerce").astype("Int64")
    genes["end"] = pd.to_numeric(genes["end"], errors="coerce").astype("Int64")
    for col, prefix in (("ID", "gene-"), ("locus_tag", ""), ("old_locus_tag", "")):
        if col in genes.columns:
            genes[col] = genes[col].astype(str).str.removeprefix(prefix)

    return genes


def load_ta_hits(hits_path: Path, column: Optional[str]) -> list[int]:
    """
    Load a list of TA insertion-site coordinates from a single-column CSV.

    ``column`` may name a header to select from (if the file has a header
    row); otherwise the first column of the file is used and any row that
    doesn't parse as an integer is skipped.
    """
    LOG.info("Loading TA insertion sites from %s", hits_path)
    if column:
        df = pd.read_csv(hits_path)
        hits = pd.to_numeric(df[column], errors="coerce").dropna().astype(int).tolist()
    else:
        df = pd.read_csv(hits_path, header=None, names=["pos"])
        hits = pd.to_numeric(df["pos"], errors="coerce").dropna().astype(int).tolist()
    LOG.info("Loaded %d TA insertion sites", len(hits))
    return hits


# --------------------------------------------------------------------------- #
# Core classification logic
# --------------------------------------------------------------------------- #

GENE_OUTPUT_COLS = ["ID", "locus_tag", "old_locus_tag", "start", "end", "strand"]


def _gene_record(row: pd.Series, suffix: str) -> dict:
    """Extract a small, uniformly-named dict of gene fields for output."""
    out = {}
    for col in GENE_OUTPUT_COLS:
        out[f"{col.lower()}{suffix}"] = row.get(col)
    return out


def classify_hit(pos: int, genes: pd.DataFrame) -> list[dict]:
    """
    Classify a single TA insertion coordinate against the gene annotation.

    Returns one or more result rows (as dicts):
      * one row per overlapping gene if the site is INTRAGENIC, or
      * a single INTERGENIC row describing the nearest upstream/downstream
        flanking genes if the site falls in no gene.

    This is the per-site unit of work that gets distributed across worker
    processes -- it only touches ``pos`` and the read-only ``genes`` table,
    so sites are fully independent of one another and trivially parallel.
    """
    overlapping = genes[(genes["start"] <= pos) & (genes["end"] >= pos)]

    if len(overlapping):
        rows = []
        for _, gene in overlapping.iterrows():
            gene_len = int(gene["end"]) - int(gene["start"])
            rel_pos = pos - int(gene["start"])
            rows.append({
                "pos": pos,
                "category": "intragenic",
                **_gene_record(gene, "1"),
                "rel_pos_on_gene1": rel_pos,
                "gene1_len": gene_len,
                "pct_of_gene1": (rel_pos / gene_len) if gene_len else None,
            })
        return rows

    # No gene contains this position -> find the flanking genes.
    upstream_candidates = genes[genes["end"] <= pos]
    downstream_candidates = genes[genes["start"] >= pos]

    upstream = (
        upstream_candidates.loc[[upstream_candidates["end"].idxmax()]]
        if len(upstream_candidates) else None
    )
    downstream = (
        downstream_candidates.loc[[downstream_candidates["start"].idxmin()]]
        if len(downstream_candidates) else None
    )

    row = {"pos": pos, "category": "intergenic"}
    if upstream is not None:
        g = upstream.iloc[0]
        row.update(_gene_record(g, "1"))
        row["dist_from_gene1_end"] = pos - int(g["end"])
    if downstream is not None:
        g = downstream.iloc[0]
        row.update(_gene_record(g, "2"))
        row["dist_to_gene2_start"] = int(g["start"]) - pos

    return [row]


def _classify_chunk(chunk: list[int], genes: pd.DataFrame) -> list[dict]:
    """Classify a batch of positions in one worker call (reduces IPC overhead)."""
    results: list[dict] = []
    for pos in chunk:
        results.extend(classify_hit(pos, genes))
    return results


def _chunked(seq: list, size: int) -> Iterable[list]:
    for i in range(0, len(seq), size):
        yield seq[i:i + size]


def classify_all_hits(
    hits: list[int],
    genes: pd.DataFrame,
    workers: int = 1,
    chunk_size: int = 200,
) -> pd.DataFrame:
    """
    Classify every TA insertion site, optionally in parallel.

    Each TA site is independent of every other, so we split the hit list
    into chunks and hand them out across a process pool. Every worker gets
    its own copy of the (read-only) ``genes`` DataFrame via the initial
    task argument -- fine for annotation tables of a few thousand genes,
    which is the typical size for a bacterial genome/plasmid.
    """
    if workers <= 1 or len(hits) < chunk_size:
        LOG.info("Classifying %d insertion sites serially", len(hits))
        records = _classify_chunk(hits, genes)
        return pd.DataFrame.from_records(records)

    chunks = list(_chunked(hits, chunk_size))
    LOG.info(
        "Classifying %d insertion sites across %d worker processes "
        "(%d chunks of <=%d sites each)",
        len(hits), workers, len(chunks), chunk_size,
    )

    records: list[dict] = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(_classify_chunk, chunk, genes): idx
                   for idx, chunk in enumerate(chunks)}
        done = 0
        for future in as_completed(futures):
            records.extend(future.result())
            done += 1
            if done % max(1, len(chunks) // 10) == 0 or done == len(chunks):
                LOG.info("  ...%d/%d chunks complete", done, len(chunks))

    return pd.DataFrame.from_records(records)


# --------------------------------------------------------------------------- #
# Output
# --------------------------------------------------------------------------- #

OUTPUT_COLUMN_ORDER = [
    "pos", "category",
    "id1", "locus_tag1", "old_locus_tag1", "start1", "end1", "strand1",
    "rel_pos_on_gene1", "gene1_len", "pct_of_gene1", "dist_from_gene1_end",
    "id2", "locus_tag2", "old_locus_tag2", "start2", "end2", "strand2",
    "dist_to_gene2_start",
]


def write_results(results: pd.DataFrame, out_path: Path) -> None:
    """Write the classified insertion-site table to CSV in a stable column order."""
    ordered_cols = [c for c in OUTPUT_COLUMN_ORDER if c in results.columns]
    remaining = [c for c in results.columns if c not in ordered_cols]
    results = results[ordered_cols + remaining].sort_values("pos")
    results.to_csv(out_path, index=False, quoting=csv.QUOTE_MINIMAL)
    LOG.info("Wrote %d rows to %s", len(results), out_path)


def write_summary(results: pd.DataFrame, out_path: Path) -> None:
    """Write a small summary of intragenic vs intergenic insertion counts."""
    counts = results["category"].value_counts()
    n_sites = results["pos"].nunique()
    with open(out_path, "w") as fh:
        fh.write("TA insertion site profiling summary\n")
        fh.write("====================================\n")
        fh.write(f"Total distinct TA sites classified: {n_sites}\n")
        for cat, n in counts.items():
            fh.write(f"  {cat}: {n} row(s)\n")
    LOG.info("Wrote summary to %s", out_path)


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ta_profiler.py",
        description=(
            "Classify Tn-seq TA insertion sites as intragenic or intergenic "
            "against a GFF gene annotation, and report the relevant gene "
            "context for each site."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
EXAMPLES
--------
  # Basic run, serial
  python3 ta_profiler.py --gff annotation.tsv --hits ta_hits.csv

  # Parallel run using 8 worker processes, custom output directory
  python3 ta_profiler.py --gff annotation.tsv --hits ta_hits.csv \\
      --outdir results/ --workers 8

  # TA hit file has a header column named 'Tnhits'
  python3 ta_profiler.py --gff annotation.tsv --hits ta_hits.csv \\
      --hits-column Tnhits

  # Also keep CDS rows (excluded by default) in the annotation
  python3 ta_profiler.py --gff annotation.tsv --hits ta_hits.csv \\
      --exclude-features region

  # Verbose logging
  python3 ta_profiler.py --gff annotation.tsv --hits ta_hits.csv -v
""",
    )

    parser.add_argument(
        "--gff", required=True, type=Path,
        help="Path to the GFF-style gene annotation TSV file "
             "(9 columns: seqname, source, feature, start, end, score, "
             "strand, frame, attribute; no header row).",
    )
    parser.add_argument(
        "--hits", required=True, type=Path,
        help="Path to a CSV file listing TA transposon insertion site "
             "coordinates (one integer position per row).",
    )
    parser.add_argument(
        "--hits-column", default=None,
        help="Name of the column to read positions from, if --hits has a "
             "header row. If omitted, the first column is used.",
    )
    parser.add_argument(
        "--outdir", type=Path, default=Path("."),
        help="Directory to write output files into (created if needed). "
             "Default: current directory.",
    )
    parser.add_argument(
        "--prefix", default="ta_profile",
        help="Filename prefix for output files. Default: 'ta_profile'.",
    )
    parser.add_argument(
        "--exclude-features", nargs="*", default=list(DEFAULT_EXCLUDED_FEATURES),
        metavar="FEATURE",
        help=f"GFF 'feature' values to drop before classification. "
             f"Default: {list(DEFAULT_EXCLUDED_FEATURES)}. Pass with no "
             f"arguments to disable filtering entirely.",
    )
    parser.add_argument(
        "--workers", type=int, default=1,
        help="Number of parallel worker processes to use for insertion-site "
             "classification. Default: 1 (serial). Use e.g. --workers 0 to "
             "auto-detect and use all available CPU cores.",
    )
    parser.add_argument(
        "--chunk-size", type=int, default=200,
        help="Number of TA sites handed to each worker process per task, "
             "when running in parallel. Larger chunks reduce inter-process "
             "overhead; smaller chunks improve load balancing. Default: 200.",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Enable debug-level logging.",
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}",
    )
    return parser


def resolve_worker_count(requested: int) -> int:
    """--workers 0 means 'use all CPUs'; otherwise use the value as-is."""
    if requested == 0:
        import os
        return os.cpu_count() or 1
    return max(1, requested)


def main(argv: Optional[list[str]] = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )

    if not args.gff.exists():
        LOG.error("GFF file not found: %s", args.gff)
        return 1
    if not args.hits.exists():
        LOG.error("TA hits file not found: %s", args.hits)
        return 1

    args.outdir.mkdir(parents=True, exist_ok=True)
    workers = resolve_worker_count(args.workers)

    start_time = time.time()

    genes = load_gff(args.gff, args.exclude_features)
    if genes.empty:
        LOG.error("No annotation rows remain after filtering -- check "
                   "--exclude-features and the GFF file contents.")
        return 1

    hits = load_ta_hits(args.hits, args.hits_column)
    if not hits:
        LOG.error("No TA insertion sites were loaded -- check --hits / "
                   "--hits-column.")
        return 1

    results = classify_all_hits(hits, genes, workers=workers, chunk_size=args.chunk_size)

    results_path = args.outdir / f"{args.prefix}.csv"
    summary_path = args.outdir / f"{args.prefix}_summary.txt"
    write_results(results, results_path)
    write_summary(results, summary_path)

    elapsed = time.time() - start_time
    LOG.info("Done in %.2fs (workers=%d)", elapsed, workers)
    return 0


if __name__ == "__main__":
    sys.exit(main())
