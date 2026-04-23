#!/usr/bin/env python3
"""
Trim GenBank records to the aligned region of each gene cluster.

This is a a custom addition to this fork of clinker.
"""

import logging
from pathlib import Path
from Bio import SeqIO

from clinker.align import Globaligner

LOG = logging.getLogger(__name__)


# =============================================================================
def get_aligned_spans(
    globaligner: Globaligner,
    padding: int = 0,
) -> dict[str, tuple[int, int]]:
    """For each locus in each cluster, find the coordinate span
    covered by genes that participate in at least one alignment link.

    ---
    Args:
        globaligner (Globaligner): completed Globaligner instance
        padding (int): bp to add on either side of the aligned span
    Returns:
        dict mapping locus.name -> (trim_start, trim_end)
    """
    linked_gene_uids = set()
    for alignment in globaligner.alignments.values():
        for link in alignment.links:
            linked_gene_uids.add(link.query.uid)
            linked_gene_uids.add(link.target.uid)

    spans = {}
    for cluster in globaligner.clusters.values():
        for locus in cluster.loci:
            linked_genes = [g for g in locus.genes if g.uid in linked_gene_uids]
            if not linked_genes:
                LOG.warning(
                    "Locus %s in cluster %s has no aligned genes, skipping trim",
                    locus.name,
                    cluster.name,
                )
                continue
            trim_start = max(0, min(g.start for g in linked_genes) - padding)
            trim_end = min(locus.end, max(g.end for g in linked_genes) + padding)
            LOG.info(
                "  %s (%s): trimming to %d-%d (padding=%d)",
                locus.name,
                cluster.name,
                trim_start,
                trim_end,
                padding,
            )
            spans[locus.name] = (trim_start, trim_end)
    return spans


# =============================================================================
def trim_genbank_file(
    input_path: str | Path,
    spans: dict[str, tuple[int, int]],
    output_path: str | Path,
    force: bool = False,
) -> None:
    """Trim records in a GenBank file to the aligned spans and write output.
    ---
    Args:
        input_path (str|Path): path to input GenBank file
        spans (dict): mapping of record.name -> (start, end)
        output_path (str|Path): path to write trimmed GenBank file
        force (bool): overwrite output if it already exists

    Returns:
        None
    """
    input_path = Path(input_path)
    output_path = Path(output_path)

    if output_path.exists() and not force:
        LOG.warning(
            "Output %s already exists, skipping (use --force to overwrite)", output_path
        )
        return

    records = []
    with input_path.open() as fh:
        for record in SeqIO.parse(fh, "genbank"):
            if record.name in spans:
                start, end = spans[record.name]
                LOG.info("  Trimming %s to %d-%d", record.name, start, end)
                records.append(record[start:end])
            else:
                records.append(record)

    with output_path.open("w") as fh:
        SeqIO.write(records, fh, "genbank")
    LOG.info("  Written: %s", output_path)


# =============================================================================
def trim_cluster_files(
    paths: list[str | Path],
    globaligner: Globaligner,
    output_dir: str | Path,
    padding: int = 0,
    suffix: str = ".trimmed",
    force: bool = False,
) -> None:
    """Trim all input GenBank files based on alignment results.
    ---
    Args:
        paths (list): input file paths (same list passed to parse_files)
        globaligner (Globaligner): completed Globaligner instance
        output_dir (str|Path): directory to write trimmed files into
        padding (int): bp padding on either side of aligned region
        suffix (str): suffix appended to each output filename stem
        force (bool): overwrite existing output files

    Returns:
        None
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    LOG.info("Calculating aligned spans (padding=%d bp)...", padding)
    spans = get_aligned_spans(globaligner, padding=padding)

    LOG.info("Writing trimmed GenBank files to: %s", output_dir)
    for path in paths:
        path = Path(path)
        output_path = output_dir / (path.stem + suffix + path.suffix)
        trim_genbank_file(path, spans, output_path, force=force)
