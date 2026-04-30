#!/usr/bin/env python3

"""
Command line interface

Original: Cameron Gilchrist
Forked: James Lingford
"""
# TODO:
# - [ ] save globaligner object as pickle file
# - [ ] add arg to read in globaligner aligner object
# - [ ] add scale bar and colour bar to svg
# - [ ] change connections b/w gene arrows to shade of grey
# - [ ] get annotations from genbank file and add them as color labels to svg
# - [ ] add better colours to svg gene arrows
# - [ ] clean up args and control flow

import argparse
import logging
import csv

from pathlib import Path
from collections import defaultdict
from typing import TextIO, Dict, List, Any
from dataclasses import dataclass

from clinker import __version__, align
from clinker.plot import plot_clusters, plot_data
from clinker.classes import find_files, parse_files
from clinker.trim import trim_cluster_files
from clinker.svg import render_svg


# =============================================================================
# Global config
# =============================================================================

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s - %(message)s",
    datefmt="%H:%M:%S",
)
LOG = logging.getLogger(__name__)


# =============================================================================
# CLI parser
# =============================================================================
@dataclass
class Args:
    files: Path
    file_list: Path
    session: Path
    ranges: str
    gene_functions: TextIO
    colour_map: TextIO
    dont_set_origin: bool
    as_separate_clusters: bool
    no_align: bool
    identity: float
    cpu: int
    json_indent: int
    force: bool
    output: Any  # TODO: ?
    plot: Any  # TODO: ?
    delimiter: str
    decimals: int
    hide_link_headers: bool
    hide_aln_headers: bool
    matrix_out: Any  # TODO: ?
    use_file_order: bool
    trim: Path  # TODO: check
    trim_padding: int
    trim_suffix: str
    svg: bool
    anchor_target_genes: bool


def get_parser() -> Args:
    """Creates an ArgumentParser object."""
    parser = argparse.ArgumentParser(
        "clinker",
        description="clinker: Automatic creation of publication-ready"
        " gene cluster comparison figures.\n\n"
        "clinker generates gene cluster comparison figures from GenBank files."
        " It performs pairwise local or global alignments between every sequence"
        " in every unique pair of clusters and generates interactive, to-scale comparison figures"
        " using the clustermap.js library.",
        # epilog="Example usage\n-------------\n"
        # "Align clusters, plot results and print scores to screen:\n"
        # "  $ clinker files/*.gbk\n\n"
        # "Only save gene-gene links when identity is over 50%:\n"
        # "  $ clinker files/*.gbk -i 0.5\n\n"
        # "Save an alignment session for later:\n"
        # "  $ clinker files/*.gbk -s session.json\n\n"
        # "Save alignments to file, in comma-delimited format, with 4 decimal places:\n"
        # '  $ clinker files/*.gbk -o alignments.csv -dl "," -dc 4\n\n'
        # "Generate visualisation:\n"
        # "  $ clinker files/*.gbk -p\n\n"
        # "Save visualisation as a static HTML document:\n"
        # "  $ clinker files/*.gbk -p plot.html\n\n"
        # "Cameron Gilchrist, 2020",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"clinker v{__version__}",
    )

    inputs = parser.add_argument_group("Input options")
    inputs.add_argument(
        "files",
        nargs="*",
        help="Gene cluster GenBank files",
    )
    inputs.add_argument(
        "-fl",
        "--file_list",
        type=Path,
        required=False,
        help="Path to text file containing list of GenBank filepaths."
        " Filepaths should be absolute paths and one per line.",
    )
    inputs.add_argument(
        "-r",
        "--ranges",
        nargs="+",
        required=False,
        help="Scaffold extraction ranges. If a range is specified, only features within"
        " the range will be extracted from the scaffold. Ranges should be formatted"
        " like: scaffold:start-end (e.g. scaffold_1:15000-40000)",
    )
    inputs.add_argument(
        "-gf",
        "--gene_functions",
        type=Path,
        required=False,
        help="2-column CSV file containing gene functions, used to build gene groups"
        " from same function instead of sequence similarity (e.g. GENE_001,PKS-NRPS).",
    )
    inputs.add_argument(
        "-cm",
        "--colour_map",
        type=Path,
        required=False,
        help="2-column CSV file containing gene functions and colours (e.g. GENE_001,#FF0000).",
    )
    inputs.add_argument(
        "-dso",
        "--dont_set_origin",
        action="store_true",
        help="Don't fix features which cross the origin in circular sequences (GenBank format only)",
    )
    inputs.add_argument(
        "-asc",
        "--as_separate_clusters",
        action="store_true",
        help="Records will be parsed into separate clusters. "
        "Enable this option when the GenBank file you downloaded from NCBI contains multiple sequences.",
    )

    alignment = parser.add_argument_group("Alignment options")
    alignment.add_argument(
        "-na",
        "--no_align",
        action="store_true",
        help="Do not align clusters",
    )
    alignment.add_argument(
        "-i",
        "--identity",
        help="Minimum alignment sequence identity [default: 0.3]",
        type=float,
        default=0.3,
    )
    alignment.add_argument(
        "-C",
        "--cpu",
        help="Number of alignments to run in parallel (0 to use the number of CPUs) [default: 0]",
        type=int,
        default=0,
    )

    output = parser.add_argument_group("Output options")
    output.add_argument(
        # "-s",
        "--session",
        help="Path to clinker session",
    )
    output.add_argument(
        "-ji",
        "--json_indent",
        type=int,
        help="Number of spaces to indent JSON [default: none]",
    )
    output.add_argument(
        "-f",
        "--force",
        action="store_true",
        help="Overwrite previous output file",
    )
    output.add_argument(
        "-o",
        "--output",
        help="Save alignments to file",
    )
    output.add_argument(
        "-p",
        "--plot",
        nargs="?",
        const=True,
        default=False,
        help="Plot cluster alignments using clustermap.js. If a path is given,"
        " clinker will generate a portable HTML file at that path. Otherwise,"
        " the plot will be served dynamically using Python's HTTP server.",
    )
    output.add_argument(
        "-dl",
        "--delimiter",
        type=str,
        default=",",
        help="Character to delimit output by [default: human readable]",
    )
    output.add_argument(
        "-dc",
        "--decimals",
        help="Number of decimal places in output [default: 2]",
        type=int,
        default=2,
    )
    output.add_argument(
        "-hl",
        "--hide_link_headers",
        help="Hide alignment column headers",
        action="store_true",
    )
    output.add_argument(
        "-ha",
        "--hide_aln_headers",
        help="Hide alignment cluster name headers",
        action="store_true",
    )
    output.add_argument(
        "-mo",
        "--matrix_out",
        help="Save cluster similarity matrix to file",
    )

    viz = parser.add_argument_group("Visualisation options")
    viz.add_argument(
        "-ufo",
        "--use_file_order",
        action="store_true",
        help="Display clusters in order of input files",
    )
    viz.add_argument(
        "-s",
        "--svg",
        type=Path,
        metavar="SVG_PATH",
        required=False,
        help="Path to write an svg image of the clinker results."
        " Requires the -p|--plot option to be passed as well.",
    )
    viz.add_argument(
        "-an",
        "--anchor_target_genes",
        action="store_true",
        help="Visually align all target genes on the same x-axis position and point in the same direction in the svg."
        " NOTE: assumes target genes can be parsed from the file list",
    )

    trim = parser.add_argument_group("Trimming options")
    trim.add_argument(
        "-t",
        "--trim",
        help="Trim GenBank records to the aligned region and write to the given directory",
        type=Path,
        default=None,
        metavar="DIR",
    )
    trim.add_argument(
        "-tp",
        "--trim_padding",
        help="Padding (bp) to add on either side of the aligned region when trimming [default: 0]",
        type=int,
        default=0,
    )
    trim.add_argument(
        "-ts",
        "--trim_suffix",
        help="Suffix to append to trimmed output filenames [default: .trimmed.gbk]",
        default=".trimmed",
    )

    args = Args(**vars(parser.parse_args()))

    # check inputs
    if not (args.files or args.file_list):
        parser.error(
            "No input files provided. Please pass files as positional args or with --file_list"
        )
    if args.files and args.file_list:
        parser.error(
            "Please provide only one form of GenBank input: either positional args or --file_list"
        )

    return args


# =============================================================================
# Input processing util funcs
# =============================================================================
def parse_range(string):
    """Extracts the scaffold name, start and end of a scaffold range string.

    Expects the format scaffold:start-stop (e.g. scaffold_1:100-3000).

    Args:
        string (str): Scaffold range string
    Returns:
        scaffold (str): Scaffold name
        start (int): Range start
        end (int): Range end
    Raises:
        ValueError: If range is in invalid format (can't split by '-' or ':')
        TypeError: Range values are not integers
    """
    try:
        scaffold, coordinates = string.split(":")
        start, end = coordinates.split("-")
    except ValueError as e:
        raise ValueError(
            "Expected format scaffold:start-stop (e.g. scaf_1:100-3000)"
        ) from e
    if not start.isdigit() or not end.isdigit():
        raise TypeError("Expected range values to be type int")
    return scaffold, int(start), int(end)


def parse_ranges(strings):
    ranges = {}
    for string in strings:
        try:
            scaffold, start, end = parse_range(string)
        except (ValueError, TypeError):
            LOG.exception("Failed to read range, please check it's in the right format")
            raise
        ranges[scaffold] = (start, end)
    return ranges


def parse_gene_functions(fp: TextIO) -> Dict[str, List[str]]:
    """Parses gene functions from a table.

    Gene        Function
    GENE_001    Cytochrome P450
    GENE_002    Methyltransferase
    ...
    """
    functions = defaultdict(list)
    for gene, function in csv.reader(fp):
        functions[function].append(gene)
    return functions


def parse_colour_map(fp: TextIO) -> Dict[str, str]:
    """Parses colours from a table.

    Function   Colour
    Type I     #FF0000
    Type II    #000000
    ...
    """
    colours = {}
    for function, colour in csv.reader(fp):
        colours[function] = colour
    return colours


def parse_target_anchor_genes(file_list: list[str]) -> dict[str, str]:
    """Parses target genes from input file list.
    Should be of the dict format:

    anchor_map = {"cluster": "gene_target"}

    where cluster is the filepath stem, and gene_target is the final field after "___"
    """
    anchor_map: dict[str, str] = {}
    for file in file_list:
        file = Path(file)
        cluster = file.stem
        # guaymas
        # gene_target = file.name.removesuffix("".join(file.suffixes)).split("___")[-1]
        # globdb renamed
        gene_target = "___".join(
            file.name.removesuffix("".join(file.suffixes)).split("___")[-2::]
        )
        LOG.info(f"Using gene_target ID: {gene_target}")
        # gene_target = file.name.removesuffix("".join(file.suffixes))
        anchor_map.update({cluster: gene_target})
    return anchor_map


# ==============================================================================
# Main functions
# ==============================================================================
def clinker(
    args: Args,
    files=None,
    file_list=None,
    session=None,
    identity=0.3,
    delimiter=None,
    decimals=2,
    plot=None,
    output=None,
    force=False,
    no_align=False,
    hide_link_headers=False,
    hide_alignment_headers=False,
    use_file_order=False,
    json_indent=None,
    jobs=None,
    ranges=None,
    matrix_out=None,
    gene_functions=None,
    colour_map=None,
    set_origin=False,
    as_separate_clusters=False,
    trim=None,  # Trimming
    trim_padding=0,  # Trimming
    trim_suffix=".trimmed",  # Trimming
):
    """Entry point for running the script."""
    LOG.info("Starting clinker")

    # load files from file list
    if file_list:
        with open(file_list) as f:
            files = [line.strip() for line in f if line.strip()]

    load_session = session and Path(session).exists()

    # Parse range strings, if any specified
    if ranges:
        ranges = parse_ranges(ranges)

    # Parse any gene functions for grouping, if specified
    if gene_functions:
        gene_functions = parse_gene_functions(gene_functions)
    if colour_map:
        colour_map = parse_colour_map(colour_map)

    if load_session:
        LOG.info("Loading session from: %s", session)
        with open(session) as fp:
            try:
                globaligner = align.Globaligner.from_json(fp)
            except Exception:
                LOG.exception(
                    "Failed to load session, is '%s' a clinker session?", session
                )
        if files:
            paths = find_files(files)
            if not paths:
                LOG.error("No files found")
                raise SystemExit
            LOG.info("Parsing files:")
            clusters = parse_files(paths, ranges=ranges, set_origin=set_origin)

            LOG.info("Adding clusters to loaded session and aligning")
            globaligner.add_clusters(*clusters)
            globaligner.align_stored_clusters(cutoff=identity, jobs=jobs)
            globaligner.build_gene_groups(functions=gene_functions, colours=colour_map)
            load_session = False
    else:
        # Parse files, generate objects
        paths = find_files(files)
        if not paths:
            # Allow no files, so that user can generate a blank clinker web app
            # and load in previously saved figure data
            if plot:
                LOG.info("Opening empty clinker web app...")
                plot_data(
                    dict(clusters=[], links=[], groups=[]),
                    output=None if plot is True else plot,
                )
            else:
                LOG.error("No files provided!")
                raise SystemExit
        LOG.info("Parsing files:")
        clusters = parse_files(
            paths,
            ranges=ranges,
            set_origin=set_origin,
            as_separate_clusters=as_separate_clusters,
        )

        # Align all clusters
        if no_align:
            globaligner = align.Globaligner()
            globaligner.add_clusters(*clusters)
            globaligner.build_gene_groups(functions=gene_functions, colours=colour_map)
        elif len(clusters) == 1:
            globaligner = align.align_clusters(clusters[0], jobs=1)
        else:
            LOG.info("Starting cluster alignments")
            globaligner = align.align_clusters(*clusters, cutoff=identity, jobs=jobs)
            globaligner.build_gene_groups(functions=gene_functions, colours=colour_map)

    if globaligner.alignments:
        LOG.info("Generating results summary...")
        summary = globaligner.format(
            delimiter=delimiter,
            decimals=decimals,
            link_headers=not hide_link_headers,
            alignment_headers=not hide_alignment_headers,
        )
        if output:
            if output and Path(output).exists() and not force:
                print(summary)
                LOG.warning(
                    "File %s already exists but --force was not specified", output
                )
            else:
                LOG.info("Writing alignments to: %s", output)
                with open(output, "w") as fp:
                    fp.write(summary)
        else:
            print(summary)
        if matrix_out:
            LOG.info("Writing synteny matrix to: %s", matrix_out)
            # matrix = globaligner.format_matrix(normalise=True, as_distance=False)
            matrix = globaligner.format_matrix(normalise=False, as_distance=True)
            with open(matrix_out, "w") as fp:
                fp.write(matrix)

    else:
        LOG.info("No alignments were generated")

    if session and not load_session:
        LOG.info("Saving session to: %s", session)
        with open(session, "w") as fp:
            globaligner.to_json(fp, indent=json_indent)

    # Generate the SVG
    if plot:
        LOG.info("Building clustermap.js visualisation")
        if isinstance(plot, str):
            LOG.info("Writing to: %s", plot)
        plot_clusters(
            globaligner,
            output=None if plot is True else plot,
            use_file_order=use_file_order,
        )
        if args.anchor_target_genes:
            LOG.info("Parsing target genes for anchoring from file list")
            anchor_map = parse_target_anchor_genes(file_list=files)
        render_svg(
            globaligner,
            output=args.svg,
            use_file_order=use_file_order,
            show_gene_labels=False,
            identity_threshold=0.3,
            anchor_map=anchor_map,
        )

    # trim the genbank files to conserved/"linked" regions
    if trim:
        if not globaligner.alignments:
            LOG.warning("No alignments were generated, cannot trim")
        else:
            LOG.info("Trimming GenBank files...")
            trim_cluster_files(
                paths,
                globaligner,
                output_dir=trim,
                padding=trim_padding,
                suffix=trim_suffix,
                force=force,
            )

    LOG.info("Done!")
    return globaligner


# =============================================================================
def main():
    # collect cli args
    args = get_parser()

    # run clinker
    clinker(
        args=args,
        files=args.files,
        file_list=args.file_list,
        session=args.session,
        json_indent=args.json_indent,
        identity=args.identity,
        delimiter=args.delimiter,
        decimals=args.decimals,
        plot=args.plot,
        output=args.output,
        force=args.force,
        no_align=args.no_align,
        hide_link_headers=args.hide_link_headers,
        hide_alignment_headers=args.hide_aln_headers,
        use_file_order=args.use_file_order,
        jobs=args.cpu if args.cpu > 0 else None,
        ranges=args.ranges,
        matrix_out=args.matrix_out,
        gene_functions=args.gene_functions,
        colour_map=args.colour_map,
        set_origin=not args.dont_set_origin,
        as_separate_clusters=args.as_separate_clusters,
        trim=args.trim,  # Trimming
        trim_padding=args.trim_padding,  # Trimming
        trim_suffix=args.trim_suffix,  # Trimming
    )


# =============================================================================
if __name__ == "__main__":
    main()
