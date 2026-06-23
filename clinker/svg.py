#!/usr/bin/env python3

"""
Warning: all vibe coded... mileage may vary
"""

from pathlib import Path
from typing import Optional
import colorsys
from collections import Counter
import numpy as np
import seaborn as sns

# ==============================================================================
# Layout constants — all in pixels
TRACK_HEIGHT = 18  # height of gene arrow body
ARROW_HEAD = 10  # width of arrowhead
ARROW_STROKE_COLOUR = "black"
ARROW_STROKE_WIDTH = 1.0
ARROW_BODY_OPACITY = 1.0
TRACK_SPACING = 80  # vertical distance between cluster tracks
LABEL_OFFSET = 14  # px above track for gene labels
CLUSTER_LABEL_X = 10  # x position of cluster name label
RIBBON_OPACITY_MIN = 0.0
RIBBON_OPACITY_MAX = 1.0
RIBBON_COLOR = "#a6adc8"
SVG_PADDING = 20  # padding around entire figure
SCALE = 0.032  # bp -> px scaling factor
LABEL_COLUMN_WIDTH = 480  # px reserved for cluster name on the left
# legend constants
LEGEND_BOTTOM_MARGIN = 80  # extra height below last track for legend
SCALEBAR_BP = 2500  # genomic length the scale bar represents
SCALEBAR_HEIGHT = 8  # tick height in px
IDENTITY_BAR_WIDTH = 160
IDENTITY_BAR_HEIGHT = 14
# color key constants
LEGEND_SWATCH_SIZE = 14
LEGEND_ROW_HEIGHT = 20
LEGEND_HEIGHT = SCALEBAR_HEIGHT + 26  # bar height + text label below it
LEGEND_COLS = 2
LEGEND_COL_WIDTH = 400
TRACK_TO_SCALEBAR_GAP = 10
SCALEBAR_TO_COLOURKEY_GAP = 40


# ==============================================================================
def build_colour_key(
    groups: list,
    group_kofam_labels: dict[str, str],
    uid_to_colour: dict[str, str],
    x: float,
    y: float,
) -> tuple[list[str], float]:
    """
    Renders a colour key for gene groups labelled by representative kofam_id.
    Only renders groups that have at least one labelled kofam_id.
    Returns (svg_elements, total_height).
    """
    elements = []

    # Collect unique (colour, label) pairs — deduplicate by colour+label
    seen: set[tuple[str, str]] = set()
    entries: list[tuple[str, str]] = []  # (colour, label)
    for group in groups:
        label = group_kofam_labels.get(group["uid"], "")
        if not label:
            continue
        # Get colour from any member gene
        colour = next(
            (uid_to_colour[uid] for uid in group["genes"] if uid in uid_to_colour),
            None,
        )
        if colour is None:
            continue
        key = (colour, label)
        if key not in seen:
            seen.add(key)
            entries.append(key)

    if not entries:
        return [], 0.0

    # Title
    elements.append(
        f'<text x="{x:.1f}" y="{y:.1f}" '
        f'font-family="sans-serif" font-size="14" font-weight="normal">Annotation</text>'
    )
    y += 6  # gap after title

    # Swatch rows, two columns
    for i, (colour, label) in enumerate(entries):
        col = i % LEGEND_COLS
        row = i // LEGEND_COLS
        sx = x + col * LEGEND_COL_WIDTH
        sy = y + row * LEGEND_ROW_HEIGHT

        # Colour swatch
        elements.append(
            f'<rect x="{sx:.1f}" y="{sy:.1f}" '
            f'width="{LEGEND_SWATCH_SIZE}" height="{LEGEND_SWATCH_SIZE}" '
            f'fill="{colour}" stroke="black" stroke-width="0.6"/>'
        )
        # Label
        elements.append(
            f'<text x="{sx + LEGEND_SWATCH_SIZE + 5:.1f}" '
            f'y="{sy + LEGEND_SWATCH_SIZE - 3:.1f}" '
            f'font-family="sans-serif" font-size="12">{label}</text>'
        )

    n_rows = (len(entries) + LEGEND_COLS - 1) // LEGEND_COLS
    total_height = 6 + n_rows * LEGEND_ROW_HEIGHT
    return elements, total_height


# ==============================================================================
def build_kofam_colour_key(
    clusters: list,
    uid_to_colour: dict[str, str],
    kofam_to_colour: dict[str, str],
    x: float,
    y: float,
) -> tuple[list[str], float]:
    """
    Colour key for no_align mode: one swatch per unique kofam_id.
    """
    # # Collect unique (kofam_id, colour) pairs
    # seen: set[str] = set()
    # entries: list[tuple[str, str]] = []  # (label, colour)
    # for cluster in clusters:
    #     for locus in cluster["loci"]:
    #         for gene in locus["genes"]:
    #             kid = gene.get("kofam_id")
    #             colour = uid_to_colour.get(gene["uid"], "#dddddd")
    #             if kid and kid not in seen:
    #                 seen.add(kid)
    #                 entries.append((kid, colour))
    #
    # if not entries:
    #     return [], 0.0
    #
    # elements = []
    # elements.append(
    #     f'<text x="{x:.1f}" y="{y:.1f}" '
    #     f'font-family="sans-serif" font-size="14" font-weight="normal">Annotation</text>'
    # )
    # y += 6
    #
    # for i, (label, colour) in enumerate(entries):
    #     col = i % LEGEND_COLS
    #     row = i // LEGEND_COLS
    #     sx = x + col * LEGEND_COL_WIDTH
    #     sy = y + row * LEGEND_ROW_HEIGHT
    #
    #     elements.append(
    #         f'<rect x="{sx:.1f}" y="{sy:.1f}" '
    #         f'width="{LEGEND_SWATCH_SIZE}" height="{LEGEND_SWATCH_SIZE}" '
    #         f'fill="{colour}" stroke="black" stroke-width="0.6"/>'
    #     )
    #     elements.append(
    #         f'<text x="{sx + LEGEND_SWATCH_SIZE + 5:.1f}" '
    #         f'y="{sy + LEGEND_SWATCH_SIZE - 3:.1f}" '
    #         f'font-family="sans-serif" font-size="12">{label}</text>'
    #     )
    #
    # n_rows = (len(entries) + LEGEND_COLS - 1) // LEGEND_COLS
    # total_height = 6 + n_rows * LEGEND_ROW_HEIGHT
    # return elements, total_height

    # Order entries by palette position (insertion order of kofam_to_colour,
    # which was built from sorted+enumerated unique_kofams)
    entries = [
        (kid, colour)
        for kid, colour in kofam_to_colour.items()
        if kid != "None"  # optionally put "None" last
    ]
    # Append "None" at the end if present
    if "None" in kofam_to_colour:
        entries.append(("None", kofam_to_colour["None"]))

    if not entries:
        return [], 0.0

    elements = []
    elements.append(
        f'<text x="{x:.1f}" y="{y:.1f}" '
        f'font-family="sans-serif" font-size="14" font-weight="normal">Annotation</text>'
    )
    y += 6
    for i, (label, colour) in enumerate(entries):
        col = i % LEGEND_COLS
        row = i // LEGEND_COLS
        sx = x + col * LEGEND_COL_WIDTH
        sy = y + row * LEGEND_ROW_HEIGHT
        elements.append(
            f'<rect x="{sx:.1f}" y="{sy:.1f}" '
            f'width="{LEGEND_SWATCH_SIZE}" height="{LEGEND_SWATCH_SIZE}" '
            f'fill="{colour}" stroke="black" stroke-width="0.6"/>'
        )
        elements.append(
            f'<text x="{sx + LEGEND_SWATCH_SIZE + 5:.1f}" '
            f'y="{sy + LEGEND_SWATCH_SIZE - 3:.1f}" '
            f'font-family="sans-serif" font-size="12">{label}</text>'
        )
    n_rows = (len(entries) + LEGEND_COLS - 1) // LEGEND_COLS
    total_height = 6 + n_rows * LEGEND_ROW_HEIGHT
    return elements, total_height


# ==============================================================================
def build_group_kofam_labels(clusters: list, groups: list) -> dict[str, str]:
    """
    For each group, finds the most frequent non-None kofam_id among
    member genes. Returns group_uid -> kofam_label (or "" if none found).
    """
    # Build gene_uid -> kofam_id lookup from cluster data
    uid_to_kofam: dict[str, str] = {}
    for cluster in clusters:
        for locus in cluster["loci"]:
            for gene in locus["genes"]:
                kid = gene.get("kofam_id")
                if kid:
                    uid_to_kofam[gene["uid"]] = kid

    group_labels: dict[str, str] = {}
    for group in groups:
        kofam_hits = [
            uid_to_kofam[uid] for uid in group["genes"] if uid in uid_to_kofam
        ]
        if kofam_hits:
            most_common, _ = Counter(kofam_hits).most_common(1)[0]
            group_labels[group["uid"]] = most_common
        else:
            group_labels[group["uid"]] = ""

    return group_labels


# ==============================================================================
def build_colour_map(groups: list) -> dict:
    """
    Builds a gene_uid -> colour dict from groups.
    Auto-generates colours for groups where colour is None.
    """
    # Count how many groups need auto-colour
    auto_groups = [g for g in groups if not g.get("colour")]
    auto_colours = generate_colours(len(auto_groups))
    auto_iter = iter(auto_colours)

    uid_to_colour = {}
    for group in groups:
        colour = group.get("colour") or next(auto_iter)
        for gene_uid in group["genes"]:
            uid_to_colour[gene_uid] = colour
    return uid_to_colour


# ==============================================================================
def build_kofam_colour_map(clusters: list) -> tuple[dict[str, str], dict[str, str]]:
    """
    Builds a gene_uid -> colour dict based on shared kofam_id annotations.
    Genes with the same kofam_id get the same colour.
    Unannotated genes get a neutral grey.
    """
    # NONE_COLOUR = "#838ba7"
    NONE_COLOUR = "#a6adc8"

    uid_to_kofam: dict[str, str] = {}
    for cluster in clusters:
        for locus in cluster["loci"]:
            for gene in locus["genes"]:
                kid = gene.get("kofam_id")
                if kid:
                    uid_to_kofam[gene["uid"]] = kid

    # Exclude "None" strings from palette generation
    unique_kofams = sorted(k for k in set(uid_to_kofam.values()) if k != "None")
    palette = generate_colours(len(unique_kofams))
    kofam_to_colour: dict[str, str] = dict(zip(unique_kofams, palette))
    kofam_to_colour["None"] = NONE_COLOUR  # explicit override

    uid_to_colour: dict[str, str] = {}
    for cluster in clusters:
        for locus in cluster["loci"]:
            for gene in locus["genes"]:
                kid = uid_to_kofam.get(gene["uid"])
                uid_to_colour[gene["uid"]] = kofam_to_colour.get(kid, "#dddddd")

    return uid_to_colour, kofam_to_colour


# ==============================================================================
def generate_colours(n: int) -> list[str]:
    """Generates n visually distinct colours using HSL spacing."""
    # colours = []
    # for i in range(n):
    #     hue = i / n
    #     # Use fixed saturation/lightness similar to clustermap.js defaults
    #     r, g, b = colorsys.hls_to_rgb(hue, 0.6, 0.7)
    #     colours.append(
    #         "#{:02x}{:02x}{:02x}".format(int(r * 255), int(g * 255), int(b * 255))
    #     )
    # return colours
    cmap = sns.color_palette("Spectral", as_cmap=True)
    colours = []
    for i in range(n):
        t = i / max(n - 1, 1)  # evenly space from 0.0 to 1.0
        r, g, b, *_ = cmap(t)
        colours.append(
            "#{:02x}{:02x}{:02x}".format(int(r * 255), int(g * 255), int(b * 255))
        )
    return colours


# ==============================================================================
def find_anchor_gene_per_cluster(
    anchor_map: dict[str, str],
    clusters: list,
) -> dict:
    """
    Builds a cluster_uid -> anchor gene_uid mapping from an explicit
    dict of {cluster_name: anchor_gene_label}.

    Args:
        anchor_map: dict mapping cluster name -> anchor gene label
                    e.g. {"NiFe_Group_4f__GB30501": "GB30501_115"}
        clusters: list of cluster dicts from globaligner.to_data()

    Returns:
        dict: cluster_uid -> anchor gene_uid
    """
    cluster_to_anchor = {}

    for cluster in clusters:
        anchor_label = anchor_map.get(cluster["name"])
        if anchor_label is None:
            continue  # cluster not in map, will be left as-is

        for locus in cluster["loci"]:
            for gene in locus["genes"]:
                if gene.get("label") == anchor_label:
                    cluster_to_anchor[cluster["uid"]] = gene["uid"]
                    break

    missing = [
        name for name in anchor_map if not any(c["name"] == name for c in clusters)
    ]
    if missing:
        import logging

        logging.getLogger(__name__).warning(
            "Anchor map contains cluster names not found in data: %s", missing
        )

    return cluster_to_anchor


# def find_anchor_gene_per_cluster(
#     anchor_label: str,
#     clusters: list,
#     groups: list,
# ) -> dict:
#     """
#     Given an anchor gene label, finds the corresponding gene UID
#     in each cluster via shared group membership.
#
#     Returns dict: cluster_uid -> gene uid of the anchor gene
#     """
#     # Build gene uid -> label lookup from cluster data
#     uid_to_label = {}
#     uid_to_cluster = {}
#     for cluster in clusters:
#         for locus in cluster["loci"]:
#             for gene in locus["genes"]:
#                 uid_to_label[gene["uid"]] = gene.get("label", "")
#                 uid_to_cluster[gene["uid"]] = cluster["uid"]
#
#     # Find the group containing the anchor gene
#     anchor_uid = None
#     for uid, label in uid_to_label.items():
#         if label == anchor_label:
#             anchor_uid = uid
#             break
#
#     if anchor_uid is None:
#         raise ValueError(f"Anchor gene '{anchor_label}' not found in any cluster")
#
#     anchor_group = None
#     for group in groups:
#         if anchor_uid in group["genes"]:
#             anchor_group = group
#             break
#
#     if anchor_group is None:
#         raise ValueError(f"Anchor gene '{anchor_label}' not found in any group")
#
#     # For each cluster, find which gene in this group belongs to it
#     cluster_to_anchor = {}
#     for gene_uid in anchor_group["genes"]:
#         cluster_uid = uid_to_cluster.get(gene_uid)
#         if cluster_uid:
#             cluster_to_anchor[cluster_uid] = gene_uid
#
#     return cluster_to_anchor


# ==============================================================================
def normalise_locus(
    locus: dict,
    anchor_uid: str,
) -> tuple[dict, float]:
    """
    Flips and translates a locus so the anchor gene is:
      - always on strand +1 (pointing right)
      - centred at x=0 in locus coordinate space

    Returns:
        normalised locus dict (deep copy, original untouched)
        anchor_centre: the original centre of the anchor gene (for debugging)
    """
    import copy

    locus = copy.deepcopy(locus)

    # Find anchor gene in this locus
    anchor = None
    for gene in locus["genes"]:
        if gene["uid"] == anchor_uid:
            anchor = gene
            break

    if anchor is None:
        # Anchor not in this locus, return as-is with no translation
        return locus, 0.0

    anchor_centre = (anchor["start"] + anchor["end"]) / 2
    anchor_strand = anchor["strand"]

    # Step 1: flip if anchor is on minus strand
    if anchor_strand == -1:
        locus_len = locus["end"] - locus["start"]
        locus_mid = locus["start"] + locus_len / 2

        for gene in locus["genes"]:
            # Mirror coordinates around locus midpoint
            new_start = 2 * locus_mid - gene["end"]
            new_end = 2 * locus_mid - gene["start"]
            gene["start"] = new_start
            gene["end"] = new_end
            gene["strand"] = -gene["strand"]

        # Recompute anchor centre after flip
        anchor_centre = (anchor["start"] + anchor["end"]) / 2

    # Step 2: translate so anchor centre is at locus["start"]
    # (the render loop will then place this at track_x0)
    offset = anchor_centre - locus["start"]
    for gene in locus["genes"]:
        gene["start"] -= offset
        gene["end"] -= offset

    # Also shift locus bounds
    locus_span = locus["end"] - locus["start"]
    locus["start"] = 0
    locus["end"] = locus_span

    return locus, anchor_centre


# ==============================================================================
def identity_to_opacity(identity: float) -> float:
    """Maps alignment identity to ribbon opacity."""
    return RIBBON_OPACITY_MIN + identity * (RIBBON_OPACITY_MAX - RIBBON_OPACITY_MIN)


# def gene_colour(gene_uid: str, groups: list) -> str:
#     """Looks up colour for a gene UID from the group list."""
#     for group in groups:
#         if gene_uid in group["genes"]:
#             return group["colour"] or "#cccccc"
#     return "#dddddd"  # ungrouped genes


# ==============================================================================
def gene_to_px(
    gene_start: int,
    gene_end: int,
    locus_start: int,
    track_x0: int,
) -> tuple:
    """Converts gene genomic coordinates to pixel x positions."""
    x1 = track_x0 + (gene_start - locus_start) * SCALE
    x2 = track_x0 + (gene_end - locus_start) * SCALE
    return x1, x2


# def gene_to_px(
#     gene_start: int,
#     gene_end: int,
#     locus_start: int,
#     track_x0: int,
#     scale: float = SCALE,
# ) -> tuple:
#     x1 = track_x0 + (gene_start - locus_start) * scale
#     x2 = track_x0 + (gene_end - locus_start) * scale
#     return x1, x2


# ==============================================================================
def arrow_polygon(
    x1: float,
    x2: float,
    y: float,
    strand: int,
    height: float = TRACK_HEIGHT,
    head: float = ARROW_HEAD,
) -> str:
    """
    Generates SVG polygon points for a gene arrow.

    Arrow body is a pentagon: flat back, pointed front.
    Strand +1 = pointing right, -1 = pointing left.
    y is the top of the arrow.
    """
    mid_y = y + height / 2
    bot_y = y + height
    width = x2 - x1

    if width <= 0:
        return ""

    # Clamp head width so very small genes still render
    hw = min(head, width)

    if strand == 1:
        # Points: back-top, shoulder-top, tip, shoulder-bot, back-bot
        pts = [
            (x1, y),
            (x2 - hw, y),
            (x2, mid_y),
            (x2 - hw, bot_y),
            (x1, bot_y),
        ]
    else:
        pts = [
            (x2, y),
            (x1 + hw, y),
            (x1, mid_y),
            (x1 + hw, bot_y),
            (x2, bot_y),
        ]

    return " ".join(f"{px:.1f},{py:.1f}" for px, py in pts)


# ==============================================================================
def ribbon_path(
    x1: float,
    y1: float,
    x2: float,
    y2: float,
    w1: float,
    w2: float,
) -> str:
    """
    Generates an SVG cubic bezier ribbon path between two gene arrows.

    (x1, y1) is the bottom-centre of the source gene arrow.
    (x2, y2) is the top-centre of the target gene arrow.
    w1, w2 are the half-widths of the source and target genes.
    """
    # Control points sit halfway between the two tracks
    cy = (y1 + y2) / 2

    # Four corners of the ribbon
    tl = (x1 - w1, y1)
    tr = (x1 + w1, y1)
    bl = (x2 - w2, y2)
    br = (x2 + w2, y2)

    # # Four corners of the ribbon
    # tl = (x1, y1)
    # tr = (x1, y1)
    # bl = (x2, y2)
    # br = (x2, y2)

    return (
        f"M {tl[0]:.1f},{tl[1]:.1f} "
        f"C {tl[0]:.1f},{cy:.1f} {bl[0]:.1f},{cy:.1f} {bl[0]:.1f},{bl[1]:.1f} "
        f"L {br[0]:.1f},{br[1]:.1f} "
        f"C {br[0]:.1f},{cy:.1f} {tr[0]:.1f},{cy:.1f} {tr[0]:.1f},{tr[1]:.1f} "
        f"Z"
    )


# ==============================================================================
def build_legend(
    svg_width: float,
    # svg_height: float,
    legend_y: float,
    scale: float,
) -> list[str]:
    """Renders scale bar and identity gradient legend."""
    elements = []

    legend_x = SVG_PADDING + LABEL_COLUMN_WIDTH
    # legend_y = svg_height - LEGEND_BOTTOM_MARGIN + 10

    # ------------------------------------------------------------------
    # Scale bar
    # ------------------------------------------------------------------
    bar_px = SCALEBAR_BP * scale
    tick_y_top = legend_y
    tick_y_bot = legend_y + SCALEBAR_HEIGHT
    mid_y = legend_y + SCALEBAR_HEIGHT / 2

    # Horizontal line + two end ticks
    elements.append(
        f'<line x1="{legend_x:.1f}" y1="{mid_y:.1f}" '
        f'x2="{legend_x + bar_px:.1f}" y2="{mid_y:.1f}" '
        f'stroke="black" stroke-width="1.5"/>'
    )
    for tx in (legend_x, legend_x + bar_px):
        elements.append(
            f'<line x1="{tx:.1f}" y1="{tick_y_top:.1f}" '
            f'x2="{tx:.1f}" y2="{tick_y_bot:.1f}" '
            f'stroke="black" stroke-width="1.5"/>'
        )

    # Label — format as kb if >= 1000 bp
    if SCALEBAR_BP >= 1000:
        label = f"{SCALEBAR_BP / 1000:g} kb"
    else:
        label = f"{SCALEBAR_BP} bp"

    label_x = legend_x + bar_px / 2
    label_y = tick_y_bot + 13
    elements.append(
        f'<text x="{label_x:.1f}" y="{label_y:.1f}" '
        f'font-family="sans-serif" font-size="11" '
        f'text-anchor="middle">{label}</text>'
    )

    # ------------------------------------------------------------------
    # Identity gradient bar
    # ------------------------------------------------------------------
    grad_x = legend_x + bar_px + 40
    grad_y = legend_y
    grad_id = "identityGrad"

    elements.append(
        f"<defs>"
        f'<linearGradient id="{grad_id}" x1="0" y1="0" x2="1" y2="0">'
        # f'<stop offset="0%"   stop-color="white"  stop-opacity="1"/>'
        # f'<stop offset="100%" stop-color="black"  stop-opacity="1"/>'
        f'<stop offset="0%"   stop-color="{RIBBON_COLOR}"  stop-opacity="0"/>'
        f'<stop offset="100%" stop-color="{RIBBON_COLOR}"  stop-opacity="1"/>'
        f"</linearGradient>"
        f"</defs>"
    )
    elements.append(
        f'<rect x="{grad_x:.1f}" y="{grad_y:.1f}" '
        f'width="{IDENTITY_BAR_WIDTH}" height="{IDENTITY_BAR_HEIGHT}" '
        f'fill="url(#{grad_id})" stroke="black" stroke-width="0.8"/>'
    )

    # Labels: 0, Identity (%), 100
    text_y = grad_y + IDENTITY_BAR_HEIGHT + 13
    elements.append(
        f'<text x="{grad_x:.1f}" y="{text_y:.1f}" '
        f'font-family="sans-serif" font-size="11" text-anchor="start">0</text>'
    )
    elements.append(
        f'<text x="{grad_x + IDENTITY_BAR_WIDTH / 2:.1f}" y="{text_y:.1f}" '
        f'font-family="sans-serif" font-size="11" text-anchor="middle">Identity (%)</text>'
    )
    elements.append(
        f'<text x="{grad_x + IDENTITY_BAR_WIDTH:.1f}" y="{text_y:.1f}" '
        f'font-family="sans-serif" font-size="11" text-anchor="end">100</text>'
    )

    return elements


# ==============================================================================
def render_svg(
    globaligner,
    output: Path,
    use_file_order: bool = False,
    scale: float = SCALE,
    show_gene_labels: bool = False,
    identity_threshold: float = 0.3,
    # anchor_label: str | None = None,
    anchor_map: dict[str, str] | None = None,
    no_align: bool = True,
) -> None:
    """
    Renders a clinker Globaligner as a static SVG file.

    Args:
        globaligner: Globaligner object after alignment
        output: Path to write SVG file
        use_file_order: Skip hierarchical clustering, use input order
        scale: bp-to-pixel scaling factor
        show_gene_labels: Annotate each gene arrow with its label
        identity_threshold: Minimum identity to draw a ribbon
    """
    # WARN: test
    # anchor_label = "GB30203_D4993_C5_H3_scaffold_103059_115"

    # -------------------------------------------------------
    # 1. Collect plot data via existing Globaligner method
    # -------------------------------------------------------
    data = globaligner.to_data(use_file_order=use_file_order)
    clusters = data["clusters"]  # list of cluster dicts, in display order
    links = data["links"]  # list of link dicts
    groups = data["groups"]  # list of group dicts (colour assignments)

    # # Build gene_uid -> group colour lookup for fast access
    # uid_to_colour = {}
    # for group in groups:
    #     for gene_uid in group["genes"]:
    #         uid_to_colour[gene_uid] = group.get("colour") or "#cccccc"

    # # uid_to_colour = build_colour_map(groups)
    # if no_align:
    #     uid_to_colour = build_kofam_colour_map(clusters)
    # else:
    #     uid_to_colour = build_colour_map(groups)
    # # NOTE: fix this
    # group_kofam_labels = build_group_kofam_labels(clusters, groups)
    # n_labelled = sum(1 for v in group_kofam_labels.values() if v)
    # n_key_rows = (n_labelled + LEGEND_COLS - 1) // LEGEND_COLS
    # colour_key_height = 20 + n_key_rows * LEGEND_ROW_HEIGHT  # title + rows
    if no_align:
        uid_to_colour, kofam_to_colour = build_kofam_colour_map(clusters)
        # Count unique annotated kofam IDs for height pre-calculation
        unique_kofams = {
            gene.get("kofam_id")
            for cluster in clusters
            for locus in cluster["loci"]
            for gene in locus["genes"]
            if gene.get("kofam_id")
        }
        n_key_rows = (len(unique_kofams) + LEGEND_COLS - 1) // LEGEND_COLS
    else:
        uid_to_colour = build_colour_map(groups)
        group_kofam_labels = build_group_kofam_labels(clusters, groups)
        n_labelled = sum(1 for v in group_kofam_labels.values() if v)
        n_key_rows = (n_labelled + LEGEND_COLS - 1) // LEGEND_COLS

    colour_key_height = (
        20 + n_key_rows * LEGEND_ROW_HEIGHT
    )  # always set before svg_height

    # collect set of target genes
    target_uids: set[str] = set()
    if anchor_map:
        for cluster in clusters:
            anchor_label = anchor_map.get(cluster["name"])
            if anchor_label is None:
                continue
            for locus in cluster["loci"]:
                for gene in locus["genes"]:
                    if gene.get("label") == anchor_label:
                        target_uids.add(gene["uid"])

    # -------------------------------------------------------
    # 1b. Anchor normalisation (optional)
    # -------------------------------------------------------
    # if anchor_label:
    #     cluster_to_anchor = find_anchor_gene_per_cluster(anchor_label, clusters, groups)
    #     # Rewrite loci in each cluster with normalised coordinates
    #     for cluster in clusters:
    #         anchor_uid = cluster_to_anchor.get(cluster["uid"])
    #         if anchor_uid is None:
    #             continue  # no ortholog in this cluster, leave as-is
    #         cluster["loci"] = [
    #             normalise_locus(locus, anchor_uid)[0] for locus in cluster["loci"]
    #         ]
    if anchor_map:
        cluster_to_anchor = find_anchor_gene_per_cluster(anchor_map, clusters)
        for cluster in clusters:
            anchor_uid = cluster_to_anchor.get(cluster["uid"])
            if anchor_uid is None:
                continue
            cluster["loci"] = [
                normalise_locus(locus, anchor_uid)[0] for locus in cluster["loci"]
            ]
        # Recompute max_locus_span after normalisation
        max_locus_span = 0
        for cluster in clusters:
            for locus in cluster["loci"]:
                span = locus["end"] - locus["start"]
                max_locus_span = max(max_locus_span, span)
        svg_width = (
            SVG_PADDING + LABEL_COLUMN_WIDTH + max_locus_span * scale + SVG_PADDING
        )
        # Find leftmost gene coordinate across all clusters
        min_gene_start = 0.0
        for cluster in clusters:
            for locus in cluster["loci"]:
                for gene in locus["genes"]:
                    min_gene_start = min(min_gene_start, gene["start"])

        # # Shift all genes right so nothing is clipped
        # if min_gene_start < 0:
        #     shift = -min_gene_start
        #     for cluster in clusters:
        #         for locus in cluster["loci"]:
        #             for gene in locus["genes"]:
        #                 gene["start"] += shift
        #                 gene["end"] += shift
        #             locus["start"] += shift
        #             locus["end"] += shift

        # After the global shift block in section 1b, add:
        if min_gene_start < 0:
            shift = -min_gene_start
            for cluster in clusters:
                for locus in cluster["loci"]:
                    for gene in locus["genes"]:
                        gene["start"] += shift
                        gene["end"] += shift
                    locus["start"] += shift
                    locus["end"] += shift

        # Force locus start to 0 so gene_to_px origin is always track_x0
        for cluster in clusters:
            for locus in cluster["loci"]:
                locus_span = locus["end"] - locus["start"]
                locus["start"] = 0
                locus["end"] = locus_span

    # -------------------------------------------------------
    # 2. Compute layout geometry
    # -------------------------------------------------------
    n_clusters = len(clusters)

    # Determine total genomic span across all loci (for SVG width)
    max_locus_span = 0
    for cluster in clusters:
        for locus in cluster["loci"]:
            span = locus["end"] - locus["start"]
            max_locus_span = max(max_locus_span, span)

    # svg_width = SVG_PADDING * 2 + max_locus_span * scale
    svg_width = (
        SVG_PADDING * 2 + LABEL_COLUMN_WIDTH + max_locus_span * scale + SVG_PADDING
    )
    # svg_height = SVG_PADDING * 2 + n_clusters * TRACK_SPACING
    # svg_height = SVG_PADDING * 2 + n_clusters * TRACK_SPACING + LEGEND_BOTTOM_MARGIN
    svg_height = (
        SVG_PADDING * 2
        + n_clusters * TRACK_SPACING
        + LEGEND_HEIGHT
        + 20
        + colour_key_height
        + SVG_PADDING
        # + LEGEND_BOTTOM_MARGIN
        # + colour_key_height
    )

    key_x = SVG_PADDING + LABEL_COLUMN_WIDTH
    # key_y = svg_height - SVG_PADDING - colour_key_height + 10

    last_track_bottom = SVG_PADDING + n_clusters * TRACK_SPACING
    # legend_y = svg_height - LEGEND_BOTTOM_MARGIN - colour_key_height + 10
    legend_y = (
        last_track_bottom + TRACK_TO_SCALEBAR_GAP
    )  # gap b/w last track and scale/identity bar
    key_y = (
        legend_y + LEGEND_HEIGHT + SCALEBAR_TO_COLOURKEY_GAP
    )  # gap b/w scale bar and colour key

    # legend_elements = build_legend(svg_width, svg_height, scale)
    legend_elements = build_legend(svg_width=svg_width, legend_y=legend_y, scale=scale)

    if no_align:
        colour_key_elements, colour_key_height = build_kofam_colour_key(
            clusters,
            uid_to_colour,
            kofam_to_colour,
            key_x,
            key_y,
        )
        # n_key_rows = (len(colour_key_elements) + LEGEND_COLS - 1) // LEGEND_COLS
    else:
        colour_key_elements, _ = build_colour_key(
            groups,
            group_kofam_labels,
            uid_to_colour,
            key_x,
            key_y,
        )
    # colour_key_elements, _ = build_colour_key(
    #     groups, group_kofam_labels, uid_to_colour, key_x, key_y
    # )

    elements = []  # accumulate SVG element strings

    # -------------------------------------------------------
    # 3. Build per-cluster geometry lookup
    #    gene_uid -> (x_centre, y_bottom, half_width)
    #    needed when drawing ribbons
    # -------------------------------------------------------
    gene_geom = {}  # uid -> (cx, y_bottom, half_w)

    track_y_top = {}  # cluster_uid -> y of top of track

    for cluster_idx, cluster in enumerate(clusters):
        track_y = SVG_PADDING + cluster_idx * TRACK_SPACING
        track_y_top[cluster["uid"]] = track_y

        # CLUSTER NAME LABEL
        # make nicer looking names for each cluster, extracted from filename
        # FIX: extract info from genbank file itself (need to update genbank contruction script)
        # This is using gbk_renamer.sh method (not ideal)
        # cluster_name = cluster["name"]
        #
        # hyd_name = cluster_name.split(".")[0].split("___")[0]
        # hyd_name = " ".join(hyd_name.split("_"))
        #
        # species_name = cluster_name.split(".")[0].split("___")[1]
        # species_name = " ".join(species_name.split("_"))
        #
        # genome_and_gene = cluster_name.split(".")[0].split("___")[2:]
        # genome_and_gene = "___".join(genome_and_gene)
        #
        # new_cluster_name = f"{hyd_name}; {species_name}; {genome_and_gene}"
        # species_and_genome = f"{species_name} ({genome_and_gene})"

        # FIX: NEW way
        # genbank files are renamed so that the name is:
        # GENOMEID___GENENUM.HydGroup.Species_name.suffix.gbk
        cluster_name = cluster["name"]

        hyd_name = cluster_name.split(".")[1]
        hyd_name = " ".join(hyd_name.split("_"))

        species_name = cluster_name.split(".")[2]
        species_name = " ".join(species_name.split("_"))

        genome_and_gene = cluster_name.split(".")[0]

        new_cluster_name = f"{hyd_name}; {species_name}; {genome_and_gene}"
        species_and_genome = f"{species_name}; {genome_and_gene}"

        # # Cluster name label — right-aligned into the label column
        # elements.append(
        #     f'<text x="{SVG_PADDING + LABEL_COLUMN_WIDTH - 8}" y="{track_y + TRACK_HEIGHT / 2 + 4:.1f}" '
        #     f'font-family="sans-serif" font-size="14" '
        #     f'text-anchor="end" '
        #     # f'dominant-baseline="middle">{cluster["name"]}</text>'
        #     f'dominant-baseline="middle">{new_cluster_name}</text>'
        # )

        # Line 1: hyd_name (normal) + species_name (italic) on same line
        # Line 2: genome_and_gene underneath, smaller font
        line1_y = track_y + TRACK_HEIGHT / 2 - 6  # slightly above centre
        line2_y = line1_y + 16  # second line below

        label_x = SVG_PADDING + LABEL_COLUMN_WIDTH - 8

        elements.append(
            f'<text x="{label_x}" '
            f'font-family="sans-serif" text-anchor="end">'
            # Line 1: hyd_name + italic species_name
            f'<tspan x="{label_x:.1f}" y="{line2_y:1f}" font-size="14" font-style="italic" fill="black">'
            # f'<tspan x="{label_x}" dy="14" font-size="11" fill="#555555">'
            f"{species_and_genome}"
            f"</tspan>"
            # f'<tspan y="{label_x:.1f}" font-size="14" font-style="italic" fill="#555555">'
            # f"{species_and_genome}"
            # f"</tspan>"
            # Line 2: genome/gene identifier, smaller and slightly dimmed
            # f'<tspan x="{label_x}" y="{line1_y:.1f}" font-size="14" font-weight="bold">'
            f'<tspan x="{label_x}" y="{line1_y:.1f}" font-size="14">'
            f"{hyd_name}"
            f"</tspan>"
            f"</text>"
        )

        # Track starts after the label column
        # track_x0 = SVG_PADDING + LABEL_COLUMN_WIDTH

        # for locus in cluster["loci"]:
        #     locus_start = locus["start"]
        #     # track_x0 = SVG_PADDING + 80  # leave room for cluster name label
        #     track_x0 = SVG_PADDING + LABEL_COLUMN_WIDTH
        #
        #     # Draw backbone line for locus
        #     locus_px_start = track_x0
        #     locus_px_end = track_x0 + (locus["end"] - locus_start) * scale
        #     mid_y = track_y + TRACK_HEIGHT / 2
        #     elements.append(
        #         f'<line x1="{locus_px_start:.1f}" y1="{mid_y:.1f}" '
        #         f'x2="{locus_px_end:.1f}" y2="{mid_y:.1f}" '
        #         f'stroke="#aaaaaa" stroke-width="1.5"/>'
        #     )
        #
        for locus in cluster["loci"]:
            locus_start = locus["start"]
            track_x0 = SVG_PADDING + LABEL_COLUMN_WIDTH

            # # Backbone starts and ends where the actual locus content is
            # locus_px_start = (
            #     track_x0 + (locus["start"] - locus_start) * scale
            # )  # = track_x0
            # locus_px_end = track_x0 + (locus["end"] - locus_start) * scale
            # mid_y = track_y + TRACK_HEIGHT / 2

            mid_y = track_y + TRACK_HEIGHT / 2

            # Compute backbone extent from actual gene coordinates
            if locus["genes"]:
                gene_starts = [g["start"] for g in locus["genes"]]
                gene_ends = [g["end"] for g in locus["genes"]]
                locus_px_start = track_x0 + (min(gene_starts) - locus_start) * scale
                locus_px_end = track_x0 + (max(gene_ends) - locus_start) * scale
            else:
                locus_px_start = track_x0 + (locus["start"] - locus_start) * scale
                locus_px_end = track_x0 + (locus["end"] - locus_start) * scale

            elements.append(
                f'<line x1="{locus_px_start:.1f}" y1="{mid_y:.1f}" '
                f'x2="{locus_px_end:.1f}" y2="{mid_y:.1f}" '
                f'stroke="#aaaaaa" stroke-width="1.5"/>'
            )

            for gene in locus["genes"]:
                x1, x2 = gene_to_px(
                    gene_start=gene["start"],
                    gene_end=gene["end"],
                    locus_start=locus_start,
                    track_x0=track_x0,
                )
                colour = uid_to_colour.get(gene["uid"], "#dddddd")
                strand = gene.get("strand", 1)

                pts = arrow_polygon(x1, x2, track_y, strand)
                if not pts:
                    continue

                # NOTE: this controls the look of gene arrows
                # NOTE: make target genes have a slightly thicker line border
                arrow_stroke_width = (
                    ARROW_STROKE_WIDTH * 2
                    if gene["uid"] in target_uids
                    else ARROW_STROKE_WIDTH
                )
                elements.append(
                    f'<polygon points="{pts}" '
                    f'fill="{colour}" stroke="{ARROW_STROKE_COLOUR}" stroke-width="{arrow_stroke_width}" '
                    f'opacity="{ARROW_BODY_OPACITY}"/>'
                )

                # Store geometry for ribbon drawing
                cx = (x1 + x2) / 2
                half_w = (x2 - x1) / 2
                # gene_geom[gene["uid"]] = (cx, track_y + TRACK_HEIGHT, half_w)
                gene_geom[gene["uid"]] = (cx, track_y + TRACK_HEIGHT / 2, half_w)

                if show_gene_labels and gene.get("label"):
                    lx = (x1 + x2) / 2
                    ly = track_y - 3
                    elements.append(
                        f'<text x="{lx:.1f}" y="{ly:.1f}" '
                        f'font-family="sans-serif" font-size="12" '
                        f'text-anchor="middle">{gene["label"]}</text>'
                    )

    # -------------------------------------------------------
    # 4. Draw ribbons between adjacent clusters
    # -------------------------------------------------------
    # Ribbons go between vertically adjacent clusters in display order.
    # A link is only drawn if both genes are in adjacent displayed clusters.

    # Build cluster_uid -> display_index for adjacency check
    cluster_display_idx = {c["uid"]: i for i, c in enumerate(clusters)}

    # Build gene_uid -> cluster_uid lookup
    gene_to_cluster = {}
    for cluster in clusters:
        for locus in cluster["loci"]:
            for gene in locus["genes"]:
                gene_to_cluster[gene["uid"]] = cluster["uid"]

    ribbon_elements = []  # drawn before arrows so arrows sit on top

    for link in links:
        q_uid = link["query"]["uid"]
        t_uid = link["target"]["uid"]
        identity = link.get("identity", 0)

        if identity < identity_threshold:
            continue

        if q_uid not in gene_geom or t_uid not in gene_geom:
            continue

        q_cluster = gene_to_cluster.get(q_uid)
        t_cluster = gene_to_cluster.get(t_uid)

        if q_cluster is None or t_cluster is None:
            continue

        q_idx = cluster_display_idx.get(q_cluster, -1)
        t_idx = cluster_display_idx.get(t_cluster, -1)

        # Only draw ribbons between adjacent tracks
        if abs(q_idx - t_idx) != 1:
            continue

        # Ensure top gene is always "query" for ribbon direction
        if q_idx > t_idx:
            q_uid, t_uid = t_uid, q_uid

        qx, qy, qw = gene_geom[q_uid]  # bottom of upper gene
        tx, ty, tw = gene_geom[t_uid]  # top of lower gene
        # ty_top = ty - TRACK_HEIGHT  # we want top of lower arrow, not bottom

        path = ribbon_path(
            x1=qx,
            y1=qy,
            x2=tx,
            y2=ty,
            w1=qw,
            w2=tw,
        )
        opacity = identity_to_opacity(identity)

        # Ribbon colour from query gene group
        # colour = uid_to_colour.get(q_uid, "#aaaaaa")

        ribbon_elements.append(
            f'<path d="{path}" fill="{RIBBON_COLOR}" opacity="{opacity:.2f}" stroke="none"/>'
        )

    # legend
    # legend_elements = build_legend(svg_width, svg_height, scale)

    # -------------------------------------------------------
    # 5. Assemble final SVG
    # -------------------------------------------------------
    svg_parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="{svg_width:.0f}" height="{svg_height:.0f}" '
        f'viewBox="0 0 {svg_width:.0f} {svg_height:.0f}">',
        '<rect width="100%" height="100%" fill="white"/>',
        "<!-- ribbons -->",
        *ribbon_elements,
        "<!-- gene arrows -->",
        *elements,
        "<!-- legend -->",
        *legend_elements,
        "<!-- colour key -->",
        *colour_key_elements,
        "</svg>",
    ]

    Path(output).write_text("\n".join(svg_parts))
