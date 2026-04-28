"""
Warning: all vibe coded... mileage may vary
"""

from pathlib import Path
from typing import Optional
import colorsys

# Layout constants — all in pixels
TRACK_HEIGHT = 18  # height of gene arrow body
ARROW_HEAD = 10  # width of arrowhead
TRACK_SPACING = 120  # vertical distance between cluster tracks
LABEL_OFFSET = 14  # px above track for gene labels
CLUSTER_LABEL_X = 10  # x position of cluster name label
RIBBON_OPACITY_MIN = 0.1
RIBBON_OPACITY_MAX = 0.8
SVG_PADDING = 40  # padding around entire figure
SCALE = 0.05  # bp -> px scaling factor


def identity_to_opacity(identity: float) -> float:
    """Maps alignment identity to ribbon opacity."""
    return RIBBON_OPACITY_MIN + identity * (RIBBON_OPACITY_MAX - RIBBON_OPACITY_MIN)


def gene_colour(gene_uid: str, groups: list) -> str:
    """Looks up colour for a gene UID from the group list."""
    for group in groups:
        if gene_uid in group["genes"]:
            return group["colour"] or "#cccccc"
    return "#dddddd"  # ungrouped genes


def gene_to_px(
    gene_start: int, gene_end: int, locus_start: int, track_x0: int
) -> tuple:
    """Converts gene genomic coordinates to pixel x positions."""
    x1 = track_x0 + (gene_start - locus_start) * SCALE
    x2 = track_x0 + (gene_end - locus_start) * SCALE
    return x1, x2


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


def ribbon_path(
    x1: float, y1: float, x2: float, y2: float, w1: float, w2: float
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

    return (
        f"M {tl[0]:.1f},{tl[1]:.1f} "
        f"C {tl[0]:.1f},{cy:.1f} {bl[0]:.1f},{cy:.1f} {bl[0]:.1f},{bl[1]:.1f} "
        f"L {br[0]:.1f},{br[1]:.1f} "
        f"C {br[0]:.1f},{cy:.1f} {tr[0]:.1f},{cy:.1f} {tr[0]:.1f},{tr[1]:.1f} "
        f"Z"
    )


def render_svg(
    globaligner,
    output: Path,
    use_file_order: bool = False,
    scale: float = SCALE,
    show_gene_labels: bool = False,
    identity_threshold: float = 0.3,
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

    # -------------------------------------------------------
    # 1. Collect plot data via existing Globaligner method
    # -------------------------------------------------------
    data = globaligner.to_data(use_file_order=use_file_order)
    clusters = data["clusters"]  # list of cluster dicts, in display order
    links = data["links"]  # list of link dicts
    groups = data["groups"]  # list of group dicts (colour assignments)

    # Build gene_uid -> group colour lookup for fast access
    uid_to_colour = {}
    for group in groups:
        for gene_uid in group["genes"]:
            uid_to_colour[gene_uid] = group.get("colour") or "#cccccc"

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

    svg_width = SVG_PADDING * 2 + max_locus_span * scale
    svg_height = SVG_PADDING * 2 + n_clusters * TRACK_SPACING

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

        # Cluster name label
        elements.append(
            f'<text x="{CLUSTER_LABEL_X}" y="{track_y + TRACK_HEIGHT / 2 + 4:.1f}" '
            f'font-family="sans-serif" font-size="12" '
            f'dominant-baseline="middle">{cluster["name"]}</text>'
        )

        for locus in cluster["loci"]:
            locus_start = locus["start"]
            track_x0 = SVG_PADDING + 80  # leave room for cluster name label

            # Draw backbone line for locus
            locus_px_start = track_x0
            locus_px_end = track_x0 + (locus["end"] - locus_start) * scale
            mid_y = track_y + TRACK_HEIGHT / 2
            elements.append(
                f'<line x1="{locus_px_start:.1f}" y1="{mid_y:.1f}" '
                f'x2="{locus_px_end:.1f}" y2="{mid_y:.1f}" '
                f'stroke="#aaaaaa" stroke-width="1.5"/>'
            )

            for gene in locus["genes"]:
                x1, x2 = gene_to_px(gene["start"], gene["end"], locus_start, track_x0)
                colour = uid_to_colour.get(gene["uid"], "#dddddd")
                strand = gene.get("strand", 1)

                pts = arrow_polygon(x1, x2, track_y, strand)
                if not pts:
                    continue

                elements.append(
                    f'<polygon points="{pts}" '
                    f'fill="{colour}" stroke="white" stroke-width="0.8" '
                    f'opacity="0.9"/>'
                )

                # Store geometry for ribbon drawing
                cx = (x1 + x2) / 2
                half_w = (x2 - x1) / 2
                gene_geom[gene["uid"]] = (cx, track_y + TRACK_HEIGHT, half_w)

                if show_gene_labels and gene.get("label"):
                    lx = (x1 + x2) / 2
                    ly = track_y - 3
                    elements.append(
                        f'<text x="{lx:.1f}" y="{ly:.1f}" '
                        f'font-family="sans-serif" font-size="8" '
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
        ty_top = ty - TRACK_HEIGHT  # we want top of lower arrow, not bottom

        path = ribbon_path(qx, qy, tx, ty_top, qw, tw)
        opacity = identity_to_opacity(identity)

        # Ribbon colour from query gene group
        colour = uid_to_colour.get(q_uid, "#aaaaaa")

        ribbon_elements.append(
            f'<path d="{path}" fill="{colour}" opacity="{opacity:.2f}" stroke="none"/>'
        )

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
        "</svg>",
    ]

    Path(output).write_text("\n".join(svg_parts))
