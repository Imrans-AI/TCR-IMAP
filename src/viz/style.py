"""Shared manuscript figure style: one place for the JCIM formatting rules, the
bootstrap-floor-aware p-value formatter, and the save convention, so every p82v1 figure
script applies them identically instead of re-declaring them per script.

Style per the project notes: 7.2 in wide, 600 dpi PNG + vector PDF, Arial 8 pt body / bold 10 pt
panel letters. No significance stars anywhere -- point/bar fill (or explicit p-text) is
the only significance encoding used in this project's figures.
"""
import matplotlib as mpl
mpl.use("Agg")

PAGE_WIDTH_IN = 7.2  # JCIM two-column figure width

MODEL_LABEL = "TCR-IMAP"          # every displayed model label in every figure
SIM_LABEL = "similarity ratio"   # indel-based (Levenshtein.ratio normalises the indel
                                 # distance); never "identity" / "sequence identity"

# =====================================================================================
# PAPER-WIDE COLOUR SCHEME -- one meaning per colour, identical in every figure.
# Base: Okabe-Ito (colour-blind safe). Figure 1 is an approved asset drawn with its own
# matplotlib source (figures/source/fig1_arial_20261102/) and a PyMOL panel; it uses Arial, as
# set below, and its own restrained palette rather than these hex codes.
# Rules:
#  * a method, molecule or region keeps its colour in every figure it appears in
#  * method hues are never reused for molecules, and vice versa
#  * significance is encoded by FILL only: filled = CI excludes zero, hollow = includes
#    zero. Hollow markers must never encode anything else (e.g. loop class).
#  * loop class is encoded by LINE STYLE: germline CDR1/2 solid, CDR3 dashed
#  * TCR chain is encoded by colour: MOLECULE_COLORS alpha / beta
#  * greys are reserved for baselines and MHC structure only
#  * protocols are a violet strictness ramp, light -> dark (random ... tcr70)
#  * TCR alpha (sky blue) and TCR-IMAP (blue) should not share a panel; re-check if they must
# =====================================================================================
from matplotlib.colors import LinearSegmentedColormap as _LSC

METHOD_COLORS = {
    MODEL_LABEL:            "#0072B2",  # blue
    "Struct-TCS":           "#009E73",  # bluish green
    "Seq-TCS":              "#CC79A7",  # reddish purple
    "TEIM-Res":             "#E69F00",  # orange
    "Position prior":       "#555555",  # dark grey: within-fold frequency baselines
    "Pair-frequency prior": "#555555",
    "NN transfer":          "#AAAAAA",  # light grey -- NN transfer ONLY; reference lines use REFERENCE_LINE_STYLE
    "AngleNet":             "#333333",  # near-black: the separate docking-angle regressor (Figure 3 C/D), never TCR-IMAP blue
}
# Reference lines that are not methods (e.g. a base rate): black, dotted. Never a method colour.
REFERENCE_LINE_STYLE = dict(color="k", ls=(0, (1, 1.3)), lw=1.3)
MOLECULE_COLORS = {
    "TCR alpha": "#56B4E9",  # sky blue
    "TCR beta":  "#D55E00",  # vermillion
    "MHC":       "#BFBFBF",  # grey
    "peptide":   "#FFD92F",  # yellow (as in v4 structural panels)
}
REGION_COLORS = {            # MHC stays grey; darker = more contact-bearing
    "alpha1 helix":         "#4D4D4D",
    "alpha2 helix":         "#969696",
    "other G-domain":       "#CCCCCC",
    "beta-sheet floor":     "#EBEBEB",
    "unresolved/uncertain": "#F5F5F5",
}
LOOP_CHAIN = {"cdr1a": "TCR alpha", "cdr2a": "TCR alpha", "cdr3a": "TCR alpha",
              "cdr1b": "TCR beta",  "cdr2b": "TCR beta",  "cdr3b": "TCR beta"}
GERMLINE_LOOPS = {"cdr1a", "cdr2a", "cdr1b", "cdr2b"}
LOOP_LINESTYLE = {"germline": "-", "cdr3": (0, (3, 2))}
PROTOCOL_COLORS = {          # light -> dark by strictness; tcr70 (primary) darkest
    "random": "#BCBDDC", "double": "#9E9AC8", "allele": "#807DBA",
    "temporal": "#6A51A3", "tcr70": "#3F007D",
}
# MD distance strata (crystal-frame minimum heavy-atom distance; Figure 7): an ordinal brown ramp, dark = near.
# Brown is used by no method or molecule; the ramp's meaning is distance from the interface and nothing else.
STRATUM_COLORS = {"4.5-6 A": "#543005", "6-8 A": "#BF812D", "> 8 A": "#DFC27D"}
# MD replicates (SI characterisation figure): an olive ramp, used for nothing else. Order carries no meaning beyond
# the replicate number (seeds 271828 / 314159 / 161803).
REPLICATE_COLORS = {1: "#3D4F1F", 2: "#7A8F3A", 3: "#B5C46E"}
CMAP_SEQ = _LSC.from_list("seq", ["#F7F7F7", "#C6DBEF", "#4292C6", "#08306B"])
CMAP_DIV = _LSC.from_list("div", ["#8C510A", "#F7F7F7", "#35978F"])


def loop_color(loop):
    return MOLECULE_COLORS[LOOP_CHAIN[loop]]


def loop_linestyle(loop):
    return LOOP_LINESTYLE["germline" if loop in GERMLINE_LOOPS else "cdr3"]


def sig_marker(color, significant):
    """Marker kwargs: filled if the CI excludes zero, hollow otherwise."""
    return dict(mec=color, mfc=color if significant else "white", mew=1.0)


def pymol_hex(hexcolor):
    """'#56B4E9' -> '0x56B4E9' for cmd.color()."""
    return "0x" + hexcolor.lstrip("#").upper()


# legacy names kept so existing scripts still import; they now point into the scheme
BLUE = METHOD_COLORS[MODEL_LABEL]
RED = METHOD_COLORS["Pair-frequency prior"]
GREY = "#AAAAAA"
COMPARATOR_COLORS = dict(METHOD_COLORS)
COMPARATOR_COLORS["PosFreq"] = METHOD_COLORS["Position prior"]
COMPARATOR_COLORS["Contact-frequency prior"] = METHOD_COLORS["Pair-frequency prior"]

_m = {v.upper() for k, v in METHOD_COLORS.items() if "prior" not in k and k != "NN transfer"}
_o = {v.upper() for v in MOLECULE_COLORS.values()}
assert not (_m & _o), f"method/molecule colour collision: {_m & _o}"
_s = {v.upper() for v in list(STRATUM_COLORS.values()) + list(REPLICATE_COLORS.values())}
assert not (_s & ({v.upper() for v in METHOD_COLORS.values()} | _o | {v.upper() for v in PROTOCOL_COLORS.values()})), \
    "stratum/replicate colour collides with a method, molecule or protocol colour"

PANEL_LABEL_STYLE = dict(fontweight="bold", fontsize=10, va="bottom", ha="left")
# One panel-letter offset for every figure, in inches, from the top-left corner of the panel's CONTENT (the union of its
# axes with their tick labels, axis labels, titles and legends, plus any figure-level text passed with the panel):
# the letter's left edge at the content's left edge, its bottom 0.04 in above the content's top.
PANEL_LETTER_OFFSET_IN = (0.0, 0.04)
PANEL_LETTER_HEIGHT_IN = 0.16          # extent of a 10 pt bold letter (0.14 in measured) plus a small margin


def panel_letters(fig, panels, offset_in=PANEL_LETTER_OFFSET_IN):
    """Place every panel letter at the same offset from its panel's content. `panels`: {letter: artist or list of
    artists} (axes, or figure-level Text belonging to the panel). Call once, after the layout is final and immediately
    before save_fig. Raises if a letter would leave the page or overlap any other visible text, so a figure can never
    be written with a clipped or colliding letter."""
    from matplotlib.transforms import Bbox
    from matplotlib.text import Text
    W, H = fig.get_size_inches()
    eng = fig.get_layout_engine()
    if eng is not None and type(eng).__name__ == "ConstrainedLayoutEngine":
        # constrained layout ignores figure-level text, so reserve the letter's height + offset at the top
        # and room for a letter in every gap between rows (h_pad applies on both sides of each gap)
        l, b, w, h = eng.get()["rect"]
        need_in = offset_in[1] + PANEL_LETTER_HEIGHT_IN
        eng.set(rect=(l, b, w, h - need_in / H), h_pad=max(eng.get()["h_pad"], need_in / 2))
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    dpi = fig.dpi
    placed = []
    for lab, arts in panels.items():
        arts = arts if isinstance(arts, (list, tuple)) else [arts]
        bb = Bbox.union([a.get_tightbbox(r) if hasattr(a, "get_tightbbox") and not isinstance(a, Text)
                         else a.get_window_extent(r) for a in arts])
        x_in, y_in = bb.x0 / dpi + offset_in[0], bb.y1 / dpi + offset_in[1]
        placed.append(fig.text(x_in / W, y_in / H, lab, **PANEL_LABEL_STYLE))
    fig.canvas.draw()
    others = [t for t in fig.findobj(Text) if t.get_visible() and t.get_text().strip() and t not in placed]
    for t in placed:
        e = t.get_window_extent(r)
        if e.x0 < 0 or e.y0 < 0 or e.x1 > W * dpi or e.y1 > H * dpi:
            raise ValueError(f"panel letter {t.get_text()} leaves the page: extent {e.bounds} (page {W * dpi:.0f} x {H * dpi:.0f} px)")
        for o in others:
            if e.overlaps(o.get_window_extent(r)):
                raise ValueError(f"panel letter {t.get_text()} overlaps text {o.get_text()[:40]!r}")
    return placed

BOOTSTRAP_FLOOR = 0.0005   # 1 / n_boot=2000, per src.eval.metrics.cluster_bootstrap_diff


def apply_style():
    """Call once at the top of every figure script, before creating any figure."""
    mpl.rcParams.update({
        # Type 42 (TrueType) required: matplotlib defaults to Type 3, which ACS
        # production rejects and which carries no selectable text layer.
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 8,
        "axes.labelsize": 8,
        "legend.fontsize": 6.5,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })


def pstr(p):
    """Bootstrap-floor-aware p-value text. Never emits a raw floored value (e.g.
    "p = 0.0005"); never emits a significance star -- this project's convention is
    p-text or point/bar fill, not asterisks."""
    return "p < 0.001" if p <= BOOTSTRAP_FLOOR else ("p = %.3f" % p)


MINUS = "−"   # typographic minus; ASCII '-' is a hyphen


def signed(x, nd=3):
    """Explicitly signed number for annotation text with a proper minus:
    signed(0.039) -> '+0.039', signed(-0.0046) -> '−0.005'."""
    return f"{x:+.{nd}f}".replace("-", MINUS)


def minus(s):
    """Replace ASCII '-' immediately before a digit with U+2212 in an already formatted
    annotation string. Use only on annotation text, never on tick labels (matplotlib already
    renders those with the proper minus) or on category labels such as '4-5'."""
    import re
    return re.sub(r"-(?=\d)", MINUS, s)


def save_fig(fig, stem, dpi=600):
    """Write <stem>.png (raster, dpi) and <stem>.pdf (vector) together, and print a
    one-line confirmation so the console log records exactly what was written."""
    png_path = f"{stem}.png"
    pdf_path = f"{stem}.pdf"
    fig.savefig(png_path, dpi=dpi, facecolor="white")
    fig.savefig(pdf_path, facecolor="white")
    print(f"wrote {png_path} / {pdf_path.split('/')[-1]}")
    return png_path, pdf_path


def print_validation_table(rows, columns):
    """Print a panel-to-source numerical validation table.

    rows: list of dicts. columns: list of (key, header, fmt) tuples, fmt a %-style
    format string applied via `fmt % row[key]` (or str() if fmt is None).
    """
    headers = [h for _, h, _ in columns]
    widths = [len(h) for h in headers]
    formatted_rows = []
    for row in rows:
        vals = []
        for key, _, fmt in columns:
            v = row.get(key, "")
            vals.append(v if fmt is None else (fmt % v))
        formatted_rows.append(vals)
    for vals in formatted_rows:
        for i, v in enumerate(vals):
            widths[i] = max(widths[i], len(str(v)))
    def fmt_line(vals):
        return "  ".join(str(v).ljust(w) for v, w in zip(vals, widths))
    print(fmt_line(headers))
    print(fmt_line(["-" * w for w in widths]))
    for vals in formatted_rows:
        print(fmt_line(vals))
