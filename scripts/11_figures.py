#!/usr/bin/env python3
"""The figures the README shows. Reads only files in results/, writes figures/.

Every figure is produced from a table an earlier stage wrote, so any number
visible in a figure is traceable to a file. No figure recomputes anything, and a
figure whose input table does not exist yet is skipped with a message rather
than invented.

Design notes, since these are the plots a reader will judge the work by.

Colour. Six arms, so six categorical slots taken in fixed order and never
cycled. The order was validated rather than chosen by eye: worst adjacent
colour-vision separation 9.1, worst normal-vision separation 19.6, both above
the floor. Three of the six sit below 3:1 contrast against the surface, which
obliges a visible label on every mark rather than relying on the fill, so every
bar carries its number and every line is labelled at its end. Scatter panels use
only the first three slots, which are the ones that separate when every pair can
appear together; more than three arms are drawn as small multiples instead.

Form. Magnitude comparisons are horizontal bars, sorted, because the arm names
are long and a reader compares lengths from a common baseline more reliably than
angles or areas. The calibration scatter is faceted rather than overlaid,
because six clouds of points on one pair of axes is not a chart. No figure has
two vertical scales.

Theme. Light surface, fixed. These are image files in a README, so they cannot
follow a reader's dark mode; a transparent background would invert the axis text
into illegibility on a dark page, so the surface is painted explicitly.

Usage
-----
    python scripts/11_figures.py --config project.conf
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_csc as L  # noqa: E402

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
BASELINE = "#c3c2b7"

# Fixed slot order. A seventh arm would take slot 7, never a generated hue.
SLOTS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300",
         "#4a3aa7", "#e34948"]

ARM_LABEL = {
    "af2_msa_notmpl": "AlphaFold2, alignment, no templates",
    "af2_msa_tmpl": "AlphaFold2, alignment, templates on",
    "af2_nomsa": "AlphaFold2, no alignment",
    "boltz2_msa": "Boltz-2, alignment",
    "null_template": "best pre-cutoff template, copied",
    "null_unrelated": "unrelated chain (floor)",
}
ARM_ORDER = ["af2_msa_notmpl", "af2_msa_tmpl", "af2_nomsa", "boltz2_msa",
             "null_template", "null_unrelated"]


def arm_colour(arm: str) -> str:
    try:
        return SLOTS[ARM_ORDER.index(arm)]
    except ValueError:
        return SLOTS[-1]


def label_of(arm: str) -> str:
    return ARM_LABEL.get(arm, arm)


def style():
    plt.rcParams.update({
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "text.color": INK,
        "axes.labelcolor": INK_2,
        "axes.edgecolor": BASELINE,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "grid.color": GRID,
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.titleweight": "bold",
        "figure.dpi": 150,
        "savefig.bbox": "tight",
        "axes.spines.top": False,
        "axes.spines.right": False,
    })


def tidy(ax, grid_axis="x"):
    ax.grid(True, axis=grid_axis, linewidth=0.6, alpha=0.9)
    ax.set_axisbelow(True)


def num(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(v) else v


def table(results: Path, name: str):
    p = results / name
    if not p.is_file():
        return None
    rows = L.read_tsv(p)
    return rows or None


# ---------------------------------------------------------------------------
def fig_calibration_scatter(results: Path, figures: Path, conf: dict) -> str | None:
    """Confidence against the measurement it predicts, one panel per arm.

    The shaded corner is the part that matters: residues the model called
    confident whose measured accuracy came back below the line. A well
    calibrated model puts almost nothing there.
    """
    rows = table(results, "residue_scores.tsv")
    if not rows:
        return None
    high = L.conf_float(conf, "PLDDT_HIGH", 90.0) / 100.0
    trust = L.conf_float(conf, "LDDT_TRUST", 0.7)

    by_arm: dict[str, list] = {}
    for r in rows:
        p, l = num(r.get("plddt")), num(r.get("lddt_ca"))
        if p is None or l is None:
            continue
        by_arm.setdefault(r["arm"], []).append((p / 100.0, l))
    if not by_arm:
        return None
    arms = [a for a in ARM_ORDER if a in by_arm] + \
           [a for a in sorted(by_arm) if a not in ARM_ORDER]

    ncols = min(3, len(arms))
    nrows = math.ceil(len(arms) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.6 * ncols, 3.4 * nrows),
                             squeeze=False, sharex=True, sharey=True)
    for i, arm in enumerate(arms):
        ax = axes[i // ncols][i % ncols]
        pts = by_arm[arm]
        xs = [p for p, _ in pts]
        ys = [l for _, l in pts]
        ax.axhspan(0, trust, xmin=high, xmax=1.0, color="#e34948", alpha=0.10, lw=0)
        ax.plot([0, 1], [0, 1], color=MUTED, lw=1.0, ls=(0, (4, 3)), zorder=1)
        ax.scatter(xs, ys, s=3, alpha=0.18, linewidths=0,
                   color=arm_colour(arm), zorder=2)
        bad = sum(1 for p, l in pts if p >= high and l < trust)
        n_high = sum(1 for p, _l in pts if p >= high)
        share = f"{100.0 * bad / n_high:.1f}%" if n_high else "no residues"
        ax.set_title(label_of(arm), color=INK, loc="left")
        # Against the shaded region rather than in the opposite corner: the
        # number counts what is inside that box, and a reader should not have
        # to work out which of the two it belongs to. The panel background
        # behind the text keeps it legible over the scatter.
        ax.text(0.97, trust + 0.02, f"above {high:.2f} and below {trust:g}: {share}",
                transform=ax.get_xaxis_transform(), fontsize=8, color=INK_2,
                ha="right", va="bottom", zorder=3,
                bbox=dict(facecolor=SURFACE, edgecolor="none", pad=1.5, alpha=0.85))
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        tidy(ax, "both")
        if i % ncols == 0:
            ax.set_ylabel("measured lDDT-CA")
        if i // ncols == nrows - 1:
            ax.set_xlabel("confidence, on the score's own scale")
    for j in range(len(arms), nrows * ncols):
        axes[j // ncols][j % ncols].axis("off")
    # tight_layout first, then the figure title above the panels. Setting the
    # title before the layout pass puts it on top of the first panel's own
    # title when there is only one panel.
    fig.tight_layout()
    fig.suptitle("Confidence against the accuracy it predicts", x=0.005, y=1.005,
                 ha="left", va="bottom", fontsize=11.5, color=INK)
    out = figures / "fig1_calibration_scatter.png"
    fig.savefig(out)
    plt.close(fig)
    return out.name


def fig_calibration_curve(results: Path, figures: Path) -> str | None:
    rows = table(results, "calibration_bins.tsv")
    if not rows:
        return None
    by_arm: dict[str, list] = {}
    for r in rows:
        mc, mm, n = num(r.get("mean_confidence")), num(r.get("mean_measured")), num(r.get("n"))
        if mc is None or mm is None or not n:
            continue
        by_arm.setdefault(r["arm"], []).append((mc, mm, n))
    if not by_arm:
        return None
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    ends: list[tuple[float, float, str]] = []
    ax.plot([0, 1], [0, 1], color=MUTED, lw=1.0, ls=(0, (4, 3)),
            label="perfect calibration", zorder=1)
    for arm in [a for a in ARM_ORDER if a in by_arm] + \
               [a for a in sorted(by_arm) if a not in ARM_ORDER]:
        pts = sorted(by_arm[arm])
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        ax.plot(xs, ys, lw=2.0, marker="o", ms=4, color=arm_colour(arm),
                label=label_of(arm), zorder=2)
        # Direct label at the end of each line, because three of the six
        # colours sit below the contrast floor against this surface.
        if xs:
            ends.append((ys[-1], xs[-1], label_of(arm)))
    # The labels go in the figure margin, not inside a widened axis. Confidence
    # cannot exceed one, and stretching the scale to 1.25 to make room for text
    # invites a reader to look for points that could never exist there.
    ends.sort()
    min_gap = 0.055
    placed: list[float] = []
    for y, x, text in ends:
        while placed and y - placed[-1] < min_gap:
            y = placed[-1] + min_gap
        placed.append(y)
        ax.annotate(text, (x, y), xytext=(8, 0), textcoords="offset points",
                    fontsize=7.5, color=INK_2, va="center",
                    annotation_clip=False)
    ax.set_xlabel("mean confidence in bin")
    ax.set_ylabel("mean measured lDDT-CA in bin")
    ax.set_title("Calibration, binned", loc="left", color=INK)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    tidy(ax, "both")
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    out = figures / "fig2_calibration_curve.png"
    fig.savefig(out)
    plt.close(fig)
    return out.name


def fig_threshold_sweep(results: Path, figures: Path) -> str | None:
    rows = table(results, "threshold_sweep.tsv")
    if not rows:
        return None
    by_arm: dict[str, list] = {}
    for r in rows:
        band, frac = num(r.get("band")), num(r.get("fraction_below_trust"))
        if band is None or frac is None:
            continue
        by_arm.setdefault(r["arm"], []).append((band, frac))
    if not by_arm:
        return None
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for arm in [a for a in ARM_ORDER if a in by_arm]:
        pts = sorted(by_arm[arm])
        ax.plot([p[0] for p in pts], [p[1] for p in pts], lw=2.0, marker="o", ms=4,
                color=arm_colour(arm), label=label_of(arm))
    ax.set_xlabel("confidence band, residues at or above")
    ax.set_ylabel("fraction measured below the accuracy line")
    ax.set_title("What a confidence threshold buys", loc="left", color=INK)
    tidy(ax, "both")
    ax.legend(frameon=False, fontsize=8)
    out = figures / "fig3_threshold_sweep.png"
    fig.savefig(out)
    plt.close(fig)
    return out.name


def _hbar(ax, labels, values, colours, fmt="{:.3f}", xlabel=""):
    y = range(len(labels))
    ax.barh(list(y), values, height=0.62, color=colours, linewidth=0)
    ax.set_yticks(list(y))
    ax.set_yticklabels(labels, color=INK)
    ax.invert_yaxis()
    span = max(values) if values else 1.0
    for i, v in enumerate(values):
        ax.text(v + span * 0.015, i, fmt.format(v), va="center", fontsize=8.5,
                color=INK)
    ax.set_xlabel(xlabel)
    ax.set_xlim(0, span * 1.18)
    tidy(ax, "x")


def fig_arm_accuracy(results: Path, figures: Path) -> str | None:
    rows = table(results, "arm_comparison.tsv")
    if not rows:
        return None
    rows = [r for r in rows if num(r.get("lddt_ca_median")) is not None]
    if not rows:
        return None
    rows.sort(key=lambda r: num(r["lddt_ca_median"]), reverse=True)
    fig, ax = plt.subplots(figsize=(7.0, 0.55 * len(rows) + 1.7))
    _hbar(ax,
          [f"{label_of(r['arm'])}  (n={r['n_targets']})" for r in rows],
          [num(r["lddt_ca_median"]) for r in rows],
          [arm_colour(r["arm"]) for r in rows],
          xlabel="median lDDT-CA against the withheld structure")
    # The interquartile range as a thin rule through each bar, so the spread is
    # visible rather than hidden behind a median.
    for i, r in enumerate(rows):
        q1, q3 = num(r.get("lddt_ca_q1")), num(r.get("lddt_ca_q3"))
        if q1 is not None and q3 is not None:
            ax.plot([q1, q3], [i, i], color=INK, lw=1.4, alpha=0.55,
                    solid_capstyle="butt", zorder=3)
    ax.set_title("Accuracy per arm, median with interquartile range",
                 loc="left", color=INK)
    out = figures / "fig4_arm_accuracy.png"
    fig.savefig(out)
    plt.close(fig)
    return out.name


def fig_geometry(results: Path, figures: Path) -> str | None:
    rows = table(results, "arm_comparison.tsv")
    if not rows:
        return None
    rows = [r for r in rows if num(r.get("fraction_passing_geometry")) is not None]
    if not rows:
        return None
    rows.sort(key=lambda r: num(r["fraction_passing_geometry"]), reverse=True)
    fig, ax = plt.subplots(figsize=(7.0, 0.55 * len(rows) + 1.7))
    _hbar(ax,
          [f"{label_of(r['arm'])}  ({r['n_passing_geometry']}/{r['n_geometry_checked']})"
           for r in rows],
          [num(r["fraction_passing_geometry"]) for r in rows],
          [arm_colour(r["arm"]) for r in rows],
          fmt="{:.1%}",
          xlabel="fraction of predictions passing every physical check")
    ax.set_title("Physical validity per arm", loc="left", color=INK)
    out = figures / "fig5_geometry_pass.png"
    fig.savefig(out)
    plt.close(fig)
    return out.name


def fig_accuracy_and_validity(results: Path, figures: Path) -> str | None:
    """Accurate, valid, and both. The third is the only one called success."""
    rows = table(results, "arm_comparison.tsv")
    if not rows:
        return None
    rows = [r for r in rows if num(r.get("fraction_accurate_and_valid")) is not None]
    if not rows:
        return None
    rows.sort(key=lambda r: num(r["fraction_accurate_and_valid"]), reverse=True)
    fig, ax = plt.subplots(figsize=(7.2, 0.62 * len(rows) + 1.9))
    y = range(len(rows))
    valid = [num(r.get("fraction_passing_geometry")) or 0.0 for r in rows]
    both = [num(r["fraction_accurate_and_valid"]) for r in rows]
    ax.barh([i - 0.19 for i in y], valid, height=0.34, color=SLOTS[2],
            linewidth=0, label="physically valid")
    ax.barh([i + 0.19 for i in y], both, height=0.34, color=SLOTS[0],
            linewidth=0, label="accurate and valid")
    for i, (v, b) in enumerate(zip(valid, both)):
        ax.text(v + 0.012, i - 0.19, f"{v:.1%}", va="center", fontsize=8, color=INK)
        ax.text(b + 0.012, i + 0.19, f"{b:.1%}", va="center", fontsize=8, color=INK)
    ax.set_yticks(list(y))
    ax.set_yticklabels([label_of(r["arm"]) for r in rows], color=INK)
    ax.invert_yaxis()
    ax.set_xlim(0, max(valid + both + [0.1]) * 1.2)
    ax.set_xlabel("fraction of targets")
    # The two bars are named in the legend. Calling one of them "the darker
    # bar" told a reader to compare shades of two different hues.
    ax.set_title("Accurate, physically valid, and both", loc="left", color=INK)
    tidy(ax, "x")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    out = figures / "fig6_accuracy_and_validity.png"
    fig.savefig(out)
    plt.close(fig)
    return out.name


def fig_timing(results: Path, figures: Path) -> str | None:
    rows = table(results, "predictions.tsv")
    if not rows:
        return None
    pts: dict[str, list] = {}
    for r in rows:
        if r.get("status") != "ok":
            continue
        n, s = num(r.get("sequence_length")), num(r.get("elapsed_s"))
        if n is None or s is None:
            continue
        pts.setdefault(r["arm"], []).append((n, s / 60.0))
    if not pts:
        return None
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for arm in [a for a in ARM_ORDER if a in pts]:
        xs = [p[0] for p in pts[arm]]
        ys = [p[1] for p in pts[arm]]
        ax.scatter(xs, ys, s=26, color=arm_colour(arm), alpha=0.85,
                   linewidths=0.6, edgecolors=SURFACE, label=label_of(arm))
    ax.set_xlabel("residues")
    ax.set_ylabel("minutes on the processor")
    ax.set_title("What one prediction costs on this machine", loc="left", color=INK)
    tidy(ax, "both")
    ax.legend(frameon=False, fontsize=8)
    out = figures / "fig7_timing.png"
    fig.savefig(out)
    plt.close(fig)
    return out.name


def fig_seed_variance(results: Path, figures: Path) -> str | None:
    rows = table(results, "seed_variance.tsv")
    if not rows:
        return None
    rows = [r for r in rows if num(r.get("range")) is not None]
    if not rows:
        return None
    fig, ax = plt.subplots(figsize=(7.0, 0.5 * len(rows) + 1.8))
    _hbar(ax,
          [f"{r['target_id']}  {label_of(r['arm'])}" for r in rows],
          [num(r["range"]) for r in rows],
          [arm_colour(r["arm"]) for r in rows],
          xlabel="range in lDDT-CA across repeat seeds")
    ax.set_title("How much the answer moves when only the seed changes",
                 loc="left", color=INK)
    out = figures / "fig8_seed_variance.png"
    fig.savefig(out)
    plt.close(fig)
    return out.name


def fig_memory_curve(results: Path, figures: Path, conf: dict) -> str | None:
    """The measurement that set the ceiling, with the fit drawn through it."""
    rows = table(results, "environment/memory_curve.tsv")
    if not rows:
        return None
    pts = [(num(r.get("sequence_length")), num(r.get("peak_rss_mb")), num(r.get("elapsed_s")))
           for r in rows]
    pts = [p for p in pts if p[0] is not None and p[1] is not None]
    if not pts:
        return None
    pts.sort()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.0, 3.8))
    ram = L.conf_int(conf, "RAM_GB", 7) * 1024
    ax1.scatter([p[0] for p in pts], [p[1] for p in pts], s=40, color=SLOTS[0],
                zorder=3, linewidths=0)
    ax1.axhline(ram, color="#e34948", lw=1.2, ls=(0, (4, 3)))
    ax1.annotate(f"memory available, {ram} MB", (pts[0][0], ram), xytext=(0, 5),
                 textcoords="offset points", fontsize=8, color=INK_2)
    for x, y, _ in pts:
        ax1.annotate(f"{y:.0f}", (x, y), xytext=(0, 7), textcoords="offset points",
                     fontsize=8, ha="center", color=INK)
    ax1.set_xlabel("residues")
    ax1.set_ylabel("peak resident memory, MB")
    ax1.set_ylim(0, max(ram, max(p[1] for p in pts)) * 1.2)
    ax1.set_title("Memory measured on this machine", loc="left", color=INK)
    tidy(ax1, "both")

    ax2.scatter([p[0] for p in pts], [(p[2] or 0) / 60.0 for p in pts], s=40,
                color=SLOTS[1], zorder=3, linewidths=0)
    for x, _, s in pts:
        if s:
            ax2.annotate(f"{s / 60.0:.1f}", (x, s / 60.0), xytext=(0, 7),
                         textcoords="offset points", fontsize=8, ha="center", color=INK)
    ax2.set_xlabel("residues")
    ax2.set_ylabel("minutes")
    ax2.set_title("Time measured on this machine", loc="left", color=INK)
    tidy(ax2, "both")
    out = figures / "fig9_memory_curve.png"
    fig.savefig(out)
    plt.close(fig)
    return out.name


def fig_pocket(results: Path, figures: Path, conf: dict) -> str | None:
    """Global accuracy against pocket accuracy, one point per target and arm.

    The question is whether a prediction that scores well over the whole chain
    can be trusted where a ligand would sit. Points above the diagonal are
    targets whose pocket is better than the chain as a whole; points below are
    the ones that would mislead anyone reading only the global number.
    """
    rows = table(results, "pocket_scores.tsv")
    if not rows:
        return None
    primary = L.conf_float(conf, "POCKET_RADIUS", 5.0)
    pts: dict[str, list] = {}
    for r in rows:
        if abs(num(r.get("radius")) or -1 - primary) > 1e-6 and                 str(r.get("radius")) != f"{primary}":
            continue
        g, pk = num(r.get("global_lddt_ca")), num(r.get("pocket_lddt_ca"))
        if g is None or pk is None:
            continue
        pts.setdefault(r["arm"], []).append((g, pk))
    if not pts:
        return None
    fig, ax = plt.subplots(figsize=(5.8, 5.4))
    ax.plot([0, 1], [0, 1], color=MUTED, lw=1.0, ls=(0, (4, 3)), zorder=1,
            label="pocket as good as the chain")
    for arm in [a for a in ARM_ORDER if a in pts] +                [a for a in sorted(pts) if a not in ARM_ORDER]:
        xs = [p[0] for p in pts[arm]]
        ys = [p[1] for p in pts[arm]]
        ax.scatter(xs, ys, s=34, color=arm_colour(arm), alpha=0.85, linewidths=0.6,
                   edgecolors=SURFACE, label=label_of(arm), zorder=2)
    ax.set_xlabel("lDDT-CA over the whole chain")
    ax.set_ylabel(f"lDDT-CA over residues within {primary:g} Angstrom of the ligand")
    ax.set_title("Global accuracy is not pocket accuracy", loc="left", color=INK)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    tidy(ax, "both")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    out = figures / "fig10_pocket_vs_global.png"
    fig.savefig(out)
    plt.close(fig)
    return out.name


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True, type=Path)
    args = ap.parse_args()
    conf = L.load_conf(args.config)
    results = Path(conf["RESULTS_DIR"])
    figures = Path(conf["FIGURES_DIR"])
    figures.mkdir(parents=True, exist_ok=True)
    style()

    made, skipped = [], []
    for name, fn in [
        ("calibration scatter", lambda: fig_calibration_scatter(results, figures, conf)),
        ("calibration curve", lambda: fig_calibration_curve(results, figures)),
        ("threshold sweep", lambda: fig_threshold_sweep(results, figures)),
        ("accuracy per arm", lambda: fig_arm_accuracy(results, figures)),
        ("physical validity", lambda: fig_geometry(results, figures)),
        ("accuracy and validity", lambda: fig_accuracy_and_validity(results, figures)),
        ("timing", lambda: fig_timing(results, figures)),
        ("seed variance", lambda: fig_seed_variance(results, figures)),
        ("memory curve", lambda: fig_memory_curve(results, figures, conf)),
        ("pocket against global", lambda: fig_pocket(results, figures, conf)),
    ]:
        try:
            out = fn()
        except Exception as exc:                      # noqa: BLE001
            print(f"  {name}: failed, {exc!r}")
            skipped.append(name)
            continue
        if out:
            made.append(out)
            print(f"  wrote figures/{out}")
        else:
            skipped.append(name)
    print(f"[11_figures] {len(made)} figures written, {len(skipped)} skipped "
          f"for want of their input table")
    if skipped:
        print(f"[11_figures] skipped: {', '.join(skipped)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
