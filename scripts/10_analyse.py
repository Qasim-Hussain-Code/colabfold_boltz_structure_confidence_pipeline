#!/usr/bin/env python3
"""Calibration, thresholds, arm comparisons, seed variance, timing.

This stage recomputes nothing about structures. It reads the tables the
scoring and geometry stages wrote and turns them into the numbers the README
quotes. Every output is a file, and every number in the README comes from one
of them.

The question this stage exists to answer
----------------------------------------
How well does the confidence value predict the accuracy it claims to predict.
Not the correlation: the correlation is one number and it hides everything that
matters. What matters is the scatter, and the behaviour at the top of the
range, because the top of the range is where people make decisions. A model
whose confidence is right on average and wrong for one residue in eight above
its highest band is a model that will mislead anyone who trusts a confident
region.

So the reported quantities are, in order of how much they should be trusted:

  the fraction of residues above the high-confidence band whose measured score
  falls below the accuracy line, which is the number that decides whether a
  confident region can be believed;

  the calibration curve, binned, with the count in each bin, so a reader can
  see which bins carry the evidence;

  the expected calibration error over those bins, which is one number and is
  reported with the binning next to it because it depends on the binning;

  the correlation, reported beside the scatter rather than in place of it.

Intervals are computed by resampling targets, never residues. Residues of one
protein are not independent observations: treating three hundred residues of
one helix bundle as three hundred pieces of evidence gives an interval several
times too narrow. The target is the unit of replication.

Usage
-----
    python scripts/10_analyse.py --config project.conf
"""

from __future__ import annotations

import argparse
import math
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_csc as L  # noqa: E402


def num(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(v) else v


def pearson(xs, ys):
    n = len(xs)
    if n < 3:
        return None
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    if sxx <= 0 or syy <= 0:
        return None
    return sxy / math.sqrt(sxx * syy)


def spearman(xs, ys):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    return pearson(rank(xs), rank(ys))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    conf = L.load_conf(args.config)
    results_dir = Path(conf["RESULTS_DIR"])
    n_bins = L.conf_int(conf, "CALIBRATION_BINS", 10)
    plddt_high = L.conf_float(conf, "PLDDT_HIGH", 90.0)
    lddt_trust = L.conf_float(conf, "LDDT_TRUST", 0.7)
    seed = L.conf_int(conf, "SEED", 0)

    res_path = results_dir / "residue_scores.tsv"
    if not res_path.is_file():
        L.eprint("[error] results/residue_scores.tsv missing; run 06_score_structures.py")
        return 1
    residues = L.read_tsv(res_path)
    # One prediction per target in the calibration pool. A seed-variance run
    # predicts one target several times, and pooling residues over all of them
    # weighs that target once per repeat: five times here, which moved the
    # calibration error of a 30-target arm from 0.0141 to 0.0208. The repeats
    # are what the spread below is computed from and they belong there, not in
    # a pool that is meant to hold one prediction per target.
    canonical = str(conf.get("SEED", "")).strip()
    seeds_per_pair: dict[tuple, set] = defaultdict(set)
    for r in residues:
        seeds_per_pair[(r.get("target_id"), r.get("arm"))].add(r.get("seed", ""))
    repeated = {k for k, v in seeds_per_pair.items() if len(v) > 1}
    if repeated:
        before = len(residues)
        residues = [r for r in residues
                    if (r.get("target_id"), r.get("arm")) not in repeated
                    or r.get("seed", "") == canonical]
        print(f"[10_analyse] {len(repeated)} target and arm pairs were predicted "
              f"more than once; the calibration keeps seed {canonical} for them "
              f"and drops {before - len(residues)} repeat residue rows")
    all_scores = L.read_tsv(results_dir / "structure_scores.tsv") \
        if (results_dir / "structure_scores.tsv").is_file() else []

    # Every per-target statistic below uses one prediction per target. A repeat
    # is a repeat of a target, not another target, and counting it again moved
    # the median accuracy of a thirty-target arm from 0.925 to 0.922 on the
    # strength of five predictions of one protein. The full set is kept under
    # its own name because the spread between repeats is computed from it.
    scores = [s for s in all_scores
              if (s.get('target_id'), s.get('arm')) not in repeated
              or s.get('seed', '') == canonical]
    geometry = L.read_tsv(results_dir / "geometry_checks.tsv") \
        if (results_dir / "geometry_checks.tsv").is_file() else []
    predictions = L.read_tsv(results_dir / "predictions.tsv") \
        if (results_dir / "predictions.tsv").is_file() else []

    # --- per arm, per residue -----------------------------------------------
    by_arm: dict[str, list] = defaultdict(list)
    for r in residues:
        p, l = num(r.get("plddt")), num(r.get("lddt_ca"))
        if p is None or l is None:
            continue
        by_arm[r["arm"]].append((r["target_id"], p, l))

    calib_rows, summary_rows, bin_rows = [], [], []
    for arm, rows in sorted(by_arm.items()):
        xs = [p for _t, p, _l in rows]
        ys = [l for _t, _p, l in rows]
        n = len(rows)
        if n == 0:
            continue

        # The headline reliability number: residues the model was most
        # confident about, whose measurement came back below the line.
        high = [(t, p, l) for t, p, l in rows if p >= plddt_high]
        high_bad = [x for x in high if x[2] < lddt_trust]
        frac_bad = len(high_bad) / len(high) if high else None

        # Grouped by target so the interval resamples proteins, not residues.
        groups: dict[str, list] = defaultdict(list)
        for t, p, l in rows:
            groups[t].append((p, l))

        def frac_bad_stat(picked):
            hi = [(p, l) for g in picked for p, l in g if p >= plddt_high]
            if not hi:
                return float("nan")
            return sum(1 for p, l in hi if l < lddt_trust) / len(hi)

        lo_ci, hi_ci = L.cluster_bootstrap(groups, frac_bad_stat, n_boot=2000, seed=seed)

        # Calibration curve. The confidence value is on a 0 to 100 scale and
        # the measurement on 0 to 1, so the confidence is divided by 100 to put
        # them on one axis. That is the only rescaling anywhere in this
        # repository and it is what the confidence value was defined to mean:
        # a prediction of the score, in the score's own units.
        edges = [i / n_bins for i in range(n_bins + 1)]
        ece, total = 0.0, 0
        for b in range(n_bins):
            lo, hi = edges[b], edges[b + 1]
            sel = [(p / 100.0, l) for _t, p, l in rows
                   if (lo <= p / 100.0 < hi or (b == n_bins - 1 and p / 100.0 == hi))]
            if not sel:
                bin_rows.append({"arm": arm, "bin_low": lo, "bin_high": hi, "n": 0,
                                 "mean_confidence": "", "mean_measured": "",
                                 "gap": "", "recorded": L.now_iso()})
                continue
            mc = sum(c for c, _ in sel) / len(sel)
            mm = sum(m for _, m in sel) / len(sel)
            ece += len(sel) * abs(mc - mm)
            total += len(sel)
            bin_rows.append({"arm": arm, "bin_low": lo, "bin_high": hi, "n": len(sel),
                             "mean_confidence": L.fmt(mc, 4),
                             "mean_measured": L.fmt(mm, 4),
                             "gap": L.fmt(mc - mm, 4), "recorded": L.now_iso()})
        ece = ece / total if total else None

        calib_rows.append({
            "arm": arm,
            "n_residues": n,
            "n_targets": len(groups),
            "pearson_r": L.fmt(pearson(xs, ys), 4),
            "spearman_rho": L.fmt(spearman(xs, ys), 4),
            "expected_calibration_error": L.fmt(ece, 4),
            "calibration_bins": n_bins,
            "mean_confidence": L.fmt(sum(xs) / n / 100.0, 4),
            "mean_measured_lddt_ca": L.fmt(sum(ys) / n, 4),
            "n_residues_above_high_band": len(high),
            "n_of_those_below_trust": len(high_bad),
            "fraction_above_band_below_trust": L.fmt(frac_bad, 4),
            "fraction_ci_low": L.fmt(lo_ci, 4),
            "fraction_ci_high": L.fmt(hi_ci, 4),
            "high_band": plddt_high,
            "trust_line": lddt_trust,
            "recorded": L.now_iso(),
        })

    L.write_tsv(results_dir / "calibration.tsv",
                ["arm", "n_residues", "n_targets", "pearson_r", "spearman_rho",
                 "expected_calibration_error", "calibration_bins", "mean_confidence",
                 "mean_measured_lddt_ca", "n_residues_above_high_band",
                 "n_of_those_below_trust", "fraction_above_band_below_trust",
                 "fraction_ci_low", "fraction_ci_high", "high_band", "trust_line",
                 "recorded"], calib_rows)
    L.write_tsv(results_dir / "calibration_bins.tsv",
                ["arm", "bin_low", "bin_high", "n", "mean_confidence",
                 "mean_measured", "gap", "recorded"], bin_rows)

    # --- per arm, per target -------------------------------------------------
    arm_rows = []
    by_arm_scores: dict[str, list] = defaultdict(list)
    for s in scores:
        if s.get("status") == "ok":
            by_arm_scores[s["arm"]].append(s)
    geom_by = {(g["target_id"], g["arm"]): g for g in geometry}
    for arm, rows in sorted(by_arm_scores.items()):
        lddt = [num(r.get("lddt_ca")) for r in rows]
        lddt = [v for v in lddt if v is not None]
        tm = [v for v in (num(r.get("tm_score")) for r in rows) if v is not None]
        # The verdict used here is whether a structure is inside the range the
        # deposited structures in this same set span, check by check. The
        # stricter verdict, none of anything, is reported beside it: under a
        # third of the deposited structures meet that one, so a rate computed
        # from it says more about the standard than about the model.
        passes = [geom_by.get((r["target_id"], arm), {}).get("within_deposited_range")
                  for r in rows]
        strict = [geom_by.get((r["target_id"], arm), {}).get("passes_all") for r in rows]
        n_geo = sum(1 for p in passes if p in ("1", 1))
        n_geo_known = sum(1 for p in passes if p not in ("", None))
        n_strict = sum(1 for p in strict if p in ("1", 1))
        # Accurate and physically valid, which is the only figure called
        # success here, following the same rule stage 1 used for poses.
        both = 0
        for r in rows:
            g = geom_by.get((r["target_id"], arm), {})
            v = num(r.get("lddt_ca"))
            if v is not None and v >= lddt_trust and                     g.get("within_deposited_range") in ("1", 1):
                both += 1
        q = L.quantiles(lddt) if lddt else [float("nan")] * 3
        arm_rows.append({
            "arm": arm, "n_targets": len(rows),
            "lddt_ca_median": L.fmt(q[1], 4),
            "lddt_ca_q1": L.fmt(q[0], 4), "lddt_ca_q3": L.fmt(q[2], 4),
            "tm_score_median": L.fmt(L.quantiles(tm)[1] if tm else None, 4),
            "n_geometry_checked": n_geo_known,
            "n_passing_geometry": n_geo,
            "fraction_passing_geometry": L.fmt(n_geo / n_geo_known if n_geo_known else None, 4),
            "n_passing_geometry_strict": n_strict,
            "fraction_passing_geometry_strict": L.fmt(
                n_strict / n_geo_known if n_geo_known else None, 4),
            "n_accurate_and_valid": both,
            "fraction_accurate_and_valid": L.fmt(both / len(rows) if rows else None, 4),
            "accuracy_line": lddt_trust,
            "recorded": L.now_iso(),
        })
    L.write_tsv(results_dir / "arm_comparison.tsv",
                ["arm", "n_targets", "lddt_ca_median", "lddt_ca_q1", "lddt_ca_q3",
                 "tm_score_median", "n_geometry_checked", "n_passing_geometry",
                 "fraction_passing_geometry", "n_passing_geometry_strict",
                 "fraction_passing_geometry_strict", "n_accurate_and_valid",
                 "fraction_accurate_and_valid", "accuracy_line", "recorded"], arm_rows)

    # --- arm against arm, on the targets both of them have -------------------
    # The arms no longer cover the same set. An alignment arm costs about
    # twelve times what the single-sequence arm costs and grows with the square
    # of the length, so each was given a budget and each bought what it could;
    # config/arm_subsets.tsv records which targets that was. Comparing the
    # medians in the table above would then be comparing two different sets of
    # proteins and calling the difference an effect of the arm.
    #
    # So every pair of arms is also compared on the targets both of them
    # predicted, as a difference per target. The interval resamples targets,
    # which is the unit of replication, and the count of shared targets is
    # reported beside it because with twelve of them it is the number that
    # decides what the comparison can support.
    lddt_by = {}
    for s in scores:
        if s.get("status") == "ok" and num(s.get("lddt_ca")) is not None:
            lddt_by[(s["arm"], s["target_id"])] = num(s["lddt_ca"])
    arms_present = sorted({a for a, _t in lddt_by})
    pair_rows = []
    for i, a in enumerate(arms_present):
        for b in arms_present[i + 1:]:
            shared = sorted({t for (arm, t) in lddt_by if arm == a}
                            & {t for (arm, t) in lddt_by if arm == b})
            if len(shared) < 3:
                continue
            diffs = {t: [lddt_by[(b, t)] - lddt_by[(a, t)]] for t in shared}

            def med_diff(picked):
                flat = [v for grp in picked for v in grp]
                return L.quantiles(flat)[1] if flat else float("nan")

            lo, hi = L.cluster_bootstrap(diffs, med_diff, n_boot=2000, seed=11)
            b_better = sum(1 for t in shared if lddt_by[(b, t)] > lddt_by[(a, t)])
            pair_rows.append({
                "arm_a": a, "arm_b": b, "n_shared_targets": len(shared),
                "median_lddt_ca_a": L.fmt(L.quantiles([lddt_by[(a, t)] for t in shared])[1], 4),
                "median_lddt_ca_b": L.fmt(L.quantiles([lddt_by[(b, t)] for t in shared])[1], 4),
                "median_paired_difference": L.fmt(
                    L.quantiles([lddt_by[(b, t)] - lddt_by[(a, t)] for t in shared])[1], 4),
                "ci_low": L.fmt(lo, 4), "ci_high": L.fmt(hi, 4),
                "n_targets_b_higher": b_better,
                "note": "b minus a, per target, over the targets both arms have; "
                        "the interval resamples targets",
                "recorded": L.now_iso(),
            })
    L.write_tsv(results_dir / "arm_pairs.tsv",
                ["arm_a", "arm_b", "n_shared_targets", "median_lddt_ca_a",
                 "median_lddt_ca_b", "median_paired_difference", "ci_low",
                 "ci_high", "n_targets_b_higher", "note", "recorded"], pair_rows)
    for r in pair_rows:
        print(f"  {r['arm_a']} against {r['arm_b']}: {r['n_shared_targets']} shared "
              f"targets, median difference {r['median_paired_difference']} "
              f"[{r['ci_low']}, {r['ci_high']}]")

    # --- the threshold, derived rather than chosen ---------------------------
    # The question a reader has is where to stop trusting a confident region.
    # It is answered by sweeping the band and reporting, at each value, the
    # fraction of residues above it that came back below the accuracy line.
    # No single value is privileged: the sweep is the answer, and the README
    # says the threshold is arbitrary unless the sweep shows a knee.
    sweep_rows = []
    for arm, rows in sorted(by_arm.items()):
        for band in range(50, 100, 5):
            hi = [(t, p, l) for t, p, l in rows if p >= band]
            if not hi:
                continue
            bad = sum(1 for _t, _p, l in hi if l < lddt_trust)
            sweep_rows.append({
                "arm": arm, "band": band, "n_residues_at_or_above": len(hi),
                "n_below_trust": bad,
                "fraction_below_trust": L.fmt(bad / len(hi), 4),
                "trust_line": lddt_trust, "recorded": L.now_iso(),
            })
    L.write_tsv(results_dir / "threshold_sweep.tsv",
                ["arm", "band", "n_residues_at_or_above", "n_below_trust",
                 "fraction_below_trust", "trust_line", "recorded"], sweep_rows)

    # --- seed variance -------------------------------------------------------
    seed_groups: dict[tuple, list] = defaultdict(list)
    for s in all_scores:
        if s.get("status") == "ok" and s.get("seed"):
            seed_groups[(s["target_id"], s["arm"])].append(s)
    seed_rows = []
    for (tid, arm), rows in sorted(seed_groups.items()):
        if len(rows) < 2:
            continue
        vals = [num(r.get("lddt_ca")) for r in rows]
        vals = [v for v in vals if v is not None]
        if len(vals) < 2:
            continue
        mean = sum(vals) / len(vals)
        sd = math.sqrt(sum((v - mean) ** 2 for v in vals) / (len(vals) - 1))
        seed_rows.append({
            "target_id": tid, "arm": arm, "n_seeds": len(vals),
            "lddt_ca_values": ", ".join(f"{v:.4f}" for v in sorted(vals)),
            "range": L.fmt(max(vals) - min(vals), 4),
            "standard_deviation": L.fmt(sd, 4),
            "recorded": L.now_iso(),
        })
    L.write_tsv(results_dir / "seed_variance.tsv",
                ["target_id", "arm", "n_seeds", "lddt_ca_values", "range",
                 "standard_deviation", "recorded"], seed_rows)

    # --- cost ----------------------------------------------------------------
    timing_rows = []
    by_arm_pred: dict[str, list] = defaultdict(list)
    for p in predictions:
        if p.get("status") == "ok":
            by_arm_pred[p["arm"]].append(p)
    for arm, rows in sorted(by_arm_pred.items()):
        secs = [v for v in (num(r.get("elapsed_s")) for r in rows) if v is not None]
        mems = [v for v in (num(r.get("peak_rss_mb")) for r in rows) if v is not None]
        lens = [v for v in (num(r.get("sequence_length")) for r in rows) if v is not None]
        qs = L.quantiles(secs, (0.5, 0.9, 1.0)) if secs else [float("nan")] * 3
        qm = L.quantiles(mems, (0.5, 0.9, 1.0)) if mems else [float("nan")] * 3
        timing_rows.append({
            "arm": arm, "n": len(rows),
            "seconds_median": L.fmt(qs[0], 1), "seconds_p90": L.fmt(qs[1], 1),
            "seconds_max": L.fmt(qs[2], 1),
            "peak_rss_mb_median": L.fmt(qm[0], 1), "peak_rss_mb_p90": L.fmt(qm[1], 1),
            "peak_rss_mb_max": L.fmt(qm[2], 1),
            "residues_median": L.fmt(L.quantiles(lens)[1] if lens else None, 0),
            "total_hours": L.fmt(sum(secs) / 3600.0 if secs else None, 2),
            "recorded": L.now_iso(),
        })
    L.write_tsv(results_dir / "timing.tsv",
                ["arm", "n", "seconds_median", "seconds_p90", "seconds_max",
                 "peak_rss_mb_median", "peak_rss_mb_p90", "peak_rss_mb_max",
                 "residues_median", "total_hours", "recorded"], timing_rows)

    # --- the headline numbers, in one file -----------------------------------
    headline = []
    for row in calib_rows:
        headline.append({
            "key": f"fraction_above_{int(plddt_high)}_below_{lddt_trust}__{row['arm']}",
            "value": row["fraction_above_band_below_trust"],
            "note": (f"{row['n_of_those_below_trust']} of {row['n_residues_above_high_band']} "
                     f"residues above confidence {int(plddt_high)} measured below "
                     f"{lddt_trust} lDDT-CA, over {row['n_targets']} targets"),
        })
    for row in arm_rows:
        headline.append({
            "key": f"lddt_ca_median__{row['arm']}", "value": row["lddt_ca_median"],
            "note": f"over {row['n_targets']} targets",
        })
    L.write_tsv(results_dir / "headline.tsv", ["key", "value", "note"], headline)

    print(f"[10_analyse] {len(residues)} residue rows over {len(by_arm)} arms")
    for row in calib_rows:
        print(f"  {row['arm']}: correlation {row['pearson_r']}, "
              f"calibration error {row['expected_calibration_error']}, "
              f"{row['n_of_those_below_trust']}/{row['n_residues_above_high_band']} "
              f"confident residues below the accuracy line")
    return 0


if __name__ == "__main__":
    sys.exit(main())
