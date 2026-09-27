#!/usr/bin/env python3
"""
=============================================================================
 02a_select_arm_subsets.py - which targets each arm can afford
=============================================================================
 The alignment arms cost about twelve times what the single-sequence arm
 costs, and the cost grows with the square of the sequence length. Measured
 on this machine, one alignment arm over the whole 150-target set is 76
 hours and the two of them together are 150. That does not fit, and the
 design says to cut the set rather than to cut an arm.

 This decides the cut, from times that were measured rather than assumed:

   t = a * L^2 * f

 where a is fitted by least squares against every completed single-sequence
 prediction, and f is the ratio measured on the alignment predictions that
 have finished. f is fitted twice, once for the targets whose alignment came
 back deep and once for those where the server found almost nothing, because
 the second group costs about what the single-sequence arm costs and
 treating them alike would misprice a fifth of the set.

 Three rules decide the membership.

 1. A target the docking handoff needs is in the arm that handoff reads.
    Docking into a predicted receptor cannot be measured for a target that
    was never predicted, and those 15 targets were chosen earlier on
    criteria of their own.

 2. A target already predicted stays. It costs nothing to keep and it adds
    to the count the interval is computed over.

 3. Everything else is chosen to span the length range, not to be cheap.
    Taking the shortest targets buys roughly three times as many of them,
    and every statement the arm supports would then be a statement about
    small single-domain proteins. The spread is worth more than the count.

 The templates arm draws from the no-templates arm, so that the comparison
 between them is paired on a target rather than across two different sets.

 What it writes:
    config/arm_subsets.tsv    target_id, arm, why, projected_s
    results/excluded.tsv      one row per target an arm cannot afford

 Usage:
     python scripts/02a_select_arm_subsets.py --hours af2_msa_notmpl=10,af2_msa_tmpl=5
     python scripts/02a_select_arm_subsets.py --hours af2_msa_notmpl=10 --dry-run
=============================================================================
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_csc as L  # noqa: E402

SUBSET_COLUMNS = ["target_id", "arm", "why", "sequence_length", "msa_depth",
                  "projected_s"]
DEEP_ALIGNMENT = 100


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="which targets each arm can afford")
    here = Path(__file__).resolve().parent.parent
    p.add_argument("--config", default=str(here / "project.conf"))
    p.add_argument("--hours", required=True,
                   help="arm=hours pairs, comma separated")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--from-results", action="store_true",
                   help="record the set each arm actually covers, rather than "
                        "the set a budget buys. A fresh run then reproduces "
                        "the arms this repository reports on.")
    return p.parse_args(argv)


def fit_cost(preds, length):
    """t = a * L^2 for the single-sequence arm, and what an alignment multiplies
    it by. Both from the runs that have finished, not from a guess."""
    num = den = 0.0
    for r in preds:
        if r.get("arm") != "af2_nomsa" or not r.get("elapsed_s"):
            continue
        n = length.get(r["target_id"])
        if not n:
            continue
        num += float(r["elapsed_s"]) * n * n
        den += (n * n) ** 2
    if den == 0:
        raise SystemExit("[02a] no completed single-sequence predictions to fit against")
    a = num / den

    deep, shallow = [], []
    for r in preds:
        if not r.get("arm", "").startswith("af2_msa") or not r.get("elapsed_s"):
            continue
        n = length.get(r["target_id"])
        if not n:
            continue
        ratio = float(r["elapsed_s"]) / (a * n * n)
        (deep if int(r.get("msa_depth") or 0) > DEEP_ALIGNMENT else shallow).append(ratio)
    f_deep = sum(deep) / len(deep) if deep else 12.0
    f_shallow = sum(shallow) / len(shallow) if shallow else 1.1
    return a, f_deep, f_shallow, len(deep), len(shallow)


def spread(candidates, k):
    """k targets spread evenly over the length-sorted candidates."""
    if k <= 0 or not candidates:
        return []
    if k >= len(candidates):
        return list(candidates)
    step = (len(candidates) - 1) / (k - 1) if k > 1 else 0
    picked, seen = [], set()
    for i in range(k):
        j = round(i * step)
        while j in seen and j < len(candidates) - 1:
            j += 1
        if j in seen:
            continue
        seen.add(j)
        picked.append(candidates[j])
    return picked


SELECTION_COST_COLUMNS = ["arm", "n_targets", "shortest_first_n_targets",
                          "ratio", "length_min", "length_max",
                          "shortest_first_length_max", "projected_hours",
                          "seconds_per_kilo_residue_squared", "n_runs_fitted",
                          "recorded"]


def record_selection_cost(conf: dict, rows: list[dict]) -> None:
    """What taking the shortest targets first would have bought instead.

    The README says the arms span the length range rather than taking the
    cheapest targets, and that spanning it costs targets. How many it costs is
    a number, and a number in the README has to come from a file. This writes
    that file rather than leaving the claim to be estimated: it fits the cost
    model on the runs that finished, then spends the same total on the
    shortest targets in the set and reports how many that buys.
    """
    paths = L.repo_paths(conf)
    results, config_dir = paths["RESULTS_DIR"], paths["CONFIG_DIR"]
    targets = {r["target_id"]: r
               for r in L.read_tsv(config_dir / "targets.tsv")
               if r.get("sequence_length")}
    preds = L.read_tsv(results / "predictions.tsv")         if (results / "predictions.tsv").is_file() else []
    out_rows = []
    for arm in sorted({r["arm"] for r in rows}):
        runs = [(int(p["sequence_length"]), float(p["elapsed_s"])) for p in preds
                if p.get("arm") == arm and p.get("status") == "ok"
                and p.get("elapsed_s") and p.get("sequence_length")]
        chosen = sorted(int(targets[r["target_id"]]["sequence_length"])
                        for r in rows
                        if r["arm"] == arm and r["target_id"] in targets)
        if len(runs) < 3 or not chosen:
            continue
        denom = sum((length / 1000.0) ** 2 for length, _s in runs)
        if denom <= 0:
            continue
        per_kres2 = sum(s for _l, s in runs) / denom

        def cost(length: int) -> float:
            return per_kres2 * (length / 1000.0) ** 2

        spent = sum(cost(x) for x in chosen)
        taken: list[int] = []
        spend = 0.0
        for length in sorted(int(r["sequence_length"])
                             for r in targets.values()):
            if spend + cost(length) > spent:
                break
            spend += cost(length)
            taken.append(length)
        if not taken:
            continue
        out_rows.append({
            "arm": arm, "n_targets": len(chosen),
            "shortest_first_n_targets": len(taken),
            "ratio": L.fmt(len(taken) / len(chosen), 2),
            "length_min": chosen[0], "length_max": chosen[-1],
            "shortest_first_length_max": taken[-1],
            "projected_hours": L.fmt(spent / 3600.0, 2),
            "seconds_per_kilo_residue_squared": L.fmt(per_kres2, 0),
            "n_runs_fitted": len(runs),
            "recorded": L.now_iso(),
        })
        print(f"[02a] {arm}: {len(chosen)} targets spanning {chosen[0]} to "
              f"{chosen[-1]} residues; the shortest first on the same spend "
              f"would be {len(taken)} targets up to {taken[-1]} residues")
    if out_rows:
        L.write_tsv(results / "arm_selection_cost.tsv",
                    SELECTION_COST_COLUMNS, out_rows)


def main(argv=None) -> int:
    args = parse_args(argv)
    conf = L.load_conf(args.config)
    paths = L.repo_paths(conf)
    config_dir, results = paths["CONFIG_DIR"], paths["RESULTS_DIR"]

    if args.from_results:
        # What the arms cover, read back from the predictions themselves.
        #
        # The budget mode answers "what can this arm afford", and the answer
        # changes as the cost fit improves: with 41 alignment predictions to
        # fit against rather than 14, the same budget buys 51 targets where it
        # bought 30. That is the right answer to that question and the wrong
        # file for a fresh clone to read, because the arms this repository
        # reports on are the ones it ran. This mode records those, so
        # 05_predict.sh reproduces them.
        preds = L.read_tsv(results / "predictions.tsv")
        targets = {t["target_id"]: t for t in L.read_tsv(config_dir / "targets.tsv")}
        depth = {r["target_id"]: r.get("depth", "")
                 for r in L.read_tsv(results / "msa_depth.tsv")}
        seen: dict[tuple, dict] = {}
        for r in preds:
            arm = r.get("arm", "")
            tid = r.get("target_id", "")
            if not arm.startswith("af2_msa") or r.get("status") != "ok":
                continue
            if tid not in targets or (tid, arm) in seen:
                continue
            seen[(tid, arm)] = {
                "target_id": tid, "arm": arm,
                "why": "predicted in the run this repository reports on",
                "sequence_length": targets[tid].get("sequence_length", ""),
                "msa_depth": depth.get(tid, ""),
                "projected_s": r.get("elapsed_s", ""),
            }
        rows = [seen[k] for k in sorted(seen, key=lambda k: (k[1], int(seen[k]["sequence_length"] or 0)))]
        by_arm: dict[str, int] = {}
        for r in rows:
            by_arm[r["arm"]] = by_arm.get(r["arm"], 0) + 1
        for arm, n in sorted(by_arm.items()):
            print(f"[02a] {arm}: {n} targets, read back from the predictions")
        if args.dry_run:
            print("[02a] dry run; nothing written")
            return 0
        out = config_dir / "arm_subsets.tsv"
        L.write_tsv(out, SUBSET_COLUMNS, rows)
        print(f"[02a] wrote {out}")
        record_selection_cost(conf, rows)
        return 0

    budgets: dict[str, float] = {}
    for part in args.hours.split(","):
        arm, _, hours = part.partition("=")
        budgets[arm.strip()] = float(hours) * 3600.0

    targets = L.read_tsv(config_dir / "targets.tsv")
    length = {t["target_id"]: int(t["sequence_length"]) for t in targets}
    preds = L.read_tsv(results / "predictions.tsv")
    depth = {r["target_id"]: int(r["depth"] or 0)
             for r in L.read_tsv(results / "msa_depth.tsv") if r.get("depth")}

    a, f_deep, f_shallow, n_deep, n_shallow = fit_cost(preds, length)
    print(f"[02a] fitted t = {a:.4e} * L^2 seconds without an alignment")
    print(f"[02a] an alignment multiplies that by {f_deep:.1f} when it is deep "
          f"(n={n_deep}) and {f_shallow:.1f} when the server found almost "
          f"nothing (n={n_shallow})")

    def cost(tid: str) -> float:
        n = length[tid]
        f = f_deep if depth.get(tid, 10_000) > DEEP_ALIGNMENT else f_shallow
        return a * n * n * f

    dock_file = config_dir / "docking_subset.tsv"
    mandatory = {r["target_id"] for r in L.read_tsv(dock_file)} if dock_file.is_file() else set()
    dock_arm = conf.get("DOCKING_RECEPTOR_ARM", "af2_msa_notmpl")

    by_length = sorted(length, key=lambda t: (length[t], t))
    rows: list[dict] = []
    chosen: dict[str, list[str]] = {}

    # The no-templates arm first. The templates arm then draws from it, so the
    # order here is not cosmetic.
    for arm in sorted(budgets, key=lambda x: (x != dock_arm, x)):
        budget = budgets[arm]
        done = {r["target_id"] for r in preds if r.get("arm") == arm}
        pool = chosen.get(dock_arm) if arm != dock_arm and dock_arm in chosen else by_length
        pool = [t for t in by_length if t in set(pool)]

        picked: dict[str, str] = {}
        spent = 0.0
        for t in pool:
            if t in done:
                picked[t] = "already predicted, so it costs nothing to keep"
        # Only the arm the handoff reads has to carry the docking targets, and
        # it carries all of them whatever the budget says, because a docking
        # measurement on a receptor that was never predicted does not exist.
        # The templates arm has no such obligation and is free to spend its
        # smaller budget on spread instead.
        if arm == dock_arm:
            for t in [t for t in pool if t in mandatory and t not in picked]:
                picked[t] = "the docking handoff reads this arm for its receptor"
                spent += cost(t)

        # Then as many of the rest as the remainder buys, spread over the
        # length range rather than taken from the cheap end.
        rest = [t for t in pool if t not in picked]
        best: list[str] = []
        for k in range(len(rest), 0, -1):
            trial = spread(rest, k)
            if spent + sum(cost(t) for t in trial) <= budget:
                best = trial
                break
        for t in best:
            picked[t] = "chosen to spread the set over the length range"
            spent += cost(t)

        chosen[arm] = [t for t in by_length if t in picked]
        lens = [length[t] for t in chosen[arm]]
        print(f"[02a] {arm}: {len(chosen[arm])} targets, lengths {min(lens)} to "
              f"{max(lens)}, {spent / 3600:.1f} h of inference still to run "
              f"against a budget of {budget / 3600:.1f} h")

        for t in chosen[arm]:
            rows.append(dict(target_id=t, arm=arm, why=picked[t],
                             sequence_length=length[t],
                             msa_depth=depth.get(t, ""),
                             projected_s=f"{cost(t):.0f}"))

    if args.dry_run:
        print("[02a] dry run; nothing written")
        return 0

    out = config_dir / "arm_subsets.tsv"
    L.write_tsv(out, SUBSET_COLUMNS, rows)
    print(f"[02a] wrote {out}")

    for arm, keep in chosen.items():
        budget_h = budgets[arm] / 3600.0
        # Same reasoning as the stage that builds the set: this arm's subset
        # has just been decided again, so its earlier rows describe a subset
        # that no longer exists.
        stale = L.clear_exclusions(results, "02a_select_arm_subsets", arm=arm)
        if stale:
            print(f"[02a] cleared {stale} exclusion rows from an earlier "
                  f"selection for {arm}")
        for t in by_length:
            if t in keep:
                continue
            L.record_exclusion(
                results, t, "02a_select_arm_subsets",
                f"beyond the {budget_h:g} hours of inference this arm was given",
                detail=f"projected {cost(t):.0f} s at {length[t]} residues",
                arm=arm)
    print("[02a] every target an arm cannot afford is recorded in excluded.tsv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
