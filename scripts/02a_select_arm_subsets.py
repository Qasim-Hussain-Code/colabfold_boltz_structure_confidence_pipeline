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


def main(argv=None) -> int:
    args = parse_args(argv)
    conf = L.load_conf(args.config)
    paths = L.repo_paths(conf)
    config_dir, results = paths["CONFIG_DIR"], paths["RESULTS_DIR"]

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
