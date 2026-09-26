#!/usr/bin/env python3
"""
=============================================================================
 lib_memory.py - how long a sequence this machine can predict, from measurement
=============================================================================
 project.conf carries a memory curve, and that curve was fitted during
 calibration on single-sequence runs. Predicting from an alignment costs far
 more memory at the same length, because the Evoformer carries the alignment
 through every block. Measured here: at 179 residues a single-sequence
 prediction peaked at 2.9 GB and the same target with an alignment peaked at
 4.1 GB, and the gap widens with length.

 So the ceiling in project.conf is a single-sequence ceiling. Applying it to an
 alignment run says 612 residues are affordable when the real answer is closer
 to 230, and the difference is not a slow run: it is the kernel killing the
 process partway through.

 This fits the curve separately for whichever arms are asked about, using the
 peak memory every finished prediction already recorded, and answers two
 questions: what a given length will cost, and what length fits a budget. It
 improves as the run proceeds, because every prediction adds a point.

 The form is the same one the configuration uses:

     peak_mb = base + quad * (length / 1000)^2

 base is taken as the smallest peak observed, which is the interpreter, the
 weights and the framework before any sequence-dependent allocation. quad is
 least squares on the rest.

 Usage:
     python scripts/lib_memory.py --config project.conf --arms af2_msa
     python scripts/lib_memory.py --config project.conf --arms af2_msa --length 350
     python scripts/lib_memory.py --config project.conf --arms af2_msa \\
         --budget-mb 5120 --quiet
=============================================================================
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_csc as L  # noqa: E402

# How much of the budget a projection is allowed to claim. The fit is made on
# lengths well below the one being asked about, so the extrapolation carries
# real uncertainty, and the failure it guards against is a kill rather than a
# slow run. Spending the whole budget on a central estimate would be the wrong
# way round.
SAFETY = 0.85


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="memory against sequence length")
    here = Path(__file__).resolve().parent.parent
    p.add_argument("--config", default=str(here / "project.conf"))
    p.add_argument("--arms", default="af2_msa",
                   help="prefix of the arms to fit against")
    p.add_argument("--length", type=int, default=0)
    p.add_argument("--budget-mb", type=int, default=0)
    p.add_argument("--quiet", action="store_true",
                   help="print only the answer, for a shell to read")
    p.add_argument("--measured", default="",
                   help="LENGTH:PEAK_MB actually observed for the thing being "
                        "sized. A measurement of the sequence in hand beats an "
                        "extrapolation from other sequences, so the longer of "
                        "the two allowances is returned.")
    return p.parse_args(argv)


def fit(conf: dict, arm_prefix: str):
    paths = L.repo_paths(conf)
    targets = L.read_tsv(paths["CONFIG_DIR"] / "targets.tsv")
    length = {t["target_id"]: int(t["sequence_length"]) for t in targets}
    pts = []
    for r in L.read_tsv(paths["RESULTS_DIR"] / "predictions.tsv"):
        if not r.get("arm", "").startswith(arm_prefix):
            continue
        peak, tid = r.get("peak_rss_mb"), r.get("target_id")
        if not peak or tid not in length:
            continue
        try:
            pts.append((length[tid], float(peak)))
        except ValueError:
            continue
    if len(pts) < 3:
        return None
    base = min(p[1] for p in pts)
    num = sum((p[1] - base) * (p[0] / 1000.0) ** 2 for p in pts)
    den = sum(((p[0] / 1000.0) ** 2) ** 2 for p in pts)
    quad = num / den if den else 0.0
    return {"base": base, "quad": quad, "n": len(pts),
            "min_length": min(p[0] for p in pts),
            "max_length": max(p[0] for p in pts),
            "largest_peak": max(p[1] for p in pts)}


def project(f: dict, length: int) -> float:
    return f["base"] + f["quad"] * (length / 1000.0) ** 2


def longest_within(f: dict, budget_mb: int) -> int:
    room = budget_mb * SAFETY - f["base"]
    if room <= 0 or f["quad"] <= 0:
        return 0
    return int(1000.0 * (room / f["quad"]) ** 0.5)


def main(argv=None) -> int:
    args = parse_args(argv)
    conf = L.load_conf(args.config)
    f = fit(conf, args.arms)
    if f is None:
        if not args.quiet:
            print(f"[memory] fewer than three finished predictions match "
                  f"'{args.arms}'; nothing to fit")
        print(0)
        return 1

    if not args.quiet:
        print(f"[memory] {f['n']} predictions matching '{args.arms}', lengths "
              f"{f['min_length']} to {f['max_length']}, largest peak "
              f"{f['largest_peak']:.0f} MB")
        print(f"[memory] peak_mb = {f['base']:.0f} + {f['quad']:.0f} * "
              f"(length/1000)^2")
        for n in (150, 200, 250, 300, 350, 400):
            print(f"[memory]   {n:4d} residues projects to {project(f, n):7.0f} MB")

    if args.length:
        if args.quiet:
            print(int(project(f, args.length)))
        else:
            print(f"[memory] {args.length} residues projects to "
                  f"{project(f, args.length):.0f} MB")
        return 0

    budget = args.budget_mb or (L.conf_int(conf, "RAM_GB", 5) * 1024)
    n = longest_within(f, budget)

    # A measurement of the sequence being sized, where one exists, beats an
    # extrapolation from other sequences. Memory depends on how deep the
    # alignment is as well as how long the sequence is, and a fit made over
    # targets whose alignments run to thousands of sequences badly overstates
    # the cost of one whose alignment holds 89. Measured here: a 201-residue
    # construct with an 89-sequence alignment peaked at 3478 MB where that fit
    # predicted about 4900, and the difference cost the construct its body.
    if args.measured:
        try:
            m_len, m_peak = args.measured.split(":")
            m_len, m_peak = int(m_len), float(m_peak)
        except ValueError:
            m_len = 0
        if m_len > 0 and m_peak > f["base"]:
            q = (m_peak - f["base"]) / (m_len / 1000.0) ** 2
            from_measured = longest_within({"base": f["base"], "quad": q}, budget)
            if not args.quiet:
                print(f"[memory] measured {m_peak:.0f} MB at {m_len} residues, which "
                      f"implies {q:.0f} MB per kres^2 and allows {from_measured} "
                      f"residues")
            n = max(n, from_measured)
    if args.quiet:
        print(n)
    else:
        print(f"[memory] against a {budget} MB budget, and claiming at most "
              f"{int(SAFETY * 100)} per cent of it because the fit is an "
              f"extrapolation, the longest sequence that fits is {n} residues")
    return 0


if __name__ == "__main__":
    sys.exit(main())
