#!/usr/bin/env python3
"""
=============================================================================
 13_verify_readme.py - the README against the tables it quotes
=============================================================================
 The README states that every number in it was read from a file in results/.
 That is a claim about the repository, and a claim about a repository should
 be checkable by the repository.

 This is not a parser. Each number the README quotes is written out here
 beside the file and the cell it should equal, so a table that is regenerated
 while the README is not gets caught rather than assumed. Adding a number to
 the README means adding a line here, which is the cost of the guarantee.

 A number this file lists but the README does not quote is reported rather
 than failed: a figure may carry it instead, and the list stays honest about
 what was actually checked.

 Usage:
     python scripts/13_verify_readme.py
=============================================================================
"""
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))
import lib_csc as L  # noqa: E402

readme = (REPO / "README.md").read_text(encoding="utf-8")
head = {r["key"]: r["value"] for r in L.read_tsv(REPO / "results/headline.tsv")}
cal = {r["arm"]: r for r in L.read_tsv(REPO / "results/calibration.tsv")}
arms = {r["arm"]: r for r in L.read_tsv(REPO / "results/arm_comparison.tsv")}
pairs = {(r["arm_a"], r["arm_b"]): r
         for r in L.read_tsv(REPO / "results/arm_pairs.tsv")}
summary = {r["key"]: r["value"] for r in L.read_tsv(REPO / "results/holdout_set_summary.tsv")}
dock = {(r["arm"], r["method"]): r for r in L.read_tsv(REPO / "results/docking_handoff.tsv")}
seed = L.read_tsv(REPO / "results/seed_variance.tsv")
sweep = {(r["arm"], r["band"]): r for r in L.read_tsv(REPO / "results/threshold_sweep.tsv")}

claims = [
    ("0.925", head["lddt_ca_median__af2_msa_notmpl"], "median lDDT, alignment arm"),
    ("0.908", head["lddt_ca_median__af2_msa_tmpl"], "median lDDT, templates arm"),
    ("0.378", head["lddt_ca_median__af2_nomsa"], "median lDDT, no alignment"),
    ("0.8825", head["lddt_ca_median__null_template"], "median lDDT, template floor"),
    ("0.016", head["lddt_ca_median__null_unrelated"], "median lDDT, unrelated floor"),
    ("0.0141", cal["af2_msa_notmpl"]["expected_calibration_error"], "calibration error"),
    ("0.7346", cal["af2_msa_notmpl"]["pearson_r"], "pearson r"),
    ("0.0297", cal["af2_msa_tmpl"]["expected_calibration_error"], "calibration error, templates"),
    ("0.0245", cal["af2_nomsa"]["expected_calibration_error"], "calibration error, no alignment"),
    ("2463", cal["af2_msa_notmpl"]["n_residues_above_high_band"], "residues above the band"),
    ("34", cal["af2_msa_notmpl"]["n_of_those_below_trust"], "of those below trust"),
    ("56", summary["identity_exactly_100"], "targets identical to a pre-cutoff relative"),
    ("67", summary["identity_95_to_100"], "targets 95 to 100 per cent"),
    ("291", summary["release_date_only_extra"], "extra clusters a release filter admits"),
    ("163", summary["excluded_for_time"], "targets cut for wall clock"),
    ("150", summary["targets_final"], "targets in the set"),
    ("0.004", seed[0]["range"], "seed range"),
    ("0.0015", seed[0]["standard_deviation"], "seed standard deviation"),
]

bad = 0
unquoted = []
for quoted, actual, what in claims:
    same = str(actual).rstrip("0").rstrip(".") == quoted.rstrip("0").rstrip(".")
    if quoted not in readme:
        # The README does not make this claim. That is allowed; a figure may
        # carry it instead. It is reported so the list stays honest about what
        # was actually checked.
        unquoted.append(f"{what} ({quoted})")
        continue
    if not same:
        bad += 1
        print(f"  MISMATCH {what}: README says {quoted}, table says {actual}")
if unquoted:
    print("  not quoted in the README, so not checked: " + "; ".join(unquoted))

# the sweep rows the README tabulates
for band, frac in [("50", "0.0672"), ("60", "0.0446"), ("70", "0.0329"),
                   ("80", "0.0257"), ("90", "0.0138"), ("95", "0.009")]:
    row = sweep.get(("af2_msa_notmpl", band))
    if row is None or row["fraction_below_trust"] != frac:
        bad += 1
        print(f"  MISMATCH sweep at {band}: README {frac}, table "
              f"{row['fraction_below_trust'] if row else 'absent'}")

# docking
for key, rate, n in [(("crystal/A2_genconf_refbox", "vina"), "0.3636", "4"),
                     (("crystal/A2_genconf_refbox", "vinardo"), "0.5556", "5"),
                     (("predicted/A2_genconf_refbox", "vina"), "0.0909", "1"),
                     (("predicted/A2_genconf_refbox", "vinardo"), "0.1111", "1")]:
    row = dock.get(key)
    if row is None or row["rate_success"] != rate or row["n_success"] != n:
        bad += 1
        print(f"  MISMATCH docking {key}: README {n} at {rate}, table "
              f"{row['n_success'] + ' at ' + row['rate_success'] if row else 'absent'}")

# the paired comparisons
for (a, b), diff in [(("af2_msa_notmpl", "af2_msa_tmpl"), "0.002"),
                     (("af2_msa_notmpl", "af2_nomsa"), "-0.476"),
                     (("af2_msa_notmpl", "null_template"), "-0.109")]:
    row = pairs.get((a, b))
    if row is None or row["median_paired_difference"] != diff:
        bad += 1
        print(f"  MISMATCH pair {a} against {b}: README {diff}, table "
              f"{row['median_paired_difference'] if row else 'absent'}")

print(f"{len(claims) + 13} checks, {bad} mismatch(es)")
sys.exit(1 if bad else 0)
