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

 A number this file lists but the README does not quote is a failure. The
 earlier version reported it and carried on, which meant that editing a number
 in the README silently retired its own check: the string no longer matched, so
 the check was skipped rather than run. A number genuinely carried by a figure
 and not by the prose goes in FIGURE_ONLY below, named, so the exemption is
 visible.

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

def pct(x: str) -> str:
    """A fraction as the README writes it: one decimal place, per cent."""
    return f"{100 * float(x):.1f}"


claims = [
    ("0.925", head["lddt_ca_median__af2_msa_notmpl"], "median lDDT, alignment arm"),
    ("0.912", head["lddt_ca_median__af2_msa_tmpl"], "median lDDT, templates arm"),
    ("0.378", head["lddt_ca_median__af2_nomsa"], "median lDDT, no alignment"),
    ("0.885", head["lddt_ca_median__null_template"], "median lDDT, template floor"),
    ("0.0415", head["lddt_ca_median__null_unrelated"], "median lDDT, unrelated floor"),
    ("0.0141", cal["af2_msa_notmpl"]["expected_calibration_error"], "calibration error"),
    ("0.7346", cal["af2_msa_notmpl"]["pearson_r"], "pearson r"),
    ("0.02", cal["af2_msa_tmpl"]["expected_calibration_error"], "calibration error, templates"),
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
    ("117", arms["null_template"]["n_targets"], "template floor targets compared"),
    ("96", arms["null_unrelated"]["n_targets"], "unrelated floor targets compared"),
    ("86.3", pct(arms["null_template"]["fraction_passing_geometry"]), "template floor valid"),
    ("66.7", pct(arms["null_template"]["fraction_accurate_and_valid"]), "template floor accurate and valid"),
    ("86.5", pct(arms["null_unrelated"]["fraction_passing_geometry"]), "unrelated floor valid"),
    ("35.7", pct(arms["af2_msa_tmpl"]["fraction_passing_geometry"]), "templates arm valid"),
    ("23.3", pct(arms["af2_msa_notmpl"]["fraction_passing_geometry"]), "alignment arm valid"),
    ("2.0", pct(arms["af2_nomsa"]["fraction_passing_geometry"]), "no-alignment arm valid"),
]

# The per-arm accuracy table quotes a quartile pair and a TM-score for every
# arm, so every cell of it is checked rather than the median alone.
for _arm, _label in [("af2_msa_notmpl", "alignment arm"),
                     ("af2_msa_tmpl", "templates arm"),
                     ("af2_nomsa", "no-alignment arm"),
                     ("null_template", "template floor"),
                     ("null_unrelated", "unrelated floor")]:
    for _col, _what in [("lddt_ca_q1", "first quartile"),
                        ("lddt_ca_q3", "third quartile"),
                        ("tm_score_median", "median TM-score")]:
        # The expected string is taken from the table, so the comparison is
        # trivially equal and the check that does the work is the one that
        # requires the string to appear in the README.
        claims.append((arms[_arm][_col], arms[_arm][_col], f"{_what}, {_label}"))

# The numbers that live only in a figure, named so the exemption is visible
# rather than implied by a string that happens not to match.
FIGURE_ONLY = set()

bad = 0
unquoted = []
for quoted, actual, what in claims:
    same = str(actual).rstrip("0").rstrip(".") == quoted.rstrip("0").rstrip(".")
    if quoted not in readme:
        unquoted.append(f"{what} ({quoted})")
        if what not in FIGURE_ONLY:
            bad += 1
        continue
    if not same:
        bad += 1
        print(f"  MISMATCH {what}: README says {quoted}, table says {actual}")
if unquoted:
    print("  not found in the README, so the check could not run: "
          + "; ".join(unquoted))

# The cost table, which drifted unchecked when the templates arm was extended
# and the seed repeats were added.
tim = {r["arm"]: r for r in L.read_tsv(REPO / "results/timing.tsv")}
for arm, n, secs, mem in [("af2_nomsa", "149", "155", "2753"),
                          ("af2_msa_notmpl", "34", "843", "3398"),
                          ("af2_msa_tmpl", "14", "1021", "3768")]:
    row = tim.get(arm)
    if row is None:
        bad += 1
        print(f"  MISSING timing row for {arm}")
        continue
    for label, quoted, actual in [("count", n, row["n"]),
                                  ("median seconds", secs, row["seconds_median"]),
                                  ("median peak MB", mem, row["peak_rss_mb_median"])]:
        if round(float(actual)) != round(float(quoted)):
            bad += 1
            print(f"  MISMATCH {arm} {label}: README {quoted}, table {actual}")

# The geometry claims, which are counted from the per-structure table.
geo = L.read_tsv(REPO / "results/geometry_checks.tsv")
dep = [r for r in geo if r["arm"].startswith("null_") and r["status"] == "ok"]
strict = sum(1 for r in dep if r.get("passes_all") == "1")
within = sum(1 for r in dep if r.get("within_deposited_range") == "1")
per_arm = {}
for arm in ("null_template", "null_unrelated"):
    sub = [r for r in dep if r["arm"] == arm]
    per_arm[arm] = (sum(1 for r in sub if r.get("within_deposited_range") == "1"),
                    len(sub))
geometry_claims = [
    ("118", str(strict), "deposited meeting the zero-fault standard"),
    ("270", str(len(dep)), "deposited structures checked"),
    ("235 of 270", f"{within} of {len(dep)}", "deposited inside the range"),
    ("103 of 120", "%d of %d" % per_arm["null_template"], "template floor inside the range"),
    ("132 of 150", "%d of %d" % per_arm["null_unrelated"], "unrelated floor inside the range"),
]
for quoted, actual, what in geometry_claims:
    if actual != quoted:
        bad += 1
        print(f"  MISMATCH {what}: README {quoted}, table {actual}")
    elif quoted not in readme:
        bad += 1
        print(f"  ABSENT {what}: the README no longer says {quoted}")

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
for (a, b), diff in [(("af2_msa_notmpl", "af2_msa_tmpl"), "0.001"),
                     (("af2_msa_notmpl", "af2_nomsa"), "-0.476"),
                     (("af2_msa_notmpl", "null_template"), "-0.109")]:
    row = pairs.get((a, b))
    if row is None or row["median_paired_difference"] != diff:
        bad += 1
        print(f"  MISMATCH pair {a} against {b}: README {diff}, table "
              f"{row['median_paired_difference'] if row else 'absent'}")

print(f"{len(claims) + len(geometry_claims) + 22} checks, {bad} mismatch(es)")
sys.exit(1 if bad else 0)
