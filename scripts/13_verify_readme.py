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
import argparse
import sys
from pathlib import Path

_ap = argparse.ArgumentParser(
    description="Check every number the README quotes against its table.",
    epilog="Takes no options. --config is accepted and ignored: this reads "
           "the tables beside it rather than a configuration.")
_ap.add_argument("--config", help=argparse.SUPPRESS)
_ap.parse_args()

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
for arm, n, dist, secs, mem in [("af2_nomsa", "149", "149", "155", "2753"),
                                ("af2_msa_notmpl", "34", "30", "843", "3398"),
                                ("af2_msa_tmpl", "14", "14", "1021", "3768")]:
    row = tim.get(arm)
    if row is None:
        bad += 1
        print(f"  MISSING timing row for {arm}")
        continue
    for label, quoted, actual in [("count", n, row["n"]),
                                  ("distinct targets", dist, row["n_distinct_targets"]),
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

# The receptor superposition, which is the quantity that joins the docking arm
# to the accuracy arms and so is quoted rather than left in the table.
align = [r for r in L.read_tsv(REPO / "results/docking_alignment.tsv")
         if r["arm"] == "predicted" and r["status"] == "ok"
         and r.get("superposition_rmsd_ca") not in ("", None)]
av = sorted(float(r["superposition_rmsd_ca"]) for r in align)
aq = L.quantiles(av)
align_claims = [
    ("13", str(len(av)), "predicted receptors prepared"),
    ("1.57", f"{aq[1]:.2f}", "median receptor superposition"),
    ("1.02", f"{aq[0]:.2f}", "first quartile receptor superposition"),
    ("2.39", f"{aq[2]:.2f}", "third quartile receptor superposition"),
    ("0.61", f"{min(av):.2f}", "closest receptor superposition"),
    ("15.65", f"{max(av):.2f}", "furthest receptor superposition"),
    ("Four of the 13", "Four of the 13" if sum(1 for x in av if x > 2.0) == 4
     else f"{sum(1 for x in av if x > 2.0)} of the {len(av)}",
     "receptors above 2 Angstroms"),
]
for quoted, actual, what in align_claims:
    if actual != quoted:
        bad += 1
        print(f"  MISMATCH {what}: README {quoted}, table {actual}")
    elif quoted not in readme:
        bad += 1
        print(f"  ABSENT {what}: the README no longer says {quoted}")

# The per-target direction of the templates effect, and the multiple of the
# seed range. Both were stated from memory once and both were wrong: the
# README said seven targets up and seven down where three went down and four
# did not move, and called the alignment effect four hundred times the spread
# where it is 119 times the range.
sc = {}
for r in L.read_tsv(REPO / "results/structure_scores.tsv"):
    if r.get("status") == "ok" and r.get("lddt_ca"):
        sc.setdefault((r["target_id"], r["arm"]), []).append(
            (r.get("seed", ""), float(r["lddt_ca"])))


def primary(target, arm):
    vals = sc.get((target, arm)) or []
    pick = [v for s, v in vals if s == "20260925"] or [v for _s, v in vals]
    return pick[0] if pick else None


shared = sorted({t for (t, a) in sc if a == "af2_msa_tmpl"}
                & {t for (t, a) in sc if a == "af2_msa_notmpl"})
diffs = [primary(t, "af2_msa_tmpl") - primary(t, "af2_msa_notmpl") for t in shared]
seedrow = L.read_tsv(REPO / "results/seed_variance.tsv")[0]
seed_range = float(seedrow["range"])
direction_claims = [
    ("fourteen shared targets", str(len(shared)), "14", "targets both arms cover"),
    ("it moved the score up for", str(sum(1 for d in diffs if d > 0)), "7",
     "targets templates raised"),
    ("down for three", str(sum(1 for d in diffs if d < 0)), "3", "targets templates lowered"),
    ("remaining four", str(sum(1 for d in diffs if d == 0)), "4", "targets templates did not move"),
    ("119 times that range", f"{0.476 / seed_range:.0f}", "119",
     "the alignment effect as a multiple of the seed range"),
    ("a quarter of the range", f"{0.001 / seed_range:.2f}", "0.25",
     "the templates effect as a fraction of the seed range"),
]
for phrase, actual, expected, what in direction_claims:
    if actual != expected:
        bad += 1
        print(f"  MISMATCH {what}: README implies {expected}, table says {actual}")
    elif phrase not in readme:
        bad += 1
        print(f"  ABSENT {what}: the README no longer says {phrase!r}")

# The two checks the README says flag no prediction. Stated as never failing
# at all once, which ten deposited structures contradict.
over = {"ca_chirality_errors": [], "cis_nonpro": []}
thr = {r["check"]: r for r in L.read_tsv(REPO / "results/geometry_thresholds.tsv")}
for r in L.read_tsv(REPO / "results/geometry_checks.tsv"):
    for c in over:
        lim, got = thr.get(c, {}).get("threshold"), r.get(c)
        if lim not in (None, "") and got not in (None, "") and float(got) > float(lim):
            over[c].append(r)
pred_flagged = sum(1 for c in over for r in over[c] if r["arm"].startswith("af2_"))
quiet_claims = [
    ("prediction at all.", str(pred_flagged), "0",
     "predictions flagged on chirality or cis"),
    ("The ten structures they do flag are all deposited",
     str(sum(len(v) for v in over.values())), "10",
     "structures flagged in total"),
    ("are all deposited, two on", str(len(over["ca_chirality_errors"])), "2",
     "structures over the chirality threshold"),
    ("eight on cis peptides", str(len(over["cis_nonpro"])), "8",
     "structures over the cis threshold"),
]
for phrase, actual, expected, what in quiet_claims:
    if actual != expected:
        bad += 1
        print(f"  MISMATCH {what}: README implies {expected}, table says {actual}")
    elif phrase not in readme:
        bad += 1
        print(f"  ABSENT {what}: the README no longer says {phrase!r}")

# The application arm, which had no check at all while this file's own
# docstring said every number the README quotes is checked against its table.
pedv_dir = REPO / "results" / "pedv"
pedv_claims = []
if (pedv_dir / "deposited_pairs.tsv").is_file():
    dep = [float(r["domain_displacement"])
           for r in L.read_tsv(pedv_dir / "deposited_pairs.tsv")
           if r.get("domain_displacement")]
    pro = L.read_tsv(pedv_dir / "protomer_pairs.tsv")         if (pedv_dir / "protomer_pairs.tsv").is_file() else []
    pred = {r["structure_b"]: r["domain_displacement"]
            for r in L.read_tsv(pedv_dir / "prediction_pairs.tsv")
            if r.get("structure_a") == "prediction/af2_msa_notmpl"}
    con = L.read_tsv(pedv_dir / "construct.tsv")[0]
    by_entry = {}
    for r in pro:
        by_entry.setdefault(r["entry"], []).append(r["domain_displacement"])
    pedv_claims = [
        ("55.4 Angstroms", f"{max(dep):.1f} Angstroms" if dep else "none",
         "the largest displacement between deposited entries"),
        ("plus 98 residues", f"plus {con['body_residues']} residues",
         "the body carried by the construct"),
        ("299 positions", f"{con['construct_length']} positions",
         "the length of the construct"),
        ("reference entry is 7W6M", f"reference entry is {con['source_entry']}",
         "the entry the construct was cut from"),
    ]
    for entry, vals in sorted(by_entry.items()):
        shown = " and ".join(f"{float(v):g}" for v in vals) + " Angstroms"
        pedv_claims.append((f"| {entry} | {shown} |", f"| {entry} | {shown} |",
                            f"the protomer displacements in {entry}"))
    for entry, val in sorted(pred.items()):
        pedv_claims.append((f"| {entry} | {float(val):.1f} Angstroms |",
                            f"| {entry} | {float(val):.1f} Angstroms |",
                            f"the prediction against {entry}"))
for quoted, actual, what in pedv_claims:
    if actual != quoted:
        bad += 1
        print(f"  MISMATCH {what}: README {quoted}, table {actual}")
    elif quoted not in readme:
        bad += 1
        print(f"  ABSENT {what}: the README does not say {quoted!r}")

# What spanning the length range cost, which the README used to state as
# "three times as many targets" with nothing behind it. It is 1.43.
sel_path = REPO / "results/arm_selection_cost.tsv"
sel_claims = []
if sel_path.is_file():
    sel = {r["arm"]: r for r in L.read_tsv(sel_path)}
    for arm, label in (("af2_msa_notmpl", "alignment arm"),
                       ("af2_msa_tmpl", "templates arm")):
        r = sel.get(arm)
        if not r:
            continue
        sel_claims += [
            (r["shortest_first_n_targets"], r["shortest_first_n_targets"],
             f"targets the shortest first would buy, {label}"),
            (r["ratio"], r["ratio"], f"the ratio, {label}"),
            (r["shortest_first_length_max"], r["shortest_first_length_max"],
             f"the longest target the shortest first reaches, {label}"),
        ]
for quoted, actual, what in sel_claims:
    if actual != quoted:
        bad += 1
        print(f"  MISMATCH {what}: README {quoted}, table {actual}")
    elif quoted not in readme:
        bad += 1
        print(f"  ABSENT {what}: the README does not say {quoted!r}")

# The machine settings the cost section quotes, which used to come only from
# project.conf and so from no file in results/ or logs/.
ver_path = REPO / "results/environment/versions.tsv"
env_claims = []
if ver_path.is_file():
    ver = {r["tool"]: r["version"] for r in L.read_tsv(ver_path)}
    env_claims = [
        (f"{ver.get('threads', '?')} threads", f"{ver.get('threads', '?')} threads",
         "the thread count"),
        (f"{ver.get('memory_budget_gb', '?')} GB memory budget",
         f"{ver.get('memory_budget_gb', '?')} GB memory budget",
         "the memory budget"),
    ]
for quoted, actual, what in env_claims:
    if actual != quoted:
        bad += 1
        print(f"  MISMATCH {what}: README {quoted}, table {actual}")
    elif quoted not in readme:
        bad += 1
        print(f"  ABSENT {what}: the README does not say {quoted!r}")

# What the configure stage projects, which the reproduction section quotes.
proj_path = REPO / "results/environment/projection.tsv"
proj_claims = []
if proj_path.is_file():
    proj = {r["key"]: r["value"] for r in L.read_tsv(proj_path)}
    proj_claims = [
        (f"{proj.get('disk_total', '?')} MB of disk",
         f"{proj.get('disk_total', '?')} MB of disk", "the projected disk"),
        (f"{proj.get('disk_peak_results', '?')} MB the results",
         f"{proj.get('disk_peak_results', '?')} MB the results",
         "the disk the results occupy"),
        (f"{proj.get('wall_clock', '?')} hours over the three arms",
         f"{proj.get('wall_clock', '?')} hours over the three arms",
         "the projected wall clock"),
    ]
for quoted, actual, what in proj_claims:
    if actual != quoted:
        bad += 1
        print(f"  MISMATCH {what}: README {quoted}, table {actual}")
    elif quoted not in readme:
        bad += 1
        print(f"  ABSENT {what}: the README does not say {quoted!r}")

# the paired comparisons
for (a, b), diff in [(("af2_msa_notmpl", "af2_msa_tmpl"), "0.001"),
                     (("af2_msa_notmpl", "af2_nomsa"), "-0.476"),
                     (("af2_msa_notmpl", "null_template"), "-0.109")]:
    row = pairs.get((a, b))
    if row is None or row["median_paired_difference"] != diff:
        bad += 1
        print(f"  MISMATCH pair {a} against {b}: README {diff}, table "
              f"{row['median_paired_difference'] if row else 'absent'}")

n_checks = (len(claims) + len(geometry_claims) + len(align_claims)
            + len(direction_claims) + len(quiet_claims) + len(pedv_claims) + len(sel_claims) + len(env_claims) + len(proj_claims)
            + 12 + 6 + 4 + 3)   # timing cells, sweep bands, docking rows, pairs
print(f"{n_checks} checks, {bad} mismatch(es)")
sys.exit(1 if bad else 0)
