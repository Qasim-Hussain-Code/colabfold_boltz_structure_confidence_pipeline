#!/usr/bin/env python3
"""
=============================================================================
 05a_reconcile_predictions.py - the prediction table against the files
=============================================================================
 A prediction leaves two traces: a confidence file under results/confidence/
 and a row in results/predictions.tsv. The row is the record. Every later
 stage reads the table, so a target with a file and no row is absent from the
 analysis while looking, on disk, as though it were done.

 That is not hypothetical. Two ways for the two traces to disagree have
 already happened here:

 1. The prediction stage uses the confidence file as its skip marker. The
    memory calibration writes confidence files for the targets it times, so a
    target used for calibration was marked done before the arm ever ran. Its
    file existed, its row did not, and the full run reported success having
    skipped it.

 2. A run killed between writing the confidence file and appending the row
    leaves the same disagreement, and a resume then skips that target
    permanently.

 Both are silent. The count at the end of the stage says "already present"
 and a reader has no reason to doubt it. This stage makes the disagreement
 loud, and with --repair moves the unaccounted files aside so the next run of
 the arm predicts the target properly rather than inheriting a file whose
 provenance nobody can state.

 Files are moved, not deleted. Something produced them, and until it is clear
 what, throwing them away would destroy the evidence.

 What it reports:
     orphan_artifact   a confidence file or structure with no row in the table
     missing_artifact  a row that claims success with no confidence file
     not_in_set        artifacts for a target that is not in config/targets.tsv
     length_mismatch   a row whose residue count disagrees with its own file

 Usage:
     python scripts/05a_reconcile_predictions.py
     python scripts/05a_reconcile_predictions.py --repair
     python scripts/05a_reconcile_predictions.py --arm af2_nomsa

 Options:
     --config PATH   project.conf (default: alongside this script's repository)
     --arm NAME      restrict to one arm
     --repair        move unaccounted artifacts into a quarantine directory
     --quiet         report counts only
=============================================================================
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_csc as L  # noqa: E402

PROBLEM_COLUMNS = [
    "target_id", "arm", "problem", "detail", "confidence_file",
    "structure_file", "action", "recorded",
]


def parse_args(argv=None):
    p = argparse.ArgumentParser(add_help=True, description=__doc__.strip()[:200])
    here = Path(__file__).resolve().parent.parent
    p.add_argument("--config", default=str(here / "project.conf"))
    p.add_argument("--arm", default="")
    p.add_argument("--repair", action="store_true")
    p.add_argument("--quiet", action="store_true")
    return p.parse_args(argv)


def confidence_key(path: Path) -> tuple[str, str] | None:
    """Split TARGET__ARM.json.gz back into its two parts.

    The seed-variance runs name themselves TARGET__seedN__ARM, and those are
    deliberately outside the main table, so they are recognised and left
    alone rather than reported as orphans.
    """
    stem = path.name
    for suffix in (".json.gz", ".json"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    parts = stem.split("__")
    if len(parts) < 2:
        return None
    return parts[0], parts[-1]


def main(argv=None) -> int:
    args = parse_args(argv)
    conf = L.load_conf(args.config)
    paths = L.repo_paths(conf)
    results, data, config_dir = paths["RESULTS_DIR"], paths["DATA_DIR"], paths["CONFIG_DIR"]

    targets_file = config_dir / "targets.tsv"
    if not targets_file.is_file():
        L.eprint("[05a] config/targets.tsv is missing; build the set first")
        return 3
    in_set = {r["target_id"] for r in L.read_tsv(targets_file)}

    pred_file = results / "predictions.tsv"
    rows = L.read_tsv(pred_file) if pred_file.is_file() else []

    # Only the arms that run inference. The null floors are copied structures
    # rather than predictions and keep no confidence file, so including them
    # would report every one of them as missing an artifact.
    def is_inference(arm: str) -> bool:
        return arm.startswith("af2_") or arm.startswith("boltz")

    recorded: dict[tuple[str, str], dict] = {}
    for r in rows:
        if not is_inference(r.get("arm", "")):
            continue
        recorded[(r["target_id"], r["arm"])] = r

    conf_dir = results / "confidence"
    seen_files: dict[tuple[str, str], Path] = {}
    seed_variance = 0
    for f in sorted(conf_dir.glob("*.json.gz")) if conf_dir.is_dir() else []:
        key = confidence_key(f)
        if key is None:
            continue
        if "__seed" in f.name:
            seed_variance += 1
            continue
        seen_files[key] = f

    def structure_for(target: str, arm: str) -> Path | None:
        for ext in (".pdb.gz", ".cif.gz"):
            p = data / "predictions" / arm / f"{target}{ext}"
            if p.is_file():
                return p
        return None

    problems: list[dict] = []
    keys = set(seen_files) | set(recorded)
    if args.arm:
        keys = {k for k in keys if k[1] == args.arm}

    quarantine = data / "quarantine"
    now = L.now_iso()

    def quarantine_it(p: Path, why: str) -> str:
        """Move a file aside, keeping the arm directory it came from."""
        if not args.repair:
            return "reported only"
        dest = quarantine / why / p.parent.name / p.name
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(p), str(dest))
        return f"moved to {dest.relative_to(data)}"

    for target, arm in sorted(keys):
        cfile = seen_files.get((target, arm))
        sfile = structure_for(target, arm)
        row = recorded.get((target, arm))

        if target not in in_set:
            action = []
            if cfile:
                action.append(quarantine_it(cfile, "not_in_set"))
            if sfile:
                action.append(quarantine_it(sfile, "not_in_set"))
            problems.append(dict(
                target_id=target, arm=arm, problem="not_in_set",
                detail="artifacts exist for a target that config/targets.tsv "
                       "does not list; the memory calibration is the usual source",
                confidence_file=cfile.name if cfile else "",
                structure_file=sfile.name if sfile else "",
                action="; ".join(action), recorded=now))
            continue

        if row is None and cfile is not None:
            action = [quarantine_it(cfile, "orphan_artifact")]
            if sfile:
                action.append(quarantine_it(sfile, "orphan_artifact"))
            problems.append(dict(
                target_id=target, arm=arm, problem="orphan_artifact",
                detail="a confidence file with no row in predictions.tsv, so the "
                       "arm skips the target and the analysis never sees it",
                confidence_file=cfile.name,
                structure_file=sfile.name if sfile else "",
                action="; ".join(action), recorded=now))
            continue

        if row is not None and cfile is None and row.get("status") == "ok":
            problems.append(dict(
                target_id=target, arm=arm, problem="missing_artifact",
                detail="the table records a successful prediction whose "
                       "confidence file is not on disk",
                confidence_file="", structure_file=sfile.name if sfile else "",
                action="reported only; nothing to move", recorded=now))
            continue

        if row is not None and cfile is not None:
            n_scored = row.get("n_residues_scored") or ""
            if n_scored.isdigit():
                try:
                    plddt = L.read_json_gz(cfile).get("plddt") or []
                except OSError:
                    plddt = []
                if plddt and len(plddt) != int(n_scored):
                    problems.append(dict(
                        target_id=target, arm=arm, problem="length_mismatch",
                        detail=f"the row says {n_scored} residues scored and the "
                               f"file holds {len(plddt)}",
                        confidence_file=cfile.name,
                        structure_file=sfile.name if sfile else "",
                        action="reported only", recorded=now))

    # Everything moved aside by this stage on any earlier run, read back off
    # the quarantine directory. The table is rewritten each run, so without
    # this a clean rerun would erase the only record that a file was ever
    # removed from results/. The directory is the evidence; this reads it
    # rather than keeping a second list that could disagree with it.
    history: list[dict] = []
    for reason_dir in sorted(quarantine.glob("*")) if quarantine.is_dir() else []:
        if not reason_dir.is_dir():
            continue
        for f in sorted(reason_dir.rglob("*.gz")):
            key = confidence_key(f) if f.parent.name == "confidence" else None
            if key is not None:
                target, arm = key
            else:
                target, arm = f.name.split(".")[0], f.parent.name
            history.append(dict(
                target_id=target, arm=arm, problem="quarantined",
                detail=f"moved out of results by an earlier run of this stage "
                       f"as {reason_dir.name}",
                confidence_file=f.name if f.parent.name == "confidence" else "",
                structure_file="" if f.parent.name == "confidence" else f.name,
                action=f"held at {f.relative_to(data)}", recorded=""))

    out = results / "reconciliation.tsv"
    L.write_tsv(out, PROBLEM_COLUMNS, problems + history)

    kinds: dict[str, int] = {}
    for p in problems:
        kinds[p["problem"]] = kinds.get(p["problem"], 0) + 1

    if not args.quiet:
        print(f"[05a] {len(recorded)} recorded predictions, {len(seen_files)} "
              f"confidence files, {seed_variance} from the seed-variance runs")
        if not problems:
            print("[05a] the table and the files agree")
        for kind, n in sorted(kinds.items()):
            print(f"[05a] {kind}: {n}")
        for p in problems:
            print(f"    {p['target_id']} {p['arm']}: {p['problem']}; {p['action']}")
        if history:
            print(f"[05a] {len(history)} file(s) held in quarantine from earlier "
                  f"runs, listed in the table")
        print(f"[05a] wrote {out}")
        if problems and not args.repair:
            print("[05a] nothing was moved; rerun with --repair to quarantine "
                  "the unaccounted files so the arm predicts those targets again")

    # A disagreement is a finding, not a crash. The exit status says whether
    # one was found so that run_all.sh can stop rather than carry on into an
    # analysis whose denominator is quietly short.
    return 1 if any(k != "not_in_set" for k in kinds) else 0


if __name__ == "__main__":
    sys.exit(main())
