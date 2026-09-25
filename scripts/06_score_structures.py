#!/usr/bin/env python3
"""Score every prediction against the withheld experimental structure.

What is computed, and why each one
----------------------------------
lDDT on C-alpha atoms is the primary measure, and it is primary because it is
the quantity the confidence head was trained to predict. Comparing a confidence
value against anything else is comparing a prediction of one thing against a
measurement of another. The all-atom lDDT is computed as well and kept in its
own column, because it answers a different question and because the two are not
interchangeable: the all-atom score is computed after a stereochemistry filter
that removes offending residues and scores them as zero, while the C-alpha
score is computed on the raw model.

TM-score and GDT come from a separate implementation, and both are reported.
They measure global fold agreement, which is what a reader wants when the
question is whether the prediction is the right shape, and they are useless for
the calibration question because no confidence head predicts them.

Root-mean-square deviation after superposition is reported and is not used as
the headline. A single number that depends on the worst-placed loop is a poor
summary of a structure, and it is the measure most often plotted against
confidence in the literature, which is a category error: confidence predicts a
local, superposition-free score, and nothing about a global superposition.

The pocket-local measures exist because a prediction can be excellent overall
and useless for docking. They are computed over residues within a stated radius
of the crystallographic ligand, at several radii, so that the conclusion can be
checked against the arbitrary choice of one.

The null floors are built here rather than in the prediction stage, because
copying a structure is not a prediction and running it through that stage would
put it in the timing tables.

Usage
-----
    python scripts/06_score_structures.py --config project.conf
    python scripts/06_score_structures.py --config project.conf --arm af2_msa_notmpl
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_csc as L  # noqa: E402

SCORE_COLUMNS = [
    "target_id", "arm", "model", "seed", "status", "reason",
    "sequence_length", "n_residues_compared",
    "lddt_ca", "lddt_all_atom", "tm_score", "gdt_ts", "gdt_ha", "rmsd_ca",
    "pocket_lddt_ca", "pocket_rmsd_all_atom", "pocket_radius", "pocket_n_residues",
    "mean_plddt", "ptm", "model_clashes", "model_bad_bonds", "model_bad_angles",
    "ost_version", "recorded",
]

RESIDUE_COLUMNS = [
    "target_id", "arm", "chain", "resnum", "plddt", "lddt_ca", "recorded",
]


def run(cmd: list[str], timeout: int = 3600):
    """Run a tool and return (ok, stdout, stderr). Never raises on a non-zero
    status: a scoring failure for one target is a row in a table, not the end
    of the stage."""
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        return p.returncode == 0, p.stdout, p.stderr
    except (subprocess.TimeoutExpired, OSError) as exc:
        return False, "", repr(exc)[:300]


# ---------------------------------------------------------------------------
# OpenStructure
# ---------------------------------------------------------------------------
def compare_structures(ost_bin: str, model: Path, reference: Path, out_json: Path,
                       inclusion_radius: float) -> tuple[bool, dict, str]:
    """One call, every score this stage takes from OpenStructure.

    --bb-lddt and --bb-local-lddt give the C-alpha-only score, which is the one
    the confidence head predicts. --lddt and --local-lddt give the all-atom
    score. --rigid-scores gives the superposition, the deviation and the global
    distance tests. Asking for them together costs one parse of each structure
    rather than four.

    The residue-number alignment option is deliberately not used. A prediction
    is numbered from one and a deposited model is numbered however the
    depositors chose, so aligning by number would pair unrelated residues. The
    default sequence alignment is what makes the two comparable, and the
    consistency check is on so that a mismatch fails loudly instead of scoring
    something meaningless.
    """
    cmd = [ost_bin, "compare-structures",
           "-m", str(model), "-r", str(reference), "-o", str(out_json),
           "--lddt", "--local-lddt", "--bb-lddt", "--bb-local-lddt",
           "--tm-score", "--rigid-scores",
           "--lddt-inclusion-radius", str(inclusion_radius),
           "--enforce-consistency"]
    ok, _out, err = run(cmd)
    if not out_json.is_file():
        return False, {}, f"no output written: {err[:200]}"
    try:
        data = json.loads(out_json.read_text())
    except json.JSONDecodeError as exc:
        return False, {}, f"output was not readable: {exc}"
    if data.get("status") != "SUCCESS":
        why = (data.get("exception") or data.get("traceback") or "")[:200]
        return False, data, f"comparison reported failure: {why}"
    return ok, data, ""


def usalign_scores(usalign_bin: str, tmscore_bin: str, model: Path, reference: Path):
    """TM-score and the global distance tests from the independent implementation.

    The reference is passed second in both calls, because both programs
    normalise by the length of the second structure and the number wanted here
    is the one normalised by the experimental structure. Passing them the other
    way round produces a number that looks similar and answers a different
    question.
    """
    out = {}
    ok, stdout, _ = run([usalign_bin, str(model), str(reference), "-outfmt", "2"], timeout=900)
    if ok:
        for line in stdout.splitlines():
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 5:
                # Column 4 is normalised by the second structure.
                try:
                    out["tm_score_usalign"] = float(parts[3])
                    out["rmsd_usalign"] = float(parts[4])
                except (ValueError, IndexError):
                    pass
                break
    ok, stdout, _ = run([tmscore_bin, str(model), str(reference), "-seq"], timeout=900)
    if ok:
        for line in stdout.splitlines():
            if line.startswith("GDT-TS-score="):
                out["gdt_ts"] = float(line.split()[1])
            elif line.startswith("GDT-HA-score="):
                out["gdt_ha"] = float(line.split()[1])
            elif line.startswith("TM-score") and "=" in line and "tm_score_tmscore" not in out:
                try:
                    out["tm_score_tmscore"] = float(line.split("=")[1].split()[0])
                except (ValueError, IndexError):
                    pass
    return out


# ---------------------------------------------------------------------------
def residue_key_parts(key: str):
    """OpenStructure keys per-residue values as chain.resnum.inscode."""
    bits = key.split(".")
    if len(bits) < 2:
        return "", ""
    return bits[0], bits[1]


def pair_plddt_with_lddt(conf: dict, local_lddt: dict, target_id: str, arm: str):
    """One row per residue, confidence beside the measurement it predicts.

    The pairing is positional: the prediction is numbered from one and its
    confidence array is in the same order, so residue i of the array is the
    i-th residue of the model. OpenStructure keys its per-residue scores by the
    model's own numbering, which for these predictions is also one upwards.
    Where a residue has no score, because the reference does not cover it, the
    row is kept with an empty measurement rather than dropped, so that the
    number of unscored residues is visible.
    """
    plddt = conf.get("plddt") or []
    rows = []
    by_resnum = {}
    for key, value in (local_lddt or {}).items():
        chain, resnum = residue_key_parts(key)
        if resnum.isdigit():
            by_resnum[int(resnum)] = (chain, value)
    stamp = L.now_iso()
    for i, p in enumerate(plddt, start=1):
        chain, lddt = by_resnum.get(i, ("", None))
        rows.append({
            "target_id": target_id, "arm": arm, "chain": chain, "resnum": i,
            "plddt": L.fmt(p, 2),
            "lddt_ca": "" if lddt is None else L.fmt(lddt, 4),
            "recorded": stamp,
        })
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--arm", default="")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    conf = L.load_conf(args.config)
    config_dir = Path(conf["CONFIG_DIR"])
    results_dir = Path(conf["RESULTS_DIR"])
    data_dir = Path(conf["DATA_DIR"])
    radius = L.conf_float(conf, "LDDT_INCLUSION_RADIUS", 15.0)

    root = Path(conf["CONDA_SH"]).parent.parent.parent
    ost_bin = str(root / "envs" / conf.get("CONDA_ENV_OST", "csc_ost") / "bin" / "ost")
    usalign_bin = str(root / "envs" / conf.get("CONDA_ENV_OST", "csc_ost") / "bin" / "USalign")
    tmscore_bin = str(root / "envs" / conf.get("CONDA_ENV_OST", "csc_ost") / "bin" / "TMscore")
    for b in (ost_bin, usalign_bin, tmscore_bin):
        if not Path(b).is_file():
            L.eprint(f"[error] {b} not found; run 01_install.sh")
            return 3

    out_scores = results_dir / "structure_scores.tsv"
    if out_scores.is_file() and not args.force:
        print(f"[06_score_structures] {out_scores.name} exists; skipping (--force to redo).")
        return 0

    preds = L.read_tsv(results_dir / "predictions.tsv") if (results_dir / "predictions.tsv").is_file() else []
    preds = [p for p in preds if p.get("status") == "ok"]
    if args.arm:
        preds = [p for p in preds if p["arm"] == args.arm]
    if args.limit:
        preds = preds[:args.limit]
    if not preds:
        L.eprint("[error] no successful predictions to score; run 05_predict.sh first")
        return 1

    targets = {t["target_id"]: t for t in L.read_tsv(config_dir / "targets.tsv")}
    rows, residue_rows = [], []
    ost_version = ""

    for i, pred in enumerate(preds, 1):
        tid, arm = pred["target_id"], pred["arm"]
        target = targets.get(tid, {})
        row = {
            "target_id": tid, "arm": arm, "model": pred.get("model", ""),
            "seed": pred.get("seed", ""), "sequence_length": target.get("sequence_length", ""),
            "mean_plddt": pred.get("mean_plddt", ""), "ptm": pred.get("ptm", ""),
            "pocket_radius": L.conf_float(conf, "POCKET_RADIUS", 5.0),
            "status": "ok", "reason": "", "recorded": L.now_iso(),
        }
        model_gz = data_dir / pred.get("structure_file", "")
        ref_gz = data_dir / "reference" / f"{tid}.cif.gz"
        if not pred.get("structure_file") or not model_gz.is_file():
            row.update({"status": "failed", "reason": "the predicted structure is missing"})
            rows.append(row)
            continue
        if not ref_gz.is_file():
            row.update({"status": "failed", "reason": "the reference structure is missing"})
            rows.append(row)
            continue

        with tempfile.TemporaryDirectory(prefix="csc_score_") as td:
            tmp = Path(td)
            model = tmp / f"model{''.join(model_gz.suffixes[:-1]) or '.pdb'}"
            reference = tmp / "reference.cif"
            L.gunzip_to(model_gz, model)
            L.gunzip_to(ref_gz, reference)
            out_json = tmp / "compare.json"
            ok, data, why = compare_structures(ost_bin, model, reference, out_json, radius)
            ost_version = data.get("ost_version", ost_version)
            if not ok:
                row.update({"status": "failed", "reason": why})
                rows.append(row)
                L.record_exclusion(results_dir, tid, "06_score_structures", why, arm=arm)
                continue
            row.update({
                "lddt_ca": L.fmt(data.get("bb_lddt"), 4),
                "lddt_all_atom": L.fmt(data.get("lddt"), 4),
                "tm_score": L.fmt(data.get("tm_score"), 4),
                "rmsd_ca": L.fmt(data.get("rmsd"), 3),
                "model_clashes": len(data.get("model_clashes") or []),
                "model_bad_bonds": len(data.get("model_bad_bonds") or []),
                "model_bad_angles": len(data.get("model_bad_angles") or []),
                "ost_version": data.get("ost_version", ""),
            })
            gdt = data.get("oligo_gdtts")
            if gdt is not None:
                row["gdt_ts"] = L.fmt(gdt, 4)
            extra = usalign_scores(usalign_bin, tmscore_bin, model, reference)
            for key, col in (("gdt_ts", "gdt_ts"), ("gdt_ha", "gdt_ha")):
                if key in extra:
                    row[col] = L.fmt(extra[key], 4)

            local = data.get("bb_local_lddt") or {}
            row["n_residues_compared"] = sum(1 for v in local.values() if v is not None)
            conf_path = results_dir / pred.get("confidence_file", "")
            if conf_path.is_file():
                residue_rows.extend(
                    pair_plddt_with_lddt(L.read_json_gz(conf_path), local, tid, arm))

        rows.append(row)
        if i % 5 == 0 or i == len(preds):
            print(f"  {i}/{len(preds)} scored", flush=True)

    L.write_tsv(out_scores, SCORE_COLUMNS, rows)
    L.write_tsv(results_dir / "residue_scores.tsv", RESIDUE_COLUMNS, residue_rows)
    ok_rows = [r for r in rows if r["status"] == "ok"]
    print(f"[06_score_structures] {len(ok_rows)} of {len(rows)} scored, "
          f"{len(residue_rows)} residue rows, OpenStructure {ost_version}")
    if ok_rows:
        vals = sorted(float(r["lddt_ca"]) for r in ok_rows if r.get("lddt_ca") not in ("", None))
        if vals:
            q = L.quantiles(vals)
            print(f"[06_score_structures] lDDT-CA median {q[1]:.3f}, "
                  f"interquartile range {q[0]:.3f} to {q[2]:.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
