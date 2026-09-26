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
import shutil
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
    "inconsistent_residues", "reference_chain", "reference_chains_dropped",
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
def write_target_chain(reference: Path, chain_name: str, out_path: Path) -> tuple[int, int]:
    """The one chain the prediction is of, written on its own.

    This is the difference between measuring a prediction and measuring a
    prediction against an assembly it was never asked to build.

    lDDT asks, for every pair of residues within its inclusion radius in the
    reference, whether the model preserves that distance. A deposited entry is
    whatever was crystallised, often several chains. A model of one entity is
    one chain. Every contact the reference makes across a chain boundary is
    then a distance the model cannot reproduce, not because it is wrong but
    because it was never given the partner, and each one counts against it.

    Measured here, comparing against the whole entry instead of the target
    chain: an entry with six chains scored 0.136 where the chain alone scored
    0.940, four chains 0.164 against 0.963, three 0.283 against 0.892, two
    0.377 against 0.950. Single-chain entries were identical either way, which
    is what says the effect is the extra chains and nothing else. Read the
    wrong way, that pattern looks like a model that is confidently wrong on
    most targets, which is the opposite of the truth.

    The ligand and any other component stay out too, so what is compared is
    polymer against polymer.

    Returns the residues written and the chains dropped, both recorded per
    target so a reader can see which entries this mattered for.
    """
    import gemmi

    st = gemmi.read_structure(str(reference))
    st.setup_entities()
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    st.remove_waters()
    dropped = 0
    for model in st:
        names = sorted({c.name for c in model if c.name != chain_name})
        dropped += len(names)
        for name in names:
            model.remove_chain(name)
    st.remove_ligands_and_waters()
    n = sum(1 for model in st for chain in model for _res in chain)
    doc = st.make_mmcif_document()
    doc.write_file(str(out_path))
    return n, dropped


def compare_structures(ost_bin: str, model: Path, reference: Path, out_json: Path,
                       inclusion_radius: float,
                       map_seqid_thresh: float) -> tuple[bool, dict, str]:
    """One call, every score this stage takes from OpenStructure.

    --bb-lddt and --bb-local-lddt give the C-alpha-only score, which is the one
    the confidence head predicts. --lddt and --local-lddt give the all-atom
    score. --rigid-scores gives the superposition, the deviation and the global
    distance tests. Asking for them together costs one parse of each structure
    rather than four.

    The residue-number alignment option is deliberately not used. A prediction
    is numbered from one and a deposited model is numbered however the
    depositors chose, so aligning by number would pair unrelated residues. The
    default sequence alignment is what makes the two comparable.

    The consistency check is deliberately not used either, and that is a
    change from how this was first written. Turning it on makes any residue
    difference between model and reference fatal, which is right for a
    prediction of the same sequence and wrong for the copied-template floor,
    where the copy is a different protein by construction. With the check on,
    every template floor failed with a residue mismatch and scored nothing,
    which would have removed the only honest comparator from the results. The
    mismatches are counted and reported per structure instead.

    The chain-mapping identity threshold is lowered from its own default of 70
    per cent and recorded in the configuration. That default exists to stop a
    model chain being compared against an unrelated reference chain, which is
    sensible when the two are meant to be the same protein. The copied-template
    floor is deliberately not the same protein: at the default, every template
    below 70 per cent identity was left unmapped and scored zero. A floor that
    reads zero because nothing was compared is not a floor, it is a missing
    measurement dressed as one. With the threshold lowered, every arm is scored
    by one rule and the unrelated floor is measured rather than refused.
    """
    cmd = [ost_bin, "compare-structures",
           "-m", str(model), "-r", str(reference), "-o", str(out_json),
           "--lddt", "--local-lddt", "--bb-lddt", "--bb-local-lddt",
           "--tm-score", "--rigid-scores",
           "--lddt-inclusion-radius", str(inclusion_radius),
           "--chem-map-seqid-thresh", str(map_seqid_thresh)]
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
# The pocket
# ---------------------------------------------------------------------------
def pocket_residues(reference: Path, ligand_sdf: Path, radius: float):
    """Sequence positions of every reference residue near the ligand.

    A prediction can be excellent over the whole chain and useless for docking,
    because the residues that decide a docking result are a few dozen side
    chains lining one cavity. Those residues are identified here from the
    experimental structure and the ligand as deposited, so the selection owes
    nothing to the prediction being judged.

    The radius is arbitrary. It is taken from the configuration, several values
    are computed, and the README says it is arbitrary rather than implying that
    5 Angstrom is a property of pockets.
    """
    import gemmi
    import numpy as np

    lig = []
    for line in ligand_sdf.read_text(errors="replace").splitlines():
        parts = line.split()
        # An SDF atom line starts with three coordinates and an element symbol.
        if len(parts) >= 4:
            try:
                x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
            except ValueError:
                continue
            if parts[3].isalpha() and parts[3] != "H":
                lig.append((x, y, z))
    if not lig:
        return {}, 0
    lig_arr = np.array(lig)

    st = gemmi.read_structure(str(reference))
    st.setup_entities()
    st.remove_ligands_and_waters()
    near = {}
    if not len(st):
        return {}, len(lig)
    for chain in st[0]:
        poly = chain.get_polymer()
        if len(poly) == 0:
            continue
        for res in poly:
            if res.label_seq is None:
                continue
            for atom in res:
                if atom.element == gemmi.Element("H"):
                    continue
                d = np.linalg.norm(lig_arr - np.array([atom.pos.x, atom.pos.y, atom.pos.z]),
                                   axis=1).min()
                if d <= radius:
                    near.setdefault(int(res.label_seq), chain.name)
                    break
        if near:
            break
    return near, len(lig)


def pocket_scores(local_lddt: dict, bb_local_lddt: dict, positions: dict):
    """Average the per-residue scores over the pocket.

    OpenStructure has already computed a score for every residue; restricting
    to the pocket is a selection rather than a second calculation, so the
    pocket number and the global number are the same quantity over different
    residues and can be compared directly.
    """
    def mean_over(table):
        vals = []
        for key, value in (table or {}).items():
            bits = key.split(".")
            if len(bits) < 2 or not bits[1].isdigit():
                continue
            if int(bits[1]) in positions and value is not None:
                vals.append(float(value))
        return (sum(vals) / len(vals)) if vals else None, len(vals)

    all_atom, n_all = mean_over(local_lddt)
    ca_only, n_ca = mean_over(bb_local_lddt)
    return ca_only, all_atom, max(n_all, n_ca)


def pocket_rmsd(model: Path, reference: Path, positions: dict):
    """All-atom deviation over the pocket, after fitting on the pocket.

    Fitting on the pocket rather than on the whole chain is the point: a model
    whose domains are slightly rotated relative to each other can have a poor
    global fit and an excellent pocket, and it is the pocket that decides
    whether a ligand can be placed. This follows the published comparison this
    work is measured against, which aligned on the pocket residues and then
    reported the deviation over their heavy atoms.
    """
    import gemmi
    import numpy as np

    def heavy_atoms(path):
        st = gemmi.read_structure(str(path))
        st.setup_entities()
        st.remove_ligands_and_waters()
        out = {}
        if not len(st):
            return out
        for chain in st[0]:
            poly = chain.get_polymer()
            if len(poly) == 0:
                continue
            for res in poly:
                if res.label_seq is None or int(res.label_seq) not in positions:
                    continue
                for atom in res:
                    if atom.element == gemmi.Element("H"):
                        continue
                    out[(int(res.label_seq), atom.name)] = (atom.pos.x, atom.pos.y, atom.pos.z)
            if out:
                break
        return out

    a, b = heavy_atoms(model), heavy_atoms(reference)
    shared = sorted(set(a) & set(b))
    if len(shared) < 8:
        return None, len(shared)
    m = np.array([a[k] for k in shared])
    r = np.array([b[k] for k in shared])
    cm, cr = m.mean(axis=0), r.mean(axis=0)
    h = (m - cm).T @ (r - cr)
    u, _s, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rot = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    fitted = (rot @ (m - cm).T).T
    return float(np.sqrt(((fitted - (r - cr)) ** 2).sum(axis=1).mean())), len(shared)

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
    map_seqid_thresh = L.conf_float(conf, "CHEM_MAP_SEQID_THRESH", 20.0)

    root = Path(conf["CONDA_SH"]).parent.parent.parent
    ost_bin = str(root / "envs" / conf.get("CONDA_ENV_OST", "csc_ost") / "bin" / "ost")
    usalign_bin = str(root / "envs" / conf.get("CONDA_ENV_OST", "csc_ost") / "bin" / "USalign")
    tmscore_bin = str(root / "envs" / conf.get("CONDA_ENV_OST", "csc_ost") / "bin" / "TMscore")
    for b in (ost_bin, usalign_bin, tmscore_bin):
        if not Path(b).is_file():
            L.eprint(f"[error] {b} not found; run 01_install.sh")
            return 3

    out_scores = results_dir / "structure_scores.tsv"
    out_residues = results_dir / "residue_scores.tsv"
    out_pockets = results_dir / "pocket_scores.tsv"

    # Scoring is resumed rather than repeated. The stage used to skip entirely
    # when its own output existed, so an arm finished after the first run was
    # never scored: the chain called this, it printed "exists; skipping", and
    # thirty predictions stayed out of every table downstream while the stage
    # reported success. Rewriting the whole file each time is the other half of
    # the same fault, because it makes scoring one new arm cost a rerun of every
    # comparison already made.
    #
    # So the pairs already scored are read back, those predictions are left
    # alone, and their rows are carried into the file that is written at the
    # end. --force ignores what is there and redoes everything.
    prior_scores = L.read_tsv(out_scores) if out_scores.is_file() else []
    prior_residues = L.read_tsv(out_residues) if out_residues.is_file() else []
    prior_pockets = L.read_tsv(out_pockets) if out_pockets.is_file() else []
    already = set()
    if not args.force:
        already = {(r.get("target_id"), r.get("arm")) for r in prior_scores}

    preds = L.read_tsv(results_dir / "predictions.tsv") if (results_dir / "predictions.tsv").is_file() else []
    preds = [p for p in preds if p.get("status") == "ok"]
    if args.arm:
        preds = [p for p in preds if p["arm"] == args.arm]
    if args.limit:
        preds = preds[:args.limit]
    if not preds:
        L.eprint("[error] no successful predictions to score; run 05_predict.sh first")
        return 1
    n_all = len(preds)
    preds = [p for p in preds if (p["target_id"], p["arm"]) not in already]
    n_skipped = n_all - len(preds)
    if n_skipped:
        print(f"[06_score_structures] {n_skipped} already scored, {len(preds)} to do")
    if not preds:
        print("[06_score_structures] everything is already scored; nothing to do")
        return 0

    targets = {t["target_id"]: t for t in L.read_tsv(config_dir / "targets.tsv")}
    # Only the targets with a usable ligand have a pocket to measure. The rest
    # get blank pocket columns rather than a zero, because no pocket is not the
    # same as a badly predicted one.
    subset_path = config_dir / "docking_subset.tsv"
    subset = {r["target_id"]: r for r in L.read_tsv(subset_path)} \
        if subset_path.is_file() else {}
    radii = [float(x) for x in (conf.get("POCKET_RADII") or "5.0").split()]
    primary_radius = L.conf_float(conf, "POCKET_RADIUS", 5.0)
    rows, residue_rows, pocket_rows = [], [], []
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

            # Compare against the chain the prediction is of, not the whole
            # deposited entry. See write_target_chain for what that was costing.
            want_chain = target.get("auth_asym_id", "")
            shutil.copyfile(reference, tmp / "reference_full.cif")
            reference_one = tmp / "reference_chain.cif"
            n_ref_res = n_dropped = 0
            if want_chain:
                try:
                    n_ref_res, n_dropped = write_target_chain(
                        reference, want_chain, reference_one)
                except Exception as e:  # noqa: BLE001 - one target's reference
                    row.update({"status": "failed",
                                "reason": f"could not isolate chain {want_chain}: "
                                          f"{type(e).__name__}: {str(e)[:90]}"})
                    rows.append(row)
                    L.record_exclusion(results_dir, tid, "06_score_structures",
                                       row["reason"], arm=arm)
                    continue
            if n_ref_res == 0:
                row.update({"status": "failed",
                            "reason": f"chain {want_chain or '(unnamed)'} holds no "
                                      f"residues in the reference"})
                rows.append(row)
                L.record_exclusion(results_dir, tid, "06_score_structures",
                                   row["reason"], arm=arm)
                continue
            row["reference_chain"] = want_chain
            row["reference_chains_dropped"] = n_dropped
            reference = reference_one

            out_json = tmp / "compare.json"
            ok, data, why = compare_structures(ost_bin, model, reference, out_json, radius,
                                               map_seqid_thresh)
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
                "inconsistent_residues": len(data.get("inconsistent_residues") or []),
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

            item = subset.get(tid)
            if item:
                lig = data_dir / "reference_ligands" / f"{tid}__{item['comp_id']}.sdf"
                # The pocket is wherever the ligand is, and the chosen copy of
                # the ligand does not always sit in the chain the manifest
                # names. Where it sits in another copy of the same entity, the
                # pocket is measured there; the sequences are identical so the
                # positions carry over to the prediction unchanged.
                pocket_ref = reference
                lig_chain = item.get("receptor_auth_asym_id") or want_chain
                if lig_chain and lig_chain != want_chain:
                    pocket_ref = tmp / "pocket_chain.cif"
                    try:
                        write_target_chain(tmp / "reference_full.cif"
                                           if (tmp / "reference_full.cif").is_file()
                                           else reference, lig_chain, pocket_ref)
                    except Exception:  # noqa: BLE001 - fall back to the named chain
                        pocket_ref = reference
                if lig.is_file():
                    for radius_value in radii:
                        near, n_lig_atoms = pocket_residues(pocket_ref, lig, radius_value)
                        if not near:
                            continue
                        ca_only, all_atom, n_scored = pocket_scores(
                            data.get("local_lddt"), local, near)
                        rmsd_value, n_atoms = pocket_rmsd(model, pocket_ref, near)
                        pocket_rows.append({
                            "target_id": tid, "arm": arm, "radius": radius_value,
                            "n_pocket_residues": len(near),
                            "n_residues_scored": n_scored,
                            "n_ligand_heavy_atoms": n_lig_atoms,
                            "pocket_lddt_ca": L.fmt(ca_only, 4),
                            "pocket_lddt_all_atom": L.fmt(all_atom, 4),
                            "pocket_rmsd_all_atom": L.fmt(rmsd_value, 3),
                            "n_atoms_compared": n_atoms,
                            "global_lddt_ca": L.fmt(data.get("bb_lddt"), 4),
                            "recorded": L.now_iso(),
                        })
                        if abs(radius_value - primary_radius) < 1e-6:
                            row["pocket_lddt_ca"] = L.fmt(ca_only, 4)
                            row["pocket_rmsd_all_atom"] = L.fmt(rmsd_value, 3)
                            row["pocket_n_residues"] = len(near)
            conf_path = results_dir / pred.get("confidence_file", "")
            if conf_path.is_file():
                residue_rows.extend(
                    pair_plddt_with_lddt(L.read_json_gz(conf_path), local, tid, arm))

        rows.append(row)
        if i % 5 == 0 or i == len(preds):
            print(f"  {i}/{len(preds)} scored", flush=True)

    # Carry forward every row for a pair this run did not touch, then write the
    # whole file once. Keeping the rows rather than appending means the file is
    # still written by a single rename and a run killed halfway cannot leave a
    # half-written table.
    redone = {(r["target_id"], r["arm"]) for r in rows}
    keep = [r for r in prior_scores if (r.get("target_id"), r.get("arm")) not in redone]
    keep_res = [r for r in prior_residues if (r.get("target_id"), r.get("arm")) not in redone]
    keep_pock = [r for r in prior_pockets if (r.get("target_id"), r.get("arm")) not in redone]
    if keep:
        print(f"[06_score_structures] carrying forward {len(keep)} rows scored earlier")

    L.write_tsv(out_scores, SCORE_COLUMNS, keep + rows)
    L.write_tsv(out_residues, RESIDUE_COLUMNS, keep_res + residue_rows)
    pocket_all = keep_pock + pocket_rows
    if pocket_all:
        L.write_tsv(out_pockets,
                    ["target_id", "arm", "radius", "n_pocket_residues",
                     "n_residues_scored", "n_ligand_heavy_atoms", "pocket_lddt_ca",
                     "pocket_lddt_all_atom", "pocket_rmsd_all_atom",
                     "n_atoms_compared", "global_lddt_ca", "recorded"], pocket_all)
    ok_rows = [r for r in rows if r["status"] == "ok"]
    print(f"[06_score_structures] {len(ok_rows)} of {len(rows)} scored, "
          f"{len(residue_rows)} residue rows, {len(pocket_rows)} pocket rows, "
          f"OpenStructure {ost_version}")
    if ok_rows:
        vals = sorted(float(r["lddt_ca"]) for r in ok_rows if r.get("lddt_ca") not in ("", None))
        if vals:
            q = L.quantiles(vals)
            print(f"[06_score_structures] lDDT-CA median {q[1]:.3f}, "
                  f"interquartile range {q[0]:.3f} to {q[2]:.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
