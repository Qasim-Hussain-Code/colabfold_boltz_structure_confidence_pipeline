#!/usr/bin/env python3
"""Prepare receptors for the docking pipeline, and collect what it returns.

Called by scripts/09_dock_into_predictions.sh. Two jobs, and the first is the
one that decides whether the comparison means anything.

Putting the prediction in the right frame
-----------------------------------------
The docking pipeline measures pose accuracy in place, without superposing, and
centres its search box on the reference ligand's own coordinates. That is the
right way to measure a pose. It also means the receptor has to be in the same
frame as that ligand. A predicted structure is not: it comes out wherever the
model put it.

So every predicted receptor is superposed onto its experimental structure
before it is handed over. The correspondence is exact rather than guessed:
the prediction was made from the deposited sequence, so residue i of the
prediction is position i of that sequence, and the deposited model records the
same position for each of its residues. Matching on that gives a residue pair
list with no alignment step and no chance of an off-by-one shift. Residues the
experiment did not resolve simply have no partner.

The superposition is the standard least-squares fit on the matched alpha
carbons, and its deviation is recorded per target. A large deviation is not a
reason to drop the target; it is the measurement this stage exists to make, and
it is what the docking result should be read against.

What the crystal arm is for
---------------------------
The same ligand docked into the experimental receptor of the same target, with
the same protocol. Without it the predicted number cannot be read: this set is
not the set the docking pipeline published on, so its published figure is not a
fair comparator for these targets. This one is.
"""

from __future__ import annotations

import argparse
import math
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_csc as L  # noqa: E402

MANIFEST_COLUMNS = ["dataset", "complex_id", "pdb_id", "ccd_id", "protein_pdb",
                    "ligand_sdf", "ligands_sdf", "start_conf_sdf"]

ALIGN_COLUMNS = ["target_id", "arm", "status", "reason", "n_matched_ca",
                 "superposition_rmsd_ca", "n_model_residues", "n_reference_residues",
                 "receptor_chain", "chains_dropped", "components_dropped",
                 "recorded"]

HANDOFF_COLUMNS = ["arm", "method", "n_attempted", "n_scored", "n_assessed",
                   "n_rmsd_within_2a", "rate_rmsd_within_2a",
                   "n_success", "rate_success", "ci_low", "ci_high",
                   "median_top1_rmsd", "note", "recorded"]


def kabsch(mob, ref):
    """Rotation and translation putting mob onto ref, least squares.

    Returns (rotation, centre of mob, centre of ref, deviation after fitting).
    """
    import numpy as np

    mob = np.asarray(mob, dtype=float)
    ref = np.asarray(ref, dtype=float)
    cm, cr = mob.mean(axis=0), ref.mean(axis=0)
    p, q = mob - cm, ref - cr
    h = p.T @ q
    u, _s, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rot = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    fitted = (rot @ p.T).T
    rmsd = float(np.sqrt(((fitted - q) ** 2).sum(axis=1).mean()))
    return rot, cm, cr, rmsd


def polymer_residues(structure, chain_name: str | None = None):
    """Residues of the first protein chain, keyed by position in the deposited
    sequence. Waters and other non-polymer components are left out.

    Returns the residues and the name of the chain they came from, because the
    receptor writer needs to keep that same chain and no other.
    """
    out = {}
    picked = ""
    if not len(structure):
        return out, picked
    for chain in structure[0]:
        if chain_name and chain.name != chain_name:
            continue
        poly = chain.get_polymer()
        if len(poly) == 0:
            continue
        # A deposited entry carries label_seq, the position in the sequence the
        # entry was built from. A prediction written as PDB carries none: the
        # format has no such field, so every label_seq reads as absent and
        # keying on it returned an empty set for every predicted receptor and
        # matched zero alpha carbons against the reference.
        #
        # The fallback is the position along the chain, which is the same
        # quantity for a prediction. The model was given the deposited sequence
        # and returns one residue per position of it in order, so residue i of
        # the prediction is position i of that sequence, which is what the
        # reference records as label_seq i. That is the correspondence this
        # module's header describes; it was simply never available through the
        # field it was being read from.
        residues = list(poly)
        has_label = all(r.label_seq is not None for r in residues)
        for index, res in enumerate(residues, start=1):
            key = int(res.label_seq) if has_label else index
            ca = res.find_atom("CA", "*")
            if ca is not None:
                out[key] = (res, (ca.pos.x, ca.pos.y, ca.pos.z))
        if out:
            picked = chain.name
            break
    return out, picked


def write_receptor_pdb(structure, path: Path, keep_chain: str = "") -> tuple[int, int, int]:
    """Write a receptor the docking pipeline can read.

    That pipeline reads fixed-column records and cannot carry a chain
    identifier longer than one character, so the chain is renamed. It strips
    waters and additives itself; nothing is removed here beyond what is needed
    to make a well formed file, so the two arms are cleaned by the same code.

    keep_chain is the reason this function takes an argument at all. The model
    predicts one chain, so a predicted receptor is one chain. A deposited entry
    is whatever was crystallised: other protein chains, cofactors, metals. If
    the crystal arm hands over all of that while the predicted arm hands over a
    lone chain, the difference between the two numbers is not the difference
    between predicted and experimental coordinates. It is that plus everything
    the prediction never had. This stage exists to measure the first, so both
    arms are cut to the same chain and what was dropped is counted rather than
    quietly discarded.

    The cost of that choice is real and runs the other way: a pocket formed
    between two chains is more open in a lone chain than it is in the entry, so
    the crystal arm here is a harder target than docking into the full assembly
    would be. The counts this returns are what the README needs to say so.

    Returns the residues written, the chains dropped, and the non-water
    components dropped along with them.
    """
    import gemmi

    st = structure.clone()
    st.setup_entities()
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    st.remove_waters()

    dropped_chains = dropped_components = 0
    if keep_chain:
        for model in st:
            # The names are collected before anything is removed. Removing a
            # chain while iterating over the model skips the next one, which
            # left three of six chains standing in an entry here and produced a
            # 438-residue receptor for a 147-residue target.
            to_drop = sorted({c.name for c in model if c.name != keep_chain})
            dropped_chains += len(to_drop)
            for chain in model:
                if chain.name in to_drop:
                    dropped_components += sum(
                        1 for res in chain if res.het_flag == "H")
            for name in to_drop:
                model.remove_chain(name)
        # A component that sits in the kept chain rather than in one of its
        # own, which is how many entries record a cofactor.
        for model in st:
            for chain in model:
                for res in list(chain):
                    if res.het_flag == "H":
                        dropped_components += 1
        st.remove_ligands_and_waters()

    for model in st:
        for chain in model:
            if len(chain.name) > 1:
                chain.name = chain.name[0]
    n = sum(1 for model in st for chain in model for _res in chain)
    st.write_pdb(str(path))
    return n, dropped_chains, dropped_components


def prepare(conf: dict, arm: str, stage1_conf: dict, limit: int) -> int:
    import gemmi

    config_dir = Path(conf["CONFIG_DIR"])
    data_dir = Path(conf["DATA_DIR"])
    results_dir = Path(conf["RESULTS_DIR"])
    s1_config = Path(stage1_conf["CONFIG_DIR"])
    s1_config.mkdir(parents=True, exist_ok=True)

    subset = L.read_tsv(config_dir / "docking_subset.tsv")
    if limit:
        subset = subset[:limit]
    targets = {t["target_id"]: t for t in L.read_tsv(config_dir / "targets.tsv")}
    predictions = {}
    pred_path = results_dir / "predictions.tsv"
    if pred_path.is_file():
        for p in L.read_tsv(pred_path):
            if p.get("status") == "ok" and p.get("arm") == conf.get("DOCK_SOURCE_ARM",
                                                                    "af2_msa_notmpl"):
                predictions[p["target_id"]] = p

    out_dir = data_dir / "docking" / arm
    out_dir.mkdir(parents=True, exist_ok=True)
    rows, align_rows = [], []

    for item in subset:
        tid = item["target_id"]
        target = targets.get(tid)
        if target is None:
            continue
        ref_gz = data_dir / "reference" / f"{tid}.cif.gz"
        lig = data_dir / "reference_ligands" / f"{tid}__{item['comp_id']}.sdf"
        align = {"target_id": tid, "arm": arm, "status": "ok", "reason": "",
                 "recorded": L.now_iso()}
        if not ref_gz.is_file() or not lig.is_file():
            align.update({"status": "failed",
                          "reason": "the reference structure or its ligand is missing"})
            align_rows.append(align)
            continue

        ref_cif = out_dir / f"{tid}_reference.cif"
        L.gunzip_to(ref_gz, ref_cif)
        ref_st = gemmi.read_structure(str(ref_cif))
        ref_st.setup_entities()
        # The chain the ligand actually sits in, which is not always the one
        # the manifest names. An entry with several copies of the protein has
        # a copy of the ligand in each; 03a_choose_ligand_instances.py picks
        # one and records which chain it belongs to, and the receptor has to
        # be that chain or the search box lands beside it.
        want_chain = item.get("receptor_auth_asym_id") or None
        ref_res, ref_chain = polymer_residues(ref_st, want_chain)
        if not ref_res and want_chain:
            ref_res, ref_chain = polymer_residues(ref_st)
        receptor = out_dir / f"{tid}.pdb"

        if arm == "crystal":
            n, dropped_chains, dropped_comps = write_receptor_pdb(
                ref_st, receptor, keep_chain=ref_chain)
            align.update({"n_matched_ca": "", "superposition_rmsd_ca": "",
                          "n_model_residues": n, "n_reference_residues": len(ref_res),
                          "receptor_chain": ref_chain,
                          "chains_dropped": dropped_chains,
                          "components_dropped": dropped_comps})
        else:
            pred = predictions.get(tid)
            if pred is None or not pred.get("structure_file"):
                align.update({"status": "failed",
                              "reason": "no prediction for this target in the chosen arm"})
                align_rows.append(align)
                ref_cif.unlink(missing_ok=True)
                continue
            model_gz = data_dir / pred["structure_file"]
            if not model_gz.is_file():
                align.update({"status": "failed", "reason": "the predicted structure is missing"})
                align_rows.append(align)
                ref_cif.unlink(missing_ok=True)
                continue
            suffix = "".join(model_gz.suffixes[:-1]) or ".pdb"
            model_path = out_dir / f"{tid}_model{suffix}"
            L.gunzip_to(model_gz, model_path)
            model_st = gemmi.read_structure(str(model_path))
            model_st.setup_entities()
            mod_res, mod_chain = polymer_residues(model_st)

            pairs = [(mod_res[i][1], ref_res[i][1]) for i in sorted(mod_res)
                     if i in ref_res]
            if len(pairs) < 4:
                align.update({"status": "failed",
                              "reason": f"only {len(pairs)} alpha carbons matched the reference"})
                align_rows.append(align)
                model_path.unlink(missing_ok=True)
                ref_cif.unlink(missing_ok=True)
                continue
            rot, cm, cr, rmsd = kabsch([p[0] for p in pairs], [p[1] for p in pairs])
            import numpy as np
            for model in model_st:
                for chain in model:
                    for res in chain:
                        for atom in res:
                            v = np.array([atom.pos.x, atom.pos.y, atom.pos.z]) - cm
                            v = rot @ v + cr
                            atom.pos = gemmi.Position(float(v[0]), float(v[1]), float(v[2]))
            n, dropped_chains, dropped_comps = write_receptor_pdb(
                model_st, receptor, keep_chain=mod_chain)
            align.update({"n_matched_ca": len(pairs),
                          "superposition_rmsd_ca": L.fmt(rmsd, 3),
                          "n_model_residues": n, "n_reference_residues": len(ref_res),
                          "receptor_chain": mod_chain,
                          "chains_dropped": dropped_chains,
                          "components_dropped": dropped_comps})
            model_path.unlink(missing_ok=True)

        ref_cif.unlink(missing_ok=True)
        align_rows.append(align)
        if align["status"] != "ok":
            continue
        rows.append({
            "dataset": arm,
            "complex_id": tid,
            "pdb_id": item["pdb_id"],
            "ccd_id": item["comp_id"],
            "protein_pdb": str(receptor),
            "ligand_sdf": str(lig),
            "ligands_sdf": "",
            "start_conf_sdf": "",
        })

    L.write_tsv(s1_config / f"dataset_{arm}.tsv", MANIFEST_COLUMNS, rows)
    # Replaced, not appended. Both tables are rewritten every time an arm is
    # rerun, and appending left three generations of rows in them with
    # nothing to say which was current.
    # Keyed on the arm alone, so a rerun that covers fewer targets does not
    # leave the ones it dropped behind. The subset shrank from fifteen to
    # thirteen and the two that went stayed in this table.
    L.replace_rows(results_dir / "docking_alignment.tsv", ALIGN_COLUMNS,
                   ("arm",), align_rows)
    ok = [a for a in align_rows if a["status"] == "ok"]
    print(f"[handoff] {len(rows)} receptors written for the {arm} arm, "
          f"{len(align_rows) - len(ok)} could not be prepared")
    if arm == "predicted":
        vals = [float(a["superposition_rmsd_ca"]) for a in ok
                if a.get("superposition_rmsd_ca") not in ("", None)]
        if vals:
            q = L.quantiles(vals)
            print(f"[handoff] superposition deviation onto the experimental structure: "
                  f"median {q[1]:.2f} Angstrom, range {min(vals):.2f} to {max(vals):.2f}")
    return 0


def collect(conf: dict, arm: str, stage1_conf: dict) -> int:
    """Read the docking pipeline's own per-run table and compute the rates.

    The rates are computed the way that pipeline computes them, from the same
    columns, so the numbers here and its own published numbers mean the same
    thing. One difference is recorded rather than silently absorbed: its
    denominator counts runs that produced a pose, so runs that failed outright
    are not in it. Both denominators are reported.
    """
    results_dir = Path(conf["RESULTS_DIR"])
    s1_results = Path(stage1_conf["RESULTS_DIR"])
    scored = s1_results / "scored" / "run_scores.tsv"
    if not scored.is_file():
        L.eprint(f"[error] {scored} not found; the docking stage wrote no scores")
        return 1
    runs = [r for r in L.read_tsv(scored) if r.get("dataset") == arm]
    if not runs:
        L.eprint(f"[error] no rows for dataset {arm} in the docking pipeline's table")
        return 1

    # Everything the pipeline attempted, including runs that produced no pose.
    attempted: dict[tuple, int] = {}
    run_dir = s1_results / "runs"
    if run_dir.is_dir():
        for tsv in run_dir.glob(f"*_{arm}.tsv"):
            for r in L.read_tsv(tsv):
                if r.get("dataset") != arm:
                    continue
                key = (r.get("arm", ""), r.get("method", ""))
                attempted[key] = attempted.get(key, 0) + 1

    # Which targets dropped, and why, rather than only how many.
    #
    # The counts below say 15 attempted and 13 scored. A reader cannot tell
    # from that whether two targets failed for a reason that matters or the
    # pipeline lost them. Both reasons here are real ligand chemistry and
    # belong in the record: an iron-bearing heme that the charge model has no
    # parameters for, and a boron-bearing ligand the docking program has no
    # atom type for. Neither is a fault in this repository and neither should
    # be invisible.
    failures: dict[str, str] = {}
    if run_dir.is_dir():
        for tsv in sorted(run_dir.glob(f"*_{arm}.tsv")):
            for r in L.read_tsv(tsv):
                if r.get("dataset") != arm or r.get("status") == "ok":
                    continue
                # That pipeline names this column "key". Reading it as
                # complex_id gave an empty string every time, so the block
                # below never recorded anything and the two targets that
                # cannot be docked stayed invisible.
                tid = r.get("key") or r.get("complex_id") or ""
                why = (r.get("reason") or "").strip()
                if tid and why and tid not in failures:
                    failures[tid] = why[:180]
    if failures:
        L.clear_exclusions(results_dir, "09_dock_into_predictions", arm=arm)
        for tid, why in sorted(failures.items()):
            L.record_exclusion(results_dir, tid, "09_dock_into_predictions",
                               "the docking pipeline could not produce a pose",
                               detail=why, arm=arm)
        print(f"[handoff] {len(failures)} targets could not be docked; each is in "
              f"excluded.tsv with the reason the pipeline gave")

    rows = []
    groups: dict[tuple, list] = {}
    for r in runs:
        groups.setdefault((r.get("arm", ""), r.get("method", "")), []).append(r)
    for (dock_arm, method), rs in sorted(groups.items()):
        def numeric(rec, key):
            try:
                return float(rec.get(key, ""))
            except (TypeError, ValueError):
                return None
        scored_rows = [r for r in rs if numeric(r, "top1_rmsd") is not None]
        assessed = [r for r in scored_rows if r.get("top1_pb_valid") not in ("", None)]
        within = [r for r in scored_rows if numeric(r, "top1_rmsd") <= 2.0]
        success = [r for r in assessed if r.get("top1_success_2a") == "1"]
        n_scored = len(scored_rows)
        n_assessed = len(assessed)
        lo, hi = L.wilson(len(success), n_assessed) if n_assessed else (float("nan"),) * 2
        med = L.quantiles([numeric(r, "top1_rmsd") for r in scored_rows])[1] \
            if scored_rows else None
        rows.append({
            "arm": f"{arm}/{dock_arm}", "method": method,
            "n_attempted": attempted.get((dock_arm, method), len(rs)),
            "n_scored": n_scored, "n_assessed": n_assessed,
            "n_rmsd_within_2a": len(within),
            "rate_rmsd_within_2a": L.fmt(len(within) / n_scored if n_scored else None, 4),
            "n_success": len(success),
            "rate_success": L.fmt(len(success) / n_assessed if n_assessed else None, 4),
            "ci_low": L.fmt(lo, 4), "ci_high": L.fmt(hi, 4),
            "median_top1_rmsd": L.fmt(med, 3),
            "note": "success requires the pose within 2 Angstrom and every physical check passed",
            "recorded": L.now_iso(),
        })
    L.replace_rows(results_dir / "docking_handoff.tsv", HANDOFF_COLUMNS,
                   ("arm", "method"), rows)
    for r in rows:
        print(f"[handoff] {r['arm']} {r['method']}: {r['n_success']}/{r['n_assessed']} "
              f"success ({r['rate_success']}), median top-1 RMSD {r['median_top1_rmsd']}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--stage1-config", required=True, type=Path)
    ap.add_argument("--arm", required=True, choices=["predicted", "crystal"])
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--collect", action="store_true")
    args = ap.parse_args()

    conf = L.load_conf(args.config)
    stage1_conf = L.load_conf(args.stage1_config)
    if args.collect:
        return collect(conf, args.arm, stage1_conf)
    return prepare(conf, args.arm, stage1_conf, args.limit)


if __name__ == "__main__":
    sys.exit(main())
