#!/usr/bin/env python3
"""The two floors every prediction is measured against.

A success rate without a floor is uninterpretable. The previous stage of this
work reported that dropping a ligand into the search box without searching at
all succeeded 0.3 per cent of the time, and that figure is what made every
other number in it readable. The same applies here, and there are two floors
worth having.

The copied template
-------------------
The highest-identity chain that was in the public archive before the cutoff,
copied verbatim, superposed onto the target and scored exactly as a prediction
is scored. No modelling, no side-chain repacking, no refinement: the
coordinates as deposited. This is what a person would have done before these
models existed, and it is the honest comparator. If a model does not beat it by
a wide margin, that is the headline however uncomfortable.

The identity of that chain is already known: the set builder recorded it while
measuring how much of this held-out set the models could have seen. So the
floor costs one download per target and no computation.

The unrelated chain
-------------------
A chain of similar length from an entry that is not related to the target. It
measures what the scoring gives you for a protein-shaped object of the right
size, which is the floor below the floor.

Neither floor produces a confidence value, because neither involves a model.
They appear in the accuracy tables and not in the calibration ones, which is
the correct place for them.

Usage
-----
    python scripts/06a_build_null_floors.py --config project.conf
"""

from __future__ import annotations

import argparse
import gzip
import random
import shutil
import sys
import tempfile
import urllib.error
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_csc as L  # noqa: E402
from lib_predict import PREDICTION_COLUMNS  # noqa: E402

FLOOR_COLUMNS = ["target_id", "arm", "source_entry", "source_chain", "identity",
                 "aligned_identity", "n_matched_ca", "superposition_rmsd_ca",
                 "status", "reason", "recorded"]


def fetch_entry(pdb_id: str, dest: Path) -> bool:
    if dest.is_file() and dest.stat().st_size > 0:
        return True
    url = f"https://files.rcsb.org/download/{pdb_id}.cif.gz"
    try:
        with urllib.request.urlopen(url, timeout=180) as r, open(dest, "wb") as fh:
            shutil.copyfileobj(r, fh)
        return dest.stat().st_size > 0
    except (urllib.error.URLError, TimeoutError, OSError):
        dest.unlink(missing_ok=True)
        return False


def chain_ca(path: Path, entity_id: str | None = None,
             match_seq: str | None = None):
    """One polymer chain: its alpha carbons in order, and its one-letter sequence.

    The coordinates are returned as a list in chain order rather than keyed by
    residue number, because the two structures paired here come from different
    entries and their numbering has nothing to do with each other. The pairing
    is made by aligning the sequences, which is the only correspondence that
    means anything between two different depositions.
    """
    import gemmi

    st = gemmi.read_structure(str(path))
    st.setup_entities()
    if not len(st):
        return None, [], "", None
    # Which chain, and why it cannot be whichever comes first.
    #
    # A template is found by searching sequences, and what comes back is an
    # entity inside an entry. That entry may hold one chain or it may hold a
    # ribosome. Taking the first polymer chain of a large assembly hands back
    # something unrelated to the target, and the alignment then pairs residues
    # that have nothing to do with each other: on this set that produced
    # superpositions of 16 to 34 Angstroms and floors that scored near zero
    # for templates 71 to 86 per cent identical to their target.
    #
    # So when the sequence being matched is known, every polymer chain is
    # aligned against it and the best match is the one returned. That is
    # self-correcting and needs no entity bookkeeping.
    found = []
    for chain in st[0]:
        poly = chain.get_polymer()
        if len(poly) < 20:
            continue
        coords, letters = [], []
        for res in poly:
            atom = res.find_atom("CA", "*")
            if atom is None:
                continue
            info = gemmi.find_tabulated_residue(res.name)
            code = info.one_letter_code.upper() if info else "X"
            letters.append(code if code.isalpha() else "X")
            coords.append((atom.pos.x, atom.pos.y, atom.pos.z))
        if coords:
            found.append((chain.name, coords, "".join(letters)))
    if not found:
        return None, [], "", None
    if match_seq and len(found) > 1:
        def identity_to(seq: str) -> float:
            res = gemmi.align_string_sequences(list(match_seq), list(seq), [],
                                               gemmi.AlignmentScoring())
            return res.calculate_identity()
        found.sort(key=lambda f: identity_to(f[2]), reverse=True)
    name, coords, letters = found[0]
    return name, coords, letters, st


def aligned_pairs(seq_a: str, coords_a: list, seq_b: str, coords_b: list):
    """Index pairs of residues that align between two chains.

    A global alignment of the two sequences, walked through its own
    description of matches and gaps. Pairing by residue number instead looks
    reasonable and is wrong whenever the two entries number their chains
    differently, which is most of the time: it produced a ten Angstrom fit
    between two chains that are 98.8 per cent identical.
    """
    import gemmi

    result = gemmi.align_string_sequences(list(seq_a), list(seq_b), [],
                                          gemmi.AlignmentScoring())
    pairs = []
    i = j = 0
    for token in _cigar_tokens(result.cigar_str()):
        n, op = token
        if op == "M":
            for k in range(n):
                if i + k < len(coords_a) and j + k < len(coords_b):
                    pairs.append((i + k, j + k))
            i += n
            j += n
        elif op == "I":
            i += n
        elif op == "D":
            j += n
    return pairs, result.calculate_identity()


def _cigar_tokens(cigar: str):
    n = ""
    for ch in cigar:
        if ch.isdigit():
            n += ch
        else:
            yield (int(n) if n else 1), ch
            n = ""


def superpose_and_write(src_st, pairs, src_coords, ref_coords, out_gz: Path,
                        keep_chain: str = ""):
    """Fit the copied chain onto the target and write it where a prediction goes.

    A copied template is only a comparator if it is put in the same frame as
    the structure it is compared against, and the deviation of that fit is
    recorded so a reader can see how well the copy lines up before reading how
    well it scores. The scores themselves do not depend on this placement: the
    local score is superposition free and the global ones fit for themselves.
    """
    import gemmi
    import numpy as np

    if len(pairs) < 8:
        return None, 0
    a = np.array([src_coords[i] for i, _ in pairs])
    b = np.array([ref_coords[j] for _, j in pairs])
    ca_, cb = a.mean(axis=0), b.mean(axis=0)
    h = (a - ca_).T @ (b - cb)
    u, _s, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rot = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    rmsd = float(np.sqrt((((rot @ (a - ca_).T).T - (b - cb)) ** 2).sum(axis=1).mean()))

    for model in src_st:
        for chain in model:
            for res in chain:
                for atom in res:
                    v = np.array([atom.pos.x, atom.pos.y, atom.pos.z]) - ca_
                    v = rot @ v + cb
                    atom.pos = gemmi.Position(float(v[0]), float(v[1]), float(v[2]))
    # Only the chain that was matched. A floor stands where a prediction would,
    # and a prediction is one chain; writing the whole source entry would hand
    # the scoring stage an assembly, which is the same mismatch that made every
    # accuracy number in this repository wrong once already.
    if keep_chain:
        for model in src_st:
            names = sorted({c.name for c in model if c.name != keep_chain})
            for name in names:
                model.remove_chain(name)
        src_st.remove_ligands_and_waters()

    out_gz.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td) / "floor.cif"
        src_st.setup_entities()
        src_st.make_mmcif_document().write_file(str(tmp))
        L.gzip_file(tmp, out_gz)
    return rmsd, len(pairs)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    conf = L.load_conf(args.config)
    config_dir = Path(conf["CONFIG_DIR"])
    results_dir = Path(conf["RESULTS_DIR"])
    data_dir = Path(conf["DATA_DIR"])
    seed = L.conf_int(conf, "SEED", 0)

    targets = L.read_tsv(config_dir / "targets.tsv")
    if args.limit:
        targets = targets[:args.limit]
    floor_dir = data_dir / "floors"
    floor_dir.mkdir(parents=True, exist_ok=True)

    rows, pred_rows = [], []
    # The unrelated floor draws from the pool of copied templates, which are
    # pre-cutoff chains of comparable length. Drawing with a recorded seed
    # keeps the floor reproducible; drawing from the same pool keeps it a
    # protein of the right size rather than an arbitrary object.
    rng = random.Random(seed)
    pool = [t for t in targets if t.get("closest_pre_cutoff_entry")]

    for i, t in enumerate(targets, 1):
        tid = t["target_id"]
        ref_gz = data_dir / "reference" / f"{tid}.cif.gz"
        if not ref_gz.is_file():
            continue
        with tempfile.TemporaryDirectory(prefix="csc_floor_") as td:
            ref_cif = Path(td) / "ref.cif"
            L.gunzip_to(ref_gz, ref_cif)
            _ref_chain, ref_coords, ref_seq, _st = chain_ca(ref_cif)
            if not ref_coords:
                continue

            for arm in ("null_template", "null_unrelated"):
                out_gz = data_dir / "predictions" / arm / f"{tid}.cif.gz"
                row = {"target_id": tid, "arm": arm, "status": "ok", "reason": "",
                       "recorded": L.now_iso()}
                if arm == "null_template":
                    source = t.get("closest_pre_cutoff_entry", "")
                    row["identity"] = t.get("closest_pre_cutoff_identity", "")
                else:
                    others = [o for o in pool if o["target_id"] != tid]
                    source = (rng.choice(others).get("closest_pre_cutoff_entry", "")
                              if others else "")
                    row["identity"] = ""
                if not source:
                    row.update({"status": "failed",
                                "reason": "no pre-cutoff chain was found for this target"})
                    rows.append(row)
                    continue
                row["source_entry"] = source
                src_pdb = source.partition("_")[0]
                src_gz = floor_dir / f"{src_pdb}.cif.gz"
                if not fetch_entry(src_pdb, src_gz):
                    row.update({"status": "failed",
                                "reason": f"{src_pdb} could not be downloaded"})
                    rows.append(row)
                    continue
                if out_gz.is_file() and not args.force:
                    row["reason"] = "already built"
                src_cif = Path(td) / f"{src_pdb}.cif"
                L.gunzip_to(src_gz, src_cif)
                src_chain, src_coords, src_seq, src_st = chain_ca(
                    src_cif, match_seq=ref_seq)
                if not src_coords:
                    row.update({"status": "failed", "reason": "no usable chain in the source"})
                    rows.append(row)
                    continue
                row["source_chain"] = src_chain
                pairs, aligned_identity = aligned_pairs(src_seq, src_coords,
                                                        ref_seq, ref_coords)
                row["aligned_identity"] = L.fmt(aligned_identity / 100.0, 4)
                fit = superpose_and_write(src_st, pairs, src_coords, ref_coords, out_gz,
                                          keep_chain=src_chain or "")
                if fit[0] is None:
                    row.update({"status": "failed",
                                "reason": "too few shared positions to place the copy"})
                    rows.append(row)
                    continue
                row["superposition_rmsd_ca"] = L.fmt(fit[0], 3)
                row["n_matched_ca"] = fit[1]
                rows.append(row)
                pred_rows.append({
                    "target_id": tid, "arm": arm, "model": "none",
                    "seed": "", "status": "ok",
                    "reason": "a copied structure, not a prediction; no confidence is produced",
                    "sequence_length": t.get("sequence_length", ""),
                    "structure_file": str(out_gz.relative_to(data_dir)),
                    "confidence_file": "", "elapsed_s": "", "peak_rss_mb": "",
                    "n_models_run": "0", "recorded": L.now_iso(),
                })
        if i % 20 == 0 or i == len(targets):
            print(f"  {i}/{len(targets)} targets", flush=True)

    L.write_tsv(results_dir / "null_floors.tsv", FLOOR_COLUMNS, rows)
    if pred_rows:
        L.append_tsv(results_dir / "predictions.tsv", PREDICTION_COLUMNS, pred_rows)
    ok = [r for r in rows if r["status"] == "ok"]
    print(f"[06a_build_null_floors] {len(ok)} floors built of {len(rows)} attempted")
    for arm in ("null_template", "null_unrelated"):
        n = sum(1 for r in ok if r["arm"] == arm)
        print(f"  {arm}: {n}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
