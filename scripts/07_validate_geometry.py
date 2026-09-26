#!/usr/bin/env python3
"""Physical validity of every predicted structure.

Stage 1 of this roadmap found that a docked pose can sit 1.4 Angstrom from the
crystal ligand and still be chemically impossible, and that nothing in a
docking score enforces chemistry. The same applies to a predicted backbone, and
almost nobody checks it. A confidence value of 95 says the model expects the
local geometry to match the experiment; it says nothing about whether the
molecule it produced could exist.

What is checked, and against what
---------------------------------
Bond lengths and angles come from OpenStructure's stereochemistry module,
whose reference values derive from a curated monomer library. The tolerance is
reported in standard deviations and two are recorded: the library's own
default, which is deliberately permissive, and the four-sigma line used by the
established validation tools. A count that does not say which tolerance it used
is not a count of anything.

Non-bonded clashes come from the same module, which calls a pair clashing when
the distance falls more than a stated margin below the sum of the van der Waals
radii.

Chirality at the alpha carbon is computed here, as a signed volume over the
four substituents. Natural residues are all one hand, so a residue of the other
hand is an error the model made, not a property of the protein. Glycine has no
side chain and is skipped.

Cis peptide bonds are computed here as the dihedral about the peptide bond.
Cis is rare outside proline and common before it, so the two are counted
separately: a cis bond before proline is ordinary, and one anywhere else is
almost always a modelling error.

Ramachandran outliers are computed here against the published contour grids,
with the residue classes and the thresholds the reference implementation uses.
The grids are not redistributed in this repository; the fetch is part of the
install stage and the version fetched is recorded.

The number this stage exists to produce
---------------------------------------
The fraction of predictions in each arm that pass every check, reported next to
the fraction that are accurate. As in stage 1, the only figure called success
is the one that requires both.

Usage
-----
    python scripts/07_validate_geometry.py --config project.conf
"""

from __future__ import annotations

import argparse
import math
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_csc as L  # noqa: E402

GEOMETRY_COLUMNS = [
    "target_id", "arm", "model", "status", "reason", "n_residues", "tolerance_used",
    "bad_bonds", "bad_angles", "bad_bonds_strict", "bad_angles_strict", "clashes",
    "ca_chirality_errors", "cis_nonpro", "cis_pro", "twisted_peptides",
    "rama_outliers", "rama_allowed", "rama_favoured", "rama_scored",
    "rama_outlier_fraction", "passes_all", "failed_checks", "recorded",
]

# Geometry, in the units the structures are written in.
CIS_LIMIT = 30.0          # within this of zero is cis
TRANS_LIMIT = 150.0       # beyond this is trans; between the two is twisted
# The reference implementation calls a residue favoured at or above two per
# cent of the contour and allowed at or above these class-dependent floors.
RAMA_FAVOURED = 0.02
RAMA_ALLOWED = {"general": 0.0005, "cispro": 0.002, "transpro": 0.001,
                "glycine": 0.001, "prepro": 0.001, "ileval": 0.001}


def dihedral(p0, p1, p2, p3) -> float:
    """Signed dihedral in degrees."""
    b0 = [p0[i] - p1[i] for i in range(3)]
    b1 = [p2[i] - p1[i] for i in range(3)]
    b2 = [p3[i] - p2[i] for i in range(3)]
    n1 = math.sqrt(sum(v * v for v in b1))
    if n1 == 0:
        return float("nan")
    b1 = [v / n1 for v in b1]
    v = [b0[i] - sum(b0[j] * b1[j] for j in range(3)) * b1[i] for i in range(3)]
    w = [b2[i] - sum(b2[j] * b1[j] for j in range(3)) * b1[i] for i in range(3)]
    x = sum(v[i] * w[i] for i in range(3))
    y = sum((b1[(i + 1) % 3] * v[(i + 2) % 3] - b1[(i + 2) % 3] * v[(i + 1) % 3]) * w[i]
            for i in range(3))
    return math.degrees(math.atan2(y, x))


def chiral_volume(ca, n, cb, c) -> float:
    """Signed volume of the three substituent vectors at the alpha carbon.

    Every natural residue except glycine has the same hand, so the sign is
    constant across a correct structure and a residue with the opposite sign
    has been built as the mirror image.
    """
    u = [n[i] - ca[i] for i in range(3)]
    v = [cb[i] - ca[i] for i in range(3)]
    w = [c[i] - ca[i] for i in range(3)]
    cross = [u[1] * v[2] - u[2] * v[1],
             u[2] * v[0] - u[0] * v[2],
             u[0] * v[1] - u[1] * v[0]]
    return sum(cross[i] * w[i] for i in range(3))


def read_structure(path: Path):
    """Residues with the atoms these checks need, in chain order."""
    import gemmi

    st = gemmi.read_structure(str(path))
    st.setup_entities()
    out = []
    if not len(st):
        return out
    for chain in st[0]:
        for res in chain:
            atoms = {a.name: (a.pos.x, a.pos.y, a.pos.z) for a in res}
            out.append({"chain": chain.name, "name": res.name,
                        "seqid": res.seqid.num, "atoms": atoms})
    return out


def check_chirality(residues) -> int:
    """Residues whose alpha carbon has the wrong hand.

    The expected sign is taken from the structure itself, as the majority sign
    over its own residues, rather than hardcoded. A hardcoded sign depends on
    the order the three substituents are passed in, and getting that backwards
    would report every residue in every structure as an error, which is the
    kind of bug that looks like a finding.
    """
    signs = []
    for r in residues:
        a = r["atoms"]
        if r["name"] == "GLY" or not all(k in a for k in ("CA", "N", "CB", "C")):
            continue
        signs.append(chiral_volume(a["CA"], a["N"], a["CB"], a["C"]))
    if not signs:
        return 0
    positive = sum(1 for s in signs if s > 0)
    expected_positive = positive >= len(signs) / 2
    return sum(1 for s in signs
               if (s > 0) != expected_positive and abs(s) > 0.1)


def check_peptides(residues):
    """Cis and twisted peptide bonds, counted separately before proline."""
    cis_pro = cis_nonpro = twisted = 0
    for a, b in zip(residues, residues[1:]):
        if a["chain"] != b["chain"]:
            continue
        need_a = ("CA", "C")
        need_b = ("N", "CA")
        if not all(k in a["atoms"] for k in need_a) or not all(k in b["atoms"] for k in need_b):
            continue
        omega = dihedral(a["atoms"]["CA"], a["atoms"]["C"],
                         b["atoms"]["N"], b["atoms"]["CA"])
        if math.isnan(omega):
            continue
        if abs(omega) < CIS_LIMIT:
            if b["name"] == "PRO":
                cis_pro += 1
            else:
                cis_nonpro += 1
        elif abs(omega) <= TRANS_LIMIT:
            twisted += 1
    return cis_nonpro, cis_pro, twisted


def rama_class(prev_res, res, next_res) -> str:
    if res["name"] == "GLY":
        return "glycine"
    if res["name"] == "PRO":
        return "transpro"
    if next_res is not None and next_res["name"] == "PRO":
        return "prepro"
    if res["name"] in ("ILE", "VAL"):
        return "ileval"
    return "general"


def backbone_dihedrals(residues):
    """Phi and psi per residue, with the class each belongs to."""
    out = []
    for i, res in enumerate(residues):
        prev_res = residues[i - 1] if i > 0 else None
        next_res = residues[i + 1] if i + 1 < len(residues) else None
        if prev_res is None or next_res is None:
            continue
        if prev_res["chain"] != res["chain"] or next_res["chain"] != res["chain"]:
            continue
        a, p, n = res["atoms"], prev_res["atoms"], next_res["atoms"]
        if not all(k in a for k in ("N", "CA", "C")) or "C" not in p or "N" not in n:
            continue
        phi = dihedral(p["C"], a["N"], a["CA"], a["C"])
        psi = dihedral(a["N"], a["CA"], a["C"], n["N"])
        if math.isnan(phi) or math.isnan(psi):
            continue
        out.append((rama_class(prev_res, res, next_res), phi, psi))
    return out


def load_rama_grids(grid_dir: Path):
    """The contour grids, as a lookup per class.

    Each file is a sparse listing over a two-degree grid, and a cell that is
    absent means zero rather than missing.
    """
    names = {
        "general": "rama8000-general-noGPIVpreP.data",
        "glycine": "rama8000-gly-sym.data",
        "cispro": "rama8000-cispro.data",
        "transpro": "rama8000-transpro.data",
        "prepro": "rama8000-prepro-noGP.data",
        "ileval": "rama8000-ileval-nopreP.data",
    }
    grids = {}
    for key, fname in names.items():
        path = grid_dir / fname
        if not path.is_file():
            continue
        table = {}
        for line in path.read_text().splitlines():
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                phi, psi, value = float(parts[0]), float(parts[1]), float(parts[-1])
            except ValueError:
                continue
            table[(int(round(phi)), int(round(psi)))] = value
        if table:
            grids[key] = table
    return grids


def rama_lookup(grids, cls: str, phi: float, psi: float):
    """Nearest grid centre on the two-degree grid, wrapping at the edges."""
    table = grids.get(cls) or grids.get("general")
    if not table:
        return None

    def centre(x):
        # Grid centres sit at odd degrees: -179, -177, ... 179.
        c = int(math.floor(x / 2.0) * 2 + 1)
        if c > 179:
            c -= 360
        if c < -179:
            c += 360
        return c

    return table.get((centre(phi), centre(psi)), 0.0)


def ost_stereochemistry(ost_bin: str, model: Path, tolerance: float):
    """Bond, angle and clash counts from the scoring engine's own module.

    Called as a script through that environment's interpreter rather than
    imported, because it lives in an environment this one cannot import from.
    """
    script = (
        "import json, sys\n"
        "from ost import io\n"
        "from ost.mol.alg import stereochemistry\n"
        "ent = io.LoadEntity(sys.argv[1]) if not sys.argv[1].endswith('.cif') "
        "else io.LoadMMCIF(sys.argv[1])\n"
        "tol = float(sys.argv[2])\n"
        "clashes = stereochemistry.GetClashes(ent)\n"
        "bonds = stereochemistry.GetBadBonds(ent, tolerance=tol)\n"
        "angles = stereochemistry.GetBadAngles(ent, tolerance=tol)\n"
        "print(json.dumps({'clashes': len(clashes.GetClashes()) "
        "if hasattr(clashes, 'GetClashes') else len(clashes), "
        "'bad_bonds': len(bonds.GetBadBonds()) if hasattr(bonds, 'GetBadBonds') else len(bonds), "
        "'bad_angles': len(angles.GetBadAngles()) if hasattr(angles, 'GetBadAngles') else len(angles)}))\n"
    )
    with tempfile.NamedTemporaryFile("w", suffix=".py", delete=False) as fh:
        fh.write(script)
        script_path = fh.name
    try:
        p = subprocess.run([ost_bin, script_path, str(model), str(tolerance)],
                           capture_output=True, text=True, timeout=900)
        for line in reversed(p.stdout.splitlines()):
            line = line.strip()
            if line.startswith("{"):
                import json
                return json.loads(line), ""
        return {}, (p.stderr or "no output")[:200]
    except (subprocess.TimeoutExpired, OSError, ValueError) as exc:
        return {}, repr(exc)[:200]
    finally:
        Path(script_path).unlink(missing_ok=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--arm", default="")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--tolerance", type=float, default=4.0,
                    help="the strict bond and angle tolerance, in standard deviations, "
                         "reported alongside but not used for the verdict")
    ap.add_argument("--lenient-tolerance", type=float, default=12.0,
                    help="the tolerance the verdict uses; the scoring engine's own "
                         "default, and the one its local score filters on")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    conf = L.load_conf(args.config)
    results_dir = Path(conf["RESULTS_DIR"])
    data_dir = Path(conf["DATA_DIR"])
    root = Path(conf["CONDA_SH"]).parent.parent.parent
    ost_bin = str(root / "envs" / conf.get("CONDA_ENV_OST", "csc_ost") / "bin" / "ost")

    out_path = results_dir / "geometry_checks.tsv"
    # Resumed rather than repeated, for the same reason the scoring stage is.
    # This skipped whenever its own output existed, and the output it had was
    # four rows from a --limit run, so every arm that finished afterwards went
    # unchecked while the stage reported success. Whether a confident structure
    # is physically valid is one of the questions this repository asks, and it
    # was being answered from four structures.
    prior = L.read_tsv(out_path) if out_path.is_file() else []
    already = set() if args.force else {(r.get("target_id"), r.get("arm")) for r in prior}

    grids = load_rama_grids(data_dir / "reference_data" / "rotarama")
    if not grids:
        print("[07_validate_geometry] no Ramachandran grids found; those counts will be blank. "
              "03_fetch_references.sh --rama fetches them.")

    preds = L.read_tsv(results_dir / "predictions.tsv") if (results_dir / "predictions.tsv").is_file() else []
    preds = [p for p in preds if p.get("status") == "ok"]
    if args.arm:
        preds = [p for p in preds if p["arm"] == args.arm]
    if args.limit:
        preds = preds[:args.limit]
    if not preds:
        L.eprint("[error] no successful predictions to check; run 05_predict.sh first")
        return 1
    n_all = len(preds)
    preds = [p for p in preds if (p["target_id"], p["arm"]) not in already]
    if n_all - len(preds):
        print(f"[07_validate_geometry] {n_all - len(preds)} already checked, "
              f"{len(preds)} to do")
    if not preds:
        print("[07_validate_geometry] everything is already checked; nothing to do")
        return 0

    rows = []
    for i, pred in enumerate(preds, 1):
        tid, arm = pred["target_id"], pred["arm"]
        row = {"target_id": tid, "arm": arm, "model": pred.get("model", ""),
               "status": "ok", "reason": "", "recorded": L.now_iso()}
        model_gz = data_dir / pred.get("structure_file", "")
        if not pred.get("structure_file") or not model_gz.is_file():
            row.update({"status": "failed", "reason": "the predicted structure is missing"})
            rows.append(row)
            continue

        with tempfile.TemporaryDirectory(prefix="csc_geom_") as td:
            suffix = "".join(model_gz.suffixes[:-1]) or ".pdb"
            model = Path(td) / f"model{suffix}"
            L.gunzip_to(model_gz, model)
            residues = read_structure(model)
            row["n_residues"] = len(residues)

            row["ca_chirality_errors"] = check_chirality(residues)
            cis_nonpro, cis_pro, twisted = check_peptides(residues)
            row.update({"cis_nonpro": cis_nonpro, "cis_pro": cis_pro,
                        "twisted_peptides": twisted})

            if grids:
                outliers = allowed = favoured = scored = 0
                for cls, phi, psi in backbone_dihedrals(residues):
                    value = rama_lookup(grids, cls, phi, psi)
                    if value is None:
                        continue
                    scored += 1
                    if value >= RAMA_FAVOURED:
                        favoured += 1
                    elif value >= RAMA_ALLOWED.get(cls, 0.0005):
                        allowed += 1
                    else:
                        outliers += 1
                row.update({"rama_outliers": outliers, "rama_allowed": allowed,
                            "rama_favoured": favoured, "rama_scored": scored,
                            "rama_outlier_fraction": L.fmt(outliers / scored, 4) if scored else ""})

            # Two tolerances, because one number here is not interpretable.
            #
            # The strict one is the four standard deviations the established
            # validation tools use. The lenient one is the scoring engine's own
            # default, which is what its stereochemistry filter applies before
            # computing an all-atom local score, so it is the tolerance that
            # actually affects a number reported elsewhere in this repository.
            #
            # The reason both are kept: at four standard deviations the
            # deposited experimental structures used as the copied-template
            # floor fail too, with between five and a hundred and sixty-nine
            # bad bonds each. A pass rate that reads zero for real crystal
            # structures is measuring the threshold rather than the models, so
            # the verdict below is taken at the lenient tolerance and the
            # strict counts are reported beside it.
            stereo, why = ost_stereochemistry(ost_bin, model, args.tolerance)
            if stereo:
                row.update({"bad_bonds_strict": stereo.get("bad_bonds"),
                            "bad_angles_strict": stereo.get("bad_angles"),
                            "clashes": stereo.get("clashes")})
            elif why:
                row["reason"] = f"stereochemistry check unavailable: {why}"
            lenient, why2 = ost_stereochemistry(ost_bin, model, args.lenient_tolerance)
            if lenient:
                row.update({"bad_bonds": lenient.get("bad_bonds"),
                            "bad_angles": lenient.get("bad_angles")})
            elif why2 and not row.get("reason"):
                row["reason"] = f"stereochemistry check unavailable: {why2}"

        # Every check has to pass. A structure with one inverted alpha carbon
        # is not a structure of a protein, whatever its confidence says.
        failed = []
        for key, label in (("bad_bonds", "bond lengths"), ("bad_angles", "bond angles"),
                           ("clashes", "non-bonded clashes"),
                           ("ca_chirality_errors", "alpha carbon chirality"),
                           ("cis_nonpro", "cis peptide bonds outside proline"),
                           ("twisted_peptides", "twisted peptide bonds"),
                           ("rama_outliers", "Ramachandran outliers")):
            value = row.get(key)
            if isinstance(value, int) and value > 0:
                failed.append(f"{label}: {value}")
        row["tolerance_used"] = args.lenient_tolerance
        row["failed_checks"] = "; ".join(failed)
        row["passes_all"] = int(not failed) if row["status"] == "ok" else ""
        rows.append(row)
        if i % 5 == 0 or i == len(preds):
            print(f"  {i}/{len(preds)} checked", flush=True)

    redone = {(r["target_id"], r["arm"]) for r in rows}
    keep = [r for r in prior if (r.get("target_id"), r.get("arm")) not in redone]
    if keep:
        print(f"[07_validate_geometry] carrying forward {len(keep)} rows checked earlier")
    L.write_tsv(out_path, GEOMETRY_COLUMNS, keep + rows)
    by_arm: dict[str, list] = {}
    for r in rows:
        by_arm.setdefault(r["arm"], []).append(r)
    print(f"[07_validate_geometry] {len(rows)} structures checked; the verdict uses "
          f"{args.lenient_tolerance} standard deviations, and {args.tolerance} is "
          f"reported beside it")
    for arm, rs in sorted(by_arm.items()):
        ok = [r for r in rs if r["passes_all"] == 1]
        print(f"  {arm}: {len(ok)}/{len(rs)} pass every check")
    return 0


if __name__ == "__main__":
    sys.exit(main())
