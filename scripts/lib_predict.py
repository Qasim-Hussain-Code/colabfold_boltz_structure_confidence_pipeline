#!/usr/bin/env python3
"""Per-target input preparation and output collection for the prediction stage.

Called by scripts/05_predict.sh, which owns the loop, the weight lifecycle and
the timing. This file owns what one target looks like going in and what is kept
coming out.

What is kept, and why it is this and not the coordinates
-------------------------------------------------------
The coordinates of a prediction are large and regenerable. The confidence
arrays are small and are the data every later stage reads: the calibration
curves, the domain question, the threshold derivation. So the per-target
artefact here is a compressed JSON holding the per-residue confidence, the
pairwise error matrix and the summary scores, and the structure is kept
separately and compressed.

Both models are normalised into one shape on the way in, because they do not
agree on anything. One reports per-residue confidence on a 0 to 100 scale, the
other on 0 to 1 in its data files and 0 to 100 in the structure it writes. One
writes a pairwise error matrix as JSON, the other as a compressed array. One
numbers residues from 1, the other from 1 but with its own chain naming. A
calibration plot built by reading whichever field happened to be present would
be comparing two different quantities, so the conversion happens once, here,
and is recorded per target with the scale it came from.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_csc as L  # noqa: E402

PREDICTION_COLUMNS = [
    "target_id", "arm", "model", "seed", "status", "reason",
    "sequence_length", "n_residues_scored", "mean_plddt", "ptm",
    "max_pae", "mean_pae", "msa_depth", "n_templates", "template_ids",
    "elapsed_s", "peak_rss_mb", "model_rank_1", "n_models_run",
    "plddt_scale_as_written", "structure_file", "confidence_file", "recorded",
]


# ---------------------------------------------------------------------------
# Input preparation
# ---------------------------------------------------------------------------
def write_fasta(target: dict, path: Path) -> None:
    """One sequence, with the target id as the header.

    The header becomes the job name, and every character outside letters,
    digits, underscore, dot and hyphen is replaced by the tool. Target ids here
    are already of the form 1ABC_1, so nothing is rewritten and the output file
    names stay predictable.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    seq = target["sequence"].strip().upper()
    path.write_text(f">{target['target_id']}\n{seq}\n")


def write_boltz_yaml(target: dict, path: Path, msa_path: str | None) -> None:
    """A single protein chain, with an alignment or explicitly without one.

    The alignment is passed as a file rather than letting the tool fetch its
    own. Both arms have to see the same alignment or the comparison is between
    two different inputs as well as two different models, and the alignment is
    the expensive part that a shared public service produced once.

    msa_path of None means single-sequence mode, which the format spells as the
    bare word empty. Leaving the field out instead makes the tool try to
    contact the alignment server, which would be a second query for an
    alignment already on disk.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    seq = target["sequence"].strip().upper()
    msa_line = f"      msa: {msa_path}" if msa_path else "      msa: empty"
    path.write_text(
        "version: 1\n"
        "sequences:\n"
        "  - protein:\n"
        "      id: A\n"
        f"      sequence: {seq}\n"
        f"{msa_line}\n"
    )


def a3m_depth(path: Path) -> int:
    """Sequences in an alignment, counting the query.

    Counting header lines rather than pairs of lines: an a3m can wrap, and a
    file written by one tool and read by another is exactly where an assumption
    about line breaks goes wrong.
    """
    if not path.is_file():
        return 0
    n = 0
    opener = L.gzip_open_text if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:  # type: ignore[operator]
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


# ---------------------------------------------------------------------------
# Output collection
# ---------------------------------------------------------------------------
def collect_af2(out_dir: Path, target_id: str):
    """Read what ColabFold wrote for one target.

    The rank 1 model is the one taken, because ranking is what the tool is
    asked to do and picking any other would be choosing with hindsight. Which
    model that was is recorded, since a five-model run that ranks model 3 first
    is a different result from one that ranks model 1 first, and the difference
    is what the model-count experiment measures.
    """
    scores = sorted(out_dir.glob(f"{target_id}_scores_rank_001_*.json"))
    if not scores:
        scores = sorted(out_dir.glob("*_scores_rank_001_*.json"))
    if not scores:
        return None, "no rank 1 score file was written"
    score_path = scores[0]
    data = json.loads(score_path.read_text())

    plddt = data.get("plddt") or []
    pae = data.get("pae") or []
    # ColabFold writes per-residue confidence on the 0 to 100 scale, which is
    # the scale the value was defined on.
    out = {
        "target_id": target_id,
        "source": "colabfold",
        "plddt": plddt,
        "plddt_scale": "0-100",
        "pae": pae,
        "pae_units": "angstrom",
        "max_pae": data.get("max_pae"),
        "ptm": data.get("ptm"),
        "iptm": data.get("iptm"),
        "rank_1_file": score_path.name,
    }
    # The file name carries which model and seed produced the top-ranked
    # prediction, and nothing else records it.
    stem = score_path.stem
    for part in stem.split("_"):
        if part.startswith("model") and part[5:].isdigit():
            out["model_rank_1"] = part
    if "_seed_" in stem:
        out["seed"] = stem.split("_seed_")[-1]
    structures = sorted(out_dir.glob(f"{target_id}_unrelaxed_rank_001_*.pdb"))
    if not structures:
        structures = sorted(out_dir.glob("*_unrelaxed_rank_001_*.pdb"))
    return (out, structures[0] if structures else None), ""


def collect_boltz(out_dir: Path, target_id: str):
    """Read what Boltz wrote for one target.

    Confidence comes back on a 0 to 1 scale in the data files and multiplied by
    100 in the structure. It is converted to the 0 to 100 scale here so that one
    number means one thing across the repository, and the scale it arrived on is
    recorded rather than assumed.
    """
    import numpy as np

    pred_dir = out_dir / "predictions" / target_id
    if not pred_dir.is_dir():
        candidates = list((out_dir / "predictions").glob("*")) if (out_dir / "predictions").is_dir() else []
        if not candidates:
            return None, "no predictions directory was written"
        pred_dir = candidates[0]

    conf_files = sorted(pred_dir.glob("confidence_*_model_0.json"))
    if not conf_files:
        return None, "no confidence file for model 0"
    conf = json.loads(conf_files[0].read_text())

    plddt_files = sorted(pred_dir.glob("plddt_*_model_0.npz"))
    plddt = []
    if plddt_files:
        with np.load(plddt_files[0]) as z:
            key = "plddt" if "plddt" in z else list(z.keys())[0]
            plddt = [round(float(v) * 100.0, 2) for v in np.asarray(z[key]).ravel()]

    pae, max_pae = [], None
    pae_files = sorted(pred_dir.glob("pae_*_model_0.npz"))
    if pae_files:
        with np.load(pae_files[0]) as z:
            key = "pae" if "pae" in z else list(z.keys())[0]
            arr = np.asarray(z[key])
            pae = [[round(float(v), 2) for v in row] for row in arr]
            max_pae = float(arr.max()) if arr.size else None

    structures = sorted(pred_dir.glob("*_model_0.cif")) or sorted(pred_dir.glob("*_model_0.pdb"))
    out = {
        "target_id": target_id,
        "source": "boltz2",
        "plddt": plddt,
        "plddt_scale": "0-1 as written, multiplied by 100 here",
        "pae": pae,
        "pae_units": "angstrom",
        "max_pae": max_pae,
        "ptm": conf.get("ptm"),
        "iptm": conf.get("iptm"),
        "complex_plddt": conf.get("complex_plddt"),
        "confidence_score": conf.get("confidence_score"),
        "rank_1_file": conf_files[0].name,
        "model_rank_1": "model_0",
    }
    return (out, structures[0] if structures else None), ""


def summarise(conf: dict) -> dict:
    plddt = conf.get("plddt") or []
    pae = conf.get("pae") or []
    mean_plddt = sum(plddt) / len(plddt) if plddt else None
    mean_pae = None
    if pae:
        total = sum(sum(row) for row in pae)
        n = sum(len(row) for row in pae)
        mean_pae = total / n if n else None
    return {
        "n_residues_scored": len(plddt),
        "mean_plddt": L.fmt(mean_plddt, 3),
        "ptm": L.fmt(conf.get("ptm"), 4),
        "max_pae": L.fmt(conf.get("max_pae"), 3),
        "mean_pae": L.fmt(mean_pae, 3),
        "plddt_scale_as_written": conf.get("plddt_scale", ""),
        "model_rank_1": conf.get("model_rank_1", ""),
    }


# ---------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_in = sub.add_parser("prepare", help="write the input file for one target")
    p_in.add_argument("--config", required=True, type=Path)
    p_in.add_argument("--target", required=True)
    p_in.add_argument("--format", required=True, choices=["fasta", "yaml"])
    p_in.add_argument("--msa", default="")
    p_in.add_argument("--out", required=True, type=Path)

    p_out = sub.add_parser("collect", help="collect one finished prediction")
    p_out.add_argument("--config", required=True, type=Path)
    p_out.add_argument("--target", required=True)
    p_out.add_argument("--arm", required=True)
    p_out.add_argument("--model", required=True)
    p_out.add_argument("--seed", default="")
    p_out.add_argument("--tag", default="",
                       help="name for the output files, which differs from the "
                            "target when one target is predicted more than once")
    p_out.add_argument("--source", required=True, choices=["colabfold", "boltz2"])
    p_out.add_argument("--out-dir", required=True, type=Path)
    p_out.add_argument("--elapsed", default="")
    p_out.add_argument("--peak-rss-kb", default="")
    p_out.add_argument("--msa-file", default="")
    p_out.add_argument("--n-models", default="1")

    args = ap.parse_args()
    conf = L.load_conf(args.config)
    config_dir = Path(conf["CONFIG_DIR"])
    results_dir = Path(conf["RESULTS_DIR"])
    targets = {t["target_id"]: t for t in L.read_tsv(config_dir / "targets.tsv")}
    target = targets.get(args.target)
    if target is None:
        L.eprint(f"[error] {args.target} is not in config/targets.tsv")
        return 1

    if args.cmd == "prepare":
        if args.format == "fasta":
            write_fasta(target, args.out)
        else:
            write_boltz_yaml(target, args.out, args.msa or None)
        return 0

    # collect
    if args.source == "colabfold":
        got, why = collect_af2(args.out_dir, args.target)
    else:
        got, why = collect_boltz(args.out_dir, args.target)

    row = {
        "target_id": args.target, "arm": args.arm, "model": args.model,
        "seed": args.seed, "sequence_length": target.get("sequence_length", ""),
        "elapsed_s": args.elapsed,
        "peak_rss_mb": (round(int(args.peak_rss_kb) / 1024, 1)
                        if str(args.peak_rss_kb).isdigit() else ""),
        "n_models_run": args.n_models,
        "recorded": L.now_iso(),
    }
    if got is None:
        row.update({"status": "failed", "reason": why})
        # Replaced, not appended. A rerun of the same target under the same
        # arm and seed is the same prediction, and a second row for it made
        # one target count twice in every per-target statistic.
        L.replace_rows(results_dir / "predictions.tsv", PREDICTION_COLUMNS,
                       ("target_id", "arm", "seed"), [row])
        L.record_exclusion(results_dir, args.target, "05_predict", why, arm=args.arm)
        print(f"[collect] {args.target} {args.arm}: {why}")
        return 1

    conf_obj, structure = got
    conf_obj["arm"] = args.arm
    conf_obj["model"] = args.model
    conf_obj["seed"] = args.seed

    msa_depth = a3m_depth(Path(args.msa_file)) if args.msa_file else 0
    # The name a repeat carries. A seed-variance run predicts one target
    # several times, and without the seed in the name each repeat overwrote
    # the one before it and the arm's own prediction along with them. The
    # prediction stage already builds this tag; it only had to be passed here.
    tag = args.tag or args.target
    conf_dir = results_dir / "confidence"
    conf_path = conf_dir / f"{tag}__{args.arm}.json.gz"
    L.write_json_gz(conf_path, conf_obj)

    kept_structure = ""
    if structure is not None and structure.is_file():
        struct_dir = Path(conf["DATA_DIR"]) / "predictions" / args.arm
        struct_dir.mkdir(parents=True, exist_ok=True)
        dest = struct_dir / f"{tag}{structure.suffix}.gz"
        L.gzip_file(structure, dest)
        kept_structure = str(dest.relative_to(Path(conf["DATA_DIR"])))

    row.update({"status": "ok", "reason": "", "msa_depth": msa_depth,
                "structure_file": kept_structure,
                "confidence_file": str(conf_path.relative_to(results_dir))})
    row.update(summarise(conf_obj))
    # Replaced, not appended. A rerun of the same target under the same arm
    # and seed is the same prediction, and stacking a second row for it made
    # one target count twice in every per-target statistic. This is the same
    # fault the floor builder and the docking tables had.
    L.replace_rows(results_dir / "predictions.tsv", PREDICTION_COLUMNS,
                   ("target_id", "arm", "seed"), [row])
    print(f"[collect] {args.target} {args.arm}: {row['n_residues_scored']} residues, "
          f"mean pLDDT {row['mean_plddt']}, MSA depth {msa_depth}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
