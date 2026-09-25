#!/usr/bin/env python3
"""Measure the deposited arrangements, build the construct, place the prediction.

Called by scripts/08_predict_pedv.sh. Three jobs, and the first two happen
before any prediction is made, because the construct has to be derived from the
coordinates rather than from a residue range quoted somewhere.

Finding the domain from the coordinates
---------------------------------------
The N-terminal domain is identified by its own contacts rather than by a
boundary taken from a paper. Contacts between residues are counted in a
deposited chain, and the first structural block, the stretch at the start of
the chain whose residues contact each other far more than they contact the
rest, is taken as the domain. The boundary chosen is written out with the
contact counts that produced it, so a reader can disagree with it and see
exactly what would change.

That is a mechanical definition and it will not match a curated one exactly.
It has one advantage that matters here: it is reproducible from the deposited
file, and it cannot be quietly adjusted to make a result come out.

Measuring how much the domain moves between deposited entries
-------------------------------------------------------------
Every pair of deposited chains is superposed on the body, the part of the
chain outside the domain, and the displacement of the domain is then measured
without any further fitting. That is the quantity the disagreement in the
literature is about. If the deposited models agree more closely than expected,
this table will say so and the framing changes.

Placing the prediction
----------------------
The same procedure, with the prediction as one of the structures. It is
superposed on the body of each deposited chain in turn, and the domain
displacement is reported against each. The arrangement it resembles is the one
with the smallest displacement, and that is all that is claimed: resemblance,
not adjudication.

The confidence side of the question
-----------------------------------
Two numbers are pulled out of the prediction's own output: the confidence
averaged over the domain, and the pairwise error between the domain and the
body. A model can be confident about two pieces of a structure and have no
idea how they sit relative to each other, and the second number is where that
shows up. If the domain is placed with high confidence and a high inter-domain
error, the model is reporting exactly the uncertainty the experimental record
displays, and that is worth knowing either way.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_csc as L  # noqa: E402

DOMAIN_COLUMNS = ["entry", "chain", "n_residues", "domain_first", "domain_last",
                  "domain_length", "internal_contacts", "external_contacts",
                  "contact_ratio", "method", "recorded"]

PAIR_COLUMNS = ["structure_a", "chain_a", "structure_b", "chain_b",
                "n_body_aligned", "body_rmsd", "n_domain_compared",
                "domain_displacement", "domain_centroid_shift", "note", "recorded"]

CONFIDENCE_COLUMNS = ["arm", "region", "n_residues", "mean_confidence",
                      "min_confidence", "max_confidence", "note", "recorded"]


def load_chain(path: Path, chain_name: str | None = None):
    """Alpha carbons of one polymer chain, keyed by position in the deposited
    sequence so that two structures of the same protein can be paired without
    an alignment step."""
    import gemmi

    st = gemmi.read_structure(str(path))
    st.setup_entities()
    if not len(st):
        return "", {}
    for chain in st[0]:
        if chain_name and chain.name != chain_name:
            continue
        poly = chain.get_polymer()
        if len(poly) < 30:
            continue
        out = {}
        for res in poly:
            if res.label_seq is None:
                continue
            ca = res.find_atom("CA", "*")
            if ca is not None:
                out[int(res.label_seq)] = (ca.pos.x, ca.pos.y, ca.pos.z)
        if out:
            return chain.name, out
    return "", {}


def find_domain(ca: dict, contact_cutoff: float = 12.0, min_len: int = 60):
    """The first structural block of the chain, by contact density.

    The boundary is swept over candidate positions and the one chosen maximises
    how much more the block contacts itself than the rest of the chain. The
    counts at the chosen boundary are returned with it.
    """
    import numpy as np

    keys = sorted(ca)
    if len(keys) < min_len * 2:
        return None
    coords = np.array([ca[k] for k in keys])
    d = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
    contact = (d < contact_cutoff)
    # Neighbours along the chain contact each other whatever the fold, so the
    # band around the diagonal is removed before anything is counted.
    n = len(keys)
    for off in range(-4, 5):
        idx = np.arange(max(0, -off), min(n, n - off))
        contact[idx, idx + off] = False

    best = None
    upper = min(n - min_len, int(n * 0.6))
    for cut in range(min_len, upper):
        internal = int(contact[:cut, :cut].sum() // 2)
        external = int(contact[:cut, cut:].sum())
        ratio = internal / (external + 1.0)
        if best is None or ratio > best[3]:
            best = (cut, internal, external, ratio)
    if best is None:
        return None
    cut, internal, external, ratio = best
    return {"first": keys[0], "last": keys[cut - 1], "cut_index": cut,
            "length": cut, "internal": internal, "external": external,
            "ratio": round(ratio, 3), "keys": keys}


def superpose_on(mob: dict, ref: dict, positions: list[int]):
    """Fit two structures using only the listed sequence positions."""
    import numpy as np

    shared = [p for p in positions if p in mob and p in ref]
    if len(shared) < 8:
        return None
    a = np.array([mob[p] for p in shared])
    b = np.array([ref[p] for p in shared])
    ca_, cb = a.mean(axis=0), b.mean(axis=0)
    h = (a - ca_).T @ (b - cb)
    u, _s, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rot = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    fitted = (rot @ (a - ca_).T).T
    rmsd = float(np.sqrt(((fitted - (b - cb)) ** 2).sum(axis=1).mean()))
    return rot, ca_, cb, rmsd, len(shared)


def displacement(mob: dict, ref: dict, fit, positions: list[int]):
    """How far the listed positions sit apart after a fit made elsewhere."""
    import numpy as np

    rot, ca_, cb, _rmsd, _n = fit
    shared = [p for p in positions if p in mob and p in ref]
    if not shared:
        return None, None, 0
    a = np.array([mob[p] for p in shared])
    b = np.array([ref[p] for p in shared])
    moved = (rot @ (a - ca_).T).T + cb
    per_atom = np.sqrt(((moved - b) ** 2).sum(axis=1))
    centroid = float(np.linalg.norm(moved.mean(axis=0) - b.mean(axis=0)))
    return float(np.sqrt((per_atom ** 2).mean())), centroid, len(shared)


def measure(conf: dict, entries: list[str]) -> int:
    data_dir = Path(conf["DATA_DIR"]) / "pedv"
    results_dir = Path(conf["RESULTS_DIR"]) / "pedv"
    results_dir.mkdir(parents=True, exist_ok=True)

    loaded, domain_rows = {}, []
    for entry in entries:
        path = data_dir / f"{entry}.cif.gz"
        if not path.is_file():
            print(f"  {entry}: not on disk, skipped")
            continue
        import gzip
        import tempfile
        with tempfile.NamedTemporaryFile("wb", suffix=".cif", delete=False) as fh:
            with gzip.open(path, "rb") as gz:
                fh.write(gz.read())
            tmp = Path(fh.name)
        try:
            chain, ca = load_chain(tmp)
        finally:
            tmp.unlink(missing_ok=True)
        if not ca:
            print(f"  {entry}: no usable chain")
            continue
        dom = find_domain(ca)
        if dom is None:
            print(f"  {entry}: chain too short to split")
            continue
        loaded[entry] = (chain, ca, dom)
        domain_rows.append({
            "entry": entry, "chain": chain, "n_residues": len(ca),
            "domain_first": dom["first"], "domain_last": dom["last"],
            "domain_length": dom["length"], "internal_contacts": dom["internal"],
            "external_contacts": dom["external"], "contact_ratio": dom["ratio"],
            "method": "first structural block by contact density, 12 Angstrom cutoff, "
                      "chain neighbours within four positions ignored",
            "recorded": L.now_iso(),
        })
        print(f"  {entry}: chain {chain}, {len(ca)} residues, domain "
              f"{dom['first']} to {dom['last']} ({dom['length']} residues), "
              f"contact ratio {dom['ratio']}")

    if not loaded:
        L.eprint("[error] no entry could be read")
        return 1
    L.write_tsv(results_dir / "domain_definition.tsv", DOMAIN_COLUMNS, domain_rows)

    # Every pair, superposed on the body, domain displacement measured after.
    pair_rows = []
    names = sorted(loaded)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            ch_a, ca_a, dom_a = loaded[a]
            ch_b, ca_b, dom_b = loaded[b]
            body = [k for k in dom_a["keys"][dom_a["cut_index"]:] if k in ca_b]
            dom_positions = [k for k in dom_a["keys"][:dom_a["cut_index"]] if k in ca_b]
            fit = superpose_on(ca_a, ca_b, body)
            if fit is None:
                continue
            rms, centroid, n_dom = displacement(ca_a, ca_b, fit, dom_positions)
            pair_rows.append({
                "structure_a": a, "chain_a": ch_a, "structure_b": b, "chain_b": ch_b,
                "n_body_aligned": fit[4], "body_rmsd": L.fmt(fit[3], 3),
                "n_domain_compared": n_dom,
                "domain_displacement": L.fmt(rms, 3),
                "domain_centroid_shift": L.fmt(centroid, 3),
                "note": "superposed on the body only; the domain was not fitted",
                "recorded": L.now_iso(),
            })
            print(f"  {a} against {b}: body fit {fit[3]:.2f} over {fit[4]} residues, "
                  f"domain displaced {rms:.2f}, centroid moved {centroid:.2f}")
    L.write_tsv(results_dir / "deposited_pairs.tsv", PAIR_COLUMNS, pair_rows)

    # The construct: the domain plus enough of the body to define where it sits.
    # A domain on its own cannot answer the question this arm asks.
    ref = names[0]
    _chain, ca_ref, dom_ref = loaded[ref]
    keys = dom_ref["keys"]
    body_wanted = min(len(keys) - dom_ref["cut_index"], dom_ref["cut_index"])
    take = keys[:dom_ref["cut_index"] + body_wanted]
    seq_path = Path(conf["DATA_DIR"]) / "pedv" / "construct.fasta"
    import gemmi
    import gzip
    import tempfile
    with tempfile.NamedTemporaryFile("wb", suffix=".cif", delete=False) as fh:
        with gzip.open(data_dir / f"{ref}.cif.gz", "rb") as gz:
            fh.write(gz.read())
        tmp = Path(fh.name)
    try:
        st = gemmi.read_structure(str(tmp))
        st.setup_entities()
        letters = {}
        for chain in st[0]:
            poly = chain.get_polymer()
            if len(poly) < 30:
                continue
            for res in poly:
                if res.label_seq is not None:
                    letters[int(res.label_seq)] = gemmi.find_tabulated_residue(
                        res.name).one_letter_code.upper()
            break
    finally:
        tmp.unlink(missing_ok=True)
    seq = "".join(letters.get(k, "") for k in take)
    seq = "".join(c for c in seq if c.isalpha())
    seq_path.parent.mkdir(parents=True, exist_ok=True)
    seq_path.write_text(f">pedv_construct\n{seq}\n")
    print(f"[pedv] construct written: {len(seq)} residues, the domain plus "
          f"{body_wanted} residues of the body, taken from {ref}")

    L.write_tsv(results_dir / "construct.tsv",
                ["source_entry", "domain_first", "domain_last", "body_residues",
                 "construct_length", "note", "recorded"],
                [{"source_entry": ref, "domain_first": dom_ref["first"],
                  "domain_last": dom_ref["last"], "body_residues": body_wanted,
                  "construct_length": len(seq),
                  "note": "the domain plus an equal length of the body, as one chain; "
                          "a single chain cannot say how three copies pack together",
                  "recorded": L.now_iso()}])
    return 0


def compare(conf: dict, entries: list[str]) -> int:
    """Where the prediction put the domain, and what it said about it."""
    data_dir = Path(conf["DATA_DIR"]) / "pedv"
    results_dir = Path(conf["RESULTS_DIR"]) / "pedv"
    results_dir.mkdir(parents=True, exist_ok=True)

    construct = L.read_tsv(results_dir / "construct.tsv")
    if not construct:
        L.eprint("[error] no construct record; run the measurement step first")
        return 1
    cut_len = int(construct[0]["domain_last"]) - int(construct[0]["domain_first"]) + 1

    rows, conf_rows = [], []
    for arm_dir in sorted(p for p in data_dir.iterdir() if p.is_dir()):
        arm = arm_dir.name
        models = sorted(arm_dir.glob("*_unrelaxed_rank_001_*.pdb"))
        scores = sorted(arm_dir.glob("*_scores_rank_001_*.json"))
        if not models:
            continue
        chain_p, ca_p = load_chain(models[0])
        if not ca_p:
            continue
        # The prediction is numbered from one over the construct, and the
        # construct was cut from the deposited numbering, so the offset is
        # known rather than guessed.
        offset = int(construct[0]["domain_first"]) - 1
        ca_pred = {k + offset: v for k, v in ca_p.items()}
        dom_positions = [k for k in ca_pred if k <= int(construct[0]["domain_last"])]
        body_positions = [k for k in ca_pred if k > int(construct[0]["domain_last"])]

        for entry in entries:
            path = data_dir / f"{entry}.cif.gz"
            if not path.is_file():
                continue
            import gzip
            import tempfile
            with tempfile.NamedTemporaryFile("wb", suffix=".cif", delete=False) as fh:
                with gzip.open(path, "rb") as gz:
                    fh.write(gz.read())
                tmp = Path(fh.name)
            try:
                chain_d, ca_d = load_chain(tmp)
            finally:
                tmp.unlink(missing_ok=True)
            if not ca_d:
                continue
            fit = superpose_on(ca_pred, ca_d, body_positions)
            if fit is None:
                continue
            rms, centroid, n_dom = displacement(ca_pred, ca_d, fit, dom_positions)
            rows.append({
                "structure_a": f"prediction/{arm}", "chain_a": chain_p,
                "structure_b": entry, "chain_b": chain_d,
                "n_body_aligned": fit[4], "body_rmsd": L.fmt(fit[3], 3),
                "n_domain_compared": n_dom,
                "domain_displacement": L.fmt(rms, 3),
                "domain_centroid_shift": L.fmt(centroid, 3),
                "note": "superposed on the body only; resemblance is not adjudication",
                "recorded": L.now_iso(),
            })

        if scores:
            data = json.loads(scores[0].read_text())
            plddt = data.get("plddt") or []
            pae = data.get("pae") or []
            n_dom_idx = min(cut_len, len(plddt))
            for region, values in (("domain", plddt[:n_dom_idx]),
                                   ("body", plddt[n_dom_idx:])):
                if not values:
                    continue
                conf_rows.append({
                    "arm": arm, "region": region, "n_residues": len(values),
                    "mean_confidence": L.fmt(sum(values) / len(values), 2),
                    "min_confidence": L.fmt(min(values), 2),
                    "max_confidence": L.fmt(max(values), 2),
                    "note": "per-residue confidence, on the score's own scale",
                    "recorded": L.now_iso(),
                })
            # The number that answers the question: how certain the model is
            # about where the domain sits relative to the body, as opposed to
            # how certain it is about each of them.
            if pae and n_dom_idx and len(pae) > n_dom_idx:
                cross = [pae[i][j] for i in range(n_dom_idx)
                         for j in range(n_dom_idx, len(pae[i]))]
                within_dom = [pae[i][j] for i in range(n_dom_idx)
                              for j in range(n_dom_idx)]
                if cross:
                    conf_rows.append({
                        "arm": arm, "region": "domain against body",
                        "n_residues": len(cross),
                        "mean_confidence": L.fmt(sum(cross) / len(cross), 2),
                        "min_confidence": L.fmt(min(cross), 2),
                        "max_confidence": L.fmt(max(cross), 2),
                        "note": "pairwise predicted error in Angstrom, not a confidence; "
                                "large here means the model does not know how the two "
                                "pieces sit relative to each other",
                        "recorded": L.now_iso(),
                    })
                if within_dom:
                    conf_rows.append({
                        "arm": arm, "region": "within the domain",
                        "n_residues": len(within_dom),
                        "mean_confidence": L.fmt(sum(within_dom) / len(within_dom), 2),
                        "min_confidence": L.fmt(min(within_dom), 2),
                        "max_confidence": L.fmt(max(within_dom), 2),
                        "note": "pairwise predicted error in Angstrom, not a confidence",
                        "recorded": L.now_iso(),
                    })

    L.append_tsv(results_dir / "prediction_pairs.tsv", PAIR_COLUMNS, rows)
    L.write_tsv(results_dir / "prediction_confidence.tsv", CONFIDENCE_COLUMNS, conf_rows)
    for r in rows:
        print(f"  {r['structure_a']} against {r['structure_b']}: "
              f"domain displaced {r['domain_displacement']} after fitting the body")
    for r in conf_rows:
        print(f"  {r['arm']} {r['region']}: mean {r['mean_confidence']} "
              f"over {r['n_residues']} values")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--entries", required=True)
    ap.add_argument("--measure", action="store_true")
    ap.add_argument("--compare", action="store_true")
    args = ap.parse_args()
    conf = L.load_conf(args.config)
    entries = [e.strip() for e in args.entries.split(",") if e.strip()]
    if args.measure:
        return measure(conf, entries)
    if args.compare:
        return compare(conf, entries)
    ap.error("one of --measure or --compare is required")
    return 2


if __name__ == "__main__":
    sys.exit(main())
