#!/usr/bin/env python3
"""
=============================================================================
 03a_choose_ligand_instances.py - which copy of the ligand, and in which chain
=============================================================================
 An entry with several copies of the protein has a copy of the ligand in each,
 and they are separate instances with separate identifiers. Which one is
 chosen is not cosmetic, because the docking pipeline centres its search box
 on the ligand's own coordinates and the receptor handed to it is one chain.
 Choose the copy bound to a different chain and the box lands beside the
 receptor rather than in it.

 That is what happened here. The subset was built by taking the instance with
 the best real-space correlation and ignoring which chain it sat on, and on
 four of fifteen targets the chosen instance was 6.6, 9.5, 10.7 and 17.1
 Angstroms from the chain being docked into. Those docking results were not
 poor, they were meaningless, and nothing in the output said so.

 Choosing well needs the coordinates, which is why this runs after the
 references are fetched rather than inside the stage that queries the archive.

 The rule, in order:

   1. The instance must sit on a chain whose sequence is the one the
      prediction is of. The same chain qualifies; so does any other chain of
      the entry with an identical sequence, because a prediction of that
      entity is a prediction of either copy. The chain it sits on is recorded
      and the docking handoff uses that chain as the receptor.
   2. Among those, the instance with the best real-space correlation.
   3. A target with no instance on any such chain is dropped, with the reason
      recorded. Its ligand binds a different protein in the same crystal, so
      there is nothing about the predicted entity to measure.

 Usage:
     python scripts/03a_choose_ligand_instances.py
     python scripts/03a_choose_ligand_instances.py --dry-run
=============================================================================
"""
from __future__ import annotations

import argparse
import gzip
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_csc as L  # noqa: E402

RECEPTOR_COLUMN = "receptor_auth_asym_id"


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="choose one ligand instance per target")
    here = Path(__file__).resolve().parent.parent
    p.add_argument("--config", default=str(here / "project.conf"))
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args(argv)


def chains_matching(cif_gz: Path, want: str) -> set[str]:
    """Every chain of the entry whose polymer sequence is the target chain's.

    Sequence rather than entity identifier, because the two agree for this
    purpose and the sequence is readable straight from the coordinates without
    a second request to the archive.
    """
    import gemmi

    with tempfile.NamedTemporaryFile("wb", suffix=".cif", delete=False) as fh:
        with gzip.open(cif_gz, "rb") as gz:
            fh.write(gz.read())
        tmp = Path(fh.name)
    try:
        st = gemmi.read_structure(str(tmp))
        st.setup_entities()
        seqs = {}
        for chain in st[0]:
            poly = chain.get_polymer()
            if len(poly) >= 20:
                seqs[chain.name] = gemmi.one_letter_code([r.name for r in poly])
    finally:
        tmp.unlink(missing_ok=True)
    ref = seqs.get(want)
    if ref is None:
        return set()
    return {c for c, s in seqs.items() if s == ref}


def main(argv=None) -> int:
    args = parse_args(argv)
    conf = L.load_conf(args.config)
    paths = L.repo_paths(conf)
    config_dir, results, data = paths["CONFIG_DIR"], paths["RESULTS_DIR"], paths["DATA_DIR"]

    targets = {r["target_id"]: r for r in L.read_tsv(config_dir / "targets.tsv")}
    cand_file = results / "ligand_candidates.tsv"
    subset_file = config_dir / "docking_subset.tsv"
    if not cand_file.is_file() or not subset_file.is_file():
        L.eprint("[03a] ligand_candidates.tsv or docking_subset.tsv is missing; "
                 "run 02_build_holdout_set.py first")
        return 3

    existing = L.read_tsv(subset_file)
    columns = list(existing[0].keys()) if existing else []
    if RECEPTOR_COLUMN not in columns:
        columns.append(RECEPTOR_COLUMN)
    wanted = {r["target_id"] for r in existing}

    cands = [r for r in L.read_tsv(cand_file)
             if r.get("status") == "ok" and r.get("target_id") in wanted]

    chosen: dict[str, dict] = {}
    dropped: list[tuple[str, str]] = []
    for tid in sorted(wanted):
        want = targets.get(tid, {}).get("auth_asym_id", "")
        ref_gz = data / "reference" / f"{tid}.cif.gz"
        if not ref_gz.is_file():
            dropped.append((tid, "the reference structure is not on disk"))
            continue
        allowed = chains_matching(ref_gz, want)
        here = [r for r in cands if r["target_id"] == tid]
        usable = [r for r in here if r.get("auth_asym_id") in allowed]
        if not usable:
            where = sorted({r.get("auth_asym_id", "?") for r in here})
            dropped.append((tid,
                            f"every copy of its ligand sits on chain "
                            f"{', '.join(where)}, and none of those is the chain "
                            f"the prediction is of ({want}) or a copy of it"))
            continue
        pick = max(usable, key=lambda r: float(r.get("rscc") or 0))
        pick = dict(pick)
        pick[RECEPTOR_COLUMN] = pick.get("auth_asym_id", want)
        chosen[tid] = pick

    print(f"[03a] {len(wanted)} targets in the subset, {len(chosen)} keep a ligand, "
          f"{len(dropped)} dropped")
    for tid, pick in sorted(chosen.items()):
        was = next((r for r in existing if r["target_id"] == tid), {})
        if was.get("label_asym_id") != pick.get("label_asym_id"):
            print(f"  {tid}: instance {was.get('label_asym_id')} on chain "
                  f"{was.get('auth_asym_id')} replaced by {pick['label_asym_id']} "
                  f"on chain {pick[RECEPTOR_COLUMN]}")
        elif pick[RECEPTOR_COLUMN] != targets.get(tid, {}).get("auth_asym_id"):
            print(f"  {tid}: kept instance {pick['label_asym_id']}, receptor is "
                  f"chain {pick[RECEPTOR_COLUMN]} rather than "
                  f"{targets.get(tid, {}).get('auth_asym_id')}")
    for tid, why in dropped:
        print(f"  {tid}: dropped, {why}")

    if args.dry_run:
        print("[03a] dry run; nothing written")
        return 0

    L.clear_exclusions(results, "03a_choose_ligand_instances")
    L.write_tsv(subset_file, columns, [chosen[t] for t in sorted(chosen)])
    for tid, why in dropped:
        L.record_exclusion(results, tid, "03a_choose_ligand_instances", why)
    print(f"[03a] wrote {subset_file}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
