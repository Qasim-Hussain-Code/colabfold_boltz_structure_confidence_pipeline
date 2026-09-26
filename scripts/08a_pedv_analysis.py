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
                  "normalised_cut", "method", "recorded"]

PAIR_COLUMNS = ["structure_a", "chain_a", "structure_b", "chain_b",
                "sequence_identity", "n_body_aligned", "body_rmsd",
                "n_domain_compared", "domain_displacement",
                "domain_centroid_shift", "note", "recorded"]

CONFIDENCE_COLUMNS = ["arm", "region", "n_residues", "mean_confidence",
                      "min_confidence", "max_confidence", "note", "recorded"]


def load_chain(path: Path, chain_name: str | None = None):
    """Alpha carbons of one polymer chain, with the sequence they belong to.

    Keyed by position in the deposited sequence, which is how one structure
    talks about itself. It is not how two structures of the same protein talk
    to each other: these entries use two different numbering conventions, one
    counting from the signal peptide and one from the mature chain, and pairing
    them by number produced twenty-five Angstrom fits between structures that
    differ by a hinge. The sequence is returned so that pairing across entries
    goes through an alignment instead.
    """
    import gemmi

    st = gemmi.read_structure(str(path))
    st.setup_entities()
    if not len(st):
        return "", {}, "", []
    for chain in st[0]:
        if chain_name and chain.name != chain_name:
            continue
        poly = chain.get_polymer()
        if len(poly) < 30:
            continue
        out, letters, keys = {}, [], []
        # A deposited entry carries label_seq. A prediction written as PDB
        # carries none, because the format has no such field, so keying on it
        # returned an empty set for every prediction and this arm compared
        # nothing while reporting success. Where it is absent the position
        # along the chain is used, which for a prediction is the same quantity:
        # the model folds the construct in order, one residue per position.
        residues = list(poly)
        has_label = all(r.label_seq is not None for r in residues)
        for index, res in enumerate(residues, start=1):
            ca = res.find_atom("CA", "*")
            if ca is None:
                continue
            key = int(res.label_seq) if has_label else index
            info = gemmi.find_tabulated_residue(res.name)
            code = info.one_letter_code.upper() if info else "X"
            out[key] = (ca.pos.x, ca.pos.y, ca.pos.z)
            letters.append(code if code.isalpha() else "X")
            keys.append(key)
        if out:
            return chain.name, out, "".join(letters), keys
    return "", {}, "", []


def _cigar_tokens(cigar: str):
    n = ""
    for ch in cigar:
        if ch.isdigit():
            n += ch
        else:
            yield (int(n) if n else 1), ch
            n = ""


def position_map(seq_a: str, keys_a: list, seq_b: str, keys_b: list):
    """Which deposited position in b corresponds to each position in a.

    A global alignment of the two sequences, walked through its own account of
    matches and gaps. This is the same correspondence the null-floor stage
    uses, and for the same reason: residue numbers are a property of an entry,
    not of a protein.
    """
    import gemmi

    result = gemmi.align_string_sequences(list(seq_a), list(seq_b), [],
                                          gemmi.AlignmentScoring())
    mapping, i, j = {}, 0, 0
    for n, op in _cigar_tokens(result.cigar_str()):
        if op == "M":
            for k in range(n):
                if i + k < len(keys_a) and j + k < len(keys_b):
                    mapping[keys_a[i + k]] = keys_b[j + k]
            i += n
            j += n
        elif op == "I":
            i += n
        elif op == "D":
            j += n
    return mapping, result.calculate_identity()


def find_domain(ca: dict, contact_cutoff: float = 12.0, min_len: int = 60,
                max_len: int = 350):
    """The N-terminal structural block of the chain, by contact density.

    The boundary minimises the normalised cut between the block and the rest:

        ncut(c) = x(c) / a(c) + x(c) / (2T - a(c))

    where x(c) is the number of contacts crossing the boundary, a(c) is the
    total degree of the residues before it, and T is the total contact count.
    Both terms are needed. A score that only rewards a block for contacting
    itself is maximised by making the block as large as the search allows, and
    the version this replaces did exactly that: it returned between 53 and 60
    per cent of every chain it was given, which is its own upper bound rather
    than a domain.

    The search is bounded to blocks of domain size. Without a bound the best
    normalised cut in a spike protein is the split between its two subunits,
    which is a real boundary and not the one this arm asks about.
    """
    import numpy as np

    keys = sorted(ca)
    n = len(keys)
    if n < min_len * 2:
        return None
    coords = np.array([ca[k] for k in keys])
    d = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
    contact = (d < contact_cutoff)
    # Neighbours along the chain contact each other whatever the fold, so the
    # band around the diagonal is removed before anything is counted.
    for off in range(-4, 5):
        idx = np.arange(max(0, -off), min(n, n - off))
        contact[idx, idx + off] = False

    deg = contact.sum(axis=1).astype(float)
    cum_deg = np.cumsum(deg)
    total2 = float(deg.sum())          # twice the number of contacts
    # Contacts inside the first c residues, accumulated rather than recomputed,
    # which keeps the sweep linear in the chain length per step.
    inside = np.zeros(n + 1)
    for c in range(1, n + 1):
        inside[c] = inside[c - 1] + float(contact[c - 1, :c - 1].sum())

    best = None
    upper = min(n - min_len, max_len)
    for cut in range(min_len, upper + 1):
        assoc_a = float(cum_deg[cut - 1])
        assoc_b = total2 - assoc_a
        if assoc_a <= 0 or assoc_b <= 0:
            continue
        crossing = assoc_a - 2.0 * inside[cut]
        ncut = crossing / assoc_a + crossing / assoc_b
        if best is None or ncut < best[3]:
            best = (cut, int(inside[cut]), int(crossing), ncut)
    if best is None:
        return None
    cut, internal, external, ncut = best
    return {"first": keys[0], "last": keys[cut - 1], "cut_index": cut,
            "length": cut, "internal": internal, "external": external,
            "ratio": round(ncut, 4), "keys": keys}


def deposited_coverage(path: Path, n_observed: int) -> float:
    """What fraction of the sequence an entry was built from is resolved in it.

    Which entry to measure everything against is decided with this rather than
    alphabetically. These are cryo-EM maps of a spike with a loosely held
    N-terminal domain, and how much of that domain is resolved differs a lot
    between them.
    """
    import gemmi

    st = gemmi.read_structure(str(path))
    st.setup_entities()
    for ent in st.entities:
        if ent.entity_type == gemmi.EntityType.Polymer and ent.full_sequence:
            return n_observed / float(len(ent.full_sequence))
    return 0.0


def rigid_core_fit(mob: dict, ref: dict, positions: list[int],
                   keep_within: float = 2.0, rounds: int = 6):
    """Fit two structures on the part of them that is rigid.

    A fit on every shared residue spreads the error of a moving domain over the
    whole molecule, so nothing looks moved and everything looks slightly wrong.
    Refitting on the residues that agree, repeatedly, converges on the part that
    did not move, and what the moving part then does is visible against it.

    Returns the fit, the positions it was made on, and the deviation of every
    shared position under that fit.
    """
    import numpy as np

    shared = [p for p in positions if p in mob and p in ref]
    if len(shared) < 20:
        return None, [], {}
    use = list(shared)
    fit = None
    for _ in range(rounds):
        fit = superpose_on(mob, ref, use)
        if fit is None:
            return None, [], {}
        rot, cm, cr, _rmsd, _n = fit
        a = np.array([mob[p] for p in shared])
        b = np.array([ref[p] for p in shared])
        dev = np.sqrt((((rot @ (a - cm).T).T + cr - b) ** 2).sum(axis=1))
        nxt = [p for p, d in zip(shared, dev) if d <= keep_within]
        if len(nxt) < 20 or nxt == use:
            use = nxt if len(nxt) >= 20 else use
            break
        use = nxt
    fit = superpose_on(mob, ref, use)
    if fit is None:
        return None, [], {}
    rot, cm, cr, _rmsd, _n = fit
    a = np.array([mob[p] for p in shared])
    b = np.array([ref[p] for p in shared])
    dev = np.sqrt((((rot @ (a - cm).T).T + cr - b) ** 2).sum(axis=1))
    return fit, use, {p: float(d) for p, d in zip(shared, dev)}


def mobile_boundary(loaded: dict, names: list, ref: str,
                    moved_at: float = 4.0, window: int = 15):
    """Where the moving part of the chain stops, in the reference numbering.

    The contact search cannot settle this on its own. Run on these entries it
    returns two answers, around residue 230 in three of them and around 460 in
    the other three, and both are real boundaries between compact blocks. Which
    one this arm is about is decided by the question rather than by the
    contacts: it asks about the domain the deposited record places differently,
    so the domain is the one that moves.

    Every other entry is fitted onto the reference on its rigid part, and the
    deviation of each residue under that fit is kept. A residue counts as
    moving if it deviates by more than moved_at Angstroms in any entry. The
    boundary is the end of the run of moving residues that starts at the N
    terminus, smoothed over a window so that one resolved loop does not end it.
    """
    _c, _ca, _dom, seq_ref, keys_ref = loaded[ref]
    ca_ref = loaded[ref][1]
    worst: dict[int, float] = {}
    for other in names:
        if other == ref:
            continue
        _c2, ca_o, _d2, seq_o, keys_o = loaded[other]
        fwd, _identity = position_map(seq_o, keys_o, seq_ref, keys_ref)
        in_ref = {fwd[k]: ca_o[k] for k in ca_o if k in fwd}
        _fit, _core, dev = rigid_core_fit(in_ref, ca_ref, sorted(in_ref))
        for pos, d in dev.items():
            if d > worst.get(pos, 0.0):
                worst[pos] = d
    if not worst:
        return None, {}, 0.0
    import statistics as stats

    ordered = sorted(worst)
    vals = [worst[p] for p in ordered]
    # The threshold comes from the background rather than from a number chosen
    # in advance, and it is an outlier rule because that is what this is: most
    # of the chain is rigid and a minority of it moves a long way. Taking the
    # median and a spread measured from the median keeps the line where it
    # belongs even though that minority is enormous. Three deviations puts it
    # above the scatter of a rigid fit and far below a domain that travels
    # tens of Angstroms.
    #
    # A fixed four Angstrom rule marks most of the chain as moving and finds no
    # boundary at all; a flat multiple of the median depends on how good the
    # maps are, and on this set it cut 22 residues off the domain it found.
    background = stats.median(vals)
    mad = stats.median([abs(x - background) for x in vals]) or 1.0
    threshold = max(moved_at, background + 3.0 * mad)
    half = window // 2
    smooth = []
    for i in range(len(vals)):
        lo, hi = max(0, i - half), min(len(vals), i + half + 1)
        smooth.append(sum(vals[lo:hi]) / (hi - lo))
    last_index = -1
    for i in range(len(ordered)):
        if smooth[i] > threshold:
            last_index = i
        elif last_index >= 0 and i - last_index >= window:
            break
    if last_index < 0:
        return None, worst, threshold
    return ordered[last_index], worst, threshold


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


def measure(conf: dict, entries: list[str], max_construct: int = 0) -> int:
    data_dir = Path(conf["DATA_DIR"]) / "pedv"
    results_dir = Path(conf["RESULTS_DIR"]) / "pedv"
    results_dir.mkdir(parents=True, exist_ok=True)

    loaded, domain_rows, coverage = {}, [], {}
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
            chain, ca, seq, keys = load_chain(tmp)
            coverage[entry] = deposited_coverage(tmp, len(ca))
        finally:
            tmp.unlink(missing_ok=True)
        if not ca:
            print(f"  {entry}: no usable chain")
            continue
        dom = find_domain(ca)
        if dom is None:
            print(f"  {entry}: chain too short to split")
            continue
        loaded[entry] = (chain, ca, dom, seq, keys)
        domain_rows.append({
            "entry": entry, "chain": chain, "n_residues": len(ca),
            "domain_first": dom["first"], "domain_last": dom["last"],
            "domain_length": dom["length"], "internal_contacts": dom["internal"],
            "external_contacts": dom["external"], "normalised_cut": dom["ratio"],
            "method": "N-terminal block minimising the normalised cut, 12 Angstrom "
                      "contact cutoff, chain neighbours within four positions "
                      "ignored, block length searched between 60 and 350 residues",
            "recorded": L.now_iso(),
        })
        print(f"  {entry}: chain {chain}, {len(ca)} residues, domain "
              f"{dom['first']} to {dom['last']} ({dom['length']} residues), "
              f"normalised cut {dom['ratio']}")

    if not loaded:
        L.eprint("[error] no entry could be read")
        return 1
    L.write_tsv(results_dir / "domain_definition.tsv", DOMAIN_COLUMNS, domain_rows)

    # Every pair, superposed on the body, domain displacement measured after.
    pair_rows = []
    names = sorted(loaded)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            ch_a, ca_a, dom_a, seq_a, keys_a = loaded[a]
            ch_b, ca_b, dom_b, seq_b, keys_b = loaded[b]
            # These entries number their chains two different ways, so which
            # residue of b answers to a residue of a is decided by aligning the
            # two sequences, never by the numbers themselves.
            amap, identity = position_map(seq_a, keys_a, seq_b, keys_b)
            ca_b_in_a = {k: ca_b[v] for k, v in amap.items() if v in ca_b}
            body = [k for k in dom_a["keys"][dom_a["cut_index"]:] if k in ca_b_in_a]
            dom_positions = [k for k in dom_a["keys"][:dom_a["cut_index"]]
                             if k in ca_b_in_a]
            fit = superpose_on(ca_a, ca_b_in_a, body)
            if fit is None:
                continue
            rms, centroid, n_dom = displacement(ca_a, ca_b_in_a, fit, dom_positions)
            pair_rows.append({
                "structure_a": a, "chain_a": ch_a, "structure_b": b, "chain_b": ch_b,
                "sequence_identity": L.fmt(identity, 1),
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

    # One definition of the domain, not six.
    #
    # The boundary is found independently in every entry and the answers differ,
    # because what is resolved differs: the same search returned 201 residues in
    # one entry and 347 in another. Comparing a pair against whichever of the
    # two happened to be first would make the measurement depend on the order of
    # the list. So each entry's boundary is carried onto the reference entry
    # through the same alignment used everywhere else, and the definition is the
    # median of those. The per-entry boundaries stay in domain_definition.tsv,
    # which is what says how stable the definition is.
    # The reference entry is the one that resolves the most of what it was
    # built from, not the one that sorts first. Taking the first by name picked
    # an entry that resolves 73 residues of this domain where others resolve
    # 201, and every number downstream would have been computed on that.
    ref = max(names, key=lambda e: coverage.get(e, 0.0))
    print(f"[pedv] reference entry {ref}, which resolves "
          f"{100.0 * coverage.get(ref, 0.0):.0f} per cent of its deposited "
          f"sequence, the most of the six")
    _chain, ca_ref, dom_ref, seq_ref, keys_ref = loaded[ref]
    mapped = []
    for other in names:
        _c, _ca, dom_o, seq_o, keys_o = loaded[other]
        if other == ref:
            mapped.append(dom_o["last"])
            continue
        back, _identity = position_map(seq_o, keys_o, seq_ref, keys_ref)
        if dom_o["last"] in back:
            mapped.append(back[dom_o["last"]])
    mapped.sort()
    moved_last, worst, move_threshold = mobile_boundary(loaded, names, ref)
    consensus_last = moved_last if moved_last else (
        mapped[len(mapped) // 2] if mapped else dom_ref["last"])
    keys = dom_ref["keys"]
    cut_index = sum(1 for k in keys if k <= consensus_last) or dom_ref["cut_index"]
    print(f"[pedv] contact boundary per entry, on {ref} numbering: "
          f"{', '.join(str(m) for m in mapped)}")
    print(f"[pedv] the part that moves ends at {moved_last}, past a threshold "
          f"of {move_threshold:.1f} Angstroms taken from the background, and "
          f"that is the boundary used")
    L.write_tsv(results_dir / "domain_motion_profile.tsv",
                ["position", "largest_deviation_A", "counts_as_moving"],
                [{"position": k, "largest_deviation_A": L.fmt(worst[k], 2),
                  "counts_as_moving": "yes" if worst[k] > move_threshold else "no"}
                 for k in sorted(worst)])

    # The construct: the domain plus enough of the body to define where it sits.
    # A domain on its own cannot answer the question this arm asks.
    #
    # The cap is wall clock, and it is stated rather than implied. An alignment
    # prediction on this machine costs about 9.4e-3 * L^2 * 11.6 seconds, so the
    # whole 1064-residue chain would be a week and the domain with an equal
    # length of body behind it would be most of a day. The README says what a
    # construct this size cannot tell you, which is anything about how three
    # copies pack against each other, and now also that the body is shorter than
    # the domain because of a time budget.
    # Counted in sequence positions rather than in resolved residues, and cut
    # from the deposited sequence rather than from the observed ones. A
    # construct spliced out of what happened to be modelled is not the protein:
    # this domain is poorly ordered in some of these maps, and joining the
    # resolved fragments would hand the model a sequence that does not exist.
    first_pos = keys[0]
    domain_positions = consensus_last - first_pos + 1
    body_wanted = domain_positions
    if max_construct and domain_positions + body_wanted > max_construct:
        body_wanted = max(0, max_construct - domain_positions)
    last_pos = consensus_last + body_wanted
    if body_wanted == 0:
        print(f"[pedv] WARNING: the cap of {max_construct} residues is at or "
              f"below the {domain_positions}-residue domain, so the construct "
              f"carries no body at all. A domain on its own has nothing to sit "
              f"against, and where the model places it relative to the body "
              f"cannot be measured from this. The prediction step decides what "
              f"to do about that; this records it.")
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
        # The deposited sequence, which includes the residues the map did not
        # resolve. label_seq counts along it from one, so it indexes it directly.
        full = []
        for ent in st.entities:
            if ent.entity_type == gemmi.EntityType.Polymer and ent.full_sequence:
                full = list(ent.full_sequence)
                break
        letters = []
        for pos in range(first_pos, last_pos + 1):
            if pos - 1 < 0 or pos - 1 >= len(full):
                continue
            name = gemmi.Entity.first_mon(full[pos - 1])
            info = gemmi.find_tabulated_residue(name)
            code = info.one_letter_code.upper() if info else "X"
            letters.append(code if code.isalpha() else "X")
    finally:
        tmp.unlink(missing_ok=True)
    seq = "".join(letters)
    seq_path.parent.mkdir(parents=True, exist_ok=True)
    seq_path.write_text(f">pedv_construct\n{seq}\n")
    print(f"[pedv] construct written: {len(seq)} residues, positions "
          f"{first_pos} to {last_pos} of {ref}, the domain plus {body_wanted} "
          f"residues of the body")

    L.write_tsv(results_dir / "construct.tsv",
                ["source_entry", "domain_first", "domain_last", "body_residues",
                 "construct_length", "note", "recorded"],
                [{"source_entry": ref, "domain_first": dom_ref["first"],
                  "domain_last": consensus_last, "body_residues": body_wanted,
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

    import gzip
    import tempfile

    rows, conf_rows = [], []
    for arm_dir in sorted(p for p in data_dir.iterdir() if p.is_dir()):
        arm = arm_dir.name
        models = sorted(arm_dir.glob("*_unrelaxed_rank_001_*.pdb"))
        scores = sorted(arm_dir.glob("*_scores_rank_001_*.json"))
        if not models:
            continue
        chain_p, ca_p, seq_p, keys_p = load_chain(models[0])
        if not ca_p:
            continue
        # The prediction is numbered from one over the construct, and the
        # construct was cut from the deposited numbering, so the offset is
        # known rather than guessed.
        offset = int(construct[0]["domain_first"]) - 1
        ca_pred = {k + offset: v for k, v in ca_p.items()}
        dom_positions = [k for k in ca_pred if k <= int(construct[0]["domain_last"])]
        body_positions = [k for k in ca_pred if k > int(construct[0]["domain_last"])]

        # The prediction is numbered into the source entry's convention. The
        # other entries use their own, and two of these six count from the
        # signal peptide where the rest count from the mature chain, so each
        # one is carried onto the source numbering by an alignment before
        # anything is measured. Pairing by number instead is what made
        # structures of the same protein fit each other at 25 Angstroms.
        src_entry = construct[0]["source_entry"]
        src_path = data_dir / f"{src_entry}.cif.gz"
        seq_src, keys_src = "", []
        if src_path.is_file():
            with tempfile.NamedTemporaryFile("wb", suffix=".cif", delete=False) as fh:
                with gzip.open(src_path, "rb") as gz:
                    fh.write(gz.read())
                tmp_src = Path(fh.name)
            try:
                _c, _ca, seq_src, keys_src = load_chain(tmp_src)
            finally:
                tmp_src.unlink(missing_ok=True)

        for entry in entries:
            path = data_dir / f"{entry}.cif.gz"
            if not path.is_file():
                continue
            with tempfile.NamedTemporaryFile("wb", suffix=".cif", delete=False) as fh:
                with gzip.open(path, "rb") as gz:
                    fh.write(gz.read())
                tmp = Path(fh.name)
            try:
                chain_d, ca_d, seq_d, keys_d = load_chain(tmp)
            finally:
                tmp.unlink(missing_ok=True)
            if not ca_d:
                continue
            identity = 100.0
            if entry != src_entry and seq_src:
                onto_src, identity = position_map(seq_d, keys_d, seq_src, keys_src)
                ca_d = {onto_src[k]: v for k, v in ca_d.items() if k in onto_src}
                if not ca_d:
                    continue
            fit = superpose_on(ca_pred, ca_d, body_positions)
            if fit is None:
                continue
            rms, centroid, n_dom = displacement(ca_pred, ca_d, fit, dom_positions)
            rows.append({
                "structure_a": f"prediction/{arm}", "chain_a": chain_p,
                "structure_b": entry, "chain_b": chain_d,
                "sequence_identity": L.fmt(identity, 1),
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
    ap.add_argument("--max-construct", type=int, default=0,
                    help="longest construct to write; the caller passes what "
                         "the measured memory curve says will fit, and the "
                         "smaller of that and PEDV_MAX_CONSTRUCT is used")
    args = ap.parse_args()
    conf = L.load_conf(args.config)
    entries = [e.strip() for e in args.entries.split(",") if e.strip()]
    if args.measure:
        from_conf = L.conf_int(conf, "PEDV_MAX_CONSTRUCT", 0)
        caps = [c for c in (from_conf, args.max_construct) if c > 0]
        return measure(conf, entries, max_construct=min(caps) if caps else 0)
    if args.compare:
        return compare(conf, entries)
    ap.error("one of --measure or --compare is required")
    return 2


if __name__ == "__main__":
    sys.exit(main())
