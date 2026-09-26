#!/usr/bin/env python3
"""Build the held-out prediction set and the docking subset from the RCSB.

The set is built here rather than downloaded, because no published held-out set
is held out from every model compared here, and because the filters are the
argument. Anyone can rerun this and get a set built by the same rules against
the archive as it stands on the day; the counts will differ and the manifest
records the date so the difference is visible.

Three filter decisions carry the argument, and each is a place where a careless
version produces a flattering number.

Deposition date, not release date
---------------------------------
A structure released in 2024 may have been deposited in 2021 and held. If a
model selected its training examples by deposition date, a release-date split
leaks: the entry was in the archive the trainers drew from, under embargo but
present. The filter here requires both dates to fall after the cutoff, which is
the stricter of the two readings, and both are recorded per target. The summary
counts how many entries a release-date-only filter would have admitted that
this one rejects, because that count is the size of the mistake and it belongs
in the README rather than in a footnote.

One cutoff for every model
--------------------------
The models were frozen at different dates, read from config/model_cutoffs.tsv
with the source of each. The common cutoff is the latest of them, so that one
set is out of training for every model in the comparison. That costs set size,
and the cost is reported: the summary says how many targets a per-model cutoff
would have allowed for the older model.

Sequence identity to the training set cannot be checked
-------------------------------------------------------
The training sets are not published. What can be checked is whether anything
resembling the target was in the public archive before the cutoff, which is a
proxy and is named as one everywhere it appears. Two proxies are used together:
the RCSB's own sequence-identity clusters, asking whether any member of the
target's cluster was released before the cutoff, and a direct sequence search
against pre-cutoff entries. Neither excludes leakage. A model trained on
predicted structures, as both of these were, may have seen the target's
sequence without any structure of it existing.

Usage
-----
    python scripts/02_build_holdout_set.py --config project.conf
    python scripts/02_build_holdout_set.py --config project.conf --limit 20
    python scripts/02_build_holdout_set.py --config project.conf --force
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_csc as L  # noqa: E402

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
SCHEMA_URL = "https://search.rcsb.org/rcsbsearch/v2/metadata/schema"
DATA_ENTRY = "https://data.rcsb.org/rest/v1/core/entry/{}"
DATA_ENTITY = "https://data.rcsb.org/rest/v1/core/polymer_entity/{}/{}"
GRAPHQL_URL = "https://data.rcsb.org/graphql"

TARGET_COLUMNS = [
    "target_id", "pdb_id", "entity_id", "auth_asym_id", "sequence_length",
    "resolution", "deposit_date", "initial_release_date", "method",
    "uniprot", "organism", "title", "cluster_id_30", "cluster_size_30",
    "pre_cutoff_cluster_members", "closest_pre_cutoff_identity",
    "closest_pre_cutoff_entry", "release_only_would_admit", "sequence",
    "recorded",
]

LIGAND_COLUMNS = [
    "target_id", "pdb_id", "comp_id", "label_asym_id", "auth_asym_id",
    "formula_weight", "heavy_atoms", "rscc", "rsr", "completeness",
    "is_subject_of_investigation", "covalent", "status", "reason", "recorded",
]


# ---------------------------------------------------------------------------
# A small client. Cached, throttled, and loud about failure.
# ---------------------------------------------------------------------------
class Rcsb:
    """Requests to the RCSB, cached on disk and rate limited.

    The cache is what makes this stage rerunnable without hammering a public
    service: a second run of the same query reads the same bytes off disk and
    makes no request at all. It is keyed by the exact request body, so changing
    a filter changes the key and refetches.

    The delay is deliberate and is not tuned for speed. The RCSB asks for a
    handful of requests per second and returns 429 when it has had enough; this
    client stays below that and backs off when told to.
    """

    def __init__(self, cache_dir: Path, delay: float = 0.34, retries: int = 5):
        self.cache = Path(cache_dir)
        self.cache.mkdir(parents=True, exist_ok=True)
        self.delay = delay
        self.retries = retries
        self.n_requests = 0
        self.n_cached = 0
        self._last = 0.0

    def _throttle(self) -> None:
        wait = self.delay - (time.time() - self._last)
        if wait > 0:
            time.sleep(wait)
        self._last = time.time()

    def _key(self, tag: str, payload: str) -> Path:
        h = hashlib.sha256(payload.encode()).hexdigest()[:24]
        return self.cache / f"{tag}_{h}.json"

    def get(self, url: str, tag: str = "get"):
        return self._request(url, None, tag)

    def post(self, url: str, body: dict, tag: str = "post"):
        return self._request(url, body, tag)

    def _request(self, url: str, body: dict | None, tag: str):
        """Returns the decoded body, {} when the service ran and matched
        nothing, or None when the request itself failed.

        The three cases have to stay distinct. A search that matches nothing
        answers 204 with an empty body, and if that is folded into the failure
        case then "this target has no pre-cutoff relative" and "the query broke"
        become the same row in the results table.
        """
        payload = url if body is None else url + json.dumps(body, sort_keys=True)
        path = self._key(tag, payload)
        if path.is_file():
            self.n_cached += 1
            txt = path.read_text()
            return {} if txt == "" else json.loads(txt)

        data = None if body is None else json.dumps(body).encode()
        headers = {"Accept": "application/json"}
        if data is not None:
            headers["Content-Type"] = "application/json"

        for attempt in range(self.retries):
            self._throttle()
            req = urllib.request.Request(url, data=data, headers=headers)
            try:
                with urllib.request.urlopen(req, timeout=120) as r:
                    self.n_requests += 1
                    # A search with no hits answers 204 with an empty body,
                    # including when counts were asked for. Treating that as an
                    # error would turn "nothing matched" into a crash, and
                    # treating it as a zero-length list would hide it.
                    if r.status == 204:
                        path.write_text("")
                        return {}
                    out = json.loads(r.read())
                path.write_text(json.dumps(out))
                return out
            except urllib.error.HTTPError as exc:
                if exc.code == 204:
                    path.write_text("")
                    return {}
                if exc.code in (429, 500, 502, 503, 504) and attempt < self.retries - 1:
                    back = 2.0 * (2 ** attempt)
                    L.eprint(f"  [rcsb] HTTP {exc.code}; waiting {back:.0f} s")
                    time.sleep(back)
                    continue
                L.eprint(f"  [rcsb] HTTP {exc.code} on {url}: {exc.read()[:200]!r}")
                return None
            except (urllib.error.URLError, TimeoutError) as exc:
                if attempt < self.retries - 1:
                    time.sleep(2.0 * (2 ** attempt))
                    continue
                L.eprint(f"  [rcsb] {exc!r}")
                return None
        return None


# ---------------------------------------------------------------------------
# Query construction
# ---------------------------------------------------------------------------
def term(attribute: str, operator: str, value):
    return {"type": "terminal", "service": "text",
            "parameters": {"attribute": attribute, "operator": operator, "value": value}}


def group(op: str, nodes: list):
    return {"type": "group", "logical_operator": op, "nodes": nodes}


def holdout_query(conf: dict, cutoff: str) -> dict:
    """Entities that no model in this comparison can have been trained on.

    Both dates are required to fall after the cutoff. The date range is written
    with explicit inclusive bounds because the range operator excludes both
    bounds by default, which is the kind of default that moves a set size by a
    few entries and is never noticed.
    """
    res_max = L.conf_float(conf, "SET_MIN_RESOLUTION", 2.0)
    len_min = L.conf_int(conf, "SET_MIN_LENGTH", 80)
    len_max = L.conf_int(conf, "SET_MAX_LENGTH", 300)
    n_prot = L.conf_int(conf, "SET_MAX_PROTEIN_ENTITIES", 1)
    nodes = [
        term("exptl.method", "exact_match", "X-RAY DIFFRACTION"),
        term("rcsb_accession_info.initial_release_date", "greater", cutoff),
        term("rcsb_accession_info.deposit_date", "greater", cutoff),
        term("rcsb_entry_info.resolution_combined", "less_or_equal", res_max),
        term("rcsb_entry_info.polymer_entity_count_protein", "equals", n_prot),
        term("rcsb_entry_info.polymer_entity_count_nucleic_acid", "equals", 0),
        term("entity_poly.rcsb_entity_polymer_type", "exact_match", "Protein"),
        term("entity_poly.rcsb_sample_sequence_length", "range",
             {"from": len_min, "to": len_max, "include_lower": True, "include_upper": True}),
    ]
    return {
        "query": group("and", nodes),
        "return_type": "polymer_entity",
        "request_options": {
            "return_all_hits": True,
            "results_verbosity": "compact",
            # One representative per 30 per cent identity cluster, so a
            # well-studied fold cannot supply twenty targets and dominate the
            # calibration. The representative is the highest-resolution member,
            # chosen by a stated rule rather than by whichever came first.
            "group_by": {
                "aggregation_method": "sequence_identity",
                "similarity_cutoff": L.conf_int(conf, "SET_CLUSTER_IDENTITY", 30),
                "ranking_criteria_type": {
                    "sort_by": "rcsb_entry_info.resolution_combined",
                    "direction": "asc",
                },
            },
            "group_by_return_type": "representatives",
        },
    }


def release_only_query(conf: dict, cutoff: str) -> dict:
    """The same filter with the deposition-date condition removed.

    The difference between this count and the count above is how much a
    release-date split would have leaked. It is the number the README quotes.
    """
    q = holdout_query(conf, cutoff)
    q["query"]["nodes"] = [n for n in q["query"]["nodes"]
                           if n["parameters"]["attribute"] != "rcsb_accession_info.deposit_date"]
    return q


def cluster_has_pre_cutoff(client: Rcsb, group_id: str, cutoff: str) -> int | None:
    """How many entities in this identity cluster were released before the cutoff.

    The two conditions on the cluster membership go in their own group. The
    attribute is one the schema marks as nested, so conditions placed in the
    same group are evaluated against the same membership record; scattered
    across the outer group they could be satisfied by two different records and
    the answer would be quietly wrong.
    """
    q = {
        "query": group("and", [
            group("and", [
                term("rcsb_polymer_entity_group_membership.aggregation_method",
                     "exact_match", "sequence_identity"),
                term("rcsb_polymer_entity_group_membership.group_id", "exact_match", group_id),
            ]),
            term("rcsb_accession_info.initial_release_date", "less", cutoff),
        ]),
        "return_type": "polymer_entity",
        "request_options": {"return_counts": True},
    }
    out = client.post(SEARCH_URL, q, tag="cluster")
    if out is None:
        return None          # the query failed; not the same as an empty answer
    return int(out.get("total_count", 0))


def closest_pre_cutoff(client: Rcsb, sequence: str, cutoff: str):
    """The most similar pre-cutoff entity, by sequence search.

    This is the second leakage proxy and it answers a different question from
    the cluster one. Clusters are built at fixed identity thresholds with a
    coverage requirement, so a target can sit in a cluster of its own and still
    have a close relative in the archive. The sequence service finds that
    relative.

    It is also the query that builds the template null floor: the highest
    identity pre-cutoff chain is exactly what a homology modeller would have
    copied, and copying it is the comparator every prediction here is measured
    against.
    """
    q = {
        "query": group("and", [
            {"type": "terminal", "service": "sequence",
             "parameters": {"evalue_cutoff": 0.1, "identity_cutoff": 0.0,
                            "sequence_type": "protein", "value": sequence}},
            term("rcsb_accession_info.initial_release_date", "less", cutoff),
        ]),
        "return_type": "polymer_entity",
        "request_options": {"return_all_hits": True, "results_verbosity": "verbose"},
    }
    out = client.post(SEARCH_URL, q, tag="seqsim")
    if out is None:
        return None, None, 0          # failed
    if not out:
        return "", 0.0, 0             # ran, and nothing in the archive matched
    best_id, best_identity, n = None, 0.0, 0
    for hit in out.get("result_set", []):
        n += 1
        for svc in hit.get("services", []):
            # The key is service_type. The documentation calls the block
            # match_context in one place and matching_context in another, and
            # names the service field "service"; the response uses
            # service_type and match_context. Reading the wrong key silently
            # yields no identities at all rather than an error, which is how
            # this first ran: every target came back with a closest pre-cutoff
            # identity of zero, and one of them has a 98.8 per cent match.
            if svc.get("service_type") != "sequence":
                continue
            for node in svc.get("nodes", []):
                for m in node.get("match_context", []):
                    ident = float(m.get("sequence_identity", 0.0) or 0.0)
                    if ident > best_identity:
                        best_identity = ident
                        best_id = hit.get("identifier")
    return best_id, best_identity, n


# ---------------------------------------------------------------------------
def entity_record(client: Rcsb, pdb_id: str, entity_id: str, identity_cutoff: int = 30):
    e = client.get(DATA_ENTRY.format(pdb_id), tag="entry")
    p = client.get(DATA_ENTITY.format(pdb_id, entity_id), tag="entity")
    if e is None or p is None:
        return None
    acc = e.get("rcsb_accession_info", {})
    info = e.get("rcsb_entry_info", {})
    res = info.get("resolution_combined") or []
    poly = p.get("entity_poly", {})
    ids = p.get("rcsb_polymer_entity_container_identifiers", {})
    src = (p.get("rcsb_entity_source_organism") or [{}])[0]
    # similarity_cutoff comes back as a float, so 30 arrives as 30.0 and a
    # lookup keyed on the string "30" finds nothing. Normalising to an integer
    # here keeps the rest of the code reading the way it is written.
    clusters = {}
    for c in (p.get("rcsb_polymer_entity_group_membership") or []):
        if c.get("aggregation_method") != "sequence_identity":
            continue
        try:
            clusters[int(float(c.get("similarity_cutoff")))] = c.get("group_id")
        except (TypeError, ValueError):
            continue
    return {
        "pdb_id": pdb_id,
        "entity_id": entity_id,
        "auth_asym_id": (ids.get("auth_asym_ids") or [""])[0],
        "sequence": poly.get("pdbx_seq_one_letter_code_can", "") or "",
        "sequence_length": poly.get("rcsb_sample_sequence_length"),
        "resolution": res[0] if res else None,
        "deposit_date": (acc.get("deposit_date") or "")[:10],
        "initial_release_date": (acc.get("initial_release_date") or "")[:10],
        "method": (e.get("exptl") or [{}])[0].get("method", ""),
        "uniprot": ",".join(ids.get("uniprot_ids") or []),
        "organism": src.get("ncbi_scientific_name", ""),
        "title": (e.get("struct", {}).get("title", "") or "")[:150],
        "cluster_id_30": clusters.get(identity_cutoff, ""),
    }


LIGAND_GQL = """
query($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    nonpolymer_entities {
      rcsb_nonpolymer_entity_container_identifiers { entity_id nonpolymer_comp_id asym_ids auth_asym_ids }
      nonpolymer_comp {
        chem_comp { id type formula_weight }
        rcsb_chem_comp_info { atom_count_heavy }
      }
      nonpolymer_entity_instances {
        rcsb_nonpolymer_entity_instance_container_identifiers { auth_asym_id asym_id comp_id }
        rcsb_nonpolymer_instance_validation_score {
          RSCC RSR completeness is_subject_of_investigation is_best_instance ranking_model_fit
        }
        rcsb_nonpolymer_instance_annotation { type }
      }
    }
  }
}
"""


def ligand_records(client: Rcsb, pdb_ids: list[str], conf: dict):
    """Ligand instances per entry, with the validation figures the filters need.

    The instance completeness is not a searchable attribute, so the selection
    cannot be done in one query however tempting that is. Each candidate entry
    is fetched and filtered here, which is slower and correct.
    """
    w_min = L.conf_float(conf, "LIG_MIN_WEIGHT", 100.0)
    w_max = L.conf_float(conf, "LIG_MAX_WEIGHT", 900.0)
    ha_min = L.conf_int(conf, "LIG_MIN_HEAVY_ATOMS", 10)
    rscc_min = L.conf_float(conf, "LIG_MIN_RSCC", 0.95)
    rsr_max = L.conf_float(conf, "LIG_MAX_RSR", 0.2)

    out = []
    for i in range(0, len(pdb_ids), 25):
        chunk = pdb_ids[i:i + 25]
        res = client.post(GRAPHQL_URL, {"query": LIGAND_GQL, "variables": {"ids": chunk}},
                          tag="ligands")
        if res is None or "data" not in res:
            continue
        for entry in (res["data"].get("entries") or []):
            pdb_id = entry["rcsb_id"]
            for ent in (entry.get("nonpolymer_entities") or []):
                comp = (ent.get("nonpolymer_comp") or {}).get("chem_comp") or {}
                cinfo = (ent.get("nonpolymer_comp") or {}).get("rcsb_chem_comp_info") or {}
                comp_id = comp.get("id", "")
                weight = comp.get("formula_weight")
                heavy = cinfo.get("atom_count_heavy")
                for inst in (ent.get("nonpolymer_entity_instances") or []):
                    ids = inst.get("rcsb_nonpolymer_entity_instance_container_identifiers") or {}
                    scores = inst.get("rcsb_nonpolymer_instance_validation_score") or []
                    ann = {a.get("type") for a in (inst.get("rcsb_nonpolymer_instance_annotation") or [])}
                    sc = scores[0] if scores else {}
                    row = {
                        "pdb_id": pdb_id, "comp_id": comp_id,
                        "label_asym_id": ids.get("asym_id", ""),
                        "auth_asym_id": ids.get("auth_asym_id", ""),
                        "formula_weight": weight, "heavy_atoms": heavy,
                        "rscc": sc.get("RSCC"), "rsr": sc.get("RSR"),
                        "completeness": sc.get("completeness"),
                        "is_subject_of_investigation": sc.get("is_subject_of_investigation", ""),
                        "covalent": int("HAS_COVALENT_LINKAGE" in ann),
                        "status": "ok", "reason": "",
                        "recorded": L.now_iso(),
                    }
                    reasons = []
                    if weight is None or not (w_min <= float(weight) <= w_max):
                        reasons.append(f"formula weight {weight} outside {w_min:g} to {w_max:g}")
                    if heavy is None or int(heavy) < ha_min:
                        reasons.append(f"{heavy} heavy atoms, below {ha_min}")
                    if "HAS_COVALENT_LINKAGE" in ann:
                        reasons.append("covalently bound to the protein")
                    if sc.get("is_subject_of_investigation") != "Y":
                        reasons.append("not flagged as the subject of investigation")
                    if sc.get("RSCC") is None or float(sc["RSCC"]) < rscc_min:
                        reasons.append(f"RSCC {sc.get('RSCC')} below {rscc_min}")
                    if sc.get("RSR") is not None and float(sc["RSR"]) > rsr_max:
                        reasons.append(f"RSR {sc.get('RSR')} above {rsr_max}")
                    if sc.get("completeness") is not None and float(sc["completeness"]) < 1.0:
                        reasons.append(f"modelled completeness {sc.get('completeness')}")
                    if reasons:
                        row["status"] = "rejected"
                        row["reason"] = "; ".join(reasons)
                    out.append(row)
    return out


# ---------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--limit", type=int, default=0,
                    help="stop after this many candidate entities (smoke test)")
    ap.add_argument("--max-candidates", type=int, default=0,
                    help="examine at most this many candidates in the first pass; "
                         "0 examines every cluster representative the filters return")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    conf = L.load_conf(args.config)
    config_dir = Path(conf["CONFIG_DIR"])
    results_dir = Path(conf["RESULTS_DIR"])
    data_dir = Path(conf["DATA_DIR"])
    out_targets = config_dir / "targets.tsv"
    if out_targets.is_file() and not args.force:
        print(f"[02_build_holdout_set] {out_targets.name} exists; skipping (--force to rebuild).")
        return 0

    cutoffs = L.read_tsv(config_dir / "model_cutoffs.tsv")
    if not cutoffs:
        L.eprint("[error] config/model_cutoffs.tsv is empty")
        return 1
    per_model = {c["model"]: c["cutoff_date"] for c in cutoffs}
    common_cutoff = max(per_model.values())
    earliest_cutoff = min(per_model.values())
    print(f"[02_build_holdout_set] per-model cutoffs: "
          + ", ".join(f"{m} {d}" for m, d in sorted(per_model.items())))
    print(f"[02_build_holdout_set] common cutoff {common_cutoff}, the latest of them, so one"
          f" set is out of training for every model")

    client = Rcsb(data_dir / "cache" / "rcsb")
    # The search schema carries its version in a comment field rather than in a
    # version field. It is recorded because attribute names and their operators
    # come from this schema, and a query written against one version can fail
    # or, worse, quietly mean something else against another.
    schema = client.get(SCHEMA_URL, tag="schema")
    schema_version = "unknown"
    comment = (schema or {}).get("$comment", "")
    if "version" in comment.lower():
        schema_version = comment.split(":")[-1].strip()
    print(f"[02_build_holdout_set] search schema version {schema_version}")

    # --- the candidate set ---------------------------------------------------
    q = holdout_query(conf, common_cutoff)
    res = client.post(SEARCH_URL, q, tag="holdout")
    if res is None:
        L.eprint("[error] the held-out query returned no hits at all; check the filters")
        return 1
    # Compact verbosity returns bare identifier strings; the other verbosities
    # return objects with an identifier field. Both shapes are accepted so that
    # changing the verbosity of a query does not silently return nothing.
    hits = [h if isinstance(h, str) else h.get("identifier")
            for h in res.get("result_set", [])]
    hits = [h for h in hits if h]
    total_before_grouping = res.get("total_count", len(hits))
    n_groups = res.get("group_by_count", len(hits))
    print(f"[02_build_holdout_set] {total_before_grouping} entities pass the filters, "
          f"in {n_groups} clusters at {L.conf_int(conf, 'SET_CLUSTER_IDENTITY', 30)} per cent identity")
    if args.limit:
        hits = hits[:args.limit]
    if args.max_candidates and len(hits) > args.max_candidates:
        # Examining every representative costs two archive requests each. When
        # the set wanted is far smaller than the pool, a bounded prefix of the
        # pool is examined instead, and the summary records both counts so a
        # reader knows the set was drawn from a prefix rather than from the
        # whole pool.
        print(f"[02_build_holdout_set] examining the first {args.max_candidates} "
              f"of {len(hits)} cluster representatives")
        hits = hits[:args.max_candidates]

    # How many a release-date-only filter would have admitted. The difference is
    # the leak, and it is a number rather than a warning.
    rel = client.post(SEARCH_URL, release_only_query(conf, common_cutoff), tag="release_only")
    n_release_only = (rel or {}).get("group_by_count", 0)

    # --- pass one: what each candidate is ------------------------------------
    # Only the two cheap record fetches here. The leakage proxies cost a
    # sequence search each, and running one per candidate would mean more than
    # a thousand searches against a public service to build a set of a
    # hundred. They are run in pass three, on the targets actually chosen.
    rows, excluded = [], []
    for i, ident in enumerate(hits, 1):
        pdb_id, _, entity_id = ident.partition("_")
        rec = entity_record(client, pdb_id, entity_id,
                            L.conf_int(conf, "SET_CLUSTER_IDENTITY", 30))
        if rec is None:
            excluded.append((ident, "the data API returned nothing for this entity"))
            continue
        if not rec["sequence"]:
            excluded.append((ident, "no one-letter sequence in the entity record"))
            continue
        # A canonical sequence can carry characters for modified residues that
        # neither model accepts. Anything outside the twenty standard letters is
        # dropped here rather than failing inside an inference run an hour later.
        odd = set(rec["sequence"].upper()) - set("ACDEFGHIKLMNPQRSTVWY")
        if odd:
            excluded.append((ident, f"sequence contains non-standard letters: {''.join(sorted(odd))}"))
            continue
        rec.update({
            "target_id": ident,
            "cluster_size_30": "",
            "pre_cutoff_cluster_members": "",
            "closest_pre_cutoff_identity": "",
            "closest_pre_cutoff_entry": "",
            "release_only_would_admit": "",
            "recorded": L.now_iso(),
        })
        rows.append(rec)
        if i % 50 == 0 or i == len(hits):
            print(f"  {i}/{len(hits)} candidates examined, {len(rows)} usable, "
                  f"{client.n_requests} requests, {client.n_cached} from cache", flush=True)

    # --- the memory ceiling --------------------------------------------------
    # Targets too long for this machine are refused here rather than at
    # prediction time, with the projected peak and the shortfall, so the
    # exclusion is a row in a table rather than a process killed by the kernel.
    max_len = L.conf_int(conf, "EFFECTIVE_MAX_SEQ_LEN", L.conf_int(conf, "MAX_SEQ_LEN", 300))
    mem_base = L.conf_float(conf, "MEM_BASE_MB", 2600.0)
    mem_quad = L.conf_float(conf, "MEM_QUAD_MB_PER_KRES2", 34000.0)
    ram_mb = L.conf_int(conf, "RAM_GB", 7) * 1024
    too_long = [r for r in rows if int(r["sequence_length"] or 0) > max_len]
    for r in too_long:
        n = int(r["sequence_length"])
        peak = mem_base + mem_quad * (n / 1000.0) ** 2
        excluded.append((r["target_id"],
                         f"{n} residues exceeds the {max_len} this machine supports; "
                         f"projected peak {peak:.0f} MB against {ram_mb} MB available, "
                         f"short by {(peak - ram_mb) / 1024:.1f} GB"))
    rows = [r for r in rows if int(r["sequence_length"] or 0) <= max_len]

    # --- the time budget -----------------------------------------------------
    # If the set is larger than the wall clock allows, the set is cut rather
    # than an arm dropped: the arms are the experiment. The targets kept are the
    # shortest, which is not a random sample, and the summary says so.
    budget = L.conf_int(conf, "TARGET_BUDGET", 0)
    cap = L.conf_int(conf, "MAX_TARGETS", 0)
    limit = min([x for x in (budget, cap) if x > 0], default=0)
    cut_for_time = []
    if limit and len(rows) > limit:
        rows.sort(key=lambda r: int(r["sequence_length"] or 0))
        cut_for_time = rows[limit:]
        rows = rows[:limit]
        for r in cut_for_time:
            excluded.append((r["target_id"],
                             f"beyond the {limit} targets the wall-clock budget allows; "
                             f"the set was cut by length, keeping the shortest"))

    # --- pass three: the leakage proxies, on the chosen set ------------------
    # A date filter says the models cannot have been trained on these entries.
    # It does not say they have never seen these proteins, because most newly
    # released structures are new determinations of something already in the
    # archive. These two queries measure how much of that there is. Neither
    # removes a target: a target with a close pre-cutoff relative is exactly
    # the case worth measuring, and dropping it would remove the evidence.
    print(f"[02_build_holdout_set] checking {len(rows)} chosen targets against the "
          f"archive as it stood before {common_cutoff}")
    for i, rec in enumerate(rows, 1):
        if rec["cluster_id_30"]:
            n_pre = cluster_has_pre_cutoff(client, rec["cluster_id_30"], common_cutoff)
            rec["pre_cutoff_cluster_members"] = "" if n_pre is None else n_pre
        best_id, best_identity, _ = closest_pre_cutoff(client, rec["sequence"], common_cutoff)
        rec["closest_pre_cutoff_identity"] = ("" if best_identity is None
                                              else L.fmt(best_identity, 4))
        rec["closest_pre_cutoff_entry"] = best_id or ""
        if i % 25 == 0 or i == len(rows):
            print(f"  {i}/{len(rows)} checked, {client.n_requests} requests, "
                  f"{client.n_cached} from cache", flush=True)

    rows.sort(key=lambda r: r["target_id"])
    L.write_tsv(out_targets, TARGET_COLUMNS, rows)

    # --- the docking subset --------------------------------------------------
    lig_rows = ligand_records(client, [r["pdb_id"] for r in rows], conf)
    by_target = {r["pdb_id"]: r["target_id"] for r in rows}
    for lr in lig_rows:
        lr["target_id"] = by_target.get(lr["pdb_id"], "")
    kept_lig = [lr for lr in lig_rows if lr["status"] == "ok"]
    # One ligand per target, the best-scoring instance, so a target with four
    # copies of the same compound does not count four times.
    best: dict[str, dict] = {}
    for lr in kept_lig:
        cur = best.get(lr["target_id"])
        if cur is None or float(lr["rscc"] or 0) > float(cur["rscc"] or 0):
            best[lr["target_id"]] = lr
    L.write_tsv(config_dir / "docking_subset.tsv", LIGAND_COLUMNS,
                sorted(best.values(), key=lambda r: r["target_id"]))
    L.write_tsv(results_dir / "ligand_candidates.tsv", LIGAND_COLUMNS,
                sorted(lig_rows, key=lambda r: (r["pdb_id"], r["comp_id"])))

    # --- exclusions and the summary -----------------------------------------
    # This stage has just rewritten config/targets.tsv, so its earlier rows
    # describe a set that no longer exists. Leaving them turned 53 targets that
    # are in the final set into targets the table said had been dropped.
    stale = L.clear_exclusions(results_dir, "02_build_holdout_set")
    if stale:
        print(f"[02_build_holdout_set] cleared {stale} exclusion rows from an "
              f"earlier build of the set")
    for tid, why in excluded:
        L.record_exclusion(results_dir, tid, "02_build_holdout_set", why)

    def band(lo, hi):
        """Targets whose closest pre-cutoff relative falls in [lo, hi).

        lo of None counts the targets where the search ran and matched nothing.
        Targets whose search failed are counted in neither and are visible as
        blanks in the manifest.
        """
        n = 0
        for r in rows:
            v = r["closest_pre_cutoff_identity"]
            if v == "":
                continue
            v = float(v)
            if lo is None:
                if v < hi:
                    n += 1
            elif lo <= v < hi:
                n += 1
        return n

    summary = [
        ("per_model_cutoffs", "; ".join(f"{m}={d}" for m, d in sorted(per_model.items())),
         "read from config/model_cutoffs.tsv with the source of each"),
        ("common_cutoff", common_cutoff,
         "the latest per-model cutoff, so one set is out of training for every model"),
        ("earliest_cutoff", earliest_cutoff,
         "a per-model set for the older model could have started here instead"),
        ("entities_passing_filters", total_before_grouping, "before clustering"),
        ("clusters_at_identity", n_groups,
         f"at {L.conf_int(conf, 'SET_CLUSTER_IDENTITY', 30)} per cent identity, one representative each"),
        ("candidates_examined", len(hits),
         "cluster representatives whose records were fetched; the rest of the pool was not examined"),
        ("release_date_only_clusters", n_release_only,
         "what the same filter admits when the deposition-date condition is removed"),
        ("release_date_only_extra", max(0, int(n_release_only) - int(n_groups)),
         "entries a release-date split would have admitted that were deposited before the cutoff"),
        ("targets_after_sequence_checks", len(rows) + len(too_long) + len(cut_for_time),
         "before the memory ceiling and the time budget"),
        ("excluded_for_length", len(too_long),
         f"longer than {max_len} residues, which this machine cannot hold"),
        ("excluded_for_time", len(cut_for_time),
         "beyond what the wall-clock budget allows; the set was cut, not the arms"),
        ("targets_final", len(rows), "the set every arm is run on"),
        # The distribution of pre-cutoff identity, which is the result this
        # stage exists to produce. A date filter selects entries the models
        # cannot have been trained on; it does not select proteins they have
        # never seen, because most newly released structures are new
        # determinations of something already in the archive.
        ("identity_none", band(None, 0.0001),
         "no pre-cutoff relative found by sequence search"),
        ("identity_under_30", band(0.0001, 0.30),
         "a pre-cutoff relative below 30 per cent identity"),
        ("identity_30_to_70", band(0.30, 0.70), "30 to 70 per cent"),
        ("identity_70_to_95", band(0.70, 0.95), "70 to 95 per cent"),
        ("identity_95_to_100", band(0.95, 1.0001),
         "95 per cent or better, which is the same protein by any practical reading"),
        ("identity_exactly_100", sum(1 for r in rows
                                     if str(r["closest_pre_cutoff_identity"]) == "1.0"),
         "an identical sequence was already in the archive before the cutoff"),
        ("docking_subset", len(best), "targets with a ligand passing every quality filter"),
        ("ligand_instances_considered", len(lig_rows), "before the ligand filters"),
        ("search_schema_version", schema_version, "the RCSB search schema this set was built against"),
        ("query_date", L.today(), "clusters are recomputed weekly, so this date matters"),
        ("api_requests", client.n_requests, "requests made; the rest came from the local cache"),
        ("api_cached", client.n_cached, "requests served from the local cache"),
    ]
    L.write_tsv(results_dir / "holdout_set_summary.tsv",
                ["key", "value", "note"],
                [{"key": k, "value": v, "note": n} for k, v, n in summary])

    # The leakage proxy, reported rather than acted on. A target whose cluster
    # has pre-cutoff members is not removed: it is exactly the case where a
    # model may have seen a relative, and dropping it would remove the evidence
    # instead of measuring it.
    with_relatives = sum(1 for r in rows if int(r["pre_cutoff_cluster_members"] or 0) > 0)
    high_identity = sum(1 for r in rows
                        if float(r["closest_pre_cutoff_identity"] or 0) >= 0.3)
    print(f"[02_build_holdout_set] {len(rows)} targets, {len(best)} with a usable ligand")
    print(f"[02_build_holdout_set] {with_relatives} targets sit in an identity cluster that "
          f"already had a member before the cutoff")
    print(f"[02_build_holdout_set] {high_identity} have a pre-cutoff relative at 30 per cent "
          f"identity or better; these are kept and flagged, not removed")
    print(f"[02_build_holdout_set] a release-date-only filter would have admitted "
          f"{max(0, int(n_release_only) - int(n_groups))} clusters that were deposited before the cutoff")
    print(f"[02_build_holdout_set] {client.n_requests} API requests, {client.n_cached} from cache")
    return 0


if __name__ == "__main__":
    sys.exit(main())
