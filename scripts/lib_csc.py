#!/usr/bin/env python3
"""Shared helpers for the Python stages. Imported, not run.

Holds what more than one stage needs: reading project.conf, writing tables
atomically, the per-target confidence files, structure I/O through gemmi, and
the interval arithmetic every rate in results/ is reported with.

The layout follows stage 1's lib_vgb.py so that a reader of both repositories
finds the same function in the same place.
"""

from __future__ import annotations

import csv
import datetime
import gzip
import json
import math
import os
import shutil
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# project.conf
# ---------------------------------------------------------------------------


def load_conf(path: str | Path) -> dict:
    """Read the shell project.conf into a dict of strings.

    project.conf is sourced by the bash stages, so it is shell syntax. Only
    KEY=VALUE lines are taken, inline comments are cut, surrounding quotes are
    stripped. Nothing in it needs shell expansion, and one file readable by
    both bash and python is worth more than a second format.
    """
    conf: dict[str, str] = {}
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, _, val = line.partition("=")
        key = key.strip()
        if not key.replace("_", "").isalnum():
            continue
        conf[key] = val.split("#", 1)[0].strip().strip('"').strip("'")
    return conf


def conf_int(conf: dict, key: str, default: int) -> int:
    try:
        return int(conf[key])
    except (KeyError, ValueError):
        return default


def conf_float(conf: dict, key: str, default: float) -> float:
    try:
        return float(conf[key])
    except (KeyError, ValueError):
        return default


def conf_list(conf: dict, key: str) -> list[str]:
    """Comma-separated value to a list, blanks dropped."""
    return [x.strip() for x in conf.get(key, "").split(",") if x.strip()]


def eprint(*a, **k) -> None:
    print(*a, file=sys.stderr, **k)


def now_iso() -> str:
    return datetime.datetime.now().astimezone().isoformat(timespec="seconds")


def today() -> str:
    return datetime.date.today().isoformat()


# ---------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------


def read_tsv(path: str | Path) -> list[dict]:
    """Rows of a TSV as dicts. Lines whose first field starts with # are skipped,
    so the config tables can carry their reasoning as comments."""
    with open(path, newline="") as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def write_tsv(path: str | Path, columns: list[str], rows) -> None:
    """Write a TSV via a temporary file and one rename.

    A stage killed halfway must not leave a truncated results table that the
    next stage would read as complete.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    shutil.move(str(tmp), str(path))


def append_tsv(path: str | Path, columns: list[str], rows) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    exists = path.is_file() and path.stat().st_size > 0
    with path.open("a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        if not exists:
            w.writeheader()
        for r in rows:
            w.writerow(r)


def fmt(x, nd: int = 4):
    """A number for a table cell: rounded, or empty when missing or NaN."""
    if x is None:
        return ""
    try:
        f = float(x)
    except (TypeError, ValueError):
        return x
    if math.isnan(f):
        return ""
    return round(f, nd)


# ---------------------------------------------------------------------------
# Exclusions
# ---------------------------------------------------------------------------
#
# Every target dropped anywhere in the pipeline gets one row here, with the
# stage that dropped it and why. Silent exclusions are how benchmark numbers get
# inflated, and the exclusion count in the README is read from this file rather
# than written by hand. Memory-limited exclusions are the long sequences, not a
# random sample, and the reason column says so rather than leaving a reader to
# work it out.

EXCLUSION_COLUMNS = ["target_id", "stage", "arm", "reason", "detail", "recorded"]


def record_exclusion(results_dir: str | Path, target_id: str, stage: str,
                     reason: str, detail: str = "", arm: str = "") -> None:
    append_tsv(Path(results_dir) / "excluded.tsv", EXCLUSION_COLUMNS, [{
        "target_id": target_id, "stage": stage, "arm": arm, "reason": reason,
        "detail": detail, "recorded": now_iso(),
    }])


def clear_exclusions(results_dir: str | Path, stage: str, arm: str = "") -> int:
    """Drop this stage's earlier rows before it records its new ones.

    The file is append-only, which is right for a stage that runs once and
    wrong for one that can be rebuilt. The set-building stage was rebuilt with
    a larger target count, and its first run's rows stayed: 53 of the 150
    targets in the final set were also listed as excluded, with a reason
    quoting a limit that no longer applied. A table whose purpose is to say
    what was dropped, saying that about targets which were kept and scored, is
    the same failure as a silent exclusion pointed the other way.

    A stage that rewrites its output owns its rows here and clears them first.
    Returns how many were removed, so the caller can say so.
    """
    path = Path(results_dir) / "excluded.tsv"
    if not path.is_file():
        return 0
    rows = read_tsv(path)
    keep = [r for r in rows
            if not (r.get("stage") == stage and (not arm or r.get("arm") == arm))]
    removed = len(rows) - len(keep)
    if removed:
        write_tsv(path, EXCLUSION_COLUMNS, keep)
    return removed


# ---------------------------------------------------------------------------
# Compressed JSON
# ---------------------------------------------------------------------------


def write_json_gz(path: str | Path, obj) -> None:
    """Gzipped JSON, written atomically, with a fixed mtime in the gzip header.

    mtime=0 makes the bytes depend only on the content, so re-running a stage
    that produces the same numbers does not show up as a changed file in git.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    raw = json.dumps(obj, separators=(",", ":"), sort_keys=True).encode()
    with open(tmp, "wb") as fh:
        with gzip.GzipFile(fileobj=fh, mode="wb", mtime=0, compresslevel=9) as gz:
            gz.write(raw)
    shutil.move(str(tmp), str(path))


def read_json_gz(path: str | Path):
    with gzip.open(path, "rt") as fh:
        return json.load(fh)


def gunzip_to(src: str | Path, dst: str | Path) -> None:
    """Decompress a file to a given path, for tools that will not read gzip."""
    dst = Path(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(src, "rb") as fin, open(dst, "wb") as fout:
        shutil.copyfileobj(fin, fout)


def gzip_open_text(path, mode="rt"):
    """gzip.open with a text mode, so a caller can treat plain and compressed
    alignments the same way."""
    return gzip.open(path, mode)


def gzip_file(src: str | Path, dst: str | Path) -> None:
    """Deterministic gzip of a text file (the a3m alignments)."""
    dst = Path(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_suffix(dst.suffix + ".tmp")
    with open(src, "rb") as fin, open(tmp, "wb") as fh:
        with gzip.GzipFile(fileobj=fh, mode="wb", mtime=0, compresslevel=9) as gz:
            shutil.copyfileobj(fin, gz)
    shutil.move(str(tmp), str(dst))


# ---------------------------------------------------------------------------
# Intervals
# ---------------------------------------------------------------------------


def wilson(k: int, n: int, z: float = 1.959963985) -> tuple[float, float]:
    """Wilson score interval for a binomial proportion, as fractions.

    Wilson rather than the normal approximation because several cells here
    are small (the docking subset especially) and the normal interval gives
    impossible bounds below zero on small counts. Same function as stage 1's
    08_analyse.py, so an interval here and one there mean the same thing.
    """
    if n <= 0:
        return (float("nan"), float("nan"))
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))


def quantiles(vals, qs=(0.25, 0.5, 0.75)) -> list[float]:
    """Linear-interpolation quantiles of a list, NaN when empty."""
    v = sorted(float(x) for x in vals if x is not None and not math.isnan(float(x)))
    if not v:
        return [float("nan")] * len(qs)
    out = []
    for q in qs:
        pos = q * (len(v) - 1)
        lo = math.floor(pos)
        hi = math.ceil(pos)
        out.append(v[lo] + (v[hi] - v[lo]) * (pos - lo))
    return out


def cluster_bootstrap(groups: dict, stat, n_boot: int = 2000, seed: int = 0,
                      alpha: float = 0.05) -> tuple[float, float]:
    """Percentile interval of stat() over targets resampled with replacement.

    groups maps target_id to that target's per-residue records. Residues are
    resampled as whole targets, never individually. Residues of one protein
    are not independent observations: resampling them one at a time treats
    300 residues of one helix bundle as 300 pieces of evidence and gives an
    interval several times too narrow. The target is the unit of replication.
    """
    import random

    rng = random.Random(seed)
    keys = sorted(groups)
    if not keys:
        return (float("nan"), float("nan"))
    vals = []
    for _ in range(n_boot):
        pick = [groups[keys[rng.randrange(len(keys))]] for _ in keys]
        v = stat(pick)
        if v is not None and not math.isnan(v):
            vals.append(v)
    if not vals:
        return (float("nan"), float("nan"))
    lo, hi = quantiles(vals, (alpha / 2, 1 - alpha / 2))
    return (lo, hi)


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------


def repo_paths(conf: dict) -> dict[str, Path]:
    keys = ("REPO_DIR", "CONFIG_DIR", "SCRIPTS_DIR", "DATA_DIR", "CACHE_DIR",
            "RESULTS_DIR", "FIGURES_DIR", "LOG_DIR")
    return {k: Path(conf[k]) for k in keys if k in conf}


def du_kb(path: str | Path) -> int:
    total = 0
    for root, _dirs, files in os.walk(path):
        for f in files:
            try:
                total += os.lstat(os.path.join(root, f)).st_size
            except OSError:
                pass
    return total // 1024
