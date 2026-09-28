#!/usr/bin/env python3
"""
=============================================================================
 14_check_table_shape.py - every results table still has the shape it claims
=============================================================================
 A tab-separated table is held together by nothing but its tabs. A tool
 message carrying a newline splits one row across two lines; a message
 carrying a tab adds a field; under the default quoting a single double quote
 opens a region that swallows everything to the next one. All three happened
 here. The second one merged two rows of results/excluded.tsv into one, which
 left a scored target listed as dropped and destroyed the record of a real
 exclusion, and no reader noticed because the merged row still had the right
 number of fields.

 Two properties catch all three. Every row has the same number of fields as
 its header, and every row whose table carries a recorded column has a
 timestamp in it. A merged row fails the second even when it passes the first,
 which is what makes the pair worth checking rather than either alone.

 Usage:
     python scripts/14_check_table_shape.py
=============================================================================
"""
import argparse
import re
import sys
from pathlib import Path

_ap = argparse.ArgumentParser(
    description="Check that every results table still has the shape it claims.",
    epilog="Takes no options. --config is accepted and ignored.")
_ap.add_argument("--config", help=argparse.SUPPRESS)
_ap.parse_args()

REPO = Path(__file__).resolve().parent.parent
STAMP = re.compile(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}")

problems = 0
checked = 0
for path in sorted((REPO / "results").rglob("*.tsv")):
    lines = path.read_text(encoding="utf-8", errors="replace").split("\n")
    lines = [ln for ln in lines if ln.strip()]
    if not lines:
        continue
    checked += 1
    header = lines[0].split("\t")
    for n, ln in enumerate(lines[1:], start=2):
        fields = ln.split("\t")
        if len(fields) != len(header):
            problems += 1
            print(f"  {path.relative_to(REPO)}:{n} has {len(fields)} fields "
                  f"where the header has {len(header)}")
            continue
        if "recorded" in header:
            got = fields[header.index("recorded")]
            if not STAMP.match(got):
                problems += 1
                print(f"  {path.relative_to(REPO)}:{n} has {got[:40]!r} in its "
                      f"recorded column, which is not a timestamp")

print(f"{checked} table(s) checked, {problems} malformed row(s)")
sys.exit(1 if problems else 0)
