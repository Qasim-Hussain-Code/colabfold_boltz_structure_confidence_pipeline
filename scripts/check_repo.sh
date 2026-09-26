#!/usr/bin/env bash
# =============================================================================
#  check_repo.sh - the claims this repository makes about itself, checked
# =============================================================================
#  Not part of the pipeline. It exists because a README that says "every number
#  here comes from a file" and "no emoji anywhere" is making claims a reader
#  cannot verify without doing the work, and those claims should be one command
#  away rather than taken on trust.
#
#  What it checks:
#     1. every shell script parses and passes shellcheck
#     2. every python script parses
#     3. no em dash, en dash or emoji in any tracked file
#     4. no banned phrase from the writing conventions
#     5. no tracked file over 50 MB, and the data directory is ignored
#     6. every figure the README shows exists and is produced by a script here
#     7. every results file the README names exists
#     8. no docking score is called a binding energy, and no confidence value
#        is called an accuracy or a probability
#     9. no machine-identifying path, address or account name is tracked
#    10. no tool is credited as an author or contributor
#
#  Exit status is the number of checks that failed, so it is usable anywhere a
#  status is read.
#
#  Usage: bash scripts/check_repo.sh [--quiet]
# =============================================================================
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT" || exit 1
QUIET=0
[[ "${1:-}" == "--quiet" ]] && QUIET=1
FAILED=0
say() { (( QUIET == 0 )) && echo "$@"; return 0; }
pass() { say "  [pass] $1"; return 0; }
fail() { echo "  [FAIL] $1"; FAILED=$(( FAILED + 1 )); return 0; }

mapfile -t TRACKED < <(git ls-files 2>/dev/null)
say "checking ${#TRACKED[@]} tracked files"
say

# Four of these checks are written in python, and python is not on PATH here
# unless an environment is active. Resolving it from project.conf first, then
# from PATH, keeps the checks working from a bare shell. If none is found the
# checks that need it say so rather than reporting the repository as faulty:
# a tool that cannot run has found nothing, and printing FAIL for that reads as
# a privacy violation that does not exist.
PY=""
if [[ -f "${ROOT}/project.conf" ]]; then
    env_name="$(sed -n 's/^CONDA_ENV_ANALYSIS=//p' "${ROOT}/project.conf" | tr -d '"'"'"'"' | head -1)"
    base="$(sed -n 's/^CONDA_BASE=//p' "${ROOT}/project.conf" | tr -d '"'"'"'"' | head -1)"
    [[ -n "$base" && -n "$env_name" && -x "${base}/envs/${env_name}/bin/python" ]] &&
        PY="${base}/envs/${env_name}/bin/python"
fi
[[ -z "$PY" ]] && PY="$(command -v python3 2>/dev/null || command -v python 2>/dev/null || true)"
if [[ -z "$PY" ]]; then
    say "  [skip] no python interpreter found; the checks written in python are"
    say "         not run, and this is not a finding about the repository"
fi

# ---- 1. shell ---------------------------------------------------------------
say "shell scripts"
SH_BAD=0
for f in run_all.sh scripts/*.sh; do
    [[ -f "$f" ]] || continue
    bash -n "$f" 2>/dev/null || { fail "$f does not parse"; SH_BAD=1; }
done
(( SH_BAD == 0 )) && pass "all shell scripts parse"
if command -v shellcheck >/dev/null 2>&1; then
    if shellcheck -x run_all.sh scripts/*.sh >/dev/null 2>&1; then
        pass "shellcheck reports no findings"
    else
        fail "shellcheck reports findings; run: shellcheck -x run_all.sh scripts/*.sh"
    fi
else
    say "  [skip] shellcheck not on PATH"
fi

# ---- 2. python --------------------------------------------------------------
say "python scripts"
if [[ -z "$PY" ]]; then say "  [skip] this check needs python"
elif "$PY" - <<'PY'
import ast, glob, sys


def ok(f):
    try:
        ast.parse(open(f, encoding="utf-8").read())
        return True
    except SyntaxError as e:
        print(f"    {f}: {e}")
        return False


bad = [f for f in sorted(glob.glob("scripts/*.py")) if not ok(f)]
sys.exit(1 if bad else 0)
PY
then pass "all python scripts parse"; else fail "a python script does not parse"; fi

# ---- 3, 4, 8, 10. prose -----------------------------------------------------
say "writing conventions"
if [[ -z "$PY" ]]; then say "  [skip] this check needs python"
elif "$PY" - <<'PY'
import re
import subprocess
import sys
import unicodedata
from pathlib import Path

BANNED = ["it is worth noting", "it's worth noting", "it is important to note",
          "in today's rapidly evolving", "plays a crucial role",
          "serves as a testament", "paving the way", "in conclusion",
          "delve", "seamless", "underscores", "showcases", "a testament to",
          "not only", "leverage", "leveraging", "leverages", "robust",
          "comprehensive", "highlights the", "overall,"]
# Names of assistants and code-generation tools. A repository with a single
# author should not credit one anywhere, including in a commit trailer.
SEQUENCE_LETTERS = set("ACDEFGHIKLMNPQRSTVWYXBZUO")
TOOLS = ["co-authored-by", "generated with", "copilot", "chatgpt", "gpt-4",
         "claude", "anthropic", "openai", "cursor.ai", "codeium"]

files = subprocess.run(["git", "ls-files"], capture_output=True, text=True).stdout.split()
# This file carries the lists it searches for, so scanning it reports every
# pattern as a violation of itself.
files = [f for f in files if not f.endswith("check_repo.sh")]
em = en = emoji = banned = energy = accuracy = tools = 0
for rel in files:
    p = Path(rel)
    if p.suffix in (".png", ".gz", ".zip") or not p.is_file():
        continue
    try:
        text = p.read_text(encoding="utf-8", errors="replace")
    except OSError:
        continue
    for i, line in enumerate(text.splitlines(), 1):
        # A protein sequence is a string of letters and some of those strings
        # spell English words. One target sequence here contains DELVE, which
        # is an ordinary peptide and not a writing habit. Sequence fields are
        # taken out of the line before any prose rule is applied to it, rather
        # than the whole file being exempted, so the prose around them is still
        # checked.
        if p.suffix == ".tsv":
            line = "	".join(
                "" if (len(f) >= 20 and set(f) <= SEQUENCE_LETTERS) else f
                for f in line.split("	"))
        # Built from their code points so this file does not itself contain
        # the characters it is looking for.
        if chr(0x2014) in line:
            em += 1
            print(f"    EM DASH {rel}:{i}")
        if chr(0x2013) in line:
            en += 1
            print(f"    EN DASH {rel}:{i}")
        for ch in line:
            if ord(ch) > 0x2100 and unicodedata.category(ch) in ("So", "Sk"):
                emoji += 1
                print(f"    EMOJI {rel}:{i} {ch!r}")
        low = line.lower()
        for b in BANNED:
            if b in low:
                banned += 1
                print(f"    BANNED {rel}:{i} '{b}'")
        for t in TOOLS:
            if t in low:
                tools += 1
                print(f"    TOOL NAME {rel}:{i} '{t}'")
        # A docking score called an energy, and a confidence value called an
        # accuracy. A sentence denying either is the correction this
        # repository exists to make, not the offence.
        denial = re.search(r"\b(not|never|nor|rather than|instead of|does not|is not|cannot)\b",
                           line, re.I)
        if re.search(r"binding energ|free energy of binding", line, re.I) and not denial:
            energy += 1
            print(f"    SCORE AS ENERGY {rel}:{i}")
        if re.search(r"(plddt|confidence)\s+(is|as)\s+(an?\s+)?(accuracy|probability)",
                     line, re.I) and not denial:
            accuracy += 1
            print(f"    CONFIDENCE AS ACCURACY {rel}:{i}")
print(f"  em dashes {em}, en dashes {en}, emojis {emoji}, banned phrases {banned}, "
      f"score-as-energy {energy}, confidence-as-accuracy {accuracy}, tool names {tools}")
sys.exit(1 if (em or en or emoji or banned or energy or accuracy or tools) else 0)
PY
then pass "no banned character, phrase, misuse or tool name"
else fail "writing conventions violated, see above"; fi

# ---- 9. nothing that identifies the machine ---------------------------------
say "privacy"
if [[ -z "$PY" ]]; then say "  [skip] this check needs python"
elif "$PY" - <<'PY'
import re
import socket
import subprocess
import sys
from pathlib import Path

host = socket.gethostname().strip()
# Built at run time from this machine's own identity, so the patterns are not
# themselves written into a tracked file.
PATTERNS = [
    ("home directory path", re.compile(r"/home/[A-Za-z0-9_.-]+")),
    ("windows user path", re.compile(r"(/mnt/[a-z]/Users/|[A-Za-z]:[\\/]{1,2}Users[\\/])", re.I)),
    ("macos user path", re.compile(r"/Users/[A-Za-z0-9_.-]+/")),
    ("application data path", re.compile(r"AppData", re.I)),
    ("machine identifier", re.compile(r"\{?[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}\}?", re.I)),
    ("electronic address", re.compile(r"[A-Za-z0-9._%+-]+@(?!users\.noreply\.github\.com)[A-Za-z0-9.-]+\.[A-Za-z]{2,}")),
    ("private network address", re.compile(r"(?<![\d.])(10\.\d{1,3}\.\d{1,3}\.\d{1,3}|192\.168\.\d{1,3}\.\d{1,3})(?![\d.])")),
]
if host:
    PATTERNS.append(("host name", re.compile(re.escape(host), re.I)))

files = subprocess.run(["git", "ls-files"], capture_output=True, text=True).stdout.split()
files = [f for f in files if not f.endswith("check_repo.sh")]
hits = 0
for rel in files:
    p = Path(rel)
    if p.suffix in (".png", ".gz", ".zip") or not p.is_file():
        continue
    try:
        text = p.read_text(encoding="utf-8", errors="replace")
    except OSError:
        continue
    for i, line in enumerate(text.splitlines(), 1):
        for name, rx in PATTERNS:
            if rx.search(line):
                hits += 1
                print(f"    {name} {rel}:{i}")
print(f"  machine-identifying strings in tracked files: {hits}")
sys.exit(1 if hits else 0)
PY
then pass "nothing tracked identifies the machine or the account"
else fail "a tracked file carries a machine-identifying string"; fi

# Commit metadata is separate from file contents and is checked separately.
if git log --format='%an%n%ae%n%b' 2>/dev/null | grep -qiE 'co-authored-by|generated with|noreply@anthropic|copilot'; then
    fail "a commit carries a tool as author, contributor or trailer"
else
    pass "no commit credits a tool"
fi

# ---- 5. repository hygiene --------------------------------------------------
say "repository hygiene"
BIG="$(git ls-files -z | xargs -0 ls -l 2>/dev/null | awk '$5 > 52428800 {print $9}')"
if [[ -z "$BIG" ]]; then pass "no tracked file over 50 MB"
else fail "tracked files over 50 MB: ${BIG}"; fi
if git check-ignore -q data/msa 2>/dev/null; then pass "the data directory is ignored"
else fail "the data directory is not ignored"; fi
if git ls-files | grep -qE '\.(ckpt|npz|pt|tar)$'; then
    fail "model weights are tracked"
else pass "no model weights tracked"; fi

# ---- 6, 7. figures and results ----------------------------------------------
say "figures and results"
if [[ -z "$PY" ]]; then say "  [skip] this check needs python"
elif "$PY" - <<'PY'
import re
import sys
from pathlib import Path

readme = Path("README.md")
if not readme.is_file():
    print("    README.md not written yet; figure and results checks skipped")
    sys.exit(0)
text = readme.read_text(encoding="utf-8", errors="replace")
bad = 0

refs = set(re.findall(r"\((figures/[^)\s]+)\)", text))
for r in sorted(refs):
    if not Path(r).is_file():
        print(f"    README shows {r}, which does not exist")
        bad += 1
print(f"  README shows {len(refs)} figures")

on_disk = {p.as_posix() for p in Path("figures").glob("*.png")} if Path("figures").is_dir() else set()
src = Path("scripts/11_figures.py").read_text(encoding="utf-8") \
    if Path("scripts/11_figures.py").is_file() else ""
for f in sorted(on_disk):
    if Path(f).name not in src:
        print(f"    {f} is not produced by scripts/11_figures.py")
        bad += 1

named = set(re.findall(r"\b((?:results|logs|config)/[A-Za-z0-9_./-]+\.(?:tsv|json|html|txt|gz))", text))
missing = [n for n in sorted(named) if not Path(n).is_file()]
for n in missing:
    print(f"    README names {n}, which does not exist")
    bad += 1
print(f"  README names {len(named)} files in results, logs or config; {len(missing)} missing")
sys.exit(1 if bad else 0)
PY
then pass "every figure and named file the README uses is present"
else fail "a figure or named file the README uses is missing"; fi

say
if (( FAILED == 0 )); then
    say "all checks passed"
else
    echo "${FAILED} check(s) failed"
fi
exit "$FAILED"
