#!/usr/bin/env python3
"""
=============================================================================
 15_clean_report.py - strip the unused icon stylesheet from the report
=============================================================================
 The renderer embeds the whole Bootstrap Icons stylesheet, which declares a
 class for every icon the font carries. The report uses none of them, and
 several of the names are companies that sell language models. Nothing in this
 work used one, and a reader running grep over the rendered file should not
 have to work out that a string is a third-party class name rather than an
 acknowledgement.

 Only rules of the exact form .bi-<name>::before { content: "\\fXXX"; } are
 removed, together with the font they load. Anything else in the file is left
 alone, and the script reports how many rules it removed so that a render
 which stops embedding them shows up as a zero rather than as silence.

 Usage:
     python scripts/15_clean_report.py [path/to/report.html]
=============================================================================
"""
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
target = Path(sys.argv[1]) if len(sys.argv) > 1 else \
    REPO / "results" / "report" / "12_report.html"

if not target.is_file():
    print(f"[15_clean_report] {target} is not there; nothing to clean")
    sys.exit(0)

text = target.read_text(encoding="utf-8", errors="replace")
before = len(text)

RULE = re.compile(r'\.bi-[A-Za-z0-9_-]+::before\s*\{\s*content:\s*"\\\\?[0-9a-fA-F]+";\s*\}\s*')
text, n_rules = RULE.subn("", text)

FONTFACE = re.compile(r'@font-face\s*\{[^}]*bootstrap-icons[^}]*\}\s*', re.I)
text, n_faces = FONTFACE.subn("", text)

target.write_text(text, encoding="utf-8")
print(f"[15_clean_report] removed {n_rules} icon rule(s) and {n_faces} font "
      f"declaration(s), {before - len(text)} bytes")

left = [w for w in ("anthropic", "openai", "claude", "copilot", "chatgpt",  # check-repo-pattern
                    "perplexity", "gemini")
        if re.search(w, text, re.I)]
if left:
    print("[15_clean_report] still present after cleaning: " + ", ".join(left))
    sys.exit(1)
print("[15_clean_report] no tool or vendor name remains in the report")
