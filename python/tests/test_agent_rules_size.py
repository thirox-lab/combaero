"""GEMINI.md, with CLAUDE.md inlined, stays under Antigravity's rule-file cap.

Antigravity (Gemini's CLI) expands `@[label](path)` includes in a rule file at
startup and truncates the result from the bottom, on line boundaries, once it
passes 24,000 bytes. GEMINI.md includes CLAUDE.md last (#1331), so growth in
either file would silently cut the tail of CLAUDE.md -- its git workflow and
Dependabot rules -- from what Gemini sees, with no error anywhere. This test
fails first, while there is still margin.
"""

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
#: Antigravity's per-file rule cap after include expansion (its rules spec).
ANTIGRAVITY_RULE_CAP_BYTES = 24_000
#: Fail this far below the cap, so a growing CLAUDE.md is caught before any
#: line is cut rather than at the first one.
MARGIN_BYTES = 2_000
_INCLUDE = re.compile(r"^@\[[^\]]*\]\(([^)]+)\)\s*$", re.MULTILINE)


def _expanded(path: Path, seen: tuple[Path, ...] = ()) -> str:
    """The file with its line-level @[label](path) includes expanded, as
    Antigravity loads it (paths relative to the including file)."""
    assert path not in seen, f"include cycle through {path}"
    text = path.read_text(encoding="utf-8")

    def sub(m: re.Match[str]) -> str:
        return _expanded((path.parent / m.group(1)).resolve(), (*seen, path))

    return _INCLUDE.sub(sub, text)


def test_gemini_md_includes_claude_md() -> None:
    gemini = (ROOT / "GEMINI.md").read_text(encoding="utf-8")
    assert [m.group(1) for m in _INCLUDE.finditer(gemini)] == ["CLAUDE.md"]


def test_expanded_gemini_md_stays_under_the_cap() -> None:
    size = len(_expanded(ROOT / "GEMINI.md").encode("utf-8"))
    limit = ANTIGRAVITY_RULE_CAP_BYTES - MARGIN_BYTES
    assert size <= limit, (
        f"GEMINI.md with CLAUDE.md inlined is {size} bytes, over {limit} "
        f"({ANTIGRAVITY_RULE_CAP_BYTES} cap minus {MARGIN_BYTES} margin): "
        "Antigravity would soon truncate CLAUDE.md's tail for Gemini. Shorten "
        "CLAUDE.md (point at the file that owns a fact instead of restating it) "
        "or split the include."
    )
