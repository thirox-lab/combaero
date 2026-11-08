"""The source register: every local reference document, what it holds and
what has been done with it.

The documents themselves live under sources/, which git ignores whatever
the file type (copyright, #846);
validation/sources.yaml records, per file, its citation, what it contains,
which questions it has been read for and with what verdict, and where the
repo uses it. One place to answer "do we already have this?", "what did we
conclude?" and "what else covers this?".

    uv run scripts/sources.py find impingement effusion   # search (all terms)
    uv run scripts/sources.py list --status unscreened    # filter
    uv run scripts/sources.py list --topic film --status implemented
    uv run scripts/sources.py list --wanted               # documents we lack
    uv run scripts/sources.py check                       # schema + local sources
    uv run scripts/sources.py check --sources <main checkout>/sources
    uv run scripts/sources.py new sources/x/new_paper.pdf # print a stub entry
    uv run scripts/sources.py hash --sources <main checkout>/sources [--write]

The file has two lists: `documents` (what is local) and `wanted` (what we do
not hold and the open question each would settle, with a priority). `find`
searches both.

`check` validates the register (schema, statuses, that every record and
code path it names exists, that every `file:` is under sources/) and,
where sources/ is present, lists local documents with no entry, entries
whose file is missing, and any file there git would track, and verifies
each entry's `sha256` against its file (#1057). The documents are backed
up outside git; the hashes in git are how a restore is proven complete and
byte-identical. A document replaced on purpose (a better scan) gets its new
hash in the same PR, so the swap is reviewable. `hash --write` fills
missing hashes; it never overwrites a mismatch.
"""

from __future__ import annotations

import argparse
import hashlib
import re
import sys
from pathlib import Path
from typing import Any

import yaml

ROOT = Path(__file__).resolve().parent.parent
REGISTER = ROOT / "validation" / "sources.yaml"
SOURCES = ROOT / "sources"
PREFIX = "sources/"
# Local files that are neither sources nor a leak.
IGNORED_NAMES = {".gitignore", ".DS_Store"}
DOC_SUFFIXES = {".pdf", ".djvu", ".epub", ".zip"}  # .zip: a paper's supplementary archive

STATUSES = {
    "implemented": "a shipped closure or kernel is built from it",
    "fidelity": "scored as fidelity on its own data",
    "accuracy": "used or scored as a cross-source accuracy set",
    "extracted": "data transcribed or digitised into the repo, not yet scored",
    "screened": "read, verdict recorded, not used yet",
    "unscreened": "nobody has read it for this project yet",
}
TOPICS = {
    "film",
    "effusion",
    "impingement",
    "pin_fin",
    "rib",
    "dimple",
    "channel",
    "jet_mixing",
    "junction",
    "tee",
    "ejector",
    "orifice",
    "combustion",
    "area_change",
    "vortex",
    "acoustics",
    "radiation",
    "transport",
    "review",
    "other",
}
REQUIRED = (
    "file",
    "citation",
    "topic",
    "content",
    "status",
    "verdict",
    "read_for",
    "issues",
    "records",
    "sha256",
)
SHA256 = re.compile(r"[0-9a-f]{64}")
OPTIONAL = ("id", "code", "notes")
WANTED_KEYS = ("citation", "checked", "issues", "priority", "question", "where")
PRIORITIES = {
    "now": "blocks or settles an open item",
    "later": "only when that work resumes",
}


def _register(path: Path) -> dict[str, Any]:
    data = yaml.safe_load(path.read_text())
    if not isinstance(data, dict) or not isinstance(data.get("documents"), list):
        raise SystemExit(f"{path}: expected 'documents' (and optionally 'wanted') lists")
    return data


def load(path: Path = REGISTER) -> list[dict[str, Any]]:
    """The local documents."""
    return _register(path)["documents"]


def load_wanted(path: Path = REGISTER) -> list[dict[str, Any]]:
    """Documents we do not hold, each with the question it would settle."""
    return _register(path).get("wanted") or []


def wanted_problems(wanted: list[dict[str, Any]]) -> list[str]:
    """Schema errors in the wanted list; empty when it is valid."""
    out: list[str] = []
    for i, w in enumerate(wanted):
        where = f"wanted {i} ({str(w.get('citation', '?'))[:40]})"
        if not isinstance(w, dict):
            out.append(f"{where}: not a mapping")
            continue
        out += [f"{where}: missing '{k}'" for k in WANTED_KEYS if k not in w]
        out += [f"{where}: unknown key '{k}'" for k in w if k not in WANTED_KEYS]
        if "file" in w:
            out.append(f"{where}: a held document belongs in 'documents'")
        if w.get("priority") not in PRIORITIES:
            out.append(f"{where}: priority '{w.get('priority')}' not one of {sorted(PRIORITIES)}")
        if not isinstance(w.get("checked"), bool):
            out.append(f"{where}: 'checked' is true or false")
        issues = w.get("issues")
        if not isinstance(issues, list) or any(not isinstance(n, int) for n in issues):
            out.append(f"{where}: issues are a list of integers")
    return out


def problems(entries: list[dict[str, Any]]) -> list[str]:
    """Schema errors in the register; empty when it is valid."""
    out: list[str] = []
    seen: set[str] = set()
    for i, e in enumerate(entries):
        where = f"entry {i} ({e.get('file', '?')})"
        if not isinstance(e, dict):
            out.append(f"{where}: not a mapping")
            continue
        for k in REQUIRED:
            if k not in e:
                out.append(f"{where}: missing '{k}'")
        for k in e:
            if k not in REQUIRED and k not in OPTIONAL:
                out.append(f"{where}: unknown key '{k}'")
        f = str(e.get("file", ""))
        if not f.startswith(PREFIX):
            out.append(f"{where}: file must be under {PREFIX} (#846)")
        if Path(f).suffix.lower() not in DOC_SUFFIXES:
            out.append(f"{where}: file must be a document ({sorted(DOC_SUFFIXES)})")
        if f in seen:
            out.append(f"{where}: duplicate file")
        seen.add(f)
        sha = e.get("sha256")
        if not isinstance(sha, str) or not SHA256.fullmatch(sha):
            out.append(f"{where}: sha256 must be 64 lowercase hex digits")
        if e.get("status") not in STATUSES:
            out.append(f"{where}: status '{e.get('status')}' not one of {sorted(STATUSES)}")
        topics = e.get("topic") or []
        if not isinstance(topics, list) or not topics or any(t not in TOPICS for t in topics):
            out.append(f"{where}: topic must be a non-empty list from {sorted(TOPICS)}")
        for k in ("read_for", "issues", "records"):
            if k in e and not isinstance(e[k], list):
                out.append(f"{where}: '{k}' must be a list")
        if any(not isinstance(n, int) for n in e.get("issues") or []):
            out.append(f"{where}: issues are integers")
        for r in e.get("records") or []:
            if not (ROOT / str(r)).exists():
                out.append(f"{where}: record '{r}' does not exist")
        for c in e.get("code") or []:
            path = str(c).split(":", 1)[0].strip()
            if not (ROOT / path).exists():
                out.append(f"{where}: code path '{path}' does not exist")
        if e.get("status") in ("implemented", "fidelity", "extracted") and not (
            e.get("records") or e.get("code")
        ):
            out.append(f"{where}: status '{e.get('status')}' needs a record or code reference")
    return out


def local_files(sources: Path = SOURCES) -> list[str]:
    """Every file under `sources` (the repo's sources/ by default; a worktree
    has none, so point it at the main checkout's), as paths starting
    sources/. Figure crops, transcriptions and data files included."""
    if not sources.is_dir():
        return []
    return sorted(
        PREFIX + p.relative_to(sources).as_posix()
        for p in sources.rglob("*")
        if p.is_file() and p.name not in IGNORED_NAMES
    )


def local_documents(sources: Path = SOURCES) -> list[str]:
    """The registrable documents (DOC_SUFFIXES) among `local_files`."""
    return [f for f in local_files(sources) if Path(f).suffix.lower() in DOC_SUFFIXES]


def tracked_documents(sources: Path, files: list[str]) -> list[str]:
    """Files under `sources` that git would NOT ignore -- copyrighted
    material that a `git add -A` would stage, or that is already tracked.
    `files` are paths starting sources/; the check runs in the repository
    that contains `sources`."""
    import subprocess

    if not files:
        return []
    repo = sources.resolve().parent
    try:
        r = subprocess.run(
            # -z: NUL-separated both ways. Without it git quotes a non-ASCII
            # path in its output and the name no longer matches. Without
            # --no-index a tracked file is reported as not ignored.
            ["git", "-C", str(repo), "check-ignore", "-z", "--stdin"],
            input="\0".join(files),
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError:
        return []
    if r.returncode not in (0, 1):  # 128: not a repository, nothing to say
        return []
    ignored = set(r.stdout.split("\0"))
    return [f for f in files if f not in ignored]


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def hash_problems(entries: list[dict[str, Any]], sources: Path) -> list[str]:
    """Entries whose local file does not match the recorded sha256, and
    files recorded under two entries (the same document twice). Entries
    whose file is not under `sources` are skipped (`check` reports them)."""
    out: list[str] = []
    by_hash: dict[str, str] = {}
    for e in entries:
        f = str(e.get("file", ""))
        path = sources / f.removeprefix(PREFIX)
        if not path.is_file():
            continue
        actual = file_sha256(path)
        want = e.get("sha256")
        if want and want != actual:
            out.append(f"{f}: sha256 {actual} does not match the register's {want}")
        if actual in by_hash:
            out.append(f"{f}: same content as {by_hash[actual]}")
        by_hash.setdefault(actual, f)
    return out


def write_hashes(text: str, hashes: dict[str, str]) -> str:
    """Insert `sha256:` after the `file:` value of each entry in `hashes`
    (file -> hex), editing the text so the register's comments and quoting
    are kept. The value is read with YAML: names may hold spaces, be
    quoted with escapes, or wrap onto continuation lines."""
    lines = text.splitlines(keepends=True)
    out: list[str] = []
    i = 0
    while i < len(lines):
        out.append(lines[i])
        if not lines[i].startswith("- file:"):
            i += 1
            continue
        j = i + 1  # continuation lines of a wrapped value are indented past the keys
        while j < len(lines) and lines[j].startswith("   "):
            j += 1
        out += lines[i + 1 : j]
        f = yaml.safe_load("".join(lines[i:j]))[0]["file"]
        if f in hashes:
            # Quoted: an all-digit hash would load as an int, one with a
            # single "e" among digits as a float.
            out.append(f"  sha256: '{hashes[f]}'\n")
        i = j
    return "".join(out)


def show(e: dict[str, Any], long: bool = False) -> str:
    head = f"[{e['status']}] {e['file']}\n    {e['citation']}"
    if not long:
        return head + f"\n    -> {e['verdict']}"
    lines = [
        head,
        f"    topic: {', '.join(e['topic'])}",
        f"    content: {e['content']}",
        f"    verdict: {e['verdict']}",
    ]
    for k in ("read_for", "issues", "records", "code"):
        if e.get(k):
            lines.append(f"    {k}: {', '.join(str(v) for v in e[k])}")
    if e.get("notes"):
        lines.append(f"    notes: {e['notes']}")
    return "\n".join(lines)


def show_wanted(w: dict[str, Any]) -> str:
    issues = ", ".join(f"#{n}" for n in w["issues"]) or "no issue"
    unchecked = "" if w["checked"] else " (citation as recalled)"
    return f"[wanted, {w['priority']}] {issues}\n    {w['citation']}{unchecked}\n    -> {w['question']}\n    where: {w['where']}"


def cmd_find(args: argparse.Namespace) -> int:
    terms = [t.lower() for t in args.terms]

    def match(e: dict[str, Any]) -> bool:
        return all(t in yaml.safe_dump(e, allow_unicode=True).lower() for t in terms)

    hits = [e for e in load() if match(e)]
    wants = [w for w in load_wanted() if match(w)]
    for e in hits:
        print(show(e, args.long))
    for w in wants:
        print(show_wanted(w))
    print(f"{len(hits)} match(es), {len(wants)} wanted")
    return 0


def cmd_list(args: argparse.Namespace) -> int:
    if args.wanted:
        rows = [w for w in load_wanted() if not args.priority or w["priority"] == args.priority]
        rows.sort(key=lambda w: list(PRIORITIES).index(w["priority"]))
        for w in rows:
            print(show_wanted(w))
        print(f"{len(rows)} wanted")
        return 0
    rows = [
        e
        for e in load()
        if (not args.status or e["status"] == args.status)
        and (not args.topic or args.topic in e["topic"])
    ]
    for e in rows:
        print(show(e, args.long))
    print(f"{len(rows)} entr{'y' if len(rows) == 1 else 'ies'}")
    return 0


def cmd_check(args: argparse.Namespace) -> int:
    entries = load()
    errs = problems(entries) + wanted_problems(load_wanted())
    for p in errs:
        print(f"error: {p}")
    sources = Path(args.sources) if args.sources else SOURCES
    files = local_files(sources)
    local = local_documents(sources)
    if files:
        listed = {e["file"] for e in entries}
        for f in local:
            if f not in listed:
                print(f"unregistered: {f}")
        for f in sorted(listed - set(local)):
            print(f"not found locally: {f}")
        bad = hash_problems(entries, sources)
        for p in bad:
            print(f"error: {p}")
        errs += bad
        exposed = tracked_documents(sources, files)
        for f in exposed:
            print(f"error: not gitignored (copyrighted, would be committed): {f}")
        errs += [f"not gitignored: {f}" for f in exposed]
    else:
        print(f"no local files under {sources}: local coverage not checked")
    old = sources.parent / "docs"
    if old.is_dir() and any(p.suffix.lower() in DOC_SUFFIXES for p in old.rglob("*")):
        print(
            f"warning: documents still under {old}; move them with "
            "scripts/migrate_docs_to_sources.py (#846)"
        )
    counts: dict[str, int] = {}
    for e in entries:
        counts[e["status"]] = counts.get(e["status"], 0) + 1
    print(f"{len(entries)} entries: " + ", ".join(f"{k} {counts.get(k, 0)}" for k in STATUSES))
    print(f"{len(load_wanted())} wanted")
    return 1 if errs else 0


def cmd_hash(args: argparse.Namespace) -> int:
    sources = Path(args.sources) if args.sources else SOURCES
    missing: dict[str, str] = {}
    absent = 0
    for e in load():
        if e.get("sha256"):
            continue
        path = sources / str(e["file"]).removeprefix(PREFIX)
        if path.is_file():
            missing[e["file"]] = file_sha256(path)
        else:
            absent += 1
            print(f"not found locally: {e['file']}")
    bad = hash_problems(load(), sources)
    for p in bad:
        print(f"error: {p}")
    if args.write and missing:
        REGISTER.write_text(write_hashes(REGISTER.read_text(), missing))
    verb = "written" if args.write else "missing (--write to fill)"
    print(f"{len(missing)} hash(es) {verb}; {absent} entr(ies) without a local file")
    return 1 if bad or absent else 0


def cmd_new(args: argparse.Namespace) -> int:
    sources = Path(args.sources) if args.sources else SOURCES
    path = sources / args.file.removeprefix(PREFIX)
    if not path.is_file():
        print(f"error: {path} not found (pass --sources <main checkout>/sources)")
        return 1
    stub = {
        "file": args.file,
        "citation": "?",
        "topic": ["other"],
        "content": "?",
        "status": "unscreened",
        "verdict": "unscreened",
        "read_for": [],
        "issues": [],
        "records": [],
        "sha256": file_sha256(path),
    }
    print(yaml.safe_dump([stub], sort_keys=False, allow_unicode=True), end="")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = ap.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("find", help="entries matching ALL terms (any field)")
    f.add_argument("terms", nargs="+")
    f.add_argument("--long", action="store_true")
    f.set_defaults(fn=cmd_find)
    ls_ = sub.add_parser("list", help="entries by status and/or topic, or the wanted list")
    ls_.add_argument("--wanted", action="store_true", help="documents we do not hold")
    ls_.add_argument("--priority", choices=sorted(PRIORITIES), help="with --wanted")
    ls_.add_argument("--status", choices=sorted(STATUSES))
    ls_.add_argument("--topic", choices=sorted(TOPICS))
    ls_.add_argument("--long", action="store_true")
    ls_.set_defaults(fn=cmd_list)
    c = sub.add_parser("check", help="validate the register; list unregistered local documents")
    c.add_argument("--sources", help="the sources/ directory to audit (default: this checkout's)")
    c.set_defaults(fn=cmd_check)
    n = sub.add_parser("new", help="print a stub entry for a new document")
    n.add_argument("file")
    n.add_argument("--sources", help="the sources/ directory holding it (default: this checkout's)")
    n.set_defaults(fn=cmd_new)
    h = sub.add_parser("hash", help="report (or --write) missing sha256 values")
    h.add_argument("--sources", help="the sources/ directory to hash (default: this checkout's)")
    h.add_argument("--write", action="store_true", help="fill missing hashes into the register")
    h.set_defaults(fn=cmd_hash)
    args = ap.parse_args()
    return args.fn(args)


if __name__ == "__main__":
    sys.exit(main())
