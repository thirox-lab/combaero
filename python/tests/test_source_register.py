"""The source register (validation/sources.yaml) stays valid.

Every local reference document has one entry: citation, content, what it was
read for, the verdict and where the repo uses it. The documents are not in
git, so CI checks the register's own consistency -- schema, statuses, and
that every record and code path it names exists. `scripts/sources.py check`
also lists local documents with no entry.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("sources", REPO / "scripts" / "sources.py")
sources = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(sources)


def test_the_register_is_valid() -> None:
    entries = sources.load()
    assert entries, "the register is empty"
    assert sources.problems(entries) == []


def test_the_schema_check_catches_what_it_should() -> None:
    """Falsified: a bad status, an unknown topic, a missing record path and
    an 'implemented' entry naming no code are each reported."""
    good = sources.load()[0]
    bad = [
        {**good, "file": "sources/a.pdf", "status": "read"},
        {**good, "file": "sources/b.pdf", "topic": ["plasma"]},
        {**good, "file": "sources/c.pdf", "records": ["validation/nowhere.md"]},
        {**good, "file": "sources/d.pdf", "status": "implemented", "records": [], "code": []},
    ]
    found = sources.problems(bad)
    assert any("status 'read'" in p for p in found)
    assert any("topic must be" in p for p in found)
    assert any("does not exist" in p for p in found)
    assert any("needs a record or code reference" in p for p in found)


def test_a_file_outside_sources_is_reported() -> None:
    """Falsified (#846): a register entry still pointing into docs/ is an
    error -- copyrighted documents live in the ignored sources/ only."""
    good = sources.load()[0]
    assert good["file"].startswith("sources/")
    found = sources.problems([{**good, "file": "docs/junction/x.pdf"}])
    assert any("must be under sources/" in p for p in found)


def test_local_files_are_listed_under_sources(tmp_path: Path) -> None:
    """Every file below the sources dir is seen (any type, nested), as a
    sources/ path; only registrable types count as documents."""
    root = tmp_path / "sources"
    (root / "junction" / "crops").mkdir(parents=True)
    (root / "junction" / "paper.pdf").write_bytes(b"%PDF")
    (root / "junction" / "crops" / "fig.png").write_bytes(b"png")
    (root / ".gitignore").write_text("*\n")
    assert sources.local_files(root) == [
        "sources/junction/crops/fig.png",
        "sources/junction/paper.pdf",
    ]
    assert sources.local_documents(root) == ["sources/junction/paper.pdf"]


def test_the_wanted_list_is_valid() -> None:
    assert sources.wanted_problems(sources.load_wanted()) == []


def test_the_wanted_check_catches_what_it_should() -> None:
    """Falsified: an unknown priority, a held file, a non-boolean 'checked'
    and string issue numbers are each reported."""
    good = sources.load_wanted()[0]
    bad = [
        {**good, "priority": "urgent"},
        {**good, "file": "sources/x.pdf"},
        {**good, "checked": "no"},
        {**good, "issues": ["468"]},
    ]
    found = sources.wanted_problems(bad)
    assert any("priority 'urgent'" in p for p in found)
    assert any("belongs in 'documents'" in p for p in found)
    assert any("'checked' is true or false" in p for p in found)
    assert any("list of integers" in p for p in found)


def test_a_missing_or_malformed_hash_is_reported() -> None:
    """Falsified (#1057): an entry without sha256, or with one that is not
    64 lowercase hex digits (or loaded as a number), is a schema error -- caught in CI, where
    sources/ is absent."""
    good = sources.load()[0]
    no_hash = {k: v for k, v in good.items() if k != "sha256"}
    found = sources.problems(
        [
            {**no_hash, "file": "sources/a.pdf"},
            {**good, "file": "sources/b.pdf", "sha256": good["sha256"].upper()},
            {**good, "file": "sources/c.pdf", "sha256": good["sha256"][:-1]},
            {**good, "file": "sources/d.pdf", "sha256": 1234},
        ]
    )
    assert sum("sha256 must be" in p for p in found) == 4


def test_the_hash_check_catches_a_changed_file(tmp_path: Path) -> None:
    """Falsified (#1057): one changed byte in a document, and the same
    document under two names, are each reported; an intact file and an
    absent one (reported by `check` as not found locally) are not."""
    root = tmp_path / "sources"
    (root / "x").mkdir(parents=True)
    (root / "x" / "a.pdf").write_bytes(b"%PDF-1.4 original")
    (root / "x" / "b.pdf").write_bytes(b"%PDF-1.4 other")
    good = sources.file_sha256(root / "x" / "a.pdf")
    entries = [
        {"file": "sources/x/a.pdf", "sha256": good},
        {"file": "sources/x/b.pdf", "sha256": sources.file_sha256(root / "x" / "b.pdf")},
        {"file": "sources/x/gone.pdf", "sha256": good},
    ]
    assert sources.hash_problems(entries, root) == []

    (root / "x" / "a.pdf").write_bytes(b"%PDF-1.4 originaL")
    found = sources.hash_problems(entries, root)
    assert len(found) == 1
    assert found[0].startswith("sources/x/a.pdf: sha256 ")

    (root / "x" / "a.pdf").write_bytes(b"%PDF-1.4 other")
    entries[0]["sha256"] = sources.file_sha256(root / "x" / "a.pdf")
    found = sources.hash_problems(entries, root)
    assert found == ["sources/x/b.pdf: same content as sources/x/a.pdf"]


def test_write_hashes_keeps_the_text_and_reads_awkward_names() -> None:
    """Names with spaces, a value wrapped onto a continuation line, and a
    quoted name with an escape all get their hash; comments and every
    other line are kept as they were."""
    text = (
        "# header comment\n"
        "documents:\n"
        "- file: sources/x/with space.pdf\n"
        "  citation: A\n"
        "- file: sources/x/long name wrapped onto\n"
        "    the next line.pdf\n"
        "  citation: B\n"
        '- file: "sources/x/Garci\\u0301a.pdf"\n'
        "  citation: C\n"
        "- file: sources/x/already.pdf\n"
        "  sha256: '" + "0" * 64 + "'\n"
        "  citation: D\n"
    )
    names = [
        "sources/x/with space.pdf",
        "sources/x/long name wrapped onto the next line.pdf",
        "sources/x/Garci" + chr(0x301) + "a.pdf",  # combining acute, as on disk
    ]
    hashes = {n: f"{i:064x}" for i, n in enumerate(names, start=1)}
    new = sources.write_hashes(text, hashes)
    docs = sources.yaml.safe_load(new)["documents"]
    assert [d["sha256"] for d in docs] == [hashes[n] for n in names] + ["0" * 64]
    assert [d["citation"] for d in docs] == ["A", "B", "C", "D"]
    assert [ln for ln in new.splitlines() if "sha256" not in ln] == text.splitlines()[:-3] + [
        "- file: sources/x/already.pdf",
        "  citation: D",
    ]
