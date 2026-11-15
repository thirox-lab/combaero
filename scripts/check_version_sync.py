"""Check that the Biome version is pinned in one place only.

Dependabot bumps gui/frontend/package.json and pnpm-lock.yaml and nothing
else, so every other consumer reads the version from there: CI's setup-biome
from the lockfile, the pre-commit hook from gui/frontend/node_modules
(scripts/biome-hook.sh), and biome.json's $schema from the installed package.
This check fails when a second pin comes back. Run in CI, in pre-commit and by
scripts/bulk_merge_safe.sh before it merges.

    uv run scripts/check_version_sync.py            # working tree
    uv run scripts/check_version_sync.py --ref REF  # any git ref, e.g. a PR head
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from collections.abc import Callable
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LOCAL_SCHEMA = "./gui/frontend/node_modules/@biomejs/biome/configuration_schema.json"


def reader(ref: str | None) -> Callable[[str], str]:
    if ref is None:
        return lambda path: (ROOT / path).read_text()
    return lambda path: (
        subprocess.run(
            ["git", "show", f"{ref}:{path}"], cwd=ROOT, check=True, capture_output=True, text=True
        ).stdout
    )


def biome_problems(read: Callable[[str], str]) -> tuple[str | None, list[str]]:
    """The package.json Biome version, and every second pin found elsewhere."""
    problems: list[str] = []

    pkg = json.loads(read("gui/frontend/package.json"))
    deps = {**pkg.get("dependencies", {}), **pkg.get("devDependencies", {})}
    version = deps.get("@biomejs/biome")
    if not version:
        problems.append("gui/frontend/package.json: no @biomejs/biome dependency")
    elif not re.fullmatch(r"[0-9]+\.[0-9]+\.[0-9]+", version):
        problems.append(f"gui/frontend/package.json: @biomejs/biome {version!r} is not an exact pin")

    schema = json.loads(read("biome.json")).get("$schema")
    if schema != LOCAL_SCHEMA:
        problems.append(f"biome.json: $schema {schema!r} is not {LOCAL_SCHEMA!r}")

    precommit = read(".pre-commit-config.yaml")
    if re.search(r"@biomejs/biome@", precommit) or "biomejs/pre-commit" in precommit:
        problems.append(".pre-commit-config.yaml: pins Biome; use scripts/biome-hook.sh")

    ci = read(".github/workflows/ci.yml")
    step_rest = r"((?:(?![ \t]*-[ \t])[ \t]+[^\n]*\n)*)"  # lines up to the next step
    for m in re.finditer(r"uses:\s*biomejs/setup-biome@[^\n]*\n" + step_rest, ci):
        if re.search(r"^\s*version:", m.group(1), re.MULTILINE):
            problems.append(".github/workflows/ci.yml: setup-biome has a version input")

    return version, problems


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", help="check this git ref instead of the working tree")
    args = ap.parse_args()
    version, problems = biome_problems(reader(args.ref))
    if not problems:
        print(f"version sync OK: biome {version}, pinned in gui/frontend/package.json only")
        return 0
    print("version sync FAILED: Biome must be pinned in gui/frontend/package.json only:")
    for p in problems:
        print(f"  {p}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
