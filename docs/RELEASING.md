# Releasing and Dependency Maintenance

How combaero and combaero-gui are released to PyPI, and how dependency updates are taken. The rules that bite
during a release are summarised in `CLAUDE.md`; this file owns the procedure.

## PyPI release process

- **Versioning is fully git-tag-driven (`setuptools_scm`); there is no version string to hand-edit.** Both
  `pyproject.toml` (combaero) and `gui/pyproject.toml` (combaero-gui) declare `dynamic = ["version"]`; the
  version comes from the nearest `v*` git tag.
- **Cutting a release:**
  1. `uv run scripts/changelog.py release x.y.z` folds `changelog.d/` and anything left under `[Unreleased]`
     into a dated `[x.y.z]` block, moves the compare links and deletes the fragments.
  2. PR, merge to `main`.
  3. `git tag -s vx.y.z <commit> -m "vx.y.z" && git push origin vx.y.z`.
- **Tag with `-s`, not `-a`.** An annotated tag carries its own signature, and GitHub reports the TAG's status
  rather than the signed commit underneath it, so an unsigned annotated tag shows as "Unverified" even on a
  signed commit (v0.4.0 through v0.5.0). `-s` fails loudly where `-a` silently produces an unsigned tag.
- **The tag push triggers two independent workflows in parallel:**
  - `publish.yml` (combaero core): wheels for Linux/macOS/Windows via `cibuildwheel`, plus the sdist.
  - `publish-gui.yml` (combaero-gui): frontend build, wheel, and a verify smoke test in a clean venv. The core
    build is slower (3 OSes), so `publish-gui.yml` waits for the matching combaero version to land on PyPI's
    `/simple/` index -- not the JSON API, which updates on a different CDN cache and can report "found" minutes
    early -- before its verify job installs.
- **Bump the GUI pin with every core minor release.** `gui/pyproject.toml`'s `combaero~=x.y` must move in
  lockstep (edit it via a shell command; the Edit tool cannot touch `pyproject.toml`). An unbumped pin resolves
  a stale combaero from PyPI that may predate classes the GUI backend imports at startup: an `ImportError` the
  moment the server starts (the v0.3.1 hotfix, a v0.4.0 near-miss).
- **The GUI verify job is the only packaging check.** It installs the built wheel in a clean venv, the one CI
  path that does not run from a repo checkout, so it alone catches a production module importing from an
  unshipped dev-only tree. `uv run pytest` passing is not evidence a release is packaging-clean.
- **The release build's compiler is the manylinux wheel's** (gcc-toolset-14, `-O3 -Werror`); CI's Linux C++ job
  builds with g++-14 to match. v0.7.0 never published because of a warning only that compiler raised.
- **Watch the publish.** After a tag push, follow `publish.yml` and check PyPI
  (`https://pypi.org/pypi/<pkg>/<version>/json`) before announcing. If the wheel build fails, nothing is
  published.
- **Never move a pushed tag.** If workflow-only files need a fix after tagging, cut the next patch version.
- **Never test publishing with `workflow_dispatch` on a branch:** `publish.yml` has no `if:` gate and would
  publish.
- **Trusted publishing is keyed on the repo owner.** A change of owner needs the publisher updated on PyPI for
  both `combaero` and `combaero-gui` before the next tag.
- **Yanking a bad release needs the PyPI web UI** (no API token here): project page, the version, Options, Yank,
  with a one-line reason. Pinned installs keep working with a warning; unpinned installs skip it.

## Dependabot

- `scripts/bulk_merge_safe.sh [--dry-run] [--no-audit] [--major] [PR ...]` first runs the security housekeeping
  (every floor against OSV, `pnpm audit` on the frontend lock; reported, not a gate), then queues the safe PRs
  in the merge queue (majors only with `--major`).
- PRs that rewrite the same lock file go in waves, one per manifest directory. The next is rebased by
  Dependabot (`@dependabot rebase` -- never `gh pr update-branch`, whose rebase is unsigned) after the previous
  merges, re-checked (CI green, version sync at its head) and queued at exactly that commit.

## Dependency floors

- Dependabot only moves the top of each range, so the nightly `dependency-floors` job (Validation Tests)
  installs every runtime dependency of combaero and combaero-gui at its declared floor
  (`scripts/dependency_floors.py`) and runs the fast suite and the examples. Every runtime requirement needs a
  floor (`test_dependency_floors.py`).
- A floor is the lowest version with no known advisory that passes the suite.
- **With every Dependabot round, check both sides:** `uv run scripts/dependency_floors.py --osv` (every floor
  against OSV; also nightly), and `pnpm audit --prod` plus `pnpm audit` in `gui/frontend` (Dependabot bumps
  direct packages only; transitive fixes come from `pnpm update --depth Infinity <pkgs>`).
- A flagged floor moves up to the first fixed release that passes the suite. Raise a floor that breaks rather
  than skipping it, and declare a transitive dependency directly when its floor is a security fix (starlette,
  h11).
