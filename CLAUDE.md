# CombAero — Claude Code Instructions

## Practical Tips
- Avoid shell variable expansions and use literal paths or printenv/env instead. Run shell expansions once per session and use the abs path after. Shell variable expansions require user approval per open bugs (#29616, #18160). Allow list workarounds are not always reliable.
- There are allow rules for sandboxed commands such as uv run, that are best executed in sandbox mode. Disabling sandbox requires user approval and should only be used when absolutely required.
- `Edit(*.pyproject.toml)`-style tools are globally denied for **any** `pyproject.toml` (not just the root one) in Claude Code's own permission settings. Edit these files via a shell command (`perl -i -pe` or similar) instead.
- `git push`/`git tag push`/`git pull` inside the sandbox reliably print `fatal: failed to store: 100001` / `could not lock config file .git/config: Operation not permitted` because writing the upstream-tracking config is a denied path. The underlying git operation still succeeds — verify with `git log`/`gh run list` rather than treating the message as a failure.

## Hard Rules
- **Virtual environment only:** use `uv run` and `uv pip`. Never install into system Python or use global `pip`.
- **Never bypass pre-commit hooks** or CI checks. All code must pass cleanly before pushing.
- **No manual edits to auto-generated files** (`docs/UNITS.md`, `include/thermo_transport_data.h`, `include/collision_integral_data.h`, `include/thermo_c2_refit_data.h`). Use the scripts below.
- **API sync is mandatory:** adding/removing/changing any function or property requires updating `include/units_data.h` and the relevant API reference ([docs/API_CPP.md](docs/API_CPP.md) or [docs/API_PYTHON.md](docs/API_PYTHON.md)) in the same commit.
- **CHANGELOG is mandatory -- as a fragment, never an edit to `CHANGELOG.md`:** every user-visible change adds `changelog.d/<id>.<category>.md` in the same commit (`<id>` = issue/PR number or a slug; `<category>` = added / changed / deprecated / removed / fixed / security / documentation; contents = [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) bullets). See [changelog.d/README.md](changelog.d/README.md). Editing `[Unreleased]` directly made every pair of open PRs conflict. `uv run scripts/changelog.py check` validates (pre-commit hook); `preview` shows the release section.
- **Documentation Hygiene:** Active references (API, build, workflow, schemas) live in `docs/`; historical feature reports and design records go to `docs/archive/`. Keep core guides high-signal — a doc points at the file that owns a fact rather than transcribing it. `docs/` is an allowlist (`.gitignore`): only our own `.md` is tracked there.
- **Copyrighted sources live in `sources/`, never in git (#846).** Every third-party document -- papers, handbooks, datasheets, scans, figure crops, supplementary archives, text transcriptions -- goes under `sources/<topic>/`; new arrivals land in `sources/new/`. `sources/` is ignored whatever the file type. Digitised coordinates are our measurements and go to `validation/`; the scans and crops stay in `sources/`. The `no-sources` pre-commit hook and CI (`scripts/check_no_sources.py`) reject any staged `sources/` path and any document/archive/image type outside its allowlist; an image we made needs its directory added there.
- **Unit sync test:** every exported Python symbol needs a `units_data.h` entry, or an `IGNORE_LIST` entry in `python/tests/test_units_sync.py`.
- **Solver (f, J) rule:** solver-facing calculations must expose a C++ PyBind11 API returning `std::tuple<double, double>` (value, derivative). Analytical derivatives via chain rule; no finite differences exposed to Python.
- **Source register:** every reference document under `sources/` has an entry in `validation/sources.yaml` -- citation, content, what it was read for, verdict, and the issues/records/code that use it. Search it before screening or fetching a paper (`scripts/sources.py find`); its `wanted` list holds documents we lack and the question each would settle (`list --wanted`), and an obtained one moves into `documents`; each entry records the file's `sha256`, which `check` verifies (a restore of the off-git backup is proven byte-identical; replacing a document updates its hash in the same PR); a PR that screens, extracts or implements a source updates its entry; `scripts/sources.py check --sources <main checkout>/sources` lists unregistered documents and any file there git would track.
- **Validation policy — fidelity, accuracy and tuning are three questions, never one number.** See [docs/VALIDATION_POLICY.md](docs/VALIDATION_POLICY.md).
  - **Fidelity** (does the implementation mirror the paper?) is judged **only on the author's own data**. A miss is our bug.
  - **Accuracy** (what will it do on another rig?) is judged **cross-source and must be labelled cross-source** in the report, not just in a doc beside it. A miss is the model's limitation. Same lab, different study counts as cross-source.
  - **Matching a specific rig is the user's job.** `Nu_multiplier`/`f_multiplier` exist for that; the harness never scores, recommends or applies a tuner.
  - A stated uncertainty band must be a **measurement** band, never derived from the model's own error on that series — a declared band correlating with measured RMS is circular.
  - Fidelity claimed on **overplotted** data is an upper bound, not an estimate; say how many cleanly-sampled points it rests on.
- **Define Once:** physics constants go as `constexpr` in the relevant public header namespace — never as magic numbers in `.cpp` files.
- **ASCII only:** no non-ASCII characters in any C++ or Python source file.
- **`python/combaero/` must never import from `validation/`, `cantera_validation_tests/`, `thermo_data_generator/`, or `gui/`.** Those trees are dev-only and excluded from the sdist/wheel (see `pyproject.toml`'s `sdist.exclude`); a production module importing from any of them installs fine from a repo checkout but raises `ModuleNotFoundError` on a real `pip install` (hit in the v0.4.0 release for `MultiPortChamberElement`/`ConstantKTeeElement` importing `validation.junction.models.mynard2010` — fixed by moving the implementation into `python/combaero/network/_mynard2010.py`). If a model implementation is needed by production code, it belongs in `python/combaero/`, not the validation tree.

## C++ Style (C++17)
- `#pragma once`, sorted minimal includes, `//` comments only (no `/* */`).
- Smart pointers, RAII. No raw `new`/`delete`.
- Use `math_constants.h` for `M_PI` (MSVC compatibility).
- **Explicit standard library includes.** Include what you use — any platform can be the permissive one, and two agreeing is not evidence. Two distinct strictnesses bite here, in opposite directions:
  - **libstdc++ (Linux, MinGW) is strict about transitive headers.** `std::max`/`std::min`/`std::clamp` need `<algorithm>`, `std::accumulate` needs `<numeric>`. macOS *and* Linux once passed a `std::max({a,b,c})` that only MinGW caught.
  - **MSVC is strict about POSIX-ish macros.** `M_PI` and friends are not defined by default — use `math_constants.h` (there is a pre-commit hook for it).

## Python Style
- `ruff` for lint + format. Run `./scripts/check-python-style.sh --fix` before committing.
- Type annotations on all signatures; `-> None` for void functions; `|` union syntax (Python 3.12).
- No non-ASCII characters, including docstrings and comments.

## Frontend Style (gui/frontend)
- Package manager: `pnpm`. Use `pnpm install`, `pnpm run dev`, etc.
- Linting and formatting: `biome`. Run `./scripts/check-gui-style.sh` before committing.
- Biome config: `biome.json` at repo root.

## Key Scripts
```bash
./scripts/check-source-style.sh             # C++ style
./scripts/check-python-style.sh --fix       # Python lint + format
./scripts/check-gui-style.sh                # Frontend lint + format (biome)
uv pip install -e .                         # rebuild C++ extension
uv run pytest python/tests/                 # Python tests
cd build && ctest                           # C++ tests
uv run scripts/generate_units_md.py         # regenerate docs/UNITS.md
uv run scripts/changelog.py check           # validate changelog.d/ fragments (release: ... release x.y.z)
./scripts/test-pypi-wheel.sh                # verify combaero-gui wheel before PyPI publish
uv run python scripts/stress/run.py --out tmp/stress/x.jsonl   # cold stress sweep; summarize.py --compare for A/B (scripts/README.md)
uv run scripts/sources.py find <terms>      # source register (validation/sources.yaml): do we have it, what was concluded
```

## Auto-Generated Files
| Target | Source | Regenerate with |
|--------|--------|-----------------|
| `docs/UNITS.md` | `include/units_data.h` | `uv run scripts/generate_units_md.py` |
| `include/thermo_transport_data.h` | `thermo_data_generator/` | see `thermo_data_generator/README.md` |
| `include/collision_integral_data.h` | `thermo_data_generator/collision_integrals.py` | `uv run python thermo_data_generator/collision_integrals.py` (repo root, main env) |
| `include/thermo_c2_refit_data.h` | `thermo_data_generator/generate_c2_refit.py` (reads `thermo_transport_data.h`) | `uv run python thermo_data_generator/generate_c2_refit.py` (repo root; rerun after regenerating the thermo header) |

## Git Workflow
- **Branch protection on `main`:** never push directly to `main`. Always work on a feature/fix branch and open a PR (`feature-branch → main`).
- **Opening PRs:** `git` and `gh` are allowlisted in `.claude/settings.local.json` and run inside the sandbox without approval — including `git push`, `gh pr view/checks/merge`, and (verified) most network operations to `github.com`. `dangerouslyDisableSandbox: true` is only needed when a command writes to a denied path (e.g., `git rebase` that touches `uv.lock`).
- **CI is required to merge:** GitHub Actions checks must pass before a PR can be merged.
- **Merge queue:** the repo lives in the `thirox-lab` organisation for its merge queue (ruleset "Iron", squash). A green PR is added to the queue; the queue tests it on top of `main` and the entries ahead of it (unfiltered, the `merge_group` trigger) and merges it. Open PRs need no rebase after another merge -- only when the queue reports a conflict. The queue's squash commit is signed by GitHub.
- **Releases publish through PyPI trusted publishing**, keyed on the repo owner: a change of owner needs the publisher updated on PyPI for both `combaero` and `combaero-gui` before the next tag.

## PyPI Release Process
- **Versioning is fully git-tag-driven (`setuptools_scm`) — there is no version string to hand-edit.** Both `pyproject.toml` (combaero) and `gui/pyproject.toml` (combaero-gui) declare `dynamic = ["version"]`; the version comes from the nearest `v*` git tag. Cutting a release is: `uv run scripts/changelog.py release x.y.z` (folds `changelog.d/` and anything left under `[Unreleased]` into a dated `[x.y.z]` block, moves the compare links, deletes the fragments) → PR → merge to `main` → `git tag -s vx.y.z <commit> -m "vx.y.z" && git push origin vx.y.z`. Use `-s`, not `-a`: an annotated tag carries its own signature, and GitHub reports the TAG's status rather than the signed commit underneath it, so an unsigned annotated tag shows as "Unverified" even on a signed commit (this is what happened to v0.4.0 through v0.5.0). `-s` fails loudly where `-a` would silently produce an unsigned tag.
- **The tag push triggers two independent workflows in parallel:** `publish.yml` (combaero core: wheels for Linux/macOS/Windows via `cibuildwheel` + sdist) and `publish-gui.yml` (combaero-gui: frontend build + wheel + verify smoke test in a clean venv). The core build is slower (3 OSes) — `publish-gui.yml` has to explicitly wait for the matching combaero version to land on PyPI's `/simple/` index (not the JSON API — it updates on a different CDN cache and can falsely report "found" minutes early) before its verify job installs.
- **Whenever a release bumps combaero's minor version, `gui/pyproject.toml`'s `combaero~=x.y` pin must move in lockstep** (edit via shell, see Practical Tips — the Edit tool can't touch it). An unbumped pin resolves a stale combaero from PyPI that may predate classes combaero-gui's backend imports at startup, causing an `ImportError`/`ModuleNotFoundError` the moment the server starts (root cause of the v0.3.1 hotfix and a v0.4.0 near-miss).
- **`combaero-gui`'s publish verify job only installs the built wheel from a clean venv** — it is the one CI path that does *not* run from a repo checkout, so it's the only place that catches a production module importing from an unshipped dev-only tree (see the `validation/` Hard Rule above). `uv run pytest` passing is not sufficient evidence a release is packaging-clean.
- **A pushed git tag should not be moved.** If a release's workflow-only files need a fix after tagging (nothing in the actual package changed), cut the next patch version rather than force-moving/retagging.
- **Yanking a bad PyPI release requires the PyPI web UI** (no API token is available in this environment): project page → the bad version → Options → Yank, with a one-line reason. Yanking doesn't delete anything — pinned installs still work with a warning; unpinned installs skip it.

## Version Sync Checklist
When changing Python version or tool versions, keep in sync:
`.python-version` · `ruff.toml` · `ci.yml` · `validation-tests.yml` · `.pre-commit-config.yaml` · `Makefile`

When changing frontend tool versions, keep in sync:
`gui/frontend/package.json` · `biome.json` (`$schema` version) · `.pre-commit-config.yaml` (biome rev + `@biomejs/biome` pin) · `ci.yml` (`setup-biome` version). `scripts/check_version_sync.py` enforces the biome set (CI + pre-commit).

Dependabot PRs: `scripts/bulk_merge_safe.sh [--dry-run] [--no-audit] [--major] [PR ...]` first runs the security housekeeping (every floor against OSV, `pnpm audit` on the frontend lock; reported, not a gate), then queues the safe ones in the merge queue (majors only with `--major`). PRs that rewrite the same lock file go in waves, one per manifest directory: the next is rebased by Dependabot (`@dependabot rebase` -- never `gh pr update-branch`, whose rebase is unsigned) after the previous merges, re-checked (CI green, version sync at its head) and queued at exactly that commit.

Dependency floors: Dependabot only moves the top of each range, so the nightly `dependency-floors` job (Validation Tests) installs every runtime dependency of combaero and combaero-gui at its declared floor (`scripts/dependency_floors.py`) and runs the fast suite and the examples. Every runtime requirement needs a floor (`test_dependency_floors.py`). A floor is the lowest version with no known advisory that passes the suite. **With every Dependabot round, check both sides:** `uv run scripts/dependency_floors.py --osv` (every floor against OSV; also nightly) and `pnpm audit --prod` plus `pnpm audit` in `gui/frontend` (Dependabot bumps direct packages only; transitive fixes come from `pnpm update --depth Infinity <pkgs>`). A flagged floor moves up to the first fixed release that passes the suite. Raise a floor that breaks rather than skipping it, and declare a transitive dependency directly when its floor is a security fix (starlette, h11).
