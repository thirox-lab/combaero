# CombAero — Claude Code Instructions

## Practical Tips
- Avoid shell variable expansions and use literal paths or printenv/env instead. Run shell expansions once per session and use the abs path after. Shell variable expansions require user approval per open bugs (#29616, #18160). Allow list workarounds are not always reliable.
- There are allow rules for sandboxed commands such as uv run, that are best executed in sandbox mode. Disabling sandbox requires user approval and should only be used when absolutely required.
- `Edit(*.pyproject.toml)`-style tools are globally denied for **any** `pyproject.toml` (not just the root one) in Claude Code's own permission settings. Edit these files via a shell command (`perl -i -pe` or similar) instead.
- `git push`/`git tag push`/`git pull` inside the sandbox reliably print `fatal: failed to store: 100001` / `could not lock config file .git/config: Operation not permitted` because writing the upstream-tracking config is a denied path. The underlying git operation still succeeds — verify with `git log`/`gh run list` rather than treating the message as a failure.
- `gh` inside the sandbox sometimes fails with TLS/x509 errors. Try it sandboxed first; on a TLS error re-run that command with the sandbox disabled. Worktree add/remove and pushes from a worktree often need the same.

## Hard Rules
- **Virtual environment only:** use `uv run` and `uv pip`. Never install into system Python or use global `pip`.
- **Never bypass pre-commit hooks** or CI checks. All code must pass cleanly before pushing.
- **No manual edits to auto-generated files** (`docs/UNITS.md`, `include/thermo_transport_data.h`, `include/collision_integral_data.h`, `include/thermo_c2_refit_data.h`). Use the scripts below.
- **API sync is mandatory:** adding/removing/changing any function or property requires updating `include/units_data.h` and the relevant API reference ([docs/API_CPP.md](docs/API_CPP.md) or [docs/API_PYTHON.md](docs/API_PYTHON.md)) in the same commit.
- **CHANGELOG is mandatory -- as a fragment, never an edit to `CHANGELOG.md`:** every user-visible change adds `changelog.d/<id>.<category>.md` in the same commit (`<id>` = issue/PR number or a slug; `<category>` = added / changed / deprecated / removed / fixed / security / documentation; contents = [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) bullets). See [changelog.d/README.md](changelog.d/README.md). An unreleased fragment is corrected in place, not contradicted by a second one. `uv run scripts/changelog.py check` validates (pre-commit hook); `preview` shows the release section.
- **Documentation Hygiene:** Active references (API, build, workflow, releasing, schemas) live in `docs/`; historical feature reports and design records go to `docs/archive/`. Keep core guides high-signal — a doc points at the file that owns a fact rather than transcribing it. `docs/` is an allowlist (`.gitignore`): only our own `.md` is tracked there.
- **Copyrighted sources live in `sources/`, never in git (#846).** Every third-party document -- papers, handbooks, datasheets, scans, figure crops, supplementary archives, text transcriptions -- goes under `sources/<topic>/`; new arrivals land in `sources/new/`. `sources/` is ignored whatever the file type. Digitised coordinates are our measurements and go to `validation/`; the scans and crops stay in `sources/`. The `no-sources` pre-commit hook and CI (`scripts/check_no_sources.py`) reject any staged `sources/` path and any document/archive/image type outside its allowlist; an image we made needs its directory added there.
- **Unit sync test:** every exported Python symbol needs a `units_data.h` entry, or an `IGNORE_LIST` entry in `python/tests/test_units_sync.py`.
- **Solver (f, J) rule:** solver-facing calculations must expose a C++ PyBind11 API returning `std::tuple<double, double>` (value, derivative). Analytical derivatives via chain rule; no finite differences exposed to Python.
- **Source register:** every reference document under `sources/` has an entry in `validation/sources.yaml` -- citation, content, what it was read for, verdict, and the issues/records/code that use it. Search it before screening or fetching a paper (`scripts/sources.py find`); its `wanted` list holds documents we lack and the question each would settle (`list --wanted`), and an obtained one moves into `documents`; each entry records the file's `sha256`, which `check` verifies (a restore of the off-git backup is proven byte-identical; replacing a document updates its hash in the same PR); a PR that screens, extracts or implements a source updates its entry; `scripts/sources.py check --sources <main checkout>/sources` lists unregistered documents and any file there git would track. Use and transcribe the sources we hold before acquiring new ones.
- **Validation policy — fidelity, accuracy and tuning are three questions, never one number.** See [docs/VALIDATION_POLICY.md](docs/VALIDATION_POLICY.md).
  - **Fidelity** (does the implementation mirror the paper?) is judged **only on the author's own data**. A miss is our bug.
  - **Accuracy** (what will it do on another rig?) is judged **cross-source and must be labelled cross-source** in the report, not just in a doc beside it. A miss is the model's limitation. Same lab, different study counts as cross-source.
  - **Matching a specific rig is the user's job.** `Nu_multiplier`/`f_multiplier` exist for that; the harness never scores, recommends or applies a tuner.
  - A stated uncertainty band must be a **measurement** band, never derived from the model's own error on that series — a declared band correlating with measured RMS is circular.
  - Fidelity claimed on **overplotted** data is an upper bound, not an estimate; say how many cleanly-sampled points it rests on.
- **Define Once:** physics constants go as `constexpr` in the relevant public header namespace — never as magic numbers in `.cpp` files.
- **ASCII only:** no non-ASCII characters in any C++ or Python source file.
- **`python/combaero/` must never import from `validation/`, `cantera_validation_tests/`, `thermo_data_generator/`, or `gui/`.** Those trees are excluded from the sdist/wheel, so such an import works from a checkout and fails on a real `pip install` (the v0.4.0 release). Production model code belongs in `python/combaero/`.

## Architecture
- **Physics in C++, glue in Python.** Kernels are deep: real physics behind a narrow (value, derivative) interface. The Python element layer composes and relays; it does not re-derive physics.
- **One composition point per correlation family.** A family's terms (e.g. hole Cd: base correlation, inclination, inlet and exit crossflow) compose in one chain that the elements consume, never inside each element.
- **A family's options travel as one spec object, like `GasRadiation`,** never as loose keyword arguments on each element. A new option is a field on the spec, so every element, builder and the GUI carry it without edits. Diagnostics follow the same rule: one prefixed group per family with a stable shape, not new keys per feature.
- **Opt-in is not an end state.** A new model option states its route to becoming the default (the evidence it needs) or to removal. A user has no evidence to choose an option with; the model should apply a term wherever it is in range.
- **Extrapolation is information, not a warning.** Prefer physics-informed forms with stable, low-slope continuation, and report leaving a source's envelope as a diagnostic naming the variable. Refuse only genuinely invalid input.
- **Before adding a parameter, diagnostic key or runner argument,** check whether an existing structure should absorb it. If none can cleanly, stop and propose the refactor before building on top.

## Working Practice
- **One worktree per task** under `.claude/worktrees/<name>`, branched from `origin/main`, with its own venv: `env -u VIRTUAL_ENV uv venv -q .venv && env -u VIRTUAL_ENV uv sync -q --all-extras --all-groups --python .venv/bin/python` (without `env -u VIRTUAL_ENV` installs leak into the main checkout's venv). Use absolute paths; never `cd` the primary checkout onto another branch.
- **Commits are signed:** write the message to a file under the worktree's `tmp/` and run `git commit -S -F <file>`; rebase with `git rebase --gpg-sign`. Never `--no-verify`, never a bare `git stash` (the stash stack is shared across worktrees), never `git reset --hard`.
- **PRs** follow the `write-pr` skill; physics, correlation and solver changes also follow `model-provenance`. A PR that changes a default, a residual or a Jacobian gets an independent review before it is queued.
- **Digitisation and transcription:** automated pass, then a private review page showing the scan beside the extracted values, then the user's approval, then the data PR. Digits read partly from the expected trend (overstruck, faint) keep a per-value confidence flag and are scored separately; printed misprints stay as printed, flagged, and are excluded by default.

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
- **CI is required to merge:** GitHub Actions checks must pass before a PR can be merged.
- **Merge queue:** the repo lives in the `thirox-lab` organisation for its merge queue (ruleset "Iron", squash). A green PR is added to the queue; the queue tests it on top of `main` and the entries ahead of it (the `merge_group` trigger) and merges it; the squash commit is signed by GitHub. A queued branch cannot be pushed to. A textual rebase is needed only when GitHub reports a conflict, but a clean merge is not proof: when another PR changed pinned values or shared tests, re-run the tests on top of current `main` before queueing.

## Releases and Dependencies
The procedure lives in [docs/RELEASING.md](docs/RELEASING.md); version-sync lists in [docs/DEVELOPMENT.md](docs/DEVELOPMENT.md). The rules that bite:
- Versions come from git tags (`setuptools_scm`); release via `uv run scripts/changelog.py release x.y.z`, PR, merge, then `git tag -s` (never `-a`) and push the tag.
- Every core minor bump moves `gui/pyproject.toml`'s `combaero~=x.y` pin in the same release.
- After a tag push, watch `publish.yml` and check PyPI before announcing; never move a pushed tag (cut the next patch instead).
- Dependabot rounds go through `scripts/bulk_merge_safe.sh`, with `@dependabot rebase` (never `gh pr update-branch`), and always check the floors (`scripts/dependency_floors.py --osv`) and `pnpm audit`.
