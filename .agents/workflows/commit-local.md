---
description: Run the checks and make a signed local commit in the task's worktree. Do not push.
---

Follow CLAUDE.md ("Working Practice", "Hard Rules"). This workflow is the command sequence only.

## Steps

### 1. Work in the task's worktree
All commands run in `.claude/worktrees/<name>`, never in the primary checkout. Use the worktree's `.venv` (see CLAUDE.md).

### 2. Format and lint
```bash
./scripts/check-python-style.sh --fix
./scripts/check-source-style.sh
```
Run `./scripts/check-gui-style.sh` too when `gui/` changed.

### 3. Units sync and changelog
```bash
uv run pytest python/tests/test_units_sync.py
uv run scripts/changelog.py check
```
A user-visible change needs a `changelog.d/<id>.<category>.md` fragment. Never edit `CHANGELOG.md` directly.

### 4. Stage the files you changed, by name
```bash
git add <file> <file> ...
```
Never `git add -A`: a worktree can hold build outputs and scratch files.

### 5. Write the message to a file under the worktree's `tmp/` and commit signed
```bash
git commit -S -F tmp/commit_msg.txt
```
The pre-commit hooks run here. Never `--no-verify`. The message ends with your agent's `Co-Authored-By:` trailer (see GEMINI.md or AGENTS.md).
