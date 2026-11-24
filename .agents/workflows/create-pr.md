---
description: Open a PR against main per the write-pr skill, label it, and hand it over for review.
---

## Steps

### 1. Write the body per `.claude/skills/write-pr/SKILL.md`
Physics, correlation, solver and Jacobian changes also follow `.claude/skills/model-provenance/SKILL.md`. Write the body to a file under the worktree's `tmp/`.

### 2. Open the PR against `main`, with labels
```bash
gh pr create --base main --title "type(scope): summary" --body-file tmp/pr_body.md \
  --label <type> --label area:<area> --label agent:<name>
```
For a stacked PR, the base is the PR below it; say so in the body. CI runs only on PRs against `main`, so a stacked PR has no CI until its base merges; say that too.

### 3. Check CI
```bash
gh pr checks <number>
```
A red check is yours to fix or explain before the hand-over.

### 4. Hand over for review
A PR that changes a default, a residual or a Jacobian gets an independent review before the owner queues it (CLAUDE.md). Gemini pings Claude as GEMINI.md describes. Never merge or queue it yourself.
