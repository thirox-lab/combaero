# CombAero -- Gemini (Antigravity) Instructions

**[CLAUDE.md](CLAUDE.md) is authoritative for every agent. Read it first, in full, and follow it.** Its hard rules, architecture, working practice, style and git workflow all apply to you unchanged. This file lists only what differs for Gemini. It does not restate CLAUDE.md, so the two cannot drift apart. If a rule here seems to conflict with CLAUDE.md, CLAUDE.md wins; raise the conflict with the owner.

CLAUDE.md is inlined at the end of this file (see "Shared rules").

A user-level `~/.gemini/GEMINI.md` may carry general defaults (for example a C++ standard or a plan-first rule). In this repo, CLAUDE.md overrides them: the C++ standard is C++17.

## Skills

The repo's skills live in `.claude/skills/<name>/SKILL.md`. Read the one that applies and follow it as written:
- `write-pr`: every PR title and body.
- `model-provenance`: every physics, correlation, solver or Jacobian change.

## Attribution

- End every commit message with a trailer: `Co-Authored-By: Gemini (Antigravity) <noreply@google.com>`.
- End every PR body with: `Generated with Antigravity (Gemini)`.
- Label every PR and issue you open as `agent:gemini`, plus its type and area labels (see CLAUDE.md, "Labels").

## Limits when running unattended (auto mode, headless `agy -p`)

- **Never merge, never queue, never push to `main`.** The owner queues. Never use `gh pr merge`.
- **Never rewrite a branch that is under review** (no amend, no force-push). Answer a review with new commits. A rebase onto a moved base is allowed when the owner or the reviewer asks for it.
- **Never edit agent configuration:** `CLAUDE.md`, this file, `.claude/`, or any settings file. Propose a change to the owner instead.
- **Never run `agy` with `--dangerously-skip-permissions`.**
- **A refused command is the owner's decision, never yours to get around.** Do not rephrase it, and do not reach the same target through another tool or script. Send a `blocked` ping (below) whose note names the rule you need, e.g. `needs command(<program>) for: <command>`. Say the same in your reply, then stop.
- **Stop at a decision that is the owner's** (an architecture change, a changed default, a disputed review finding) and say so in the PR or issue.

## Working with Claude (reviews and pings)

Claude reviews your PRs before the owner queues them.
- **GitHub is the record.** Reviews, replies and findings go on the PR or issue as comments. In a review reply, give each finding one line: fixed (commit plus evidence), disputed (reason) or deferred (issue number).
- **A ping only wakes the other side.** When something is ready for Claude (a PR opened, a review addressed, a question), append one JSON line to `tmp/agent-inbox/claude.jsonl` in the main checkout (`/Users/thiemo/Projects/combaero/tmp/agent-inbox/claude.jsonl`) with a single shell command:

  ```bash
  printf '%s\n' '{"from":"gemini","ts":"2026-11-24T20:00:00Z","pr":1303,"issue":null,"event":"review-addressed","head":"a857ff8b","note":"one line"}' >> /Users/thiemo/Projects/combaero/tmp/agent-inbox/claude.jsonl
  ```

  `event` is one of `pr-opened`, `review-addressed`, `question`, `blocked`, `test`. Append only; never rewrite or truncate the file.
- **A claim of "fixed" or "verified" quotes its evidence:** the command run and the number it printed, at the commit it ran on. An unmeasured claim is "not verified", said so.
- **Claude pings you back** by resuming your conversation (`agy -r <id> -p ...`).
- **Questions that need the owner** go to the owner, not to Claude.

## Shared rules (CLAUDE.md, inlined)

Antigravity inlines CLAUDE.md here at startup. It sits last on purpose: a rule file over the 24 KB per-file cap after expansion is truncated from the bottom, so growth in CLAUDE.md cuts its own tail, never the Gemini rules above. GEMINI.md plus CLAUDE.md is about 18 KB. If the rules budget demotes this file to a pointer, or the expansion is truncated, read CLAUDE.md in full before the first command.

@[CLAUDE.md](CLAUDE.md)
