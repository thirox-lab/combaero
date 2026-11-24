# CombAero -- Agent Instructions

**[CLAUDE.md](CLAUDE.md) is authoritative for every coding agent (Codex, Cursor, Gemini, Claude, any other).** Read it first, in full, and follow it. This file does not restate it, so the two cannot drift apart.

What differs per agent lives in that agent's own file:
- Gemini (Antigravity): [GEMINI.md](GEMINI.md).

Any other agent:
- **Skills:** read the skills in `.claude/skills/<name>/SKILL.md`. `write-pr` applies to every PR; `model-provenance` to every physics, solver or Jacobian change.
- **Attribution:** add a `Co-Authored-By:` trailer naming the agent to every commit. Label every PR and issue you open with `agent:<name>` (see CLAUDE.md, "Labels").
- **Limits:** never merge, queue or push to `main`. Never edit agent configuration (`CLAUDE.md`, `AGENTS.md`, `GEMINI.md`, `.claude/`, `.agents/`). Stop at decisions that are the owner's.
