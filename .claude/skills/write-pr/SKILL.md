---
name: write-pr
description: Use when writing or editing a pull request title, body or comment, including an independent review's verdict posted on the PR.
---

### Title Rules
- Conventional Commits, under 72 chars, scope in parens. The type list is open:
  `feat`, `fix`, `refactor`, `perf`, `chore`, `docs`, `test`, `ci`, plus the
  API-lifecycle pair this repo uses -- `remove(junction):`, `rename(junction):`.

### Tone & General Rules
- Be direct and high-signal. No conversational fluff, boilerplate templates, or
  statements like "All tests pass".
- Write for the final squash commit -- ignore WIP commits and review
  iterations. Branch protection means the PR body IS the permanent record.
- Size the body to blast radius, not to commit type. A one-line `docs:` change
  to the release process earns more explanation than a large `feat:` nothing
  depends on yet. Anything touching release, publish, or packaging machinery is
  never "just a chore".
- No pointers into local, untracked paths (`tmp/...`, a worktree, a scratchpad)
  in a PR body or PR comment: no reader can open them, and they are gone once
  the worktree is archived. If what is there is key, summarise it in the text
  (the numbers, the reproducer's inputs); if not, omit it. The same holds for
  review-verdict comments.
- End with the required attribution footer.

### Standard PR Body Structure
1. **Summary:** 1-3 bullet points on *what* changed and *why*.
2. **Context & decisions:** edge cases, non-obvious implementation details, and
   dead ends worth not repeating.
3. **Not changed:** what you deliberately left alone, and why -- a tag not
   moved, a name not renamed, a deprecation kept. Omit when there is nothing.
4. **Artifacts (include only when relevant):**
   - **Measurement / evidence:** the numbers that justify the change --
     validation scorecard, MAE, analytic-vs-FD agreement, a falsification
     result, coverage instantiation counts. A physics, solver or correlation
     change without one is unfinished; see the `model-provenance` skill.
   - **Architecture / flow:** minimal Mermaid diagram for new state machines or
     inter-service flows.
   - **UI / visual:** table with `Before` and `After` screenshot columns.
   - **Performance:** table comparing target-branch baseline against this branch.
   - **API / contract:** short snippet showing breaking changes or public
     signature updates.

### What not to claim
- State what you measured. State separately what you did NOT verify, and why.
  An untested path named in the PR is cheap; one discovered during a release is
  not.
- A negative result, reported, is a good outcome. Never widen a tolerance, drop
  a case, or tune against the score to make a PR look clean.

### Review-verdict comments
One comment per independent review round, posted on the PR (the residual,
Jacobian or default rule in CLAUDE.md), so a later audit of agent mistakes
finds them on GitHub.
- **First line:** round, the head SHA reviewed, the verdict (ready / not ready /
  ship), and how many bugs were confirmed.
- **Confirmed bugs:** each with its size and a reproducer a reader can rerun
  (inputs and the call, not a script path), and the fix with its measured result.
- **Retracted:** claims the reviewer, or the PR body, made and then disproved,
  and why. They show where a reading went wrong.
- **Notes, not regressions:** known limits found on the way, with the issue
  filed for each.
- **Last line:** the tests run on the head, with their counts.
- Factual, no @-mentions, no pointers into local paths.

### High-Risk / Architecture Changes (Exception)
Expand into an RFC-style write-up -- problem statement, alternatives considered,
failure modes, rollout and rollback -- when the change removes or renames public
API, alters a physics closure or a residual form, changes release or packaging
behaviour, or carries high regression risk.
