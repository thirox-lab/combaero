# Bulk Merge Safe Updates Workflow

**Only when the owner asks for a Dependabot round.** It merges PRs, which agents otherwise never do (CLAUDE.md, "Releases and Dependencies").

This workflow automates the consolidation of verified Dependabot updates into the `main` branch. It ensures that only "green" and "safe" (minor/patch/grouped) updates are merged, significantly reducing manual overhead.

## Prerequisites
- GitHub CLI (`gh`) installed and authenticated.
- Local `main` branch is clean and synchronized.

## Steps

### 1. Discovery & LLM Review
The Agent (LLM) must first perform a dry-run to identify candidates and verify their safety.
```bash
./scripts/bulk_merge_safe.sh --dry-run
```

### 2. Safety Audit
For each PR identified by the script, the Agent must:
- [ ] Confirm the CI status is truly **SUCCESS**.
- [ ] For `minor` bumps, quickly check for known breaking changes in project-critical libraries.
- [ ] Group the candidates into a summary table for the user.

### 3. User Approval
Present the candidate table to the user. Do **NOT** proceed until the user provides explicit approval (e.g., "Proceed with merge").

### 4. Execution
Once approved, run the merger:
```bash
./scripts/bulk_merge_safe.sh
```

### 5. Final Sync
The script handles the final `git pull` from `main`, but the Agent should verify the local state is up-to-date and all checks remain green on `main`.
