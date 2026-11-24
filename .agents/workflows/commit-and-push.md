---
description: Signed local commit (commit-local), then push the task branch. Never main.
---

## Steps

### 1. Follow commit-local
Run every step of [commit-local](commit-local.md).

### 2. Push the task branch
```bash
git push -u origin HEAD
```
Never push to `main` (branch protection). Never force-push a branch that is under review: answer a review with new commits.
