#!/bin/bash
# scripts/bulk_merge_safe.sh
# Queue "safe" Dependabot PRs (grouped safe-updates, or same-major bumps) in
# the merge queue (ruleset "Iron"), after a security housekeeping pass.
#
# 0. Housekeeping (skip with --no-audit), reported, never a gate: every
#    dependency floor against OSV (dependency_floors.py --osv) and the
#    frontend lock with pnpm audit. Dependabot moves only the top of each
#    range and only direct dependencies; these cover the rest.
# 1. A PR is eligible when its author is app/dependabot and the bump is safe
#    (same major, or a grouped safe set); --major also admits major bumps,
#    which the user has decided to take.
# 2. Waves: PRs that rewrite the same lock file (same manifest directory)
#    cannot share the queue -- the second one's lock was generated against
#    the old main and bounces. Each wave queues at most one PR per manifest
#    directory, then waits until each is merged or has left the queue.
# 3. A PR whose manifest directory changed on main after its branch point
#    (an earlier merge, in this run or before it) is first rebased by
#    Dependabot ("@dependabot rebase"; never gh pr update-branch, whose
#    rebase is unsigned), then waits for CI. Git may merge such a lock file
#    without a conflict, yet pnpm's frozen install rejects it (#976 after
#    #975: "Broken lockfile: no entry for 'eslint@10.13.0'").
# 4. Gate before queueing: all checks finished, none failed (SKIPPED/NEUTRAL
#    are fine), and Biome is pinned in package.json only at the PR head
#    (scripts/check_version_sync.py). Queued at exactly that head
#    (--match-head-commit); the queue picks the merge method (squash).
# Anything else is reported and skipped; nothing is forced.
#
# Usage: scripts/bulk_merge_safe.sh [--dry-run] [--no-audit] [--major] [PR ...]

set -euo pipefail

GH_BIN="$(command -v gh)"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
REPO="thirox-lab/combaero"
DRY_RUN=false
AUDIT=true
MAJOR=false
ONLY=()

for arg in "$@"; do
    case $arg in
        --dry-run) DRY_RUN=true ;;
        --no-audit) AUDIT=false ;;
        --major) MAJOR=true ;;
        [0-9]*) ONLY+=("$arg") ;;
        *) echo "Usage: $0 [--dry-run] [--no-audit] [--major] [PR ...]"; exit 1 ;;
    esac
done

is_safe_title() {
    local title="$1"
    [[ "$title" == *"safe-updates"* || "$title" == *"actions group"* ]] && return 0
    if [[ "$title" =~ from[[:space:]]v?([0-9]+)\.([0-9]+)[^[:space:]]*[[:space:]]to[[:space:]]v?([0-9]+)\. ]]; then
        [[ "${BASH_REMATCH[1]}" == "${BASH_REMATCH[3]}" ]] && return 0  # same major
    fi
    return 1
}

pr_field() { $GH_BIN pr view "$1" -R "$REPO" --json "$2" --jq ".$2"; }

manifest_dir() {
    # The directory of the PR's manifest or lock file: PRs sharing it share
    # a lock file. GitHub Actions bumps have none and share ".github".
    local d
    d=$($GH_BIN pr view "$1" -R "$REPO" --json files --jq '
      [.files[].path | select(test("(^|/)(package\\.json|pnpm-lock\\.yaml|pyproject\\.toml|uv\\.lock)$"))
       | (split("/") | .[:-1] | join("/"))] | first // ""')
    [[ -z "$d" ]] && d=".github"
    echo "$d"
}

in_queue() {
    # true while the PR has a merge-queue entry (scripts/merge_queue_status.sh);
    # unknown when the query fails, which wait_for_queue waits out.
    local out
    out=$("$ROOT/scripts/merge_queue_status.sh" "$1") || { echo unknown; return; }
    if [[ "$out" == *" queued: "* ]]; then echo true; else echo false; fi
}

checks_state() {
    # PASS when every check finished and none failed; PENDING; or FAIL.
    $GH_BIN pr view "$1" -R "$REPO" --json statusCheckRollup --jq '
      [.statusCheckRollup[] | (.conclusion // .state // "")] as $c
      | if ($c | length) == 0 then "PENDING"
        elif any($c[]; . == "" or . == "PENDING" or . == "QUEUED" or . == "IN_PROGRESS" or . == "EXPECTED") then "PENDING"
        elif all($c[]; . == "SUCCESS" or . == "SKIPPED" or . == "NEUTRAL") then "PASS"
        else "FAIL" end'
}

wait_for_new_head() {
    # Echo the new head oid once it differs from $2, or TIMEOUT (~20 min).
    local n="$1" old="$2" i=0 h
    while :; do
        sleep 30
        h=$(pr_field "$n" headRefOid)
        [[ "$h" != "$old" ]] && { echo "$h"; return; }
        (( i++ > 40 )) && { echo "TIMEOUT"; return; }
    done
}

wait_for_checks() {
    local n="$1" i=0 s
    sleep 20  # let GitHub register the checks of a new head
    while :; do
        s=$(checks_state "$n")
        [[ "$s" != "PENDING" ]] && { echo "$s"; return; }
        (( i++ > 120 )) && { echo "TIMEOUT"; return; }  # ~60 min
        sleep 30
    done
}

wait_for_queue() {
    # MERGED, or DEQUEUED once the PR is open and out of the queue (a failed
    # or conflicting entry), or TIMEOUT (~90 min; a full queue run is ~17).
    local n="$1" i=0
    while :; do
        sleep 60
        [[ "$(pr_field "$n" state)" == "MERGED" ]] && { echo "MERGED"; return; }
        (( i > 1 )) && [[ "$(in_queue "$n")" == "false" ]] && { echo "DEQUEUED"; return; }
        (( i++ > 90 )) && { echo "TIMEOUT"; return; }
    done
}

housekeeping() {
    echo "=== housekeeping (reported, not a gate)"
    uv run --no-project --with packaging "$ROOT/scripts/dependency_floors.py" --osv \
        || echo "!!! a dependency floor has an advisory: raise it (see CLAUDE.md, Dependency floors)"
    (cd "$ROOT/gui/frontend" && pnpm audit --prod) \
        || echo "!!! the frontend production dependencies have an advisory"
    (cd "$ROOT/gui/frontend" && pnpm audit 2>&1 | tail -1) || true
}

stale_base() {
    # true when main changed the PR's manifest directory after the PR's
    # branch point, so its lock file was generated against an older main.
    local n="$1" d="$2" head mb
    head=$(pr_field "$n" headRefOid)
    git -C "$ROOT" fetch -q origin main "$head"
    mb=$(git -C "$ROOT" merge-base "$head" origin/main)
    if git -C "$ROOT" diff --quiet "$mb" origin/main -- "$d"; then echo false; else echo true; fi
}

prepare_one() {
    # Gate one PR and queue it. $2 = "rebase" when its lock file predates
    # main's. Returns 0 when queued (or would be, in a dry run).
    local n="$1" rebase="$2" head
    echo "--- #$n  $(pr_field "$n" title)  [$(manifest_dir "$n")]"
    if [[ "$(in_queue "$n")" == "true" ]]; then echo "    already in the queue; waiting for it"; return 0; fi
    if [[ "$rebase" == "rebase" ]]; then
        if $DRY_RUN; then echo "    would ask Dependabot to rebase (its lock file predates main's), wait for CI"; return 0; fi
        echo "    asking Dependabot to rebase: its lock file predates main's"
        local old; old=$(pr_field "$n" headRefOid)
        $GH_BIN pr comment "$n" -R "$REPO" --body "@dependabot rebase" >/dev/null
        if [[ "$(wait_for_new_head "$n" "$old")" == "TIMEOUT" ]]; then
            echo "    skip: Dependabot did not rebase (if it says the branch was edited, comment '@dependabot recreate')"
            return 1
        fi
    fi
    if [[ "$(checks_state "$n")" == "PENDING" ]] && ! $DRY_RUN; then
        echo "    CI running; waiting..."
        wait_for_checks "$n" >/dev/null
    fi
    [[ "$(checks_state "$n")" == "PASS" ]] || { echo "    skip: checks not all passed"; return 1; }

    head=$(pr_field "$n" headRefOid)
    git -C "$ROOT" fetch -q origin "$head"
    if ! uv run --no-project "$ROOT/scripts/check_version_sync.py" --ref "$head"; then
        echo "    skip: pinned versions disagree at the PR head"; return 1
    fi
    if $DRY_RUN; then echo "    would queue $head"; return 0; fi
    $GH_BIN pr merge "$n" -R "$REPO" --squash --match-head-commit "$head"
    echo "    queued $head"
}

$AUDIT && housekeeping

if [[ ${#ONLY[@]} -eq 0 ]]; then
    # bash 3.2 (macOS): no mapfile, no associative arrays
    while read -r n; do ONLY+=("$n"); done < <($GH_BIN pr list -R "$REPO" --author "app/dependabot" --json number --jq '.[].number' | sort -n)
fi
[[ ${#ONLY[@]} -eq 0 ]] && { echo "No open Dependabot PRs."; exit 0; }

# Eligibility first; then pending holds "number:directory" entries.
pending=()
skipped=0
for n in "${ONLY[@]}"; do
    if [[ "$($GH_BIN pr view "$n" -R "$REPO" --json author --jq .author.login)" != "app/dependabot" ]]; then
        echo "--- #$n  skip: not a Dependabot PR"; skipped=$((skipped + 1)); continue
    fi
    if ! $MAJOR && ! is_safe_title "$(pr_field "$n" title)"; then
        echo "--- #$n  skip: not a safe (grouped/minor/patch) bump: $(pr_field "$n" title)"; skipped=$((skipped + 1)); continue
    fi
    pending+=("$n:$(manifest_dir "$n")")
done

merged=0
dirty=" "   # directories that had a merge in this run, space-delimited
wave=0
while [[ ${#pending[@]} -gt 0 ]]; do
    wave=$((wave + 1))
    echo "=== wave $wave"
    used=" "; this=(); rest=()
    for e in "${pending[@]}"; do
        d="${e#*:}"
        if [[ "$used" == *" $d "* ]]; then rest+=("$e"); else this+=("$e"); used="$used$d "; fi
    done
    queued=()
    for e in "${this[@]}"; do
        n="${e%%:*}"; d="${e#*:}"
        # dirty: a merge earlier in this run (all a dry run can know of);
        # stale_base: main moved the directory since the PR branched.
        r=no
        if [[ "$dirty" == *" $d "* || "$(stale_base "$n" "$d")" == "true" ]]; then r=rebase; fi
        if prepare_one "$n" "$r"; then queued+=("$e"); else skipped=$((skipped + 1)); fi
    done
    if ! $DRY_RUN; then
        for e in ${queued[@]+"${queued[@]}"}; do
            n="${e%%:*}"; d="${e#*:}"
            s=$(wait_for_queue "$n")
            echo "    #$n: $s"
            if [[ "$s" == "MERGED" ]]; then merged=$((merged + 1)); dirty="$dirty$d "; else skipped=$((skipped + 1)); fi
        done
    else
        # A dry run assumes every queued PR merges, to show the later waves.
        for e in ${queued[@]+"${queued[@]}"}; do merged=$((merged + 1)); dirty="$dirty${e#*:} "; done
    fi
    pending=(${rest[@]+"${rest[@]}"})
done
if $DRY_RUN; then echo "=== dry run: $merged would be queued, $skipped skipped"; else echo "=== done: $merged merged, $skipped skipped"; fi
