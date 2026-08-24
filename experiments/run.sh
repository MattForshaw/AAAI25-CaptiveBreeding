#!/usr/bin/env bash
#
# Run solver configurations listed in experiments/manifest.csv.
#
# All paths are resolved relative to the repository root, so this works from any
# working directory and on any machine. Output goes to results/<run_id>-r<N>.out.
#
# Usage:
#   experiments/run.sh                        # every family
#   experiments/run.sh a10-leaveout wolves    # selected families
#   experiments/run.sh --dry-run scaling      # print commands, run nothing
#   experiments/run.sh --force a10-leaveout   # re-run configs that already have output
#
# Each run writes results/<id>.out (the ###COBREEDER-* records and group listings) and
# results/<id>.log (developer traces from stderr). Repeat N runs with --seed N, so each
# repeat is a distinct but individually reproducible sample.
#
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MANIFEST="$ROOT/experiments/manifest.csv"
RESULTS="$ROOT/results"
SOLVER="$ROOT/cobreeder/solver.py"
PYTHON="${PYTHON:-python3}"

DRY_RUN=0
FORCE=0
FAMILIES=()

while [ $# -gt 0 ]; do
    case "$1" in
        --dry-run) DRY_RUN=1 ;;
        --force)   FORCE=1 ;;
        -h|--help) sed -n '2,13p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        -*)        echo "run.sh: unknown option $1" >&2; exit 2 ;;
        *)         FAMILIES+=("$1") ;;
    esac
    shift
done

[ -f "$MANIFEST" ] || { echo "run.sh: missing $MANIFEST" >&2; exit 1; }
[ -f "$SOLVER" ]   || { echo "run.sh: missing $SOLVER" >&2; exit 1; }

wanted() {
    [ ${#FAMILIES[@]} -eq 0 ] && return 0
    local f
    for f in "${FAMILIES[@]}"; do [ "$f" = "$1" ] && return 0; done
    return 1
}

mkdir -p "$RESULTS"

# Record the environment once per invocation: CP-SAT's default worker count follows
# the core count, so results are not comparable across machines (see README).
if [ "$DRY_RUN" -eq 0 ]; then
    {
        echo "date            $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
        echo "host            $(hostname)"
        echo "os              $(uname -srm)"
        echo "cores           $(getconf _NPROCESSORS_ONLN 2>/dev/null || echo unknown)"
        echo "python          $($PYTHON --version 2>&1 | tr '\n' ' ' | sed 's/ *$//')"
        echo "ortools         $($PYTHON -c 'import ortools; print(ortools.__version__)' 2>/dev/null || echo unknown)"
        echo "commit          $(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null || echo 'not a git checkout')"
        echo "dirty           $(git -C "$ROOT" status --porcelain 2>/dev/null | wc -l | tr -d ' ') uncommitted paths"
    } > "$RESULTS/environment.txt"
fi

total=0; ran=0; skipped=0; failed=0

# Field order must match the manifest header.
while IFS=, read -r family run_id individuals pr corrals objective wa wb total_individuals repeats; do
    [ "$family" = "family" ] && continue          # header
    [ -z "${family:-}" ] && continue              # blank line
    case "$family" in \#*) continue ;; esac       # comment
    wanted "$family" || continue

    for path in "$individuals" "$pr" "$corrals"; do
        if [ ! -f "$ROOT/$path" ]; then
            echo "MISSING INPUT  $run_id: $path" >&2
            failed=$((failed + 1))
            continue 2
        fi
    done

    for rep in $(seq 1 "$repeats"); do
        total=$((total + 1))
        out="$RESULTS/${run_id}-r${rep}.out"
        log="$RESULTS/${run_id}-r${rep}.log"
        tag="${run_id}-r${rep}"

        if [ -s "$out" ] && [ "$FORCE" -eq 0 ]; then
            skipped=$((skipped + 1))
            continue
        fi

        if [ "$DRY_RUN" -eq 1 ]; then
            echo "$PYTHON $SOLVER run $individuals $pr $corrals $objective $tag $wa $wb $total_individuals --seed $rep > results/${tag}.out"
            continue
        fi

        printf '%-34s ' "$tag"
        set +e
        # Results to .out, developer traces to .log. Keeping them apart is what lets
        # extract.sh read the .out without filtering, and keeps the .out small.
        #
        # The seed is the repeat number. The solver now defaults to one worker with a
        # fixed seed, so repeats of a configuration would otherwise be byte-identical
        # and measure nothing. Varying the seed makes each repeat a different, and
        # individually reproducible, sample of the search.
        ( cd "$ROOT" && "$PYTHON" "$SOLVER" run \
              "$individuals" "$pr" "$corrals" "$objective" \
              "$tag" "$wa" "$wb" "$total_individuals" \
              --seed "$rep" > "$out" 2> "$log" )
        rc=$?
        set -e

        # The solver exits 0 for any definitive answer, including INFEASIBLE, and
        # non-zero when it produced no usable result: 3 MODEL_INVALID, 4 UNKNOWN,
        # 1 bad input, 2 command-line misuse. Report the status either way.
        # || true on both: with `set -o pipefail` a grep that matches nothing fails the
        # whole substitution, and under `set -e` that aborted the sweep instead of
        # reporting the run -- so the NO COMPLETION RECORD branch below was unreachable
        # and one crashed configuration took the remaining ones with it.
        status=$(grep -m1 '^###COBREEDER-COMPLETION' "$out" | cut -d, -f2 || true)
        case "${status:-}" in
            4) label="OPTIMAL" ;;
            2) label="FEASIBLE" ;;
            3) label="INFEASIBLE" ;;
            1) label="MODEL_INVALID" ;;
            0) label="UNKNOWN" ;;
            *) label="NO COMPLETION RECORD" ;;
        esac
        secs=$(grep -m1 '^###COBREEDER-COMPLETION' "$out" | cut -d, -f5 || true)
        sols=$(grep -c '^###COBREEDER-SOLUTION' "$out" || true)

        if [ "$rc" -eq 0 ]; then
            printf '%-14s %8ss  %3s solutions\n' "$label" "${secs:-?}" "$sols"
            ran=$((ran + 1))
        else
            printf '%-14s %8ss  exit %s -- see results/%s.log\n' "$label" "${secs:-?}" "$rc" "$tag"
            failed=$((failed + 1))
        fi
    done
done < <(tr -d '\r' < "$MANIFEST")   # tolerate CRLF, e.g. a manifest edited in Excel

if [ "$DRY_RUN" -eq 1 ]; then
    echo "# dry run: $total invocations would be issued"
else
    echo
    echo "$total invocations: $ran ran, $skipped skipped (output present), $failed failed"
    echo "results in $RESULTS  (environment.txt records host, cores and versions)"
    echo "next: experiments/extract.sh"
fi
