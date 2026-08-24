#!/usr/bin/env bash
#
# Collect the machine-readable records from results/*.out into six CSVs with
# named column headers.
#
# Every record already carries its run_id, so filenames are suppressed (grep -h)
# and the output is directly loadable:
#
#     read_csv("results/completions.csv")
#
# Usage:
#   experiments/extract.sh              # read results/
#   experiments/extract.sh <directory>  # read somewhere else
#   experiments/extract.sh --legacy     # format AAAICR-Thesis.R expects
#
# --legacy writes no header row and keeps the ###COBREEDER-* tag as the first
# field, which is the shape the existing R analysis reads via col_names=:
#
#     completions_col_names <- c("Filename","Status","NumConflicts", ...)
#
set -euo pipefail

LEGACY=0
POSITIONAL=()
while [ $# -gt 0 ]; do
    case "$1" in
        --legacy)  LEGACY=1 ;;
        -h|--help) sed -n '2,21p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        -*)        echo "extract.sh: unknown option $1" >&2; exit 2 ;;
        *)         POSITIONAL+=("$1") ;;
    esac
    shift
done

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RESULTS="${POSITIONAL[0]:-$ROOT/results}"

[ -d "$RESULTS" ] || { echo "extract.sh: no such directory: $RESULTS" >&2; exit 1; }

shopt -s nullglob
outs=("$RESULTS"/*.out)
if [ ${#outs[@]} -eq 0 ]; then
    echo "extract.sh: no .out files in $RESULTS -- run experiments/run.sh first" >&2
    exit 1
fi

# Column names follow the format strings in cobreeder/solver.py.
collect() {
    local tag="$1" header="$2" dest="$3" rows expected mismatched
    if [ "$LEGACY" -eq 1 ]; then
        # Keep the tag as field 1 and emit no header: the R script supplies its own
        # names and would otherwise read the header as a data row.
        grep -h "^###COBREEDER-${tag}," "${outs[@]}" > "$dest" || true
        rows=$(wc -l < "$dest")
        expected=$(( $(printf '%s' "$header" | awk -F, '{print NF}') + 1 ))
    else
        printf '%s\n' "$header" > "$dest"
        grep -h "^###COBREEDER-${tag}," "${outs[@]}" | sed "s/^###COBREEDER-${tag},//" >> "$dest" || true
        rows=$(( $(wc -l < "$dest") - 1 ))
        expected=$(printf '%s' "$header" | awk -F, '{print NF}')
    fi
    printf '  %-28s %6d rows\n' "$(basename "$dest")" "$rows"

    # A row whose field count differs from the header means .out files from two solver
    # versions are being read together. The ARGS record lost its unused 'subst' field, so
    # a results/ directory spanning that change mixes four- and five-field rows -- which
    # loads without error and silently shifts every column.
    mismatched=$(awk -F, -v want="$expected" -v skip="$LEGACY" \
        'NR==1 && skip==0 {next} NF!=want {n++} END {print n+0}' "$dest")
    if [ "$mismatched" -gt 0 ]; then
        echo "  !! $mismatched of $rows rows in $(basename "$dest") have a field count" >&2
        echo "     other than $expected. Your .out files span a change to the record" >&2
        echo "     layout. Re-run with experiments/run.sh --force, or delete the stale" >&2
        echo "     .out files, rather than analysing the mixture." >&2
        MISMATCH=1
    fi
}

MISMATCH=0

if [ "$LEGACY" -eq 1 ]; then
    echo "reading ${#outs[@]} .out files from $RESULTS  (legacy: no header, tag retained)"
else
    echo "reading ${#outs[@]} .out files from $RESULTS"
fi

collect ARGS \
    'pr_file,corral_file,objective_function,run_id' \
    "$RESULTS/args.csv"

collect SOLUTION \
    'solution_number,wallclock,seconds_elapsed,objective,num_corrals,objective_id,num_individuals,run_id' \
    "$RESULTS/solutions.csv"

collect ALLOCATION \
    'solution_number,individual_index,corral_index,run_id' \
    "$RESULTS/allocations.csv"

collect COMPLETION \
    'status,conflicts,branches,wall_time,num_solutions,num_corrals,objective_id,num_individuals,run_id' \
    "$RESULTS/completions.csv"

collect SEARCH \
    'seed,workers,time_limit,deterministic_time,run_id' \
    "$RESULTS/search.csv"

collect BOUNDS \
    'component,lo,hi,weight,scale,run_id' \
    "$RESULTS/bounds.csv"

echo
echo "status codes: 0 unknown, 1 model invalid, 2 feasible, 3 infeasible, 4 optimal"
[ "$LEGACY" -eq 1 ] && echo "legacy layout: field 1 is the record tag; supply col_names yourself"
[ "$MISMATCH" -eq 1 ] && exit 1
exit 0
