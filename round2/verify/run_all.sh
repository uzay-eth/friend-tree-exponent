#!/usr/bin/env bash
# Run the round-two modest verifications, sequentially, in clean scratch space.
set -euo pipefail

usage() {
    printf '%s\n' 'Usage: bash round2/verify/run_all.sh [--lane astra|astra-late|neon|neon2|fable | --slow]' \
        'With no arguments, run all five modest suites. --slow runs only the Neon D10/D12 controls.' \
        'See round2/verify/README.md for larger replays.'
}

lanes=(astra astra-late neon neon2 fable)
lane_args=()
case $# in
    0) ;;
    1) case $1 in
           --help|-h) usage; exit 0 ;;
           --slow) lanes=(neon); lane_args=(--slow) ;;
           *) usage >&2; exit 2 ;;
       esac ;;
    2) if [[ $1 != --lane ]]; then usage >&2; exit 2; fi
       case $2 in
           astra|astra-late|neon|neon2|fable) lanes=("$2") ;;
           *) usage >&2; exit 2 ;;
       esac ;;
    *) usage >&2; exit 2 ;;
esac

verify_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
for lane in "${lanes[@]}"; do
    if [[ ! -f "$verify_dir/$lane/run.sh" ]]; then
        printf 'Missing suite: %s/run.sh\n' "$lane" >&2
        exit 1
    fi
done

# Some unchanged verifiers use Python assertions. Ignore inherited Python
# startup overrides, and always run against a new, isolated CPython 3.12 venv.
unset BASH_ENV ENV PYTHONHOME PYTHONPATH PYTHONOPTIMIZE PYTHONPYCACHEPREFIX
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 LC_ALL=C
if ! "${PYTHON:-python3.12}" -IB -c 'import sys; sys.exit(0 if sys.implementation.name == "cpython" and sys.version_info[:2] == (3,12) and __debug__ else 1)'; then
    printf 'CPython 3.12 with assertions enabled is required; set PYTHON to its executable.\n' >&2
    exit 1
fi

# Keep foreign BLAS/OpenMP settings (and parallel suite launches) from
# consuming the host. The quick suites do not need any downloaded inputs.
export OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OMP_DYNAMIC=FALSE OMP_WAIT_POLICY=PASSIVE
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 NEON2_THREADS=1

scratch_root=$(cd -- "${TMPDIR:-/tmp}" && pwd)
work=$(mktemp -d "$scratch_root/ft2-round2.XXXXXXXX")
trap 'rm -rf -- "$work"' EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
"${PYTHON:-python3.12}" -IB -m venv --without-pip "$work/venv"
export PYTHON="$work/venv/bin/python"
export PATH="$work/venv/bin:$PATH"
"$PYTHON" -IB --version

passes=()
declare -A seen=()
started=$SECONDS
for lane in "${lanes[@]}"; do
    mkdir -- "$work/$lane"
    printf '\n=== %s ===\n' "$lane"
    before=$SECONDS
    # A success diagnostic cannot hide a failed child. Also keep stderr out
    # of the result protocol: only structured stdout records count.
    log="$work/$lane.out"
    if (cd -- "$work/$lane" || exit; export TMPDIR="$PWD"; bash "$verify_dir/$lane/run.sh" "${lane_args[@]}") \
        | tee "$log" | sed -u "s/^/[$lane] /"; then
        count=0
        while IFS= read -r line; do
            [[ $line == 'RESULT '* ]] || continue
            if [[ ! $line =~ ^RESULT\ (${lane}-[a-z0-9]+(-[a-z0-9]+)*):\ ([[:graph:]].*)$ ]]; then
                printf 'FAIL  %s: malformed result record: %s\n' "$lane" "$line" >&2
                exit 1
            fi
            id=${BASH_REMATCH[1]}
            scope=${BASH_REMATCH[3]}
            if [[ ${seen[$id]+present} ]]; then
                printf 'FAIL  %s: duplicate result ID: %s\n' "$lane" "$id" >&2
                exit 1
            fi
            seen[$id]=1
            passes+=("PASS  $id: $scope")
            count=$((count + 1))
        done < "$log"
        if ((count == 0)); then
            printf 'FAIL  %s: no result records\n' "$lane" >&2
            exit 1
        fi
        printf '=== %s: %s results in %ss ===\n' "$lane" "$count" "$((SECONDS - before))"
    else
        status=$?
        printf '\nFAIL  %s (exit %s)\n' "$lane" "$status" >&2
        exit "$status"
    fi
done

printf '\n%s\n' 'ROUND-TWO VERIFICATION SUMMARY'
printf '%s\n' "${passes[@]}"
printf 'ALL %s SELECTED RESULTS PASSED (%ss)\n' "${#passes[@]}" "$((SECONDS - started))"
