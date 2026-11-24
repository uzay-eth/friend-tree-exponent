#!/usr/bin/env bash
set -euo pipefail

# Run as `bash round2/verify/gold/rebuild.sh NEW_DIRECTORY [gold1|gold2]`.
usage() { printf '%s\n' 'Usage: rebuild.sh NEW_DIRECTORY [gold1|gold2]' >&2; exit 2; }
[[ $# -ge 1 && $# -le 2 ]] || usage
presets=(gold1 gold2)
case ${2:-all} in
    all) ;;
    gold1|gold2) presets=("$2") ;;
    *) usage ;;
esac

here=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
src=$here/../../proofs/neon2/quintet-cert
cxx=${CXX:-g++}
flags=(-std=c++17 -O3 -fopenmp -Wall -Wextra -Werror -pedantic)
runtime=(LC_ALL=C OMP_WAIT_POLICY=PASSIVE OMP_DYNAMIC=FALSE
         OMP_NUM_THREADS=10 OMP_THREAD_LIMIT=10 OMP_PROC_BIND=FALSE)

export LC_ALL=C
# Authenticate the release sources before compiling. The historical receipt is
# a comparison target; the builder itself only reads its embedded inputs.
(cd -- "$src" && sha256sum -c "$here/inputs.sha256")
mkdir -- "$1"
work=$(CDPATH= cd -- "$1" && pwd)
{
    printf '%s\n' 'source_commit=edd297185538b5ec399c04cb7fc80ebe2c07310a'
    printf 'compiler_flags=%s\n' "${flags[*]}"
    "$cxx" --version
    "$cxx" -dumpmachine
    uname -srmo
    lscpu
    /usr/bin/time --version
    printf '%s\n' "${runtime[@]}"
} > "$work/environment.txt"
"$cxx" "${flags[@]}" "$src/build_quintet.cpp" -o "$work/build_quintet" > "$work/compile.log" 2>&1 || {
    cat "$work/compile.log" >&2
    exit 1
}
ldd "$work/build_quintet" > "$work/libraries.txt"
printf 'preset\tside\texit\twall_s\tuser_s\tsystem_s\tpeak_rss_kib\tbytes\tsha256\thistorical_sha256\tstatus\n' > "$work/results.tsv"
failed=0
for preset in "${presets[@]}"; do
    for side in lower upper; do
        tag=$preset-$side
        dir=$work/$tag
        mkdir -- "$dir"
        file=quintet_$side.bin
        expected=$(awk -v file="$tag.bin" '$2 == file {print $1}' "$here/expected.sha256")
        [[ $expected =~ ^[a-f0-9]{64}$ ]] || { printf 'Missing expected digest for %s\n' "$tag" >&2; exit 2; }
        printf 'Starting %s\n' "$tag" >&2
        rc=0
        (cd -- "$dir" && env -i "${runtime[@]}" /usr/bin/time -f '%e\t%U\t%S\t%M' -q -o time.tsv \
            ../build_quintet "$preset" "$side" "$file" > build.log 2>&1) || rc=$?
        timing=$'-\t-\t-\t-'
        if [[ -s $dir/time.tsv ]]; then timing=$(<"$dir/time.tsv"); fi
        size=- actual=-
        if [[ -f $dir/$file ]]; then
            size=$(stat --format=%s -- "$dir/$file")
            if digest=$(sha256sum -- "$dir/$file"); then actual=${digest%% *}; fi
        fi
        if [[ $rc != 0 ]]; then status=build-failed
        elif [[ $actual == - ]]; then status=unreadable-output
        elif [[ $actual != "$expected" ]]; then status=digest-mismatch
        else status=match
        fi
        if [[ $status != match ]]; then failed=1; fi
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$preset" "$side" "$rc" "$timing" "$size" "$actual" "$expected" "$status" >> "$work/results.tsv"
        printf '%s: %s (exit %s, SHA-256 %s)\n' "$tag" "$status" "$rc" "$actual" >&2
    done
done
exit "$failed"
