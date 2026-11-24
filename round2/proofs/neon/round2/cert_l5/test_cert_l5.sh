#!/usr/bin/env bash
# Scratch-only, sequential reproducibility harness. A failed stage is fatal.
set -euo pipefail
cd -- "$(dirname -- "$0")"
export OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OPENBLAS_NUM_THREADS=1 TMPDIR="${TMPDIR:-$PWD}"
export ASAN_OPTIONS=halt_on_error=1:detect_leaks=1 UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1
ulimit -c 0
tmp=$(mktemp -d "$TMPDIR/neon-cert-l5.XXXXXX")
if [[ ${L5_KEEP_SCRATCH:-0} != 1 ]]; then trap 'rm -rf -- "$tmp"' EXIT; fi
printf 'Scratch: %s\n' "$tmp"
cc=${CC:-gcc}
flags=(-std=c11 -Wall -Wextra -Werror -fopenmp)
san=(-O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined -no-pie)
"$cc" "${flags[@]}" -O3 -ffp-contract=off cert_b3f_l5_build.c -lm -o "$tmp/build"
"$cc" "${flags[@]}" -O3 cert_b3f_l5_verify_gmp.c -lgmp -o "$tmp/verify"
"$cc" "${flags[@]}" "${san[@]}" cert_b3f_l5_verify_gmp.c -lgmp -o "$tmp/san"
profile10=(--D0 10 --KT 9 --LF 4000 --KF 80 --DP4 8 --KP4 6 --DYT 7 --KYT 5 --DYB 8 --KYB 6 --DS3 6 --KS3 5 --Dc 6)
profile12=(--D0 12 --KT 11 --LF 6000 --KF 100 --DP4 10 --KP4 8 --DYT 9 --KYT 7 --DYB 10 --KYB 8 --DS3 8 --KS3 6 --Dc 7 --nit 18)
common=(--buckets 2 --l5mask 3 --Q 10000000 --cushion 5e-7 --threads 1)
timed() { /usr/bin/time -f 'wall=%e rss_kib=%M' "$@"; }
bounded() { (ulimit -v 1950000; timed "$@"); }

reject() {
    local rc=0
    "$@" >"$tmp/reject.log" 2>&1 || rc=$?
    if ((rc != 1 && rc != 2)) || grep -Eq 'AddressSanitizer|LeakSanitizer|runtime error:' "$tmp/reject.log"; then
        cat "$tmp/reject.log"; printf 'Expected clean rejection, got exit %d: %s\n' "$rc" "$*" >&2; return 1
    fi
    grep -E 'FAIL|Invalid|Unknown' "$tmp/reject.log"
}

if (($# == 0)); then set -- 10 12; fi
for d0 in "$@"; do
        case $d0 in 10) profile=("${profile10[@]}");; 12) profile=("${profile12[@]}");; *) exit 2;; esac
        for mode in 2 1; do
            cert="$tmp/d$d0-m$mode.bin"
            bounded "$tmp/build" "${profile[@]}" "${common[@]}" --mode "$mode" --out "$cert"
            bounded "$tmp/verify" "$cert" --threads 1 --tree-check 30 | tee "$tmp/verify.log"
            read -r c0 C0 < <(awk '/\[COMPARISON\]/{gsub(",", ""); print $4,$7}' "$tmp/verify.log")
            bounded python3 -B tree_identity_l5_check.py "$cert" 30
            p=$(python3 -B -c 'import struct,sys; f=open(sys.argv[1],"rb"); f.seek(16); print(struct.unpack("<q",f.read(8))[0])' "$cert")
            if ((mode == 1)); then tight=$((p-20)); else tight=$((p+20)); fi
            neg="$tmp/verify"
            if ((d0 == 10)); then
                timed "$tmp/san" "$cert" --threads 1 --tree-check 30
                neg="$tmp/san"
            fi
            reject "$neg" "$cert" --threads 1 --tree-check 0 --override-pq "$tight" 10000000
        done
        if ((d0 == 10)); then
            for invalid in 'Dc 3' 'buckets 3' 'l5mask 0' 'Dc 4294967296' 'KS3 4'; do
                read -r k value <<<"$invalid"
                reject "$tmp/build" "${profile[@]}" "${common[@]}" --size "--$k" "$value"
            done
            for variant in entry reserved truncated extra dc buckets mask rows ks3; do
                python3 -B - "$cert" "$tmp/bad.bin" "$variant" <<'PY'
import os, shutil, struct, sys
src, dst, variant = sys.argv[1:]
shutil.copyfile(src,dst)
with open(dst,'r+b') as f:
    if variant == 'entry':
        f.seek(-8,os.SEEK_END); f.write(struct.pack('<q',10**15))
    elif variant == 'truncated':
        f.truncate(os.stat(dst).st_size-1)
    elif variant == 'extra':
        f.seek(0,os.SEEK_END); f.write(b'\0')
    else:
        offset = {'dc':96,'buckets':100,'mask':104,'reserved':108,'rows':112,'ks3':92}[variant]
        f.seek(offset); f.write(struct.pack('<I',{'ks3':4,'reserved':1}.get(variant,0)))
PY
                reject "$tmp/san" "$tmp/bad.bin" --threads 1 --tree-check 0
            done
            # Dropping a direction-selected adverse crossing can only make
            # the inequality easier. Its independent term audit below is
            # the negative control, not the fitted row sign.
            for defect in BUCKET SOURCE NEWBORN; do
                "$cc" "${flags[@]}" "${san[@]}" "-DL5_PLANT_$defect" cert_b3f_l5_verify_gmp.c -lgmp -o "$tmp/plant"
                reject "$tmp/plant" "$cert" --threads 1 --tree-check 30
            done
            for defect in crossing newborn; do
                if python3 -B tree_identity_l5_check.py "$cert" 0 --scaled-bounds "$c0" "$C0" --plant "$defect" >"$tmp/plant-python.log" 2>&1; then
                    printf 'Python planted %s survived\n' "$defect" >&2; exit 1
                fi
                if [[ $defect == crossing ]]; then grep 'AssertionError: boundary crossing' "$tmp/plant-python.log";
                else grep 'AssertionError: newborn ownership' "$tmp/plant-python.log"; fi
            done
        fi
        printf 'D0=%d CERTIFICATION AND CONTROLS PASS\n' "$d0"
done

# Count only: --threads 3 models per-worker storage, with no team spawned.
"$tmp/build" --size --threads 3 --Dc 7 --buckets 2 --D0 20 --KT 19 --LF 16000 --KF 200 \
    --DP4 12 --KP4 9 --DYT 11 --KYT 8 --DYB 12 --KYB 9 --DS3 10 --KS3 8
"$tmp/build" --size --threads 3 --Dc 7 --buckets 2 --D0 24 --KT 23 --LF 18000 --KF 220 \
    --DP4 12 --KP4 10 --DYT 11 --KYT 9 --DYB 12 --KYB 10 --DS3 10 --KS3 8
