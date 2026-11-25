#!/usr/bin/env bash
# Full, fail-closed replay of the refined lower witness (paper, Theorem 2 and Appendix F)
# with the undefined-behaviour sanitizer in non-recovering mode.
# Usage: bash verify/refined-lower/replay_full.sh WITNESS NEW_OUTPUT_DIR [THREADS]
#   WITNESS  the 13,046,351,000-byte lower table (SHA-256 80e3b80a...afd4ad28); see README.md for regeneration
#   NEW_OUTPUT_DIR must not exist yet; all builds and outputs go there.
# Needs: bash, g++ (C++17, OpenMP, libubsan), patch, sha256sum, cmp, GNU time at /usr/bin/time, python3 >= 3.10.
# Measured (16-core x86-64, g++ 13.3, 6 threads): about 6 minutes, 13.4 GB peak RSS (README.md).
set -euo pipefail
[ $# -ge 2 ] || { echo "usage: $0 WITNESS NEW_OUTPUT_DIR [THREADS]" >&2; exit 2; }
PAYLOAD=$(readlink -f -- "$1"); OUT=$2; T=${3:-6}
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
QC=$(cd -- "$HERE/../../round2/proofs/neon2/quintet-cert" && pwd)
PAYLOAD_SHA=80e3b80a9f962b458854828abab5b00e2806a6c1ce709835b5e884b0afd4ad28
RATIO_SRC_SHA=f85a172ab45195d0c9d49d809f0c903c146f0bbc570cb170b5887b493954f067
die(){ echo "REPLAY FAILED: $*" >&2; exit 1; }
[ ! -e "$OUT" ] || die "output directory $OUT already exists"
mkdir -p -- "$OUT"; OUT=$(readlink -f -- "$OUT")
export LC_ALL=C OMP_WAIT_POLICY=PASSIVE UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1
FLAGS=(-std=c++17 -O2 -g -fopenmp -fsanitize=undefined -fno-sanitize-recover=all)

echo "== 1. package manifest"; (cd "$HERE" && sha256sum --quiet -c MANIFEST.sha256) || die "package files differ from MANIFEST.sha256"
echo "== 2. payload digest and LOWER header"
[ "$(sha256sum < "$PAYLOAD" | cut -d' ' -f1)" = "$PAYLOAD_SHA" ] || die "payload digest is not $PAYLOAD_SHA"
python3 -I -c 'import struct,sys; h=open(sys.argv[1],"rb").read(512); m=struct.unpack_from("<q",h,8)[0]; sys.exit(0 if m==2 else 1)' "$PAYLOAD" || die "payload header mode is not 2 (LOWER)"
echo "== 3. sources"; mkdir "$OUT/src"; cp "$QC"/{verify_quintet.cpp,quintet_format.hpp,cert_check.py,verify_tail_bernstein.py} "$HERE/ubsan_probe.cpp" "$OUT/src/"
cp "$OUT/src/verify_quintet.cpp" "$OUT/src/verify_ratio.cpp"; patch -s "$OUT/src/verify_ratio.cpp" "$HERE/verify_ratio.patch"
[ "$(sha256sum < "$OUT/src/verify_ratio.cpp" | cut -d' ' -f1)" = "$RATIO_SRC_SHA" ] || die "patched verifier digest mismatch"
echo "== 4. UBSan builds and planted-overflow control"
cd "$OUT/src"
g++ "${FLAGS[@]}" ubsan_probe.cpp -o ubsan_probe
if ./ubsan_probe 1 > probe.out 2> probe.err; then die "planted overflow was not caught: sanitizer inactive"; fi
grep -q "runtime error: signed integer overflow" probe.err || die "planted overflow control did not report signed overflow"
g++ "${FLAGS[@]}" verify_quintet.cpp -o verify_quintet_ubsan
g++ "${FLAGS[@]}" verify_ratio.cpp -o verify_ratio_ubsan
echo "== 5. instrumented verifier (per-class dump) on the witness, $T threads"
rc=0; AUX_OUT="$OUT/aux.txt" /usr/bin/time -v -o "$OUT/ratio.time" ./verify_ratio_ubsan "$PAYLOAD" --threads "$T" > "$OUT/ratio_verify.json" 2> "$OUT/ratio.err" || rc=$?
[ "$rc" = 0 ] || die "instrumented verifier exit $rc (see $OUT/ratio.err)"
[ ! -s "$OUT/ratio.err" ] || die "instrumented verifier wrote to stderr (sanitizer or error)"
python3 -I - "$OUT/ratio_verify.json" "$HERE/receipt.json" <<'PY' || die "instrumented verify JSON differs from the receipt"
import json, sys
a = json.load(open(sys.argv[1])); r = json.load(open(sys.argv[2]))
sys.exit(0 if a == r["lower"]["verify"] and a["mode"] == 2 else 1)
PY
cmp -s "$OUT/aux.txt" "$HERE/classes.txt" || die "fresh per-class dump differs from the pinned classes.txt"
python3 -I "$HERE/check_ratio.py" --aux "$OUT/aux.txt" > "$OUT/check_ratio.out" || die "check_ratio.py on the fresh per-class dump failed"
echo "== 6. unmodified checker on the unmodified verifier (UBSan), LOWER side only, plus the 47-tree oracle"
rc=0; /usr/bin/time -v -o "$OUT/check.time" python3 -B - "$PAYLOAD" "$OUT/src/verify_quintet_ubsan" "$T" "$PAYLOAD_SHA" "$HERE/receipt.json" > "$OUT/lower_side.json" 2> "$OUT/check.err" <<'PY' || rc=$?
import json, sys
sys.path.insert(0, ".")
import cert_check as cc
pay, vbin, T, sha, rec = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], json.load(open(sys.argv[5]))
lo = cc.combine_side(pay, vbin, threads=T, expected_mode=2, return_raw=True, expected_sha256=sha)
if not lo.get("ok"):
    print(json.dumps({k: v for k, v in lo.items() if k != "_raw"}), file=sys.stderr); sys.exit(1)
raw = lo.pop("_raw"); got_sha = lo.pop("sha256")
oracle = cc.run_tree_oracle(pay, lo, raw=raw)
print(json.dumps({"sha256": got_sha, "lower": lo, "tree_oracle_lower": oracle}, indent=1))
ok = got_sha == rec["lower_sha256"] == sha and lo == rec["lower"] and oracle == rec["tree_oracle_lower"]
sys.exit(0 if ok else 3)
PY
[ "$rc" = 0 ] || die "unmodified checker LOWER side did not reproduce the receipt (exit $rc, see $OUT/check.err)"
[ ! -s "$OUT/check.err" ] || die "checker wrote to stderr"
echo "== 7. summary"
grep -E "Elapsed|Maximum resident" "$OUT/ratio.time" "$OUT/check.time"
sha256sum "$OUT/aux.txt" "$OUT/ratio_verify.json"
cat "$OUT/check_ratio.out"
echo "FULL REPLAY PASSED"
