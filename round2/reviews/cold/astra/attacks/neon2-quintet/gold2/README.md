# Astra replay of the second Neon-2 five-vertex enclosure

**Verdict on original gold2.** The verifier at b435b4c824e2ba5835f58b85e67a27de7d6a0d49 proves
```math
 35083855339718511/62017003637500000\le\mu\le
 127993535256830043/226242120050000000 .
```
The outward decimals are **[0.565713486333, 0.565736986679]**. The report includes the continuous-time and discrete-size transfer for the original edge seed. The common exponent uses the earlier existence theorem. The upper fraction is below $`142/251`$ by exactly $`3697635659207/56786772132550000000`$. A denominator scan through 425 leaves $`99/175`$ in this original enclosure. The audit concerns the posted five-vertex certificate. Production receipt (`production.json`) · [Dated board verdict](../../../../../attacks/astra-on-neon2-6.md)

## Clean production run

Run from the Astra repository root. Clone into a new private directory and copy both original snapshots. The two published hashes were checked on the copies before compiling the verifier or running it.

~~~sh
w=$(mktemp -d "$PWD/round2/.work/neon2-gold2-XXXXXXXX")
git clone --no-local --no-checkout --single-branch --branch exact/neon2 [local path] "$w/source"
git -C "$w/source" checkout --detach b435b4c824e2ba5835f58b85e67a27de7d6a0d49
cp --reflink=auto round2/outputs/quintet-gold2/quintet_{lower,upper}.bin "$w/"
chmod a-w "$w/quintet_lower.bin" "$w/quintet_upper.bin"
(
  set -e
  cd "$w"
  sha256sum --strict --check - <<'PINS'
41565d49c10809eb1f732794822359410af3ed804577476563bef17f27946413  quintet_lower.bin
f1dcd7ad3e6df82a985eb3dd6fd793b5a3e38c78368d94e4b135ab174ab92876  quintet_upper.bin
PINS
  src=source/round2/neon2/quintet-cert
  g++ -std=c++17 -O3 -fopenmp -Wall -Wextra -Werror -pedantic "$src/verify_quintet.cpp" -o verify_quintet
  OMP_WAIT_POLICY=PASSIVE OMP_NUM_THREADS=2 /usr/bin/time -v -o gold.time \
    nice -n 18 python3 -IBu "$src/cert_check.py" --verify-bin ./verify_quintet \
    --lower quintet_lower.bin --upper quintet_upper.bin --threads 2 >gold.json
  cmp gold.json "$src/certificate_gold2.json"
  sha256sum gold.json
  git -C source status --short
)
python3 -IBu round2/attacks/neon2-quintet/gold2/size_law.py "$w"
~~~

The compiler issued no warnings. Both files contained 6,344,902,168 bytes. The complete production command exited zero in **5m22.14s**, with 324.19s user CPU, 36.25s system CPU and 12,414,080 KiB maximum RSS. The owner used two niced threads on CPUs 9 and 10. Both hash checks printed OK. The byte comparison exited zero, the source checkout remained clean and the receipt SHA-256 was 5e11f9c0b6eb3bd43e880be608c37d6f0b1fdbcc632b91bb5e8037db52e13284.

| Reproduced item per witness | Lower | Upper |
|---|---:|---:|
| Finite drift rows | 767,588,844 | 767,588,844 |
| Envelope inequalities | 106,953,150 | 106,953,150 |
| Analytic singleton cap, rounded up in ticks | 0 | −1524 |
| Total outward allowance, ticks per core | 19,217 | 20,062 |
| Minimum potential, ticks per core | 99,227,205,820 | 158,369,484,035 |
| Largest physical seam cap, rounded up in ticks | 0 | 9 |

The compiler input, its format header, the rational tail checker, the wrapping combiner and the additional independent audit script are byte-identical to the gold1 verification files. The updated builder emits larger tables. In the proof, a saturated arity makes the promotion weight change zero while its physical promotion rate can be positive. The revised comparison charges each higher motif to its actual host. Starting from the edge, the first jump at rate two produces the three-vertex tree. Integrating the exponent bounds from that time preserves both exponential rates. At deterministic size the continuous generator gives the one-step conditional mean directly.

## Independent physical and adverse checks

Astra's sparse evaluator checked 24 labelled trees per witness, including all motif families, joined high degrees, and the singleton seam. It checked 46,928 directed birth choices, 23,488 recipient updates, 700 signed knot checks and 32 seam rows per sign. A separate rational calculation reconstructed all 4,000 possible seam arities and the 16 Bernstein tail panels at both receive-rate endpoints. The three new low knots $`1/252,1/126,1/84`$ replace $`1/180,1/90`$, for a total of 20. [Physical trees and knots](physical/README.md) · [Tail and seed transfer](tail/README.md)

The size-law replay enumerated all 985 free trees of sizes 3 through 12. For each sign it checked 434 one-birth conditional means, all 985 potential comparisons and the nine exact products against the evolving discrete law. Uniformly choosing the recipient in a deliberately altered model violated all 434 upper inequalities. Size-law replay (`size_law.py`)

A private lower-snapshot mutation exchanged the two zero-radius projection centres at bytes 6,141,599,992 and 6,141,600,000. Both values remained in range and the complete coefficient multiset was preserved. Independent arithmetic found that both replacements violate their env23 inequalities. On this new file with its own measured digest, the stock verifier stopped at env23 with exit 2, empty stdout and a numerical error. Each inequality exceeded its radius by exactly 82,239,727,528 ticks. [Mutation, digest and production output](controls/README.md)
