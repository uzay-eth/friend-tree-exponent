# Replaying Neon-2 Claim 6

**HOLDS.** The two exact endpoints are $`2825181729593748861/4994114995300000000`$ and $`439124781075521909/776183857175000000`$. Their outward enclosure is **[0.565702177913, 0.565748407438]**. Astra reviewed the original unoriented five-vertex potential at `ft2-neon2@b04cc50c1e391e0420eb7b006ac781d484c68543`. The official board report records the verdict and attribution. [Board verdict](../../../../attacks/astra-on-neon2-6.md)

## Clean production replay

The source checkout, compiler outputs and private copies of the original binary snapshots were kept under Astra’s ignored work directory. The following commands use a fresh directory and two verifier threads. Run them from the Astra repository root.

```sh
astra_root=$PWD
claim6_work=$(mktemp -d "$astra_root/round2/.work/neon2-claim6-XXXXXXXX")
git clone --no-local --no-checkout --single-branch --branch exact/neon2 [local path] "$claim6_work/source"
git -C "$claim6_work/source" checkout --detach b04cc50c1e391e0420eb7b006ac781d484c68543
cp --reflink=auto round2/outputs/quintet-gold1/quintet_{lower,upper}.bin "$claim6_work/"
cd "$claim6_work"
sha256sum --strict --check - <<'PINS'
81dc1f250d49083716e0998dcb63b615bd55401c4e35ea388ce098d29b4a3986  quintet_lower.bin
8a36c3afdd2a8f4bbc03461a306af04594c46ada6ef037646261da4d364d9383  quintet_upper.bin
PINS
quintet_src=source/round2/neon2/quintet-cert
g++ -std=c++17 -O3 -fopenmp -Wall -Wextra -Werror -pedantic "$quintet_src/verify_quintet.cpp" -o verify_quintet
OMP_WAIT_POLICY=PASSIVE OMP_NUM_THREADS=2 /usr/bin/time -v -o gold.time nice -n 14 \
  python3 -IBu "$quintet_src/cert_check.py" --verify-bin ./verify_quintet \
  --lower quintet_lower.bin --upper quintet_upper.bin --threads 2 > gold.json
cmp gold.json "$quintet_src/certificate_gold1.json"
python3 -IBu "$astra_root/round2/attacks/neon2-quintet/size_epoch.py" \
  source quintet_lower.bin quintet_upper.bin
for part in physical tail transfer
do
  ASTRA_QUINTET_AUDIT="$claim6_work" python3 -IBu "$astra_root/round2/attacks/neon2-quintet/$part/check.py"
done
```

Both hash checks printed `OK`. The release compiler emitted zero warnings. The verifier/combiner and byte comparison each exited zero. The receipt SHA-256 was `d26d22a31d8f18790f517c451708a35b2721cbb4809613d1316c05df5071faeb`. Elapsed time was 116.88 s with 84.07 s user CPU, 21.96 s system CPU and 3,143,492 KiB maximum RSS.

| Reproduced entry | Lower | Upper |
|---|---:|---:|
| Nominal discount | 0.56570237 | 0.56574828 |
| Finite drift rows | 191,375,747 | 191,375,747 |
| Projection inequalities | 29,879,919 | 29,879,919 |
| Analytic singleton cap, ticks | 0 | −6509 |
| Total signed error cap per core, ticks | 19,186 | 19,783 |
| Positive comparison floor per core, ticks | 99,882,299,906 | 155,236,771,435 |
| Comparison ceiling per core, ticks | 1,808,571,640,916 | 1,808,494,729,937 |

The wrapper recomputes the analytic tail and requires the C++ digest to match the bytes used by the Python calculation. Its separately implemented physical oracle passed its 47 trees on each side. Our replay started from the published snapshots. The full proposal builder is outside the trust chain.

## Proof checks

Let $`D`$ be total degree and $`c`$ core degree. A gain changes $`D`$ at rate $`D-c+\sum_{w\sim u, \,w\ {\rm core}}1/D_w`$. A promotion changes $`c`$ at rate $`(D-c)/D`$, inserting one core with state $`(2,1)`$. Deleting the unique newborn identifies the lower-order source of each newly created motif. Existing supported embeddings can leave through an internal ceiling as the degrees increase. The eight families count unordered embeddings, including the separate choices of equal-looking branches. Our physical evaluator checks changes before centering or redistributing any generator terms. [Labelled checks and counterexamples to faulty multiplicities](physical/README.md)

An endpoint with $`D\ge d`$ has mass $`d/D`$ on its capped-arity knot and the rest at infinity. Averaging its reciprocal gives exactly $`1/D`$. The projection to a coarser knot composes these probabilities. At an exact hub outside the three-path support, a transfer of $`A(1/D_v)-A(1/D_w)`$ cancels when the ordered incident pair is reversed. At a high-degree endpoint $`A`$ is affine on the first segment. On the other side of the transfer the maximum of the 19 signed affine functions is convex in the edge potential change. This gives the corner bound for either witness sign. We reconstructed the projection inequalities and multiplicities for the path, broom and four-star. [Transfer proof and fixtures](transfer/README.md)

Above total degree 2000 the singleton has an explicit polynomial in $`q=(D-c)/D`$ and $`1/D`$. The adjacent newborn edge uses the post-promotion arity. Multiplying the residual by $`(1+1/D)^2`$ yields a polynomial in two variables. Bernstein controls at both receive-rate endpoints enclose it on the full continuous rectangle. Negative caps are divided by the largest denominator before upward rounding. All 2,000 arities at $`D=1999`$ and all 16 tail panels passed an independent rational evaluation. [Analytic transfer and seam checks](tail/README.md)

For comparability, root the finite core tree at any vertex and give each nonroot its parent edge. The root can take the zero edge value. A four-path contributes half its value to each internal endpoint. The five-path belongs to its middle, and the broom and four-star to their hubs. The identity $`\sum(c-2)=-2`$ permits the lower witness’s degree shift, with an additive $`2\alpha S`$ left over. A tree with one core has its singleton value directly. The finite comparison floors and analytic tail bounds then give $`m_{\min}C\le S\Phi\le m_{\max}C`$. The finite verifier’s signed cap is the maximum over every supported local state, including outward rounding at intermediate divisions. The wrapper combines the cap with the comparison floor.

The continuous-time chain has total jump rate $`N`$, its total vertex count. On stopping at vertex count $`M`$, the generator identity holds on a finite state space. The finite-horizon Yule second moment bounds $`N(2N+1)`$, allowing removal of the stop and the usual differential inequalities. At deterministic size $`n`$, the conditional expectation satisfies
```math
\mathbb E[\Phi(T_{n+1})\mid T_n]=\Phi(T_n)+(\mathcal L\Phi)(T_n)/n.
```
Put $`B_n(a)=\prod_{j=3}^{n-1}(1+a/j)`$. With $`\Phi_j`$ expressed in ticks, the verified rates imply
```math
\frac{\Phi_{\rm lo}(T_3)}{m_{\max,{\rm lo}}}\,B_n(\lambda_{\rm lo}^*)
\ \le\ \mathbb E C_n\ \le\ \frac{\Phi_{\rm up}(T_3)}{m_{\min,{\rm up}}}\,B_n(\lambda_{\rm up}^*).
```
Since $`\log B_n(a)=a\log n+O(1)`$, these are unconditional bounds for the discrete liminf and limsup. To name their common exponent the pinned project cites its round-one existence theorem. That file credits Astra’s argument. Existence proof (`proofs/existence.md`)

Our fixed-size replay enumerated every free tree through size 12. For the small-size product check it uses the pinned potential reader and independently sums the physical directed-arrow probabilities and the exact tree law. On each side it passed 434 conditional inequalities, 985 comparisons and nine finite-product tests. The uniform-recipient replacement violated all 434 upper conditional inequalities. The separate physical checker uses its own sparse binary reader and matched both endpoint inequalities on 22 labelled trees per side, with 20,994 directed birth slots. Size-law replay (`size_epoch.py`) · [Independent physical reader](physical/README.md)

## Negative controls

The unmodified production executable rejected two altered, rehashed copies at the edge-to-three-path envelope. On a raw edge coefficient change the cleared discrepancy exceeded its zero radius by 17,713,784,569,643. A changed projection midpoint gave 35,280,594,571,876 against zero. Both processes exited 2 at the numerical inequality. A separate header calculation reproduced all nine drift-row and three projection counts and both corrected fractions. [Commands, hashes and diagnostics](controls/README.md)

We also built the small preset and reran all 46 published controls. The 19 tests per side exercise specific error paths and increases in the certified residual. Eight additional tests check the header, side selection and input binding. The suite exited zero in 14.50 s.

```sh
g++ -std=c++17 -O3 -fopenmp -Wall -Wextra -pedantic "$quintet_src/build_quintet.cpp" -o build_quintet
OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 ./build_quintet tiny lower tiny-lower.bin
OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 ./build_quintet tiny upper tiny-upper.bin
OMP_WAIT_POLICY=PASSIVE python3 -IBu "$quintet_src/cert_check.py" --verify-bin ./verify_quintet \
  --lower tiny-lower.bin --upper tiny-upper.bin --threads 1 --run-mutations > tiny.json
```
