# Gold2 independent tail and transfer audit

**PASS**, conditional on the separately checked finite rows. The detached source is `b435b4c824e2ba5835f58b85e67a27de7d6a0d49`, clean when inspected. Opened only the private binary copies after `hashes.log` reported both `OK`. The 6,344,902,168-byte payloads have receipt hashes `41565d49c10809eb1f732794822359410af3ed804577476563bef17f27946413` and `f1dcd7ad3e6df82a985eb3dd6fd793b5a3e38c78368d94e4b135ab174ab92876`. No full table verifier rerun here.

## Signed seam and unbounded domain

For `D>=4000`, let `q=L/D`, `x=1/D`. The true degree gain sends the arguments to `((q+x)/(1+x),x/(1+x))`; bud sends them to `(q-x,x)` at rate `q`, hence no bud at zero leaves. The gain divided by `D` lies between `q` and `(1+q)/2` because each of the `c=D-L` physical core neighbors sends at most `1/2`. No physical arity is capped in that statement. The newborn charge uses `f(J)+g_inf+(252/D)*(g_T-g_inf)` at `T=min(c+1,23)`. There is no newborn higher motif: the promoting vertex would have to be internal and Gold2's largest internal cutoff is 80. Existing motifs can still have a high-degree endpoint, and their charges stay in their own rows. The signed upper envelope takes the unsigned minimum newborn edge for lower (class 11, slope `189238964292`), maximum for upper (class 23, slope `487806245640`).

The independent rational oracle forms the cleared residual by direct evaluation of the gain and bud jumps and inverse Bernstein evaluation matrices, not by importing production coefficients. The bidegrees are bounded by `(2,3)` upper and `(3,4)` lower after multiplication by `(1+x)^2`; therefore 16 closed `q` panels at both gain-rate endpoints enclose every real `D>=4000` under the physical rate bounds, including the `x=0` limit. A negative polynomial bound is divided by the maximum `(1+1/4000)^2` before integer ceiling. Exact off-node identities, concrete states up to `D=10^30`, saturation probes, and endpoint domain controls pass. Lower zero is at relaxed no-leaf, low-receive corner; this does not require a finite no-leaf hub with zero receive.

| Side | Analytic cap | Upward tick cap | Maximum exact seam residual | Upward seam cap |
|---|---:|---:|---:|---:|
| Lower | `0` | `0` | `0` at `c=3999` | `0` |
| Upper | `-122011194591/80040005` | `-1524` | `158899595482883/17768890000000` at `c=3918` | `9` |

Seam checks covered exactly 4,000 arities **per side**, both gain extrema, and independent rounding of all signed terms to sub-ticks before integer ceiling. The gain destination is analytic at `D=4000`; the bud destination stays tabular at `D=3999`. Seam caps fit the receipt's tabular caps `55/2359`. Comparison independently encloses the parent edge over 127,489 `(tail class,parent corner)` combinations per side, including the degree-affine segment from the cut to infinity and the root's zero edge. The resulting edge intervals are `[−175650701124,0]` and `[−175646937559,24506]`; shifted analytic minima are `99227205820` and `979444199938`. The upper tail maximum must include **positive 24506** parent-edge ticks.

## FY change and proof diff

`Q_K=10080`. The old low positive knots `56,112` in numerator units are replaced by `40,80,120` (one additional knot overall); the common rest begins at 168. The first interval stays exactly `[0,1/d2]`, now `[0,1/252]`. Therefore interpolation of a real tail reciprocal yields exactly the finite-corner/INF average on that first segment. The signed dual maximum over all 20 knots bounds each interpolation segment; the exact ordered-pair cancellation is unchanged when the same observed donor is used on both sides. Tested actual FY rows `(80,2),(126,2),(251,22)` at exact and interior knots and up to `D=10^30`. The full finite FY inequalities remain the production verifier's responsibility.

`git diff b04cc50 HEAD -- round2/neon2/quintet-cert/PROOF.md` has three substantive clarifications, all valid:

- A saturated endpoint bud has **zero weight increment**, not necessarily zero hazard. The real class stays clamped at `K+1` during that transition. Exact-class gain and promotion bounds remain monotone in `D`; no positive uniform lower bound for overflow gain is used. Cross-donor sends on a double-tail edge stay within the full gain factors.
- High-degree vertices do not **host** higher motifs under the comparison allocation. Finite internal hosts still account for incident motifs containing arbitrarily large endpoints. This is the correct reason the tail comparison uses only singleton plus assigned parent edge.
- The first jump from the two-vertex seed is `Exp(2)` and yields the unique three-vertex tree. Convolution preserves upper growth with multiplier at most `2/(2+lambda_up*)` (positive denominator here), and restriction of the lower integral to `[0,1]` yields a positive constant. The discrete recursion follows directly from jump rate `n`, conditional expectation and positive comparison; its products grow as `n^lambda*`. For continuous time, stopping at size `M` is removable by the Yule second moment and `|generator|<=const*N*(2N+1)`, before Grönwall. No finite-degree cutoff is inserted in the generator.

## Exact receipt transfer

With the separate finite premises, `cap_f=55/2359`, total errors `19217/20062`, and the original `m_min` produce `[35083855339718511/62017003637500000, 127993535256830043/226242120050000000]`. Outward at 12 decimals: `[0.565713486333,0.565736986679]`. The strict gap below `142/251` is `3697635659207/56786772132550000000`. The target `99/175` is **inside**: it lies `347026470423/434119025462500000` above lower (about `7.994e−7`) and `35951399810301/1583694840350000000` below upper.

Replay from `[local path]`: `taskset -c 8 nice -n 17 timeout 40s python3 -Bu round2/attacks/neon2-quintet/gold2/tail/check.py`. It reuses the independent rational algorithm without opening Gold1 data. Actual final line: `PASS Gold2 tail, seam, FY, comparisons and exact receipt transfer`.
