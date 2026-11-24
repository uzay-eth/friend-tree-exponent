**A-lower (`e0708ef`): HOLDS. B (`c488e00`): HOLDS. Denominator statement: HOLDS.** If the original exponent is rational, its reduced denominator is at least **426**. This was an adversarial check **inside the GPT family**. I found no surviving use of the withdrawn H3 price in either proof.

### Audit

| Item | Finding |
|---|---|
| **Singleton — 207,990 ticks/core** | Recovers `55 − q_f` from the original outward singleton row, including the once-charged newborn `f+g`. Centering changes the charge, not the physical gain rate. The conditioned solve tracks the live mother and latest child, exits to the broad value at the reached census, and mixes the roots under the actual birth law. The replay gave `207990.215716345`. A replaces this price with the marked singleton calculation; it does not add both. |
| **Small-payer FY — 395,261** | Recovers the gap between the original knot maximum and the value at a distinct donor’s reciprocal, at an exact G2 host outside M3. The absorbing donor band retains the physical incidence and bounds both reward and receive. Actual promotion orders and directed Bellman bounds give the stationary floor. |
| **Wide-payer FY — 52,341** | Recovers the same ordered term on a disjoint payer-degree range: `6…34` rather than `≤5`. The donor is distinct and has degree `2…5`. Both FY streams require a paying-edge endpoint above 18, so A’s full D18 edge rebate does not overlap them. |
| **Actual sibling H3 — 251,635** | Uses the original transferred mother–middle–child row. An extra maternal child is registered in exactly one physical birth order. Its robust reward is rounded down **before** clock scaling and uses unshifted gain continuations. For $`F(z)=[z-B]_+-[-B]_+`$, convexity and $`F(0)=0`$ give $`\sum_i F(z_i)\le F(\sum_i z_i)\le[-q_3]_+`$. Thus simultaneous siblings, and A’s distinct upstream donor, cannot overspend the H3 row. A’s corrected mother-`80…100` seam is disjoint from the compact `19…79` payments. |
| **Stationary drift and ceiling** | The price uses actual newborn law, upper discount and lower birth normalizer, once. Bounded local rewards converge in mean; the unbounded FY donor enters through a uniformly bounded reciprocal score. Stopped Dynkin and the Yule second moment justify the integrated drift. Since $`e^{-\mu t}\mathbb E C_t\to c_E>0`$ and comparability bounds $`Q(t)=e^{-\mu t}\mathbb E P_t`$, the divided boundary term vanishes. The positive net first forces $`\mu>\lambda`$; only then may one use the upper bound on the mean potential. |
| **Local allocation and lifetime recurrence** | The unshifted upper allocation counts one edge at its younger endpoint, allows an absent incoming edge, and assigns higher motifs to their actual hosts, splitting P4 in half. For $`B=M_*-u`$, the stopped inequality $`sV\le B+g(V_{\rm gain}-V)+p(V_{\rm promo}-V)`$ is affine in the physical gain. Both endpoints $`D-c`$ and $`D-c+c/2`$ suffice. Here $`p=(D-c)/D`$, and degree never returns after exit. Stationary balance gives $`\pi B\ge\mu V(2,1)+(s-\mu)\pi V`$. Thus $`M=\lceil M_*-\ell V(2,1)\rceil`$ bounds the **stationary mean** potential per core, giving $`\mu\ge\lambda+888010/M`$. |
| **A’s exact lower trial** | The marked singleton is priced as $`\pi(bV_{\rm mark})`$; the host-cap saving retains the caps that still cover the paid H3 rows. The true physical occupation lifts into the LP, and every column is repriced rationally, with the signed conversion using `min(LO*z,HI*z)`. The positive restricted margin rules out $`\mu\le s`$, including equality; it is not divided to claim a further endpoint. The cold helper parses the historical H3 amount but discards it. |

### Replay and planted errors

All commands and raw outputs are under `round2/children/a28/`; the line-by-line audit is `audit.md` (`round2/children/a28/audit.md`). Both original 6,344,902,168-byte binaries matched their full SHA-256 hashes before replay and were read in place.

- **B:** Both README invocations of `rationals/repaired/check.py` exited 0. `--rebuild-host` printed `PASS rebuilt 299 local ceilings from the authenticated Gold2 lower payload`. The cold replay independently matched all 299 cells, tested 2,792 literal fixtures and 32,580 arrows, and passed 40/340/598 Bellman corners. The singleton, sibling and both FY producer replays matched their published amounts.
- **A:** `host/frontier/check.py --live-sibling --replay-singleton` regenerated the lower marked roots and reproduced byte-identical receipt SHA-256 `7a6ddda08b83574f6b745f11b3e5445ba80a60f47d7ef75c625b093a56cad20f`. Its fixed lower margin was `27.022222885461634` ticks/core at `0.565714467904`. The cold frontier independently reproduced the compact floors and the restricted `99/175` margin of about `247577.388356574` ticks/core. The corrected lower H3 seam independently replayed at `14132.789764062283` ticks/core.
- **Plants:** A 10% small-FY increase changed its floor from `395261` to `434787` and failed the floor-vector assertion even after bypassing the hash. Reproducing the withdrawn extra rate-denominator division changed the exact adverse value from `1/2` to `3/4`; the physical residual became `−1/2`. Applying a double divisor to the 20-state saving instead changed the target margin to `−8871372887299/87500000`, losing the strict comparison.

### Independent exact recomputation

`taskset -c 10 nice -n 5 timeout 890s python3 -IB independent.py` read the frozen data but imported no Astra solver:

| States | $`V(2,1)`$ | $`M`$ | $`888010-(99/175-\lambda)M`$ |
|---:|---:|---:|---:|
| 20 | 619269501408 | 1458274085072 | `3681682957/781250` |
| 170 | 680528936614 | 1423618796411 | `321296970031/12500000` |
| 299 | 686326891945 | 1420338814887 | `2422917810989/87500000` |

The last row gives the claimed B lower endpoint `10043824822832044427/17754235186087500000`. The original Gold2 upper is below `142/251` by `3697635659207/56786772132550000000`. For any fraction $`a/b`$ strictly between the two Farey neighbours, $`175a-99b`$ and $`142b-251a`$ are positive integers. The determinant is 1, so $`b=251(175a-99b)+175(142b-251a)\ge426`$. The exact scan found only `241/426` through denominator 426, and it is inside the certified enclosure. This says nothing about whether $`\mu`$ is rational.

**Limits:** I inherited the original Gold2 production verification and the earlier physical stable-age theorem; I did not re-prove them or audit A’s upper proof. One unmodified cold repair run briefly exceeded the memory budget at 4,360,628 KiB RSS from its read-only mapping. A repeat changing only mapping eviction reproduced every result at 60,004 KiB.
