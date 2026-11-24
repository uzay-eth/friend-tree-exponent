# Gold2 exclusion of reduced rationals through denominator 425

**Pinned cold-review verdict: HOLDS**, conditional on the original pinned Gold2 certificate and the physical stable-age law. I cold-read the writer's `round2/occupation-quintet-gold2/rationals/README.md` and `check.py`, verifying their Git blobs `c9914db75b756b1d081c258b5108093da5bf4229` and `d8832e99ec690020507ca61860a4fba4a9c31f05`. Only the original Gold2 potential and its original per-core allocation are used. The physical floor replay and local-cap check are corroborated by `gold2-h3-sol.md` (d4d3).

## The three credits are additive

Set `S=2^40`, and use the original lower sign `σ=−1`. In the physical signed-row sum, FY terms in the non-M3 exact-G2 host case were bounded separately for each ordered `(host, payer, distinct donor)` incidence by a maximum over the 20 knots. The difference between that maximum and the actual signed physical term is nonnegative even for a donor reciprocal between knots, by piecewise affinity. The small-payer floor takes `D_payer≤5` and `18<D_host≤251`, host arity at most 17, with the donor's degree-12 absorbing continuation. The older-payer floor takes `6≤D_payer≤34`, payer arity at most 7, `D_donor≤5`, `D_host≤150` and host arity at most 10, with `max(D_host,D_payer)>18`. Both require exact G2 outside M3 at reward time. Opposite orientations of the same physical edge are separate FY incidences, not duplicate copies of a credit; the payer split separates the actual ordered atoms.

H3 instead uses the original signed physical H3 row after transfers, including the outgoing terms where higher support is present. Its middle is in M3, with degree at most 5 and arity at most 4; one arm is its actual nonfounder mother of degree 19–79 and arity at most 8, the other an actual child of degree at most 18 and arity at most 3. The roles select each credited unordered pair of arms at most once. The solver's `--kind threshold` call leaves `reserve=False`: it pays the negative part of the outward row itself, not the whole unused family cap. For the original **local nonnegative middle cap** `K=max(0,cap_h3_m[m])`, `h+max(0,−h)=max(0,h)≤K`. This preserves the original aggregate `E` without replacing a local cap by the global maximum 408. No other H3 reserve or host relaxation is added. Repeated appearances of a vertex across motifs do not merge the accounts: the transfer into H3 is made at M3 hosts and the FY rebates are from non-M3 hosts.

Per-core tick floors are based on these exact normalized rewards, converted by multiplying by `S`:

| Lane | Exact potential units/core | Conservative ticks/core |
|---|---:|---:|
| H3 | `1084423030952180046412413909/1704710415481443450880000000000000` | 699000 |
| Small payer | `39526146756920663734311/109951162777600000000000000000` | 395000 |
| Older payer | `1046834122639903380009/21990232555520000000000000000` | 52000 |

The stopped Bellman comparisons use the original Gold2 enclosure: upper discount for nonnegative rewards, and lower birth multiplier once at the true newborn. The H3 replay used the exact endpoints and both FY replays used the outward decimals of the same endpoints. All three replays reproduced the exact fractions. This is a stationary improvement; it is not a new uniform pointwise exponent bound on every finite tree.

## Convert the net with the original allocation

The writer's Cesàro proof does not need a stationary limit for the whole potential allocation. With expectation denoted by `𝔼`, set `F(t)=𝔼[SΦ(T_t)]` and `K(t)=𝔼[C(T_t)]`. Stopped Dynkin applies since total event rate is the Yule vertex count and the potential and its successor are bounded by a constant times that count plus one; the generator integrand has a finite-time quadratic majorant. With `Q(t)=e^(−μt)F(t)`, the physical inequality gives `Q′+(μ−λ)Q≥e^(−μt)(𝔼[R]−E K(t))`. The paid scores are bounded oriented local scores: the host caps bound multiplicity and the reward only needs degrees and port incidences in a fixed neighbourhood. Oriented empirical convergence and uniform integrability of `e^(−μt)C(T_t)` (bounded in L² by the comparable core eigenmartingale) yield `e^(−μt)𝔼[R]→k πR≥kr`, where `k=lim e^(−μt)K(t)>0`. The same conclusion follows from the physical spatial mean result for these local scores. Thus the averaged right side is at least `k(r−E)+o(1)`. Since `Q` is bounded, its endpoint difference divided by time vanishes. Since `Q` is positive and `r−E>0`, this first forces `μ>λ`; then `limsup (T^−1∫Q)≤m_max k` gives the claimed denominator. Here the limsup is of the integral average, before multiplying by `μ−λ`. This establishes the desired signed comparison without a separate termwise generator balance. The original two-vertex seed reaches the three-vertex path at its first birth; for continuous means the independent `Exp(2)` delay just multiplies the limiting `e^(−μt)K(t)` coefficient by `2/(2+μ)`, and the discrete chain from size three is identical.

The original receipt has `λ=7071421/12500000`, `m_max=1808603193693`, and `E=55+792+18370=19217`. The sum of conservative floors is `1146000`, so the **net** is `1126783>0`. Consequently `μ≥λ+1126783/m_max`. Dividing this positive improvement by `m_min` would overstate it; the original outward correction instead had a negative net and used `m_min`. At `99/175`, the exact target shift from `λ` is `53/87500000`, so the margin in ticks per core is `2737543234271/87500000>31286`. The exact lower endpoint is `12789408689335247753/22607539921162500000`, strictly greater than `99/175`. For the upper certificate, `28286843/50000000+20062/158369484035` reproduces `127993535256830043/226242120050000000`, strictly below `142/251` by `3697635659207/56786772132550000000`. No upper-side improvement is required.

## All reduced rationals, not just the target

In increasing order, `155/274`, `99/175`, `142/251` are consecutive in Farey order 425. The left fraction is the mediant of `56/99` and `99/175`; using the wider left interval directly would not suffice because its denominator sum is only 274. For either adjacent pair `a/b<c/d`, `bc−ad=1`; for a fraction `p/q` strictly between, the positive integers `bp−aq` and `cq−dp` satisfy `q=d(bp−aq)+b(cq−dp)≥b+d`. The two sums are 449 and 426, respectively. The original lower endpoint is strictly above `155/274`, so the original bracket contained only `99/175` through order 425. The strengthened lower endpoint excludes it; the unchanged upper remains below the right neighbor. Exact cross-multiplication and independent closed-interval enumeration verify both claims. The denominator-426 mediant `241/426` is still inside the strengthened bracket and is a boundary control on what these endpoints prove.

## Replays

From the Astra root:

```text
nice -n 15 timeout 120s python3 -IBu round2/occupation-quintet-gold2/post-fy/replay.py --kind threshold
  H3 normalized 1084423030952180046412413909/1704710415481443450880000000000000; 271350 states; 187105 paid rows
nice -n 15 timeout 105s python3 -IBu round2/occupation-quintet-gold2/edges/fy.py --host 251 --stock 17 --donor 12 --absorb
  PASS exact original Gold2 knot and physical Bellman replay; 3585120 states
nice -n 15 timeout 115s python3 -IBu round2/occupation-quintet-gold2/payer/price.py --host 150 --stock 10 --receiver 34 --rstock 7
  PASS pinned lower bytes and exact Bellman inequalities; 3946320 states; 516312 exact reward labels
python3 -IB round2/checks/gold2-rational-sol.py
  PASS all reduced rationals with denominator <= 425 excluded
```

The checker authenticates the original JSON receipt by SHA-256 before comparing the two signed endpoints. These are exact arithmetic and read-only Bellman replays on the existing data; no refit or LP was used.

The final writer checker, `python3 -IB round2/occupation-quintet-gold2/rationals/check.py`, passed at the pinned blobs. It authenticates the recorded receipts and recomposes the exact result; it does not redo the physical Bellman proof, which is separately replayed above. Independent unrounded recomposition gives lower endpoint `126905728962161948969119589499374909399/224328301684105796116575000000000000000` and target margin `28066074690162200126495624365793/868238050925000000000000000 > 32325` ticks/core. The exact enclosure contains `241/426` strictly, and this is its only reduced fraction with denominator at most 426. Sharpness is for this enclosure, not evidence that the exponent has denominator 426 or is rational.
