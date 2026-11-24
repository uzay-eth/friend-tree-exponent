> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits and Claude replay of earlier freezes; original Neon-2 Gold1 witnesses and inherited stationary results. See the [release index](../../../../README.md#claim-25) for exact scope and citations.

# Physical singleton lifetime, original Neon-2 quintet

**Both ends improve.** The exact-rational lifetime calculation gives the outward interval **`[0.565702523889859, 0.565747980742416]`**, inward by `3.45975945e−7 / 4.26695175e−7` from the exact original endpoints. A separate all-degree, all-arity upper control shows that reaching `5e−6` on either end is impossible using only this original singleton slack and the required original denominator rule. No other family or older witness enters the calculation.

## Charged row and transfer

Use the pinned accepted `b04cc50c1e391e0420eb7b006ac781d484c68543` witness, the two authenticated binaries under `round2/.work/replay-inputs`, and the sibling sparse [`Snapshot`](../census/census.py) reader. This uses the original signed row allocation: a centered singleton uses `c*r1` in its *charge* (the omitted actual-minus-centered term is assigned to incident edges by the original proof); the physical descriptor always uses the real gain rate. The singleton promotion owns its own weight jump and the newborn `f(J)+g((D,c+1),J)` at actual post-promotion arity and degree. The newborn edge is never charged a second time. The signed score in ticks is `F=cap_f+A+B*rho` with `A+B*rho` the inward signed trial residual from the reader. The original `cap_f` is **0 lower, −1053 upper**. Both corners of the exact score are nonnegative on every examined census; all unexamined scores are nonnegative by the accepted original verifier and analytic audit.

Every physical nonfounder starts at `(D,c)=(2,1)`, with the actual mother counted as one core neighbor. At `(D,c)`, `L=D-c`, the transition to `(D+1,c)` has rate `L+rho`, and the transition to `(D,c+1)` has rate `b=L/D`. The actual receive `rho` may adapt to the full history; it always belongs to `[0,c/2]`. The solve uses **the same** rate in the instantaneous reward, holding time and gain continuation. In grid sub-ticks `10^-12`, round both affine coefficients downward for the lower comparison (upward for the upper comparison), preserving order because `rho>=0`. The two backward quotients at the endpoints have positive denominators. After rounding the selected value, test the cross-multiplied inequalities with integers at both receive corners. Selected cases, including both cutoff joins, are also checked directly using rational quotients. Rows beyond the arity cap have zero masked reward and value; no vertex returns from overflow. Rows beyond the degree cut either have zero continuation or pass one certified analytic continuation once, at the actual degree exit without changing arity. The current code makes no assumption about independence of physical receive, arity or lifetime.

Under the [stable-age law](../../../../../astra/source/round2/degree-law.md#birth-times), the expected per-core score is `mu` times the expected discount-`mu` lifetime over the genuine physical birth law. The universal newborn lower bound is calculated at the **exact original upper endpoint** and multiplied by the **exact original lower endpoint**. Conversely the upper calculation uses the lower discount and upper normalization. The founder is outside the paid score; it need not be removed from the original drift bound. More explicitly, sum the original capped drift allocations and subtract this nonnegative nonfounder reserve once. In normalized Laplace occupation as the discount decreases to `mu`, [integration by parts](../../../../../astra/source/round2/occupation/physical-test.md) gives

```
sign*(mu-lambda)*mean_potential_ticks >= stationary_F_floor - err_per_core
```

with `sign=+1` lower and `-1` upper. Here `err_per_core=19186/19783` by the original cap aggregation. Use original `m_min=99882299906/155236771435` for a negative net and original `m_max=1808571640916/1808494729937` for a positive net. This is an improvement of the original certified interval, not a shift from a D18/D20 interval or an addition of a second reserve over the same rows.

## Unbounded continuation and ceiling

For fixed `c` and analytic degree `d>=D` write `x=1/d`. At each receive corner, the exact `P(x)=(1+x)^2 F` is a polynomial from the analytic gain and bud increment, with the **actual fixed arity** clamped only for the newborn edge lookup. The replay checks `P` against the sparse row at `D,D+1,10^6` and uses exact one-variable Bernstein coefficients for the entire interval `0<=x<=1/D`. Since `P(0)>=0`, the lower envelope of `P/x` is bounded by the Bernstein minimum for `(P-P(0))/x + D*P(0)`. With nonnegativity of the original score this yields `F>=k[c]/d` after dividing by the maximum `(1+1/D)^2`. For the lower continuation `h[c]/d`, gain costs at most `h[c]/d` in absolute drift, using real gain rate `<=d`, and arity increments occur at a rate in `[1-c/D,1]`. Thus the inequalities for `h` are `(HI+1+b)h[c] <= k[c]+b*h[c+1]` at both promotion corners. The suffix is bounded and this proves a lower comparison for arbitrarily large finite degrees.

For the upper control, the same polynomial yields `F <= k0[c]+k1[c]/d`: use `k0=P(0)` and the nonnegative Bernstein maximum of `(P-P(0))/x`. Compare with `a[c]+h[c]/d`. Gain has no constant-part cost and decreases `h[c]/d` at relative rate at least `(D-c)/(D+1)`, since the physical gain is at least `d-c`. Two outward Bellman tests at the promotion corners give the two parts of the supersolution. At `c>64`, skipped stationary reward is at most `65000*S*(1+LO)^(-64)` ticks. Indeed a vertex must have had at least 64 accepted promotions; the rate-one proposal count at exponential age with rate `mu` has geometric tail at most `(1+LO)^(-64)`. The crude score bound `65000*S` uses `r<=d`, finite weights bounded by `16*S` (so gain before `d=2000` is bounded by `1999*32*S`), one `16*S` newborn edge, and the original discount/weight bounds. In the analytic part `|A|,|B|,|C|,|E|<=16*S`, and the polynomial derivative bound on `q` and `x`, with `|Delta q|<=1/(d+1)` and `|Delta x|<=1/(d(d+1))`, keeps the absolute gain contribution bounded by `100*S`, well within the same crude score bound. This overflow is **added only to the ceiling**, never used as lower credit. The opposite comparison is independent of any improvement in the lower solver; its ceiling concerns the entire physical original singleton account.

## Evidence and replay

For an isolated replay, set `ASTRA_QUINTET_AUDIT` to the absolute directory containing the two original snapshots, `gold.json`, and the detached `source/` at the pin above, laid out as in the [primary replay](../../../../../../reviews/cold/astra/attacks/neon2-quintet/README.md). Use the existing authenticated receipt (or a copy of the pinned `certificate_gold1.json` as `gold.json`); no exhaustive verifier rerun is needed. The default is the original Astra scratch directory. The producer checks the pinned source receipt, and the independent audit hashes its selected snapshot before reading it.

From `.` on one niced core, run the following (drop `--preauthenticated` to rehash both input blobs):

```sh
taskset -c 8 nice -n 17 python3 -Bu round2/occupation-quintet/lifetime/check.py
taskset -c 8 nice -n 17 python3 -Bu round2/occupation-quintet/lifetime/check.py --degree 1999 --arity 28 --tail --preauthenticated
taskset -c 8 nice -n 17 python3 -Bu round2/occupation-quintet/lifetime/check.py --degree 6000 --arity 43 --tail
taskset -c 8 nice -n 17 python3 -Bu round2/occupation-quintet/lifetime/check.py --degree 12000 --arity 30 --tail
taskset -c 8 nice -n 17 python3 -Bu round2/occupation-quintet/lifetime/check.py --degree 6000 --arity 64 --tail --masked-upper --preauthenticated
```

| Lower solve | States per side | Stationary lower reward, ticks (approx.) | Outward interval |
|---|---:|---|---|
| `D<=64`, every arity, no tail | 2,079 | 1960.944 / 2811.051 | `[0.565702197546461, 0.565748389329441]` |
| `D<=1999,c<=28`, unbounded suffix | 55,593 | 256383.952 / 551990.409 | `[0.565702501152090, 0.565747985718021]` |
| `D<=6000,c<=43`, unbounded suffix | 257,096 | 294828.569 / 560430.819 | `[0.565702522408985, 0.565747981050930]` |
| `D<=12000,c<=30`, unbounded suffix | 359,564 | 297506.836 / 560988.765 | `[0.565702523889859, 0.565747980742416]` |

For the last lower solve the discounted newborn values, before multiplying by the exact original `LO`, are `262953589366589949/500000000000` and `991668031155825523/1000000000000` ticks, respectively. The exact endpoints are:

```text
2554773340753153512830889130454314672798089/4516107375986461313847400000000000000000000
5109739952474889043287333305354815204020697/9031830649699395524296100000000000000000000
```

The scan shows why the degree cut mattered: sampled tabular worst-receive slack is usually a few thousand ticks; at the degree-2000 analytic join it reaches `98362299` ticks at `(2000,1)` lower and `227661576` upper. Widening arity from 28 to 43 at degree 6000 changed reward by less than `0.014` tick on either side. The degree-6000 upper control covers 381,983 masked states per sign and gives all-stock ceilings **`4617552 / 5134900` ticks** after including `24701.721409` overflow ticks. Under the original denominator rule these imply inward ceilings **`0.000002734626008 / 0.000002955820730`** relative to the original exact endpoints. That rules out the requested `5e−6` scale for this singleton account by itself even if the lower price is further sharpened.
