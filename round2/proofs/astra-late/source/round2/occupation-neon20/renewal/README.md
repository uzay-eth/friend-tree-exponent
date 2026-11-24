> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../README.md#claim-21) for exact scope and citations.

# Physical occupation repricing of the original Neon D20 witness

**Outward bracket: `[0.565681716887067, 0.565766354357632]`.** The movement from the original D20 `[0.565680,0.565768]` is at least `1.716887067e−6` and `1.645642368e−6`. All charges, stock floors and denominator inputs are for the original `ft2-neon@51f4cfa` D20 payloads. No D18 numeric prices or potential bounds are used.

## Matched occupation and ownership

For each sign use the D20 physical signed charged-row allocation, with all eight families globally nonnegative and the original trial unchanged. The finite support may omit rows due to nonnegativity. First-body `f` pays the original charged singleton row (including its assigned newborn-edge term) with the true total receive `r` in `[0,(k+1)/2]`. The same `r` changes gain intensity from `L` to `L+r`. The two endpoint quotients of the discounted Bellman equation therefore suffice even when the unobserved exterior rate is chosen adaptively. Promotion rate is `L/D`; the update increases `k` at fixed `D`. At `D=65` use the D20 all-stock `k=0,…,4` floors from the authenticated rows exports. On that invariant degree tail the subsolution is independent of degree, and the promotion rate is bounded between `1−(k+1)/65` and one. At the next stock class use zero. Thus the continuation pays all accessible tail time on its admitted arity support without reusing the finite reward.

For the other seven families solve the genuine finite-fringe problem with whole planted size ≤16 and at most four cores. At each state sum the actual occurrences whose apex is the tracked root. A child donor omitted from one motif but present in the fringe contributes its exact reciprocal on the proper receiver; only the apex's exterior parent is unknown. In the joint price, sum intercepts and slopes before taking the minimum over that common exterior rate at 0 and 1/2. Its gain transition uses the same rate. The physical gain/promotion events increase size, so the finite recursion is acyclic. Crossings are zero for this nonnegative lower payment. `g` here is the full separate charged **edge row**, not the newborn term already assigned to `f`; all transfers are priced on the receiving rows from the D20 decomposition. The apex ownership counts each physical higher occurrence only once, with unordered sibling pairs/triples in `hC`, `Yt`, `S3` and distinct roles in `Yb`.

For the physical stable-age occupation use the HI discount and multiply uniform newborn bounds by LO. The newborn's actual parent reciprocal is at most 1/2, independent of the permanent birth mark; all tracked sends use actual degrees. The upper potential denominator is from the [D20 rows exports](../rows/), via the exact zero-sum recharge and the separately certified physical occupation of its deficit, not via a D18 estimate. With only original GMP `C0`, the smaller `D≤48`, size-12 finite replay already gives `[0.565681188107710,0.565766933110960]`.

We use the normalized signed generator identity `π(δ_lo)=(μ−λ_lo)π(Φ_lo)` and `π(δ_hi)=(λ_hi−μ)π(Φ_hi)`. In the stable-age limit apply normalized Dynkin to unit-time intervals before taking the mean limit; no differentiation of a large-time asymptotic is needed. The founder is outside these payments. For an upper endpoint divide the certified positive reward by the upper potential ceiling and subtract; lower is analogous with addition.

## Which D20 rows pay

Rounded down, per physical core:

| Payment | Lower witness | Upper witness |
|---|---:|---:|
| Joint first body, `D≤64` | `2.225776442e−6` | `2.177191995e−6` |
| First body with all-stock continuation | `2.272208662e−6` | `2.221826313e−6` |
| Joint finite edge and higher | `5.7189520e−8` | `2.4708513e−8` |
| Total | `2.329398182e−6` | `2.246534826e−6` |
| D20 stationary potential ceiling (rounded up) | `1.356756788233927` | `1.365141582117995` |

The degree-2–4 subobjective alone pays `1.222275367e−6 / 1.194362348e−6`. The 20–64 subobjective pays `1.93142609e−7 / 1.90921998e−7`; it is not a difference of joint prices and is not added to them. On the finite support the separate lower-witness leaders are `hC` (`2.5714509e−8`) and `hV` (`1.9890954e−8`); upper is led by `hV` (`1.2467599e−8`), then `hC` (`4.683428e−9`) and `g` (`4.205722e−9`). Four-slot prices individually stay under `5e−12` here. At size 12 the joint edge/higher price was only `3.80e−9` on either side, almost all edge. These are controlled support diagnostics, not claims about the full stationary contribution on the omitted crowns.

Against the earlier D18 old-flow bracket `[0.565676636772,0.565771871520]`, the new D20 interval is inward by `5.080115067e−6` and `5.517162368e−6`. This is an endpoint comparison only; the D18 rewards cannot be added to a D20 ledger.

## Replay

From the Astra root, for each `side` in `lower upper`:

```sh
taskset -c 8 nice -n 17 timeout 60s env OPENBLAS_NUM_THREADS=1 \
  ./round2/.venv/bin/python -Bu round2/occupation-neon20/renewal/replay.py \
  b3f20_lo.bin b3f20_hi.bin --side "$side" --degree 64 --cap 16 \
  --rows-data "round2/occupation-neon20/rows/$side.json" --preauthenticated
```

Both sign-specific runs passed and evaluated 3,914 fringe states apiece, checking both relevant Bellman corners and all 2,079 exported singleton rows per sign against direct extraction. The `--preauthenticated` invocation uses the same local payloads already hashed by the rows producer, with source and payload metadata checked and the file checked unchanged during replay. Omit it for fresh authentication. Several earlier aggregate attempts stopped at their time limits; no incomplete output is needed for this result.
