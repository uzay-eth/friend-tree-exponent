> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits and Claude replay of earlier freezes; original Neon-2 Gold1 witnesses and inherited stationary results. See the [release index](../../../../../README.md#claim-25) for exact scope and citations.

# FY continuation with an old donor and small paying receiver

Best completed lower-side floor: **`7.355843186702752e−7/core`**. It does not reach `1e−6`. The smaller `(70,5,24)` degree mask with core-arity caps `(5,5,4)` is rigorously below the target even with optimized box hazards: lower-witness ceiling **`8.611930317702526e−7/core`**. This ceiling excludes no larger mask.

## Reward and occupation

Pay exactly one ordered-pair FY difference for a physical host `h∈G2^exact\M3`, `D_h>18`, `c_h≥2`, paying neighbor `r` with `D_r≤5`, and distinct donor `t` with `D_t>5`. If `z=Δ_h^L g(h,r)` and `F(x)=σ(xz−A_h(x))`, the reward is `max_k F(ξ_k)−F(1/D_t)`. It is disjoint from the existing `175/5/8` small/small payment, from 2f74's high-receiver/small-donor extension, and from the original `D≤18` true edge reward. It may overlap the parent's two-vertex exchange reward; **compare only, never add**. Original higher and singleton terms are not used.

The new replay subclasses `fy/extended/replay.py` for its directed integer Bellman and birth registration. Retain exact degree and arity for all three roles in every transition; only the reward is binned in the **donor degree**, taking the extremal exact discrete FY value within the bin for that fixed host and paying label. Bins start at degree 6; degree 5 cannot collect via a straddling bin. `--exact` chooses singleton bins. Only two successive host degree layers are needed, with arity descending within a layer.

The old neighbor evolves from the time its edge forms. A mother can therefore already be old when her promotion creates the host, and a first child can already be old when the second relevant child is born. In either case the receiver and donor roles are separate; both sibling orientations are counted once. A donor evolves below degree 6 with zero reward until it first enters the paying mask. Prefix and terminal values are not cut by age. Physical adjacent reciprocals are exact; remaining neighbors enter the original safe gain-rate boxes. This is the existing all-age registration proof in `../fy/README.md` with the roles transposed in the reward gate, not a mean-age substitution. Exact integer corner inequalities are checked on each row.

Raw lower fixture at `h=(38,2),r=(2,1)`: old/new edge words `−71038717472,−71180026520`; donor degrees `6,18,90` yield `2816730`, `39895397/3`, `286280377/15` ticks. The upper witness for the same paying neighbor falls from `156300613/20` to `776383/60` ticks over degree 6→90. `--fixtures` checks both witnesses, both paying labels `(2,1),(4,1)`, all 19 knot values, and exact piecewise-linear donor interpolation.

## Executed receipts

From Astra root, `R=round2/occupation-quintet/edges/old-donor/replay.py`; run `taskset -c 11 nice -n 19 /usr/bin/time -f 'user=%U sys=%S' python3 -Bu "$R" ...`.

| Flags | Lower-witness interval | Upper-witness interval | Path rows per pass |
|---|---:|---:|---:|
| `--host 70 --stock 5 --donor 24 --tstock 4` (`--side lower`, then `--ceiling-only`) | `[1.7448809757e−7,8.6119303178e−7]` | ceiling `2.6975834125e−7` | 336,420 |
| `--host 100 --stock 7 --donor 36 --tstock 6 --ceiling` | `[4.4797938598e−7,2.2292258620e−6]` | `[1.3916076252e−7,6.5166054561e−7]` | 1,621,200 |
| `--host 170 --stock 12 --donor 70 --tstock 7 --side lower` | floor `7.3558431867e−7` | not run | 11,819,808 |

The widest completed run printed component floors `2.8871493355e−8` (mother pays), `1.3965207806e−7` (mother donates), and `5.6706074725e−7` (sibling directions); exact total `80878351158765429162857/109951162777600000000000000000`. It used `136.06+0.50` CPU s, checked `204820` band rows, and `1328678/1331330` positive pre-bin rewards. The four `100/5/36` passes used `73.89+0.16` s total. The first `70/5/24` lower floor used `3.96+0.01` s; both ceiling results printed before the tool capture ended. Larger `175/5/90` and `175/5/105` Python runs did not finish; no result is inferred. `--fixtures` exited zero. No production scans; the optional `--hash` mode was not used on these already pinned copies.
