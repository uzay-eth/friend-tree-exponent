# Quintet stationary host cold audit

**Holds for the pinned host-only recipe.** The independent limit-24 replay reproduces both exact endpoint fractions in the producer log, giving the outward interval **`[0.565702368490575, 0.565748281318382]`**. Both signs use zero extra singleton credit. The original signed singleton caps are retained (`0` lower, `−1053` upper); the lifetime-dual exports are outside this verdict.

| Witness | Stationary denominator / `2^40`, outward | Residual higher cap, ceiling ticks | Joint margin, ticks (downward) |
|---|---|---:|---:|
| Lower, mode 2 | `[0.974761052809, 1.237126745262]` | 1277 | 76.170884607 |
| Upper, mode 1 | `[0.992265019071, 1.237905153303]` | 1601 | 68.539526844 |

## Host accounting and normalization

Use the accepted production receipt for the pointwise signed row caps and global comparability. The audit rebuilds the host extrema from the physical binary layout with a separate reduction; the producer functions are called for parity only. Endpoint interpolation is a convex combination of physical knot corners, with the child's parent edge allowed to contribute zero at the root. The potential has exactly one `f` at each core and one `g` per nonfounder. The host for H3 is its middle, for S3/S4 its hub, for P5 its center, and for B5 its hub. P4 contributes one half at each internal vertex. For a fixed internal P4 host, the other internal vertex has `c_y−1` extensions; the `c(c−1)` ordered neighbour choices divided by two give `C(c,2)`. P5 has `C(c,2)` branch pairs and `(c_b−1)(c_d−1)` extensions; B5 has `3*C(c,3)` choices of an extended arm and two short arms, then `c_y−1` extensions. The five-spoke labelled fixture saturates all six resulting hub counts, including 20 incident P4 embeddings of half weight each, 30 B5 and 5 S4. Caps are allocated by the same counting; replacing negative family caps by zero and capping the sum by the original uniform higher cap are both outward.

On the lower witness only, use `W(u)+alpha*(c_u−2)` for the floor, with `alpha=2^38` ticks. Its tree sum is the original sum minus `2*alpha`, restored by a constant. The constant vanishes in normalized occupation; the upper comparison always uses the unrecharged `W(u)`. Therefore a lower and upper number in *one* cell need not be ordered: they bound different allocations. Each bounds the same original total after summation and normalization. The root's individual host comparisons are uniform (`O(1)` in core count), and the single-core state uses its original singleton bound. Under the physical Laplace occupation limit, omission of the founder and the constant recharge do not change the stationary denominator.

For nonfounders `c=k+1`, with the mother included. Degrees `2..24` retain every exact arity; above 24 each degree goes to `(25,min(k,7))`. The replay merges all finite degrees through 1999 before merging the analytic envelope starting at 2000. Analytic hosts carry no higher family, though an analytic vertex can occur as an endpoint on a motif owned elsewhere. For fixed `c`, `q` lies in `[max(0,1−c/DF),1]`; the saturated cell uses `[0,1]` and all reachable G2 endpoint classes. In mode 2 the lower concave-quadratic envelope subtracts `f2/DF` and the ceiling drops the nonnegative penalty; in mode 1 monotonicity gives both endpoints. The edge interpolation attains its extrema over the enlarged domain at `DF` or infinity, including zero. On the saturated lower envelope `alpha>=0` permits using its minimum arity. Intersecting with the original host floors/ceilings is valid by the production host comparison, not merely by the summed potential inequality.

## Sided LP and signed composition

Use the pinned physical birth-law builder (`occupation-neon/pressure-memory/birth-law/replay.py`) with the original exact quintet `LO,HI`, no Palm kernel. On its true lift, the occupation and edge classes have mass `1/mu`; the excess stubs and overflow upward-receive classes have mass at most `1/mu` by the arity and total-receive identities. The physical reciprocal ports of a shared edge are lifted on the same edge column. The tail preserves incoming degree flow, lower arity promotions, and the saturated arity and upward-receive identities. All audit objective scores are nonnegative and cell-constant, with zero on auxiliary columns. The independent certificate takes candidate multipliers from HiGHS, rounds to a `10^−12` grid, scans all 5,542 rational reduced costs and subtracts the maximum nonnegative violation in each class. On the per-core scale the total repair is at most the sum of these defects; the birth multiplier is multiplied by `LO` when nonnegative (if negative the nonnegative objective already gives zero). No floating primal or floating objective is accepted as evidence.

Write `x=sign*(mu−lambda)` with `sign=+1` lower and `−1` upper, and let `M` be the mean original potential in ticks. The signed physical inequality gives `x*M >= credit−cap_f−cap_g−pi(H)`, where `H` is the nonnegative higher host cap and `credit=0` here. The minus sign on the upper witness's `cap_f` is retained; this is the original singleton account. Given a net lower bound, divide by `M_lo` when negative and `M_hi` otherwise. For a trial `t`, bound `(x−t)*M` directly. For `t<0`, the joint nonnegative cell objective is `H_global−H+(-t)*(W_lo−global_lo)`; for `t>=0` it is `H_global−H+t*(global_hi−W_hi)`. Anchor the constant term at `−t*global_lo` or `−t*global_hi`, respectively. A positive joint margin can then increase the trial by `margin/M_hi` on either witness. Both production trials are negative and both joint margins positive. The positive-trial branch was also exercised with a formal intercept, with no physical reward assigned to that control.

If composed with a later verified singleton reward, `F` must be the slack of the same signed allocation (including `cap_f` once); its lifetime floor and the nonnegative masked remainder would jointly lower-bound `pi(F)`, not constitute two independent uses of `F`. This replay does not certify that remainder or its price. It did not rerun the production drift verifier or alter the board.

## Replay evidence

Both snapshots and the receipt were authenticated under `round2/.work/neon2-claim6-7M6Ys3rJ`; the detached production source was clean at `b04cc50c1e391e0420eb7b006ac781d484c68543`. The audit script (`quintet-stationary-sol.py`) pins the host-only producer source, both snapshots, the receipt, the shared builder and `astra-quintet-host24.log`. It compares exact denominator bounds, higher caps, trial offsets, margins and endpoints to the log as fractions. Run from `[local path]`:

```sh
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 nice -n 15 round2/.venv/bin/python -Bu round2/checks/quintet-stationary-sol.py
```

Output (exit 0):

```text
LP limit 24 saturation 7 rows 601 columns 5542
lower finite marks 1998999 cells 307 host comparisons 388 higher counts {'h3': 86843988, 'p4': 29674890, 's3': 17410176, 'p5': 15108750, 'b5': 26374140, 's4': 3795012}
lower denominator 0.974761052809 1.237126745262 cap ceiling 1277 margin 76.170884607 endpoint 0.565702368490575
upper finite marks 1998999 cells 307 host comparisons 388 higher counts {'h3': 86843988, 'p4': 29674890, 's3': 17410176, 'p5': 15108750, 'b5': 26374140, 's4': 3795012}
upper denominator 0.992265019071 1.237905153303 cap ceiling 1601 margin 68.539526844 endpoint 0.565748281318382
```
