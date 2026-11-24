# Cold check: Claim 6 original `g` and FY household credit

**Verdict: pass on both signs.** I reproduced the degree-18 `g` floors and all three FY roots at `48/4/5` with independent exact solvers. The strongest finished FY values (`175/5/8`, host/tip/arity caps) are bound to the pinned `stock8.log`: all six exact components and their sums check, and 96 selected full-mask path/edge values per sign match an independent corner solver. The full `175/5/8` root optimization was **not** rerun. No floor is booked here.

| Quantity | Lower | Upper |
|---|---:|---:|
| `g`, D18, ticks/core (independently solved) | `1568526522607404053651093/10000000000000000000000` | `1687123111883472663645069/5000000000000000000000` |
| FY `175/5/8`, ticks/core (receipt-bound) | `266189292981106865228871/500000000000000000` | `1873363823474418533418231/1000000000000000000` |
| FY `175/5/8`, normalized, approximately | `4.841955032699968e-7` | `1.703814471943059e-6` |

## Pointwise accounting

- Use the original trial and outward sign (`−` lower, `+` upper). At FY host `h` and paying neighbor `r`, the receiver degree fixes the jump `δ=Δ_h^L g(h,r)` and the primal input `A_h(1/D_r)`. An actual distinct donor `d` contributes `σ(z_d δ−A_h(z_d))`. With receiver and jump fixed, this is affine on each donor knot segment, so the maximum is exactly `M_h(δ)=max_k σ(ξ_k δ−A_h(ξ_k))`. This uses the sign **inside** the max. The ordered atom is `M_h(δ)−σ(z_d δ−A_h(z_d))`. Different donors can spend separate copies of the max; they do not each spend the entire edge residual.
- The signed `g` row includes both physical endpoints, the same donor reciprocals in reward and gain rates, the singleton centering transfer, and the old edge’s newborn `H3` at **post-promotion** host arity. The newborn edge itself belongs to `f`. FY primal and subtracted `A` cancel on a host when summing over reversed ordered pairs. `M3` donors are routed to `H3` rather than to FY. In the old `g` Bellman solve, reward and continuation use the same corner donor; extremizing them separately is unnecessary. Each outside donor enters affinely within a knot segment, so taking all outside donors at an adverse knot gives the endpoint corner even when their physical values differ.
- In a mixed edge, compute `δ` from the *physical* interpolated receiver. Convexity gives `M_h(E δ) ≤ E M_h(δ)`; the primal commutes with interpolation because the first interval is `[0,1/180]`. The physical max therefore consumes no more than the projected max already in the original cap. The audit tested 87 mixed FY cases per sign, including three strict inequalities upper; no Jensen difference is credited. A tail donor still uses its exact reciprocal, whether or not the receiver is a tail. A double-tail edge has no FY host; both cross-receive rates remain in its physical `g` row. The 264 allocated-edge cases per sign include the 37/38, 47/48, 179/180 seams, inner overflow, and mixed and double tails; total raw row plus all ordered FY gaps stayed below the original cap.
- Gate each atom on `max(D_h,D_r)>18` at the **time of payment**. The donor degree and the birth-time degrees are irrelevant to overlap. Summed at an eligible edge, the atoms spend at most its original `cap_g−G`; the old `g` credit spends this residual only at edges with both endpoint degrees at most 18. The same physical path history can contribute to both calculations at different ages without double-spending a live row.

## Registration and lifetime

On a promotion at host `h`, pair its new child with each old core neighbor of `h`. Each new two-edge path is born exactly this way, from its unique surviving old edge. If the old neighbor is the mother, the two paying directions are distinct modes on the **same** history. If it is a child, the old edge itself originated when `h` promoted that older child; the younger child completes the path and the reward sums **both** directions, even if both children have identical current census. The path value may then be symmetric in the two child labels, but the edge-to-path birth registration is not divided by two. The edge roots attach to the promoter of the first birth in either family. Exhaustive enumeration over 873 increasing labelled trees through seven core vertices matched 11,744 terminal ordered atoms against direct adjacency permutations, checked every prefix, and checked 29,722 live masked atoms across degree assignments.

Each stopped state uses known reciprocals from its tracked pair/path and boxes only untracked neighbors. Degrees and core arities never decrease, so exiting a cap cannot lead back into the stopped set. The independent solver explicitly enumerates every low/high corner of these adjustable rates, rather than using the production continuation sort. For the three FY stages, stationary balance telescopes from path reward to old-edge birth input to promoter birth input to `μV(J)`. At each step the extra `(HI−μ)` times a nonnegative potential can be discarded. Using `LO≤μ≤HI` thus gives `LO·V(J)` once per physical core birth; the founder supplies no stable-age birth mass. The same argument in two stages gives the old `g` floor. A root for each mother direction and one joint two-child root are added because they are different ordered atoms, not because their histories are independent.

## Reproduction scope

```sh
taskset -c 5 nice -n 19 timeout 150s /usr/bin/time -f 'wall=%e rss=%M exit=%x' \
  python3 -Bu round2/checks/quintet-fy-sol.py
```

Result: `PASS both signs`, `wall=71.38 rss=101672 exit=0`. D18 independently solves 14,535 unordered edge states and 170 promoter states per sign. `48/4/5` recomputes all stages and exactly matches all six receipt roots. The large spot checks include boundary, near-boundary and stopped host states, all three modes, and child reversal; the owner solver is used only as a spot-check oracle. Strongest receipt binding checks SHA-256, parses exact components, checks their root grid, and independently sums and normalizes them. It is not a substitute for a fresh full solve if the large receipt changes. The timed-out `lower26.log` was not used.

The script pins the three logs and uses read-only sparse payload access through the existing census reader with file identity/header checks and the pre-authenticated gold receipt; it intentionally skips a new full payload digest and verifier/table scan. The detached source HEAD was `b04cc50c1e391e0420eb7b006ac781d484c68543`. No owner source, receipt, payload, or board edit was made.
