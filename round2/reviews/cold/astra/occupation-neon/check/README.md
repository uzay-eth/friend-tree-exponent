# Independent physical pricing of Neon slack

**Result using the previously verified full eight-motif Neon certificates:**

```
Neon lower   .565673       → .565674338794
Neon upper   .565775       → .565773985336
```

These are outward rational endpoints, with strictly positive improvement on both certificates and no cancellation allowance. The improved Neon lower bound also exceeds Astra's `.565673570322` by `0.000000768472`; both endpoints improve that intersection.

## Signed price and stationary averaging

Use the original `B3-FullH-18` payloads from [Neon's pinned proof](neon-proof-pin.md), with the recorded GMP production verifications in the tracked receipt archive (`attacks/neon1/outputs/neon1-production-rehearsal.tar.gz`), members `outputs/production-final/{hi,lo}.log`. Write `σ=+1` for upper, `−1` for lower, and let `Fσ(L,k,ρ)` be the signed one-body slack, in normalized physical potential units. The published transfer decomposition exactly sums the signed slacks of all eight motif classes to `σ(λ Φ−ℒΦ)`, and the already verified row inequalities make every charge nonnegative. In particular we may throw out **all** other motifs and all unused one-body slack. A positive uncharged slack on a small standalone tree would not justify this step.

`crowns.py` verifies the complete payload hashes and independently reads the sparse required integer entries via mmap. It uses the centered one-body formula: the original gain `L+ρ` is replaced for inner vertices by `L+barρ`; for tail vertices the net gain part is `(L+ρ) Δf+(barρ−ρ)ψ` with `ψ=c2/D²+c3/D³`. The physical promotion at rate `L/D` pays the new edge with the exact Bernstein-3 weights. The upper comparison constants from the production logs are `1.95820675182` for the upper witness and `1.84124847315` for the lower witness.

The newborn measure is `β=π B/μ`, not an arbitrary fixed exterior. Conditional on its complete birth environment, follow the same physical old root. The stable-age identity gives `πF=μ E_β ∫ exp(−μt) F(root_t) dt` for any nonnegative root test. To transfer the charge to the finite-tree drift, take normalized Laplace occupation at `s>μ`. For every noninitial birth, its discounted expected price is at least the uniform baseline `V_2(s)` of the recurrence with discount `s`. A birth in class `d` gets a nonnegative bonus `V_d(s)−V_2(s)`. The discounted number of births in that class is at least `a_d(s)` times the discounted number of all noninitial births, by the single selected-birth construction below, applied once to each such newborn. Summing uses Tonelli and disjoint birth classes. The finite recurrences **without** rounding and the `a_d(s)` are continuous as `s↓μ`; the rounded replay at `HIGH≥μ` bounds these values at `μ` from below. Equivalently one can keep the replay's computed baseline and per-bin lower values and use the maximum of the baseline and bin value on each class; this avoids subtracting an upper estimate for the true baseline. Normalize by the integrated mean core count as in `../../occupation/physical-test.md`. The total discounted birth count ratio tends to `μ`. The normalized signed integrated drift is `σ(λ−s)[Φ]_s` up to a vanishing boundary by finite-time integration by parts and `Φ≤M C`. Thus `σ(λ−μ)≥ LOW*v/M`. This argument requires no integrability assumption for a particular per-core generator charging under `π`.

## Exact crowns and actual exterior bins

The complete planted crown has a silent slot only for drawing its intrinsic adjacency graph; to evolve it add the actual exterior-to-root rate `p=1/D_parent`. For each labelled attachment target the script sums reciprocal degrees of its sending neighbours, attaches a leaf, and re-encodes the resulting crown; it checks every successor size and total intrinsic rate `n−1/D`. At crown size at most 12 a reverse-size Bellman solve minimizes over `p=0,1/2`. The fractional-linear first-event expression is valid at all intermediate rates because its denominator is positive; the resulting generator inequality holds for arbitrarily time-varying predictable exterior input. It pays `Fσ` only at the focal core, not again at any of its children. With size 12 this gives lifetime `2.6621514479e−6` upper and `3.304642903e−6` lower. Size two is exactly the virgin holding calculation with denominator `HIGH+2` up to flooring.

For explicit conditional occupation, choose any exact birth degree `d≥2`. A typical newborn can make consecutive plain gains from its original leaves up to degree `d`, then its first promotion, using discount `HIGH` and the upper outside hazard `1/2`. The lower discounted probability of the selected birth is `a_d=∏(wanted rate)/(HIGH+l+l/(l+1)+1/2)` over `l=1,…,d−1`; wanted rate is `l` until `l=d−1`, then `(d−1)/d`. Hence `β(D_parent at birth=d)≥a_d`. This follows by counting the selected physical births (at most one per old vertex); under the spatial semigroup normalize the stationary input count by the birth intensity. In finite time it is ordinary birth compensation followed by the normalized Laplace limit. The exact-degree birth classes for distinct `d` are disjoint. In particular overlapping **current** degree bins in the printed examples below are separate lower bounds, not jointly additive events.

Given parent birth degree in `[d,M]`, construct the crown with one local leaf and one virgin child. For time `h=1/(2(M+6))`, require one gain in the first half, one promotion in the second, and no other crown events or parent degree gains. Until an unwanted event, the parent degree stays in its true bin; the exterior input is at least `1/M`. The wanted rates are at least `1+1/M` and `2/3`, and the total hazard of monitored events is at most `M+5`: a parent's degree gain rate is at most its current degree, and the complete crown rate is at most `4+1/2` on the scripted path. Integrate the two ordered times and use `exp(−(M+5)h)≥1−(M+5)h`. Hold for another `h` with the same hazard bound and pay the stable-age discount. The exact replay returns conditional masses before multiplying by `LOW*a_d`. For `[5,6]`, the unconditional `π` floor is `3.86851867241e−8`; slack there is `2.05480952014e−6` upper / `2.55008675601e−6` lower. The exterior bound clamps an actual rate during the monitored event, without conditioning a proposal clock to be independent of success. These small explicit floor products are examples and are **not** added to the stronger root solve.

## Stronger endpoint from arbitrary progress

The best bound here follows all root count progress up to `D=64`, including arbitrarily large descendant subtrees while in the root cut. The rewards are still only the one-body row. Let `V_d(L,k)` use parent degree at birth at least `d`, so its inverse stays at most `1/d`; the child contributions total at most `k/2`. The actual physical `ρ` lies in `[0,1/d+k/2]`. Minimize

```
( Fσ(L,k,ρ) + (L+ρ) V_d(L+1,k) + (L/D) V_d(L−1,k+1) )
 / ( HIGH + L+ρ + L/D )
```

at the two endpoints, in decreasing `D` and decreasing `k`, taking zero beyond 64 and flooring all values to `10^-16`. Root gain increases `D`; promotion at fixed `D` increases `k`. Descendant-only changes affect the adaptive rate but not the value function, so no dependence on that exterior evolution is assumed. At exact parent birth degrees `3,…,8`, replace the baseline value by `V_d(1,0)` and add `a_d*(V_d−V_2)` for that class only. Both baseline and bonus already account for the tail `ρ` terms in the same row. Sample upper tail rows benefit from positive lower `ρ`, but a positive reciprocal bound cannot be extended to all ages: it is used only inside the short monitored bin, not throughout this unbounded-time solve.

| Witness | Baseline lifetime | Disjoint birth-bin bonus | Certified improvement (display rounded) |
|---|---:|---:|---:|
| Upper | 3.50937224910e−6 | 3.11695848867e−9 | 1.01466421200e−6 |
| Lower | 4.35189425120e−6 | 5.83660406493e−9 | 1.33879441441e−6 |

In this method the bulk gain comes from arbitrary physical progress, not from the very small short-bin occupation or its positive lower rate. The latter is independently quantified above. The calculation uses one noncancelling row class of the eight; additional verified rows could be priced with their own physical support and multiplicities.

## Replay

Supply the original two binary payloads, one per sign. The checker is standard-library only and verifies both full SHA-256 hashes before doing endpoint arithmetic; a header match cannot substitute for this. Use unchanged files for the duration of the replay. Independent hash verification is also available from the directory containing both files:

```sh
sha256sum --check --strict <<'EOF'
1c4dec226208b02e2126dac5dcd911cfd714635bdaf322eac01e5ef0735fda41  b3f18_hi.bin
9921eee052c489f07e8eba63d65c1dac5e85e4e495d8427e3d67a1a0642025ec  b3f18_lo.bin
EOF
```

From the repository root, with the original payloads at the paths shown:

```
timeout 60s taskset -c 13 nice -n 10 python3 -Bu round2/occupation-neon/check/crowns.py b3f18_hi.bin b3f18_lo.bin
```

The exact final checks require outward upper `≤565773985336/10^12` and outward lower `≥565674338794/10^12`. The default also replays crown sizes 2, 4, 7, 12, four conditional exterior bins, and the signed tail `ρ` samples.
