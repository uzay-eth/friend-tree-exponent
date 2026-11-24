> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../README.md#claim-21) for exact scope and citations.

# High101 all-arm continuation

The size-seven freeze is unchanged. [`frozen.sha256`](frozen.sha256) authenticates the four original files, [`exact_receipt.py`](exact_receipt.py) and the [exact receipt](high101-exact.txt). The receipt records `hV P4 Yt sum`, in that order, on both original D18 witnesses. Here `D≥101`, apex `k≤12`, and arm size are **present-time** predicates. Arity and size exits are irreversible; they are not age cutoffs. `Old` and `young` additionally split first-arm births by whether the actual hub had reached 101 at their promotion. The seam `T=101` in the computation changes the Bellman representation, not the time or support of the reward.

## Nonzero continuation through size seven

[`continuation.py`](continuation.py) sets current reward to zero at sizes ≤7, but carries the actual attached shape and all directed transitions through that boundary. It retains the nonzero continuation when an intrinsic arrow or a hub send produces a size-eight arm. The reward is the exact original charge at sizes 8 and 9, summed over all matching occurrences and both chronological classes. It is disjoint from the frozen size-≤7 charge. The suffix at `D≥101` is solved over the entire future with the shared arity and physical send-rate comparisons; arms born below 101 still traverse the hub crossing at their true continuation state. The present-time arity mask is unchanged in this calculation.

| Original witness | Additional floor (sizes 8–9) | Pre-entrance arm floor within it | Combined with frozen ≤7 floor (downward) |
|---|---:|---:|---:|
| Lower | `2.6364088e−7` | `1.6628448e−7` | `7.8691527e−7` |
| Upper | `3.0083983e−7` | `1.8974706e−7` | `4.6267486e−7` |

The remaining certified gap to `5e−6` is about `4.2131e−6` / `4.5373e−6` per core. The shell is an actual increase on the same high101 hub support, but the gaps are not upper bounds on recoverable charge. The size-eight-only whole-mask control finished on both sides (199 states; totals `[6.6252665e−7,1.89020066e−6]` and `[3.0548998e−7,1.17221731e−6]`). The size-8–9 shell has 485 states. Both lower-price computations printed all three exact results. The upper witness's ceiling completed with `7.9388086e−7` on the shell; the lower witness ceiling did not finish under the wall limits.

## Exact suffix identity with no arm-size kill

Let `τ` be the first hit of 101 in the actual birth-Palm hub process. For one family let `q(x,d,k,u)` be the sum of its actual charged occurrences with first child `x`, and write `V_s(H,x)` for its true conditional discounted future from attached state `x` and full physical history `H`. It includes actual descendant donors and sends from the hub. It can also enforce a present-time arity predicate. Conditional on the full history at `τ`, the suffix decomposes into

```
sum_{x present at τ} V_s(Hτ,x)
+ Eτ ∫₀∞ exp(−s t) bτ+t E[V_s(Hpost,J) | pre-state and its promotion] dt,
b = 1 − (k+1)/D.
```

The second term is physical promotion compensation: each accepted promotion creates exactly one first-arm registration. The post-promotion `k` belongs inside its value. The stationary formula multiplies by `μ exp(−μτ) 1{τ<∞}` and averages over the physical birth law. This formula is for all arms and has no terminal reset. It is not an assertion that a conditional `D≥101` population is in equilibrium for the immigration process. The stationary post-entrance suffix cannot replace the pre-entrance term; the latter is already over half of the shell floor in both witnesses.

For the uncapped inert model, the same formula uses the actual existing arm states in the inert resolvent, and the virgin immigration term becomes `U_s(J)/s`. Only the first few generations are used by a running charge (including outside donors and newborn extensions). Each core has promotion intensity ≤1. For a fixed attached entrance and bounded initial census in these generations, the expected growth of any fixed-depth generation is polynomial in elapsed time: count paths from those existing vertices and integrate one unit intensity for each missing birth in chronological order. For products, use the same argument on the union of the two (or more) paths, matching equal incidences before bounding intensities. An arity factor adds a child incidence and therefore at most one more generation. A send from the hub affects only the addressed arm, and degree growth enters the capped-table charge through uniform gain-difference bounds, not through its unbounded leaf stock. These polynomial bounds give integrable discounted tails for fixed entrance data; a finite-window coupling can be extended to the all-future law after controlling the required census moments. They do not supply a uniform material-sized error at a random physical entrance.

## What the all-arm bound still needs

[`majorants.py`](majorants.py) computes a coarse all-arm ceiling using inherited D18 table maxima, sparse center/transfer extrema, re-rooted arity moments, and a bound on full planted arm size. Per occurrence with maximal tagged arity `j` the charge is bounded by `a+b j`; the linear factor follows from incoming core count. Each gain from ordinary leaves is uniformly bounded by `34` times its table norm (adjacent-control cubic bound at high degree and the inner bound at the seam). Routing terms use the difference of a physical and centered reciprocal times a bounded potential difference or the bounded transferred `g` derivative. Re-root at the largest-arity role to control the arity tail; if all tagged arities are ≤K, re-root at the unique first child to bound the arm-size tail. Its planted size plus one is dominated by unit Yule from three at the common physical age. The degree event is alternatively bounded at the apex by the high101 gain-max chain. These are bounds on the whole omitted support, not a product of independent tail probabilities.

This strategy fails quantitatively: even at `K=60`, size cap `10000`, the `Yt` omitted-arity ceiling alone is `6.0901138918e−4` and its omitted-size ceiling exceeds `581`. The latter is driven by charging every large arm with a worst-case descendant incidence count. Thus simple control maxima plus size/arity tails cannot decide whether the `5e−6` target is physically attainable. The needed replacement is an occupancy-weighted projection that retains the original motif slots and live exterior donor arcs as the rest of an old arm grows. The projection must retain only bounded-depth information and permit leaf-driven degree and arity overflow without discarding the tag; the suffix identity then averages its price at actual entrance attachments and at actual promotions. Neither the all-arm quantitative enclosure nor the target is claimed here.

## Replay evidence

From `.`, using the same pins and NumPy environment as the frozen producer:

```sh
P=./round2/.venv/bin/python
S=round2/occupation-neon/high-hub-limit/all-arms
sha256sum --strict -c "$S/frozen.sha256"
"$P" -Bu "$S/exact_receipt.py" b3f18_lo.bin
"$P" -Bu "$S/exact_receipt.py" b3f18_hi.bin
"$P" -Bu "$S/majorants.py" b3f18_lo.bin b3f18_hi.bin
# Wrap individual shell replays in nice/taskset; wall time varies with contention.
"$P" -Bu "$S/continuation.py" b3f18_lo.bin --cap 9 --side floor
"$P" -Bu "$S/continuation.py" b3f18_hi.bin --cap 9
```

The original both-side shell invocation on the lower witness printed its exact floors but timed out before the ceiling at 350 s and again at 385 s. The upper witness printed both sides and finished in 199.04 s. The separate lower-witness ceiling attempt also timed out at 220 s. Exact high101 runs for the freeze exited zero as did the manifest and majorant checks.
