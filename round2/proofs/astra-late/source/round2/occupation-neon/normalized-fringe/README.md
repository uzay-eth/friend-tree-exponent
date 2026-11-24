> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../README.md#claim-21) for exact scope and citations.

# Normalized critical forcing: what the quotient removes

**Result.** Positive real-discount normalization gives an exact paired quotient and a way to use a certified sign without separately enclosing two continuations of `Φ`. It does not, by itself, bound the remaining signed continuation. On two reachable homogeneous chain families the *complete*, potential-normalized paired forcing has induced roots separated by more than `6.3203e−5` for **every** nonnegative mixing multiplier. A common scalar subtraction therefore cannot achieve an eight-digit uniform error on the physical bulk. No new exponent endpoint is claimed.

**Neon‑2 credit.** Neon‑2 Track 3, pinned at `ft2-neon2@ac55a7b9fa9236e73a305c32fbf509ae3bf028a3` (`round2/neon2/claim14-proofs/README.md`, Theorem 1 and Corollary 2.1), explicitly constructs `𝓜_s(ψ)=s W_{ψ,s}(J)`. For bounded downward tests it is the normalized inert occupation for real `s>μ` and the physical downward-fringe marginal at `μ`. Its bounded-domain coboundary identity and exact free-tree lifts are credited there; they hold per the two reviews. Our application uses the real identities. Root ownership and any unbounded local multiplicity still need justification when transferring complete physical motif charges; the finite-support tests in [the Poisson pilot](poisson/README.md) fall directly in the bounded domain.

## 1. Exact algebra at criticality and at real trial discounts

Write `a=.565673`, `b=.565775`, `P=Φ−`, `Q=Φ+`, and the two positive signed physical forcing totals as `D=(ℒ−a)P`, `E=(b−ℒ)Q`. For a nonnegative rational `t`, use the same physical measure and set

```
H_t = P+tQ,
N_t = D−tE+aP+t bQ,
F_{s,t} = N_t−sH_t = (ℒ−s)H_t.
```

These are additive observables on the complete tree (or their equivalent allocated scores). This separates the large common potential from the small signed test *before* taking any enclosure. At the common physical critical functional `ν`, applied to per-core additive scores,

```
ν(D)=(μ−a)ν(P),      ν(E)=(b−μ)ν(Q),
μ = ν(N_t)/ν(H_t) = s+ν(F_{s,t})/ν(H_t).                (1)
```

In (1) the denominator is strictly positive. It need not be accurately evaluated to decide a sign of the **one** numerator `ν(F_{s,t})`. To bound a numerical quotient, any certified `h≤ν(H_t)` gives `|μ−s|≤|ν(F_{s,t})|/h`. For a target width `1e−8` from one center, an absolute centered ratio error of `5e−9` suffices; equivalently certify opposite numerator signs at two trial discounts separated by at most `1e−8`. In particular a small difference `ν(D−tE)` by itself is insufficient: the bulk term `(s−a)P−t(b−s)Q` must use the **same** measure. Dropping it would amount to using the unproved potential ratio `ν(Q)/ν(P)`.

One derivation uses normalized physical Laplace averages only for `s>μ`. With `M_s(X)=∫₀∞e^(−su) E[X(T_u)] du` from a fixed physical seed, the positive mean growth bound and discounted Dynkin give

```
M_s(F_{s,t}) = −H_t(seed).
M_s(N_t)/M_s(H_t) = s − H_t(seed)/M_s(H_t).             (2)
```

Use the [physical rank-one mean limit](../../../../../astra/source/round2/fringe/spatial-limit.md) and the uniform exponential moment for local core census; the bounded per-incidence charges with bounded-depth ownership are integrable as in [lineage](../lineage/README.md). Dividing by the common positive `M_s(C)` and sending `s↓μ` gives (1), including the normalized generator identity. Multiplication of each `M_s` by `s` cancels from the ratio in (2). The relation `(s−μ) M_s(H_t)→` a positive amplitude explains why the normalized residual vanishes at criticality. It is not a formula for the physical integral below `μ`, and observing the isolated-seed term `−H_t(seed)` does not determine that amplitude or give the missing digits.

## 2. A real two-sided trial sign without arbitrary continuation of the potential

There is an exact positive-discount alternative using the [original degree-two trigger](../../../../../astra/source/round2/finite/gauge/physical-trigger.md). Denote its resolvent with trigger births removed by `Ĝ_s`, its virgin trigger state by `j`, its trigger rate by `h`, and `f_s=(Ĝ_s h)(j)`. The scores `X` here must be allocated to downward fringes; if they read one outside ancestor, lift to that ancestor with the old owner distinguished. In deriving the following formula on finite trees, include any founder term in the scalar `A_X(s)`; bounded-depth ownership implies the founder term has negligible critical normalized mass. On the real interval with the inherited positive contraction, the untriggered functional `ω_s(X)=(Ĝ_s X)(j)` is well-defined and positive on positive `X`. At `μ`, it induces the normalized `ν` from (1).

One subtlety: **substituting `s` for `μ` in the critical law is not a sign test by itself.** If a finite-seed transform is justified for a chosen allocation, write its standard rank-one form for `s>μ` as

```
M_s(X) = A_X(s) + I(s) ω_s(X)/(1−f_s).
```

Here `I(s)>0` is the common discounted trigger entrance, while `A_X(s)` accounts for seed/top and nontriggered paths. From (2) one gets the *cleared* identity

```
I(s) ω_s(F_{s,t}) = −(1−f_s) [ H_t(seed) + A_{F_{s,t}}(s) ].     (3)
```

For a two-sided sign from the off-critical `ω_s(F_{s,t})`, one must justify this cleared identity at the trial discount **and** certify that the bracket in (3) has known positive sign there. Neither follows from positive comparison for `H_t` alone; `A_{F_{s,t}}` is signed. No extension of the divergent full physical integral is required for the direct positive test `f_s−1`: it already changes sign exactly at `μ`. For the *different*, bounded downward formulation, Neon‑2 Theorem 1(4) gives the exact replacement: `𝓜_s((L⁰+CΔ+b·ev_J−s)h)=s(U_s(J)−1)h(J)`. In particular a zero-virgin correction integrates to zero even off critical, and a generator score for which `h(J)>0` has the known positive prefactor. This identity alone does not give a certified numerical approximation to that score. Another safe use of (1) is to bound `ω_u(F_{s,t})` uniformly over the **prior** range of `u` (including `μ`), by signed residual integration with all error included. Fix `s` in the score while varying `u` in the occupation.

This form avoids an arbitrary continuation of two potentials. For a trial function one can enclose the single paired source via `ω_u(F)=U(j)+ω_u(F−(u−R̂)U)`. The residual must be bounded on all physical states. In particular, descendant boundary terms have to be included *before* taking absolute values or putting bounds on exits. An exit-mass normalization controls constants exactly but does not give `ω_u` of a nonconstant source that varies on omitted states.

## 3. Exact counterexample to uniform cancellation from the normalization

Use two homogeneous physical chains of cores, putting respectively 14,001 and 1,000,000 leaves on each core. These are reachable by finite attachments. Compute complete signed `D,E` and both physical `P,Q` using the existing paired evaluator from `ratio-transfer`, with separate whole-potential differences enabled. There is no new table or motif scan: all-contrast (`round2/occupation-neon/all-contrast/README.md`, historical reference) owns the detailed paired raw rows, whereas this test asks what a complete paired *normalized* score can do. For chain lengths 8, 9 and 10 the increments across adjacent lengths agree exactly. A motif spans at most three core edges, and an outside sender extends its gain stencil by at most one edge; a newborn term deletes to its unique smaller old stencil. The two ends are already separated. These increments are therefore the limiting bulk totals per core, independent of founder terms.

| Leaves per core | Bulk root using only `P` | Bulk root using only `Q` |
|---:|---:|---:|
| 14,001 | `0.5656833075923935` | `0.5656921354552727` |
| 1,000,000 | `0.5657553392551306` | `0.5657641731082980` |

The bulk root for every `t≥0` is a convex combination of the two roots on that row, weighted by the same physical `P,t Q`. The exact intervals therefore stay separated by at least `6.3203799857e−5`. Normalizing by `P+tQ` and subtracting any *one* common constant leaves a uniform radius of at least half the gap (`3.16018e−5`), even if `t` and the centering are optimized jointly. For `t=1` the measured half-gap is `3.6017329097e−5`. This is an obstruction to a uniform bulk error claim, including one relying only on normalization plus a constant centering; it is **not** a lower bound on the error after averaging under the actual critical law or after applying a state-dependent corrector. It does not assume either train has substantial stationary mass. A proof at eight digits must still identify and witness a small stationary signed defect, or give a certified absolute pressure sign by the direct trigger route.

## Replay

From this directory, with the previously verified D18 payloads:

```sh
timeout 15s taskset -c 7 nice -n 10 env OPENBLAS_NUM_THREADS=1 \
  ./round2/.venv/bin/python -Bu check.py \
  b3f18_lo.bin b3f18_hi.bin
```

Passed the six complete-tree paired drift checks, both equal-bulk-increment controls, the exact strict separation for every positive multiplier, and five explicit multiplier cases. The decimal output is display only; all assertions use rational arithmetic. The source commit and any Neon‑2 claims are outside this replay.
