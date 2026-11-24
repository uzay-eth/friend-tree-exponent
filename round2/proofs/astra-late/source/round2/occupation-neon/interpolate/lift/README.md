> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../README.md#claim-21) for exact scope and citations.

# Pointwise lift of higher-motif occupation

**Internal replay after the floor correction.** A small top-owned `P4/Yt` lift gives a pointwise upper test at `0.565773198615`, an `8e−8` improvement over the recomputed first-body-only pointwise upper `0.565773278615`. It does not improve the existing occupation upper endpoint. The lower lift is positive but too small to use at the existing reporting precision. Neither side reaches a `1e−5` shift. The earlier upper `0.565770220167`, its mixed descendant, and the initial lift composed on that upper are **withdrawn**. This branch does not publish a new combined bracket.

## Triangular lift

Write `ε=+1` for lower, `−1` for upper. Let `c_M` be a nonnegative **unspent** signed charge on an existing physical embedding `M`, in original witness normalization. Keep every rate physical. Let `A_M` evolve the *same* embedding and the auxiliary state its correction reads, with no new copies of the correction credited as part of `A_M`. Off-embedding gains and promotions still belong to `A_M` when they change that auxiliary state. Let `z_M` be its proposed correction, and `B_{M,N}z_N` the sum over genuinely new embeddings `N` born at promotions, whose old tagged set is exactly `M`, of their actual birth rates times their post-event values. A newly promoted vertex is a terminal tip, so its deletion leaves a unique old connected set. Unordered siblings must have the same multiplicity convention in `B` and in the additive potential. If only the old motif shape changes and the tags are unchanged, that is an `A_M` event, not a new embedding. Domain entry within an extant embedding is also an `A_M` event.

With `W=Σ_M z_M` there is an exact finite-tree identity:

```
(ℒ−s)W = Σ_existing M [ (A_M−s) z_M + Σ_N B_{M,N} z_N ]
           + Σ_promotions b · z_singleton(J_post).
```

Sum only selected owners, consistently in both terms. The birth terms are triangular in the number of tags, including when a correction reads additional untagged donors. Starting at the largest family, prove `(s−A_M)z_M ≤ c_M+Σ_N B_{M,N}z_N`. Move each `B` value into the next smaller problem **once**, using the shared physical rates in both the current charge and the continuation. A singleton bound `z_1(J_post)≥a` uniform over every physical birth environment can be normalized by subtracting `a` at each nonfounder: the promotion contribution is then nonnegative and the old-vertex row pays `s a(C−1)`. If this birth value is identical in all environments the promotion contribution of the normalized potential is zero. On the founder use a separately paid correction or the cap-only constant from the first-body proof. The actual potential is `Φ+ε W_normalized`; subtracting the sign here instead would consume the birth payment. Check positive potential comparison and the shift at the same `s`. A birth source priced at an occupation discount cannot be installed at a different `s` by assertion.

The equivalent top-owned form is to solve a single inequality `(s−A_apex)z≤Σ_{M apex-owned} c_M` in a live planted crown. Apex ownership ensures each new apex is a virgin, and a promotion within the old crown changes its `z` by its true post-event value. This is direct back-substitution through shaping events; it also allows pre-birth states with all higher charges zero. An independently verified first-body margin can be added because the crown proof spends only the chosen higher rows. To combine with an edge lift, the edge ledger must spend only `g` (including its GMP-owned original birth terms), and the singleton and target-shift rows must be paid only once; do not add whole endpoint shifts obtained at different trial discounts.

### Scope of a general conversion

The lift is systematic for **uniform physical birth-start Bellman certificates** of this type. An averaged stable-occupation lower bound alone is not such a certificate. E.g. on two persistent environments with local test `c=0,1`, the stationary mixture can have positive mean, but in the closed zero-reward environment the local inequality forces `s z≤0`. No positive uniform birth payment follows. More generally any admissible closed class of the physical-rate relaxation with insufficient occupation is a dual falsifier of the proposed uniform inequality: integrate it on that class. This objection applies to conversion by local subsolutions, not to existence of better physical witnesses of genuinely larger support. In the actual model the `P4/Yt` face occupations depend sharply on specific incoming child rates; a bound based on the wrong receiver, or on the two independent worst rates for reward and continuation, is not the shared-rate inequality. There is no universal numerical `1e−5` no-go theorem here.

## All-degree donor continuation and exact pilot

Use `P4` with all four protected vertices bare of degree two. The fourth vertex has a child outside the motif, whose reciprocal is `u`. Use `Yt` with apex bare degree two, middle bare degree three, one tip bare degree two and the other a virgin. The bare tip carries the outside donor. The two tips are selected physically by their distinct census, and ownership is once at the apex. In either face the parent of the apex sends at `r∈[0,1/2]`. The charge is exactly `a+p r+w u`; the replay obtains the coefficients from the pinned affine physical charging formula. The protected-slot exit rates are `3+r+u` for `P4` and `4+r+u` for `Yt`. The latter includes three bidirectional inside edges at total `5/2`, the virgin's leaf gain at one, and its promotion at `1/2`. Neither face stops at a gain or promotion of the donor, and the donor's complete descendants can grow without restriction. Its degree growth intensity is `g∈[0,D]`; any transition that leaves `D` unchanged leaves the test `z=A+B/D` unchanged.

With `H=s+fixed+r` and `u=1/D`, the verified subsolution inequalities on a live face are

```
(H+u)(A+B u) + e B u/(1+u) ≤ a+p r+w u,   e=g/D.
```

Check both `e` and `r` endpoints jointly and multiply by `1+u`; check all degree-three Bernstein controls on `u∈[0,1/2]`. The same exact actual send `u` funds the charge and the protected vertex gain. Positivity of `z` and the two-coordinate rational optimisation for the largest proved `z(1/2)` are checked exactly. Solving backwards on complete planted crowns below size seven, with exact physical arrows and the shared exterior `r`, passes these face values back to the apex newborn. At other size-seven exits use zero. This extension by zero is safe outside the active domain since all successor values and selected rewards are nonnegative. Once the protected face breaks it cannot return: some protected vertex's degree has grown or the last virgin became nonvirgin. Even if a successor enters a *different* positive face, killing that transition only lowers the Bellman right-hand side. The finite prefix really is a crown check, not an independence approximation at birth. A full whole-tree generator calculation independently checks the composed lift on physical fixtures.

At a shift of `1e−5` from the original witnesses, the extra per-nonfounder drift floor is about `3.610e−14` lower and `1.6482e−7` upper. After the corrected first-body payment the combined margins are `−1.06269e−5` and `−1.11016e−5`, respectively. The higher rows can be sizeable after entry (`P4` upper face at a virgin donor: `z≈.00108366`), but the selected entrance is rare. A matching backwards *upper* comparison takes the larger shared-rate endpoint throughout the prefix and at a face pays `max(c)/(s+fixed)` until its first protected break. Even that credits less than `4.682e−7` lower and `5.732e−7` upper per nonfounder. Thus solving the very same selected-entry problem more accurately cannot bridge the first-body margin deficit at `1e−5`. This ceiling does not include later entries at larger crowns or charge classes outside these two stopped faces. To continue, admit other protected rows and retain the donor that actually drives each positive coefficient. A degree seam or unbounded arity in *another* protected slot requires its own pointwise terminal proof; the live-donor tail here does not certify that extension.

## Withdrawal, authenticated inputs and replay

The historic upper used the `k=1` floor `861071198001/10^15` instead of `861071198001/10^18`. The first-body continuation at `s=.565770220167` built with that sole mutation has **negative** upper-signed whole-tree drift `−.00346313849618…` on a chain of seven cores with 14,001 leaves on each core. The negative control computes drift by exact complete potential differences, including all seven original higher charges, and asserts the negative sign before using the repaired rows. The corrector calls in **this directory** take their ten all-stock floors from the hash-verified integer slabs by exact independent recomputation; every floor is cross-checked against the separately hash-verified `bands.log`. No decimal floor or copied `renewal/check.py` literal is trusted. The pointwise endpoints above are recomputed, not offsets applied to any withdrawn upper result. Other branches' repairs, mixed witnesses, and compositions are outside this replay.

From this directory:

```sh
timeout 35s taskset -c 7 nice -n 10 env OPENBLAS_NUM_THREADS=1 \
  ./round2/.venv/bin/python -Bu check.py \
  b3f18_lo.bin b3f18_hi.bin
```

Passed the ten floor comparisons, 45 whole-tree lift checks per side and the seven-core negative control. At the new internal upper test the first-body-only margin is about `−1.01102e−7`; with the `P4/Yt` lift the exact margin is positive (`≈6.37156e−8`). The unchanged lower pointwise test is `0.565675148075`. The replay consumes the already verified D18 payloads and does not rehash the full files.
