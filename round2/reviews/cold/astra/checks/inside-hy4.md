CORRECT

# Independent referee audit: first physical occupation calculation

Scope. This audits `inside.md` (`occupation/inside.md`), `physical-test.md` (`occupation/physical-test.md`),
`five_reserve.cpp` (`occupation/five_reserve.cpp`) and `lifetime.py` (`occupation/lifetime.py`). The raw signed
inequality it consumes, from `proof.md` (`motif-cert/proof.md`) (line 84) and
`higher-rows.md` (`motif-high/higher-rows.md`) (lines 108-113), has its own referee and is taken as given. No author
file was changed; production was not rerun. This is independent of the sister review.

## 1. Conditional input consistency

The conditional statement is σS(𝓛Φ−λΦ)+Q₄+Q₁ ≤ (e_low+e_high)C+E_founder, i.e. after dividing by S,
σ(λΦ−𝓛Φ) ≥ F₄+F₁ − eC − E, exactly the form used in `inside.md` (1). `certificate.json` (`motif-cert/certificate.json`)
gives e_high = 627 (upper) and 9,081 (lower); the logged total caps are e = 155,608 and 115,075
(`occupation-verified.txt` (`outputs/occupation-verified.txt`)), so e_low = 154,981 and 105,994, both positive. The
trials in the certificate (56581165/10⁸, 5656453/10⁷) match `inside.md`. The identity ∑(k_x−1) = (C−1)−C = −1 quoted in
`inside.md:17` is correct. No mismatch found at this interface.

## 2. First/last set accounting and founder

- **First set.** Every promotion creates exactly one new noninitial core, hence one singleton of value v=V(2,1), and
  this includes promotions by the founder (`inside.md:29`). The dump's `(2,1)` row equals the reported newborn reserve.
  My small-tree test verifies the rate identity ∑_events rate·(#new cores) = B = ∑_promotions rate on all 249 states,
  with founder promotions contributing positively on 145 of them.
- **Last set.** At the maximum arity 5 the enlarged-set score is dropped. Dropping a nonnegative term from the right of
  uV ≤ f+… makes the row check stronger, hence conservative, so the true inequality still holds; the dump contains no
  size-6 table. The aggregated identity is verified with T truncated at size 5.
- **Founder.** The founder is excluded from T (no founder table) but its promotions are included in B. A planted error
  that removes the founder's newborn from T(next) while still counting it in B is rejected on 74 (upper) / 79 (lower)
  trees, so the founder contribution is load-bearing and correctly signed.

## 3. Laplace normalization and core-degree moments

With I(s)=∫e^{−st}𝔼C dt and [Z]_s=∫e^{−st}𝔼Z dt/I(s): [C]_s=1, and because d𝔼C/dt=𝔼B,
[B]_s = s − C(0)/I(s) → μ as s↓μ (I(s)→∞ by the mean asymptotic). For a score Z with d𝔼Z/dt=𝔼[𝓛Z],
integration by parts gives [𝓛Z]_s = s[Z]_s − Z(0)/I(s) → μ[Z]_s along a convergent subsequence. The boundary terms
vanish: for s>μ, e^{−st}𝔼|Φ(t)|→0 under the polynomial-growth bound, and Φ(0),T(0) are killed by I(s)→∞.

The moment inputs in `physical-test.md` §1 are the right ones: 𝔼∑_v(1+c_v)^j ≤ K_j e^{μt}, connected r-sets bounded by
A_r∑_v(1+c_v)^{r−1}, so a bounded per-set score gives 𝔼T(t)=O(e^{μt}) and [T]_s bounded; the lower mean comparison
𝔼C(t)≥k e^{μt} supplies I(s)→∞ and the denominator. Omitted founder-centre terms have polynomial expected size and
vanish after normalization; for a leaf-gain drift the required increment bound is supplied by the affine 1/D tail.

Normalizing (1) then gives σ(λ−μ)m ≥ μv−e; normalizing (6) gives m ≤ μy since l−μ≤0; with m>0 and v>e/l,
σ(λ−μ) ≥ (μv−e)/(μy) = (v−e/μ)/y ≥ (v−e/l)/y. The μ in [B]_s and the μ in m≤μy are both necessary; dropping the
cancellation is one of the planted arithmetic faults below.

## 4. Verified v, e, y and both exact endpoints

| | lower | upper |
|---|---:|---:|
| v = v₄+v₁ᴾ+v₁ᶜ | 7,402,100+11,253,928+301,604 = 18,957,632 | 7,709,975+8,148,713+3,226,994 = 19,085,682 |
| e (high / low) | 115,075 (9,081 / 105,994) | 155,608 (627 / 154,981) |
| y from frozen profile | 2,654,889,165,321 | 2,655,887,184,018 |

Exact `Fraction` recomputation of λ ± (v−e/l)/y reproduces both inward endpoints of `inside.md` and both outward
ticks (lower floor 5656523639, upper ceiling 5658045676). I also reimplemented the reverse recursion of eq (5)
independently in `inside_hy4.py` and ran it on the frozen profile/comparison files
(`mc-mbar-*.i64` (`outputs/mc-mbar-lower.i64`), `mc-replay-clean/*-comparison.json` (`outputs/mc-replay-clean/lower-comparison.json`));
it returns y = 2,654,889,165,321 and 2,655,887,184,018, matching the logged values. Scope: this recomputes y from
already frozen rows; it does not rebuild the profile or witness.

## 5. Independent code, small trees and mutations

`inside_hy4.py` uses my own small-grid five_reserve dump (D=10,7,6,5,4 / C=7,5,4,4,4, different from
the shipped cutoffs), my own interpolation and my own derivation of F ≥ uT−𝓖+vB with u=1133/2000. Results, from the
repository root on one CPU:

```
python3 -B round2/checks/inside_hy4.py
[A] v/e/y, both exact endpoints and both outward ticks: PASS (lower and upper)
[A] independent lifetime recursion reproduces y exactly (lower and upper)
[B] 249 real small trees: F >= uT-G+vB holds, min slack 0 (upper and lower)
[B] founder-active trees 145; drop-founder-newborn rejected on 74/79, inflated-newborn on 221/225 trees
INDEPENDENT INSIDE AUDIT: PASS
```

Planted errors rejected here: dropping one v component, doubling e, and omitting the μ-cancellation break the
endpoint arithmetic; dropping the founder newborn and inflating the newborn value break the small-tree identity on
many trees. The row rate-sign mutation is covered in the companion tail review. Commands used:

```
taskset -c 0 python3 -B round2/checks/inside_hy4.py
```

## 6. Limitations

The raw signed inequality, the witness rows and the profile generation are separately refereed and were not redone.
Production tables, the full replay and the deeper `small_fringe.py` illustration (whose weaker endpoints are not used
in the final numbers) were not rerun. Within this scope no gap or error was found.
