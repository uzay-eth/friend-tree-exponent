# Astra cold attack on Fable Claim 4

**VERDICT: HOLDS.** Scope: board `fable.md` Claim 4(a,b) only, pinned to `ft2-fable@3675c3874e37c6ba40e0b9b396ddb629a5bc67f2`, `round2/newground/q5-low-degree-forest/REPORT.md` §§1–2. Part (a) depends on positive Claim 6 convergence; (b) does not. Parent-degree frequency, low-degree forest statements and simulations outside these two board subparts receive no verdict.

**(a) Age-frequency transfer.** Since a core never ceases to be a core, the conditional empirical age survival at physical time `t` is exactly `C(t−a)/C(t)` for fixed `0≤a≤t` (with consistent endpoint conventions). The Yule mass `e^(−t)N_t→W>0` and Claim 6 give `e^(−μt)C(t)→W^μ X>0` on a single event, giving the ratio `e^(−μa)` and hence a.s. weak convergence, e.g. via rational ages. In size, use the **post**-promotion size `m`: the empirical CDF at `x∈(0,1]` is exactly `C_floor(nx)/C_n` once the index is in range. It tends to `x^μ`, so the Beta law holds quenched and, by bounded convergence, also for the annealed uniform core. No independence between an empirical core and the observation clock is needed, and the assertion concerns promotion age rather than birth age.

**(b) Two separation bounds.** At the promotion stopping time, condition on the complete history. The sole edge from the newborn fringe to the outside is its parent edge. The root-to-slot arrow has no internal effect, and the slot stays in the degree. The only discrepancy comes from a parent-to-newborn send, with intensity `1/D_u`, not `1/D_v`. Use a predictable first-ranked uniform mark at the parent's rate-one clock. A parent degree initially `d` never decreases, so the first incoming time dominates the first point of a rate-`1/d` Poisson thinning. Thus the whole-path disagreement probability up to age `a` is at most `1−exp(−a/d)` for deterministic or promotion-history-measurable `a`. No future-size or eventual-hub conditioning enters this claim.

For the all-age bound I re-derived the comparison. Write `l` for **post**-promotion parent stock. At equality between the physical stock and a separate lower birth–death count `B`, the physical upward rate is at least `B` and its downward rate is at most `B/d`; synchronize the smaller appropriate rates and add residual jumps, preserving the physical generator. Unit jumps preserve order when the counts differ. The lower walk hits three with probability `d^(−(l−3))` by finite gambler's ruin and monotone passage to arbitrarily high barriers. For `B≥4` the killed generator applied to `f(k)=1/(k−1)` is at most `−f/2`: the unkilled expression factors as `f(k)(−1+k/[d(k−2)])`, and killing at three decreases it. Therefore its expected stopped `f` integral is at most `2/(l−1)`. Before the lower walk hits three the incoming rate is at most `1/B≤f(B)`. Compensation in the joint filtration bounds the probability of any incoming arrow before that stop, and adding the hitting probability gives `2/(l−1)+d^(−(l−3))` (optionally capped at one). This uses a union bound at an adapted stop, not independence between separation and the lower walk. Using pre-promotion stock would incorrectly improve the asserted budget.

**Independent evidence.** From Astra, `taskset -c 8 nice -n 10 timeout 8s python3 -I -B round2/attacks/fable4/check.py`:

```text
PASS: 3721 killed-generator and physical-rate rows including k=4
PASS: 96 finite ruin recurrences and infinite-bound comparisons
CONTROL: using pre-promotion stock incorrectly changes 5/6 to 19/36
```

The fixture checks finite arithmetic; the infinite-time coupling and age transfer are analytical. No Fable files edited.
