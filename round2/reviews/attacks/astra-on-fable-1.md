VERDICT: HOLDS

**Scope: final audit of Fable's 15:51 UTC Claim 1, as a finite-window statistical statement and conditional bias calculation.** The genuine final code and data are pinned at `ft2-fable@470ee9f45cbc679bd257cbb1599eee890be3166e`. The independent read-only replay agrees with the final estimates, leave-file jackknives, registered controls, and conditional algebra. It does not establish (C), (C′), (C+) or (C′+), or certify any digits of the limiting exponent beyond the mathematical prior.

**Final evidence.** The primary (19,000 hold-outs) is `0.5657236605826` with run SE `6.15109e−9` and exact leave-file SE `5.94571e−9`; the late cohort (2,840 hold-outs) is `0.5657236688045` with run SE `5.98171e−9` and leave-file SE `6.26079e−9`. Late minus primary is `+8.22191e−9 ± 8.58002e−9`; the inverse-variance descriptive summary is `0.5657236648083 ± 4.28834e−9`. The two registered ANPQ controls are `−0.16` and `+1.18` SE. The original/fresh residual-block permutation checks give `T=.431880, p=.0010` and `T=1.016630, p=.6040`, respectively. Under the respective premises, the separate approximate primary bias bounds are `1.73e−9` and `1.42e−8`; the late bounds are `8.83e−11` and `3.37e−9` using the final archive's outward rounding. The `8.8e−11` board shorthand is not a strict outward rounding.

The bounded verification used **567 raw files and 752,840 distinct stream/run keys**, compared all inspected chunk moments to the committed archive, checked derived RNG pre-seeds and phase splits, and did not run a simulation or the production reducer. [Final audit and replay](../cold/astra/attacks/fable1/final.md) · Manifest (`ft2-astra/round2/attacks/fable1/final.json`) (SHA-256 `8f752cbe407758fc5648648d9e760b15eb03620aa9d12b6be4e2ace3e18ad48a`).

The derivation and earlier interim evidence follow for context; their snapshot-specific numbers are historical.

## Independent derivation

For a physical tree each directed event donor `v` → recipient `x` adds a leaf to `x` at rate `1/deg(v)`. Aggregating by recipient, a core vertex `u` gains a leaf at rate `L_u + Σ_{core v~u}1/D_v`; its leaf recipients together convert to core at rate `L_u/D_u`. Conversion keeps the parent's total degree fixed, changes its `(L,c)` by `(-1,+1)`, and appends a new core with `(1,1)`. The physical total event rate is the total number of vertices, a nonexplosive Yule process. With the stated fixed finite feature layout the observable and absolute generator integrand are polynomially bounded in vertex count, even for signed weights. All finite-time Yule moments exist; localize by size, then apply dominated convergence to Dynkin's formula. Thus for `F(t)=E den(t)`, `H(t)=E num(t)` we have `F′=H`. Defining `β=H−μF` gives `R−μ=Δβ/ΔF`. A constant component of `F` cancels; this by itself gives no shrink factor for the remaining ratio bias.

Write `B=R(16,24)−μ`, `b_j=R(j,j+4)−μ`, and `D_j=D(j,j+4)`. Under the stated positive population denominator conditions, telescoping both increments gives

```
w = D_16/(D_16+D_20),  1−w = D_20/(D_16+D_20),
B = w b_16 + (1−w)b_20.
```

These are denominator weights, not equal weights or inverse-variance weights. The reported sample ratio `D_20/D_16 ≈ 9.61072` implies `w ≈ 0.0942443`; deriving `w` from exponent matching would assume an exact denominator growth law the claim does not need. `[12,16]` and `[16,24]` only meet at the time-16 endpoint, which cancels in a telescoped `[12,24]` increment. If instead combining two windows with positive-length overlap, subtract the overlap once from **both** numerator and denominator increments, e.g. `[12,20]+[16,24]−[16,20]=[12,24]`; averaging the ratios cannot repair double counting. Reusing the time-16 endpoint on the *same paths* induces covariance, even for intervals with disjoint interiors. The nominated `[12,16]` rung uses different paths/stream, so temporal adjacency does not induce sampling covariance there. The two nominated contrasts against the primary and late ratios *both* use that same independent old rung: their covariance is `Var(R̂_12,16)`, correlation 0.801 in the live freeze. They are not independent draws of evidence.

For either factor `r = exp[-4(μ−s)]`, `s∈{0,0.3}`, let `f=w+(1−w)r`. The triangle inequality plus one contraction gives `|B|≤f|b_16|`; with two contractions `|B|≤fr|b_12|`. For any target bias `x` and reference bias `y`, if `|x|≤q|y|`, `0≤q<1`, then `|x|≤q/(1−q)|y−x|` by `|y|≤|y−x|+|x|`. This yields the one- and two-step primary multipliers `f/(1−f)` against `R(16,20)−R(16,24)` and `fr/(1−fr)` against `R(12,16)−R(16,24)`. For late `b_24`, use `q=r` against `b_20` or `q=r³` against `b_12`, yielding `r/(1−r)` or `r³/(1−r³)`. The late-vs-primary agreement only estimates `b_24−B`; it does not supply such an inequality without premises.

Monotonicity: `dr/dμ=−4r`; `∂(fr)/∂r=w+2(1−w)r>0`; `∂(fr)/∂w=r(1−r)>0`. Also `∂f/∂r=1−w>0`, `∂f/∂w=1−r>0`, and `q↦q/(1−q)` is increasing. So use the **lower** old-prior endpoint `μ=0.5650`, not `0.5665`, and the upper limit for `w` (where applicable). At `w=0.094249`:

| shift | `r` at lower μ | K1/K2 primary | K3/K4 late |
|---:|---:|---:|---:|
| 0 | 0.104350484754765 | 0.020093467464298 ≤ 0.020094 | 0.001137567487701 ≤ 0.001138 |
| 0.3 | 0.346455810330057 | 0.164648557544511 ≤ 0.164649 | 0.043390059156962 ≤ 0.043391 |

At the upper μ endpoint these drop respectively to `(0.01990986, 0.00111725)` and `(0.16262809, 0.04258306)`. The population sign assumptions are standing conditions; positive sample denominators are a check, not a proof of them. The primary's approximate 95% bound uses Bonferroni on the one-sided 97.5% `w` limit and two-sided 97.5% contrast limit, with normal quantiles 1.95996 and 2.24140. No `w`/contrast independence needed. The late bound alone needs only the contrast limit; a joint primary/late claim or a joint confidence interval for μ after adding sampling error needs its own coverage calculation. An SE and a conditional systematic bound are different quantities.

## Read-only raw reconstruction

Read endpoint differences directly for each trajectory, without Fable's moment reducer or main simulator. With per-run `(u_i,d_i)`, `R̂=Σu/Σd` and influence `ψ_i=(u_i−R̂ d_i)/Σd`; the linearized delete-run SE is `sqrt[n/(n−1) Σ(ψ_i−mean ψ)²]`. The reported cluster check is the **exact delete-one-file** jackknife: use `(Σu−U_c)/(Σd−D_c)` for each chunk, square deviations from their leave-one mean, and multiply the sum by `(G−1)/G`. For same-path differences, subtract the paired paths' influences or the matching chunk leave-outs *before* taking SE; for disjoint independent streams use quadrature. The leading finite-MC ratio bias plug-in for the primary, `−Cov(d,u−R̂d)/(n mean(d)²)`, is `+2.70e−13`; it is a delta-method diagnostic, not a finite-age bias bound.

| Quantity | Posted raw reconstruction: ratio / run SE / cluster SE | 08:09 snapshot: ratio / run SE / cluster SE | live paths / chunks |
|---|---|---|---:|
| Full 16→24, main | 0.5657236568235 / 7.274e−9 / 6.971e−9 | 0.5657236536283 / 6.817e−9 / 6.452e−9 | 15,400 / 154 |
| Full 12→16, independent rung | 0.5657237030344 / 1.836e−8 / 1.268e−8 | same | 260,000 / 52 |
| Full 16→20, independent rung | 0.5657236423605 / 1.881e−8 / 1.861e−8 | same | 25,500 / 51 |
| Full 24→28, late | 0.5657236702555 / 2.301e−8 / 2.559e−8 | 0.5657236744159 / 1.116e−8 / 9.990e−9 | 860 / 43 |
| `w`, main | 0.0942445232134 / 1.857e−6 / 1.980e−6 | 0.0942443061286 / 1.737e−6 / 1.878e−6 | 15,400 / 154 |

At the posted 13,500-main / 180-late counts, reconstructed using the completed files by a historical mtime cutoff and matching S1–S4 to posted precision, K1/K2 give `1.82e−9` / `1.49e−8`, K3/K4 give `1.13e−10` / `4.29e−9`. At 08:09, `R̂_old−R̂_primary = +4.94061e−8 ±1.95881e−8` and `R̂_old−R̂_late = +2.86185e−8 ±2.14884e−8`; corresponding **conditional** values are `1.88e−9` / `1.54e−8` and `8.74e−11` / `3.34e−9`, outward rounded. Same-path `R̂_old−R̂_primary` on main alone instead gives `+4.75357e−8 ±7.70425e−8` and is not the nominated comparison.

Freeze: `round2/attacks/fable1/current.json`, SHA-256 `b21a7a3ccfc7a11b6a80506402c88fdf472cb5b60ec96aab56612f730b1bc3f6`, pins **301,760 paths / 300 files** and records path, stream, first run, count, checkpoints, size and SHA-256 per completed file, plus all five weight and executable hashes. `posted.json` SHA-256 `423fe32c48a2cd7a2f1d6f11459ee8c2699dac4e001862b6e2408332deebc11d` pins **299,180 paths / 247 files**. Live path eligibility matches the claimed post-freeze chunk thresholds; checks require full layouts, sequential IDs, no duplicate trajectory key within stream, and stable file reads. Raw rows themselves lack embedded seed/binary/weight metadata; provenance depends on the published launch mapping and correct simulation. The executable/weight hashes check the available artifacts, not stochastic independence or every historical launch.

## Stress tests and interpretation

- Independent exact fractions on all 2,956 chronological labelled trees reachable from the three-vertex path at sizes 3–8, plus three cap-exercising structured trees. Enumerated full directed physical events and reconstructed the weighted feature value from graph structure with the stored weights represented exactly as rationals. The rival fast production generator, including its specialized R path, agrees on count drift and frozen full value/drift to maximum absolute discrepancy `2.66e−15`. This fixture checks the generator on supplied trees, not late-time simulation correctness. Under the exact Yule size law, a rational mixture through `n=8` checks the derivative with the explicit outgoing boundary term and three checkpoint increments; this is a truncated fixture, not an exact evaluation of any large-time window.
- Planted donor/recipient inversion, equal-weight or summed ratios, double-counted overlapping windows, and discarded paired covariance all fail. For example `F(t)=e^{0.56572t}+e^{0.56573t}` is a positive increasing mean with `H=F′`, positive windows and eventual rate `0.56573`. It gives old-vs-primary contrast `1.91e−10` with true bias `5e−6` and violates both finite-window contraction premises. It is a logical counterexample to inferring a premise from agreement, not asserted to be the friend-tree mean. Even two individually fast corrections `0.1^j−0.05^j` start at zero and then separate: a spectral envelope does not prove contraction of the sum.
- Snapshot evidence only: L3 is `+2.08e−8 ±1.31e−8` (1.59 SE). Exploratory all-main full−ANPQ at 20→24 is `+4.59e−8 ±1.73e−8` (2.65 SE); disjoint A-test-1 chunks ≥550 give `−2.91e−8 ±4.94e−8` (−0.59 SE), and late L4/L5 are 0.50/−0.05 SE. The current hold-outs do not reproduce the exploratory anomaly at the nominal thresholds, but these are provisional, low-power comparisons, and ANPQ is a control, not a premise of Claim 1. Do not read these as validating (C) for full.

Checked replay (15.79 s wall, one core CPU13, niceness 10, exit 0):

```sh
cd round2/attacks/fable1
taskset -c 13 nice -n 10 env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /usr/bin/time -f 'elapsed=%es' round2/.venv/bin/python -B audit.py current --replay
```

```text
exact physical generator: 2959 trees, max abs C++ discrepancy=2.66e-15; Yule n<=8 rational windows R=3.406169371; rejected donor inversion, equal/sum pooling and double overlap
negative controls: dropped covariance rejected; fake (C)/(C') rejected (contrast=1.91e-10, true bias=5e-06); fast-mode cancellation rejected
P: 0.5657236536283 runSE=6.8166882e-09 clusterSE=6.4522499e-09 Dmean=874275.9746 n=15400 chunks=154
conditional shift=0 primary: k<=0.020094, bias<=1.88e-9
conditional shift=.3 primary: k<=0.164649, bias<=1.54e-8
elapsed=15.79s
```

The separate official-final manifest and replay are linked above.
