# Cold audit: two-sided Neon D18 EDGE correction

**CORRECT in scope.** The outward interval is **`[0.565677980686114, 0.565770513279953]`**, conditional on the inherited full D18 GMP B verifications, physical mean/stable-age identities and incoming interval `[0.565673570322, 0.565775]`. The comparison is to the repaired all-seven release at `80fe02660b9f0a0d4e7fdc8240b0815615ff5896`, `[0.565676930574, 0.565771634792]`. Both trial margins are negative; it is the **net stationary numerator** that is positive. This audit does not enter a claim on the board.

## Physical generator and signed rows

- A counted edge joins a nonfounder parent to its actual core child. The incoming arrows from the opposite endpoint contribute `1/D_w` and `1/D_p`; the remaining sums `u,v` are in `[0,k_p/2]` and `[0,k_w/2]`. The parent has `k_p−1` other children and one parent, all of degree at least two. Each promotion is at `L_i/D_i`, retains the existing edge and changes its endpoint from `(D_i,k_i)` to `(D_i,k_i+1)`. Gains change `(D_i,k_i)` to `(D_i+1,k_i)`. The cuts are absorbing because neither coordinate decreases.
- The `g` charge at the original trial includes both incoming `f→g` donor transfers, outgoing recharge and promotion/newborn motif terms; it is affine in the **same** live `u,v` as the projected generator. Cross-multiplying the Bellman inequality leaves an affine residual on either side of the hinge, so rectangle corners and hinge intersections suffice. On actual rows the residual charge after reserve subtraction is nonnegative: the booked parts are lower bounds, and the other rows inherit nonnegativity from B. Taking the hinge on the relaxed box therefore does not invent physical reward. Downward rounding gives the certificate; the upward recurrence only challenges the fixed robust box/cut class.
- A nonfounder promotion creates a separate edge with value `V((D,k+1),J)` at `J=(2,0)`, booked once in the parent's census equation. Existing incident edges update separately. The newborn vertex itself pays zero because its centered value is `z(J)−z(J)`. The founder has no census value and creates no counted edge. The global signed constant `−δM/s` pays the founder's part of the pointwise cap. The cap keeps the repaired degree-64 census and clamped arity tail, and does not serve as the stationary denominator.

## Reserve and division

The unit is the current signed **row at its physical embedding**, not an underlying arrow. The size-20 replay's `g` price is the uniform inner floor on actual both-inner edges; the stationary high-parent joint replay's `g` price is the mixed virgin-edge floor at degrees 18–191, arities 1–10. These supports are disjoint. The EDGE rules reserve at least these amounts on the whole intersection with the correction's cut, including the `17/18` seam. Outside that cut, and on founder edges, the correction spends no `g`. The six higher-body families are untouched, including the high-parent `hC` pair. The repaired size-12, size-20 joint and high-parent joint occupation receipts thus still apply, once each. No other stationary `f` or competing edge credit is added.

With `δ=ε(s−λ)`, the edge and census comparisons sum to `ε(ℒ−s)Φnew ≥ m(C−1) + Qreserved`, where `m=s z(J)−δM`. The bounded-per-core correction preserves strict positivity of the original potential. Normalized Laplace integration gives `ε(μ−s)πΦnew ≥ m + πQreserved`: `[C−1]→1` even when `m<0`; the single constant in the potential vanishes after normalization. The upper sign can drop `−Σ V` in a *denominator ceiling* because `V≥0`. On the lower sign there is at most one counted edge per core, so add `sup V` to the maximum signed census correction. In both cases use the old all-degree stationary ceiling plus that one-sided correction bound. **Add the negative margin to the reserve first.** Only the resulting positive numerator can be divided by a denominator upper bound in the claimed direction.

| Quantity | Lower | Upper |
|---|---:|---:|
| Trial | `.565677930574` | `.565770634792` |
| Local margin `m` | `−2.3505248367831583e−6` | `−2.0812586098702505e−6` |
| Conservative net numerator, approximate | `6.80268861557417e−8` | `1.66128843901049e−7` |
| Exact denominator ceiling | `1357493815758092773/1000000000000000000` | `273436004214812983/200000000000000000` |
| Exact endpoint before outward rounding | `191976090122981726641370633087/339373453939523193250000000000` | `24172191945911554486185483421/42724375658564528593750000000` |

## Evidence

From `[local path]` with `ulimit -v 10485760` and `OPENBLAS_NUM_THREADS=1`, using the same NumPy-enabled Python with assertions on:

```sh
P=round2/.venv/bin/python
R=round2/occupation-neon/interpolate/edge
A=round2/checks/neon-edge-sol.py
sha256sum -c round2/occupation-neon/interpolate/compose/inputs.txt
nice -n 15 taskset -c 10 timeout 150s "$P" -Bu "$R/replay.py" b3f18_lo.bin --degree 40 --arity 8 --challenge --mask disjoint --s .565677930574
nice -n 15 taskset -c 11 timeout 150s "$P" -Bu "$R/replay.py" b3f18_hi.bin --degree 32 --arity 7 --challenge --mask disjoint --s .565770634792
nice -n 15 taskset -c 11 timeout 25s "$P" -Bu "$R/regression.py" b3f18_hi.bin
nice -n 15 taskset -c 10 timeout 140s "$P" -Bu "$A" --full b3f18_lo.bin
nice -n 15 taskset -c 11 timeout 140s "$P" -Bu "$A" --full b3f18_hi.bin
```

- All seven parent manifest entries passed. The three EDGE source hashes matched the freeze. Producer runs authenticated both whole payloads, evaluated `91732` / `44492` Bellman states and passed `119` full-tree comparisons per sign with `391` nonzero edge births per sign. Both printed the outward endpoints above.
- The companion (`neon-edge-sol.py`) uses the separate physical coefficient oracle (`occupation-neon/all-contrast/check/row_oracle.py`), which must be available alongside the frozen inputs. It checked `2079` first-body rows, five stock floors, `43776` inner charged edge rows and `91` edge port comparisons across inner/mixed/tail cases per sign. It independently contracted all `184` / `105` mixed reserve rows in the respective cuts. An omitted `f→g` source control produced `3` / `88` negative sampled rows. This control is distinct from the old signed-leaf surrogate incident; the latter's profile and toll are not dependencies of this certificate. The companion rebuilds whole-tree gain and promotion differences on labelled cases, including unequal exterior donors, old and new edges, and pre/post-promotion and founder controls. It reuses the exact Bellman solver for the final table and independently composes the result from the frozen reserve fractions.
- The obligatory upper tail-floor mutation `861071198001/10^15` in place of `861071198001/10^18` gives a misleading positive local margin at `.565770220167`. Both the producer regression and the independent physical difference recover the exact signed drift `−902875978288208312083656815181293852064686619898589/1097886886904005745377459200000000000000000000000000` on parents `(-1,0,1,2,1,0,5)`, leaf stocks `(62,63,0,1,4,65,1)`. The authenticated floor gives a negative census-only margin at that withdrawn trial. The same tree passes the full reserve-inclusive inequality at both new trials.
