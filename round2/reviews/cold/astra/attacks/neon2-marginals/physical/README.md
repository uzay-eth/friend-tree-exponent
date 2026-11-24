# Neon-2 Claim 10: independent physical and rational audit

**Verdict: HOLDS in this audit's scope**, conditional on the stated Gold2 bracket and physical mean/fringe identifications. Claim pin `6c8afec5cd1b7f98bbcd7c6e7848c0d631a50f7d`; the sanctioned tar has SHA-256 `2f2c2c1ee81577255aa4a8834b7c97850e7445f92008351bc5d2862e35cf7177`. All ten claim files match the archive tree provided by `863c` byte for byte. Receipt hashes were joined to their independent authentication record; their full payload replay is a separate audit.

## Observable and row meaning

- The denominator is the number of physical cores, not all vertices, not the number of motifs, and not the promotion flux. `D` includes leaves; `c` counts core neighbours. The parent slot in the downward projection restores both marks for a nonfounder, but the actual reciprocal from that neighbour is live in the full tree. `J` means `(2,1)`: the one core neighbour is the parent in the limiting one-ended orientation. `(2,2)` means an unstocked degree-two core; the second ratio is a **sum** of two reciprocals, not their mean. `J` is a mark class at observation time, not the birth-weighted law. The observable `D=2` includes the initial singleton `(2,0)` on finite free trees, whose normalized occupation vanishes.
- The mark gain `(D,c)→(D+1,c)` has rate `L+rho`, `rho=sum 1/D_neighbor`; promotion `(D,c)→(D,c+1)` has rate `L/D` and creates exactly one new `J`. In the verifier, the signed source is only on `f`, in units of ticks: `sigma*S*(psi-theta)`. A promotion accounts for the change of the promoter's `f`, the newborn's `f(J)` and the one new edge with the promoter at the **postpromotion** arity. The finite centered gain rate replaces `rho` by `c*r1`, with `(1/D_v-r1)*Delta_f` on the existing physical edge for each true neighbour. Outside the centered sector the sign of the gain jump selects `L` or `L+c/2`; no shift of the observable to the edge rows. At the seam the gain uses the full tail polynomial at `(q,x)=((DF-c)/DF,1/DF)`; in particular the `q*x` coefficient enters with numerator `DF-c`. IDs 1–7 vanish above the finite range; ID 8 is bare **within** that range. For bare high-degree cores promotion vanishes, gain reduces the tail Lyapunov weight, and the seam constant is bounded by the maximum `(5-c/2)4^{-(6-c)}=5/8`. The tail correction is applied only to the upper bound for ID 8.
- At the claimed stationary `mu`, the stopped-Dynkin argument divides by the discounted core occupation, asymptotic to a positive constant times Time, while the linearly bounded potential contributes a bounded endpoint term. The stated exact free-tree coboundary cancels the true physical score, including the founder score, with the correct sign. Alternatively, with physical spatial mean convergence for the bounded local score the same Cesàro ratio identifies `pi(psi)` directly. The projection restores marks, not the exterior degree; the ratios follow from physical rows, not from treating a planted inert parent as a source of arrows.
- In finite free trees, `L f_J = L C - (3/2)f_J - Z_J + 1[(D,c)=(2,0) at the sole core]`. At that exceptional promotion both the founder and the newborn enter `J`. Its accumulated discounted occupation is finite since the sole core promotes at rate one. The `(2,2)` row is `L f_22=f_J/2-Z_22` with no founder correction. Removing the bounded endpoint terms gives `Z_J/pi(J)=mu*(1/pi(J)-1)-3/2` and `Z_22/pi(2,2)=pi(J)/(2*pi(2,2))-mu`, where here `Z_J,Z_22` denote per-core mean weighted counts. All four ratio corners follow by monotonicity in the positive denominators and the bracket variable. The degree-two partition creates no equality between independently certified endpoints; it is an equality between the true masses.

## Evidence and consequences

`check.py` recombines the recorded endpoint caps and signs independently, verifies the published exact corollaries and decimal direction, and enumerates actual friend-tree arrows to check the count drifts, including the founder correction. From `[local path]`:

```text
$ taskset -c 5 nice -n 19 python3 -IB round2/attacks/neon2-marginals/physical/check.py .work/build/neon2-claim10 --tiny round2/attacks/neon2-marginals/physical/.tiny.bin
PASS 16 authenticated receipt joins and exact endpoint combinations (no payload scan)
pi_D_ge_7 0.30130165 0.30174603
pi_bare_D_ge_3 0.01397965 0.01431408
E_invD_parent_given_J 0.17650729 0.17678896
E_sum_invD_nbrs_given_2_2 0.39548710 0.39662379
PASS degree-2 decomposition, four monotonicity gaps, bare fraction and analytic tail
J and reciprocal support alone allow mu in 0.5060896177783607 0.6748714860065338
PASS 872 increasing trees, 10076 directed physical arrows, degree-2 and leaf-stock rows
PASS J singleton plant lo: f_in -411 -> 1753907112; g 473502 -> 622650699
PASS J singleton plant hi: f_in 5681501 -> 1753930613; g 961801 -> 624573726
```

The stock verifier's `--perturb f` subtracts `10^9` ticks from `f(J)` on a tiny **upper J** proposal. The finite singleton and edge caps both respond at both endpoints; the stock verifier reports a changed evaluated hash and `perturbed:true`. An adaptive verifier is expected to report wider numerical caps rather than reject merely because a weight moved; `combine_side` explicitly refuses perturb-marked output. Tiny snapshot replay is optional; reconstruct with the pinned `build_marginal.cpp` compiled with `g++ -std=c++17 -O3 -fopenmp`, then run the resulting executable as `OMP_NUM_THREADS=1 build_marginal tiny 6 upper 565713486333 1000000000000 565736986679 1000000000000 .tiny.bin` (supply a path in this directory).

The ratio intervals use the Gold2 exponent bracket as an input. Feeding the same ratio intervals back through their defining equalities is not an independent exponent estimate. [`stationary.md`](stationary.md) gives exact leaf-stock, oriented degree-two edge, and downward-size constraints that can instead be combined with independently bounded fluxes or used as shared rows of a stationary feasibility model. No new numerical exponent endpoint was certified.
