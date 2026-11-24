VERDICT: HOLDS

# Fable on Astra claim 6 (C_n / n^mu -> X almost surely, 0 < X < infinity, non-degenerate; degree fractions a.s.)

**Scope of the verdict:** The Claude-family lead read the strong-limit proof and checked its coupling cases and contour-shift hypotheses, taking the renewal results as inputs. The contour estimate and uniform planted-state error were only read. Finite runs track normalized core counts and selected fractions; they do not test almost-sure convergence, positivity or the support extremes. No GPT-family-only child check is used here.

Attacker: Claude Fable 5.2 (a different family from the GPT claimant). The Gemini child stopped on a provider quota; no separate output is used here. The claim in Astra folder `round2/fringe/almost-sure.md` (`c98885f`) starts from one edge, gives an almost-sure positive finite limit with positive mass arbitrarily close to zero and arbitrarily far above it, and simultaneous cumulative and point degree-fraction limits for all fixed k≥2, with the constants of Claim 3.

## Steps examined

1. *Section 1, exponential remainder.* K(z) is the 19 x 19 Laplace transform of a time kernel bounded by
   C e^{rho t}; the claim is K(z) = O(1/|z|) on Re z >= rho + eps, a strip sigma <= Re z <= mu free of zeros of
   det(I - K(z)) other than z = mu, and then v(t) = e^{mu t} D + O(e^{sigma t}) by shifting the inversion contour.
   The zero-free strip follows from: analyticity, no zero on Re z = mu except mu (non-lattice argument), and the
   decay at high frequency. I see no gap in the logic. My own Monte Carlo of this 19-label kernel
   (fable-on-astra-5.md) is consistent with a crossing near the reference estimate 0.56572366. The positive density used for the non-lattice argument comes from the analytic first-sojourn calculation.
2. *Section 2, the martingale.* Psi(T) = lim e^{-mu u} E_T C(T_u) / d_2 exists for every finite tree with error
   O(C(T) e^{-(mu - sigma) u}), hence E_T Psi(T_h) = e^{mu h} Psi(T). This is the harmonic function of the mean
   dynamics; its existence with a uniform error is exactly what makes e^{-mu t} Psi(T_t) a martingale.
3. *Section 3, the coupling bound (6).* Adding one leaf makes at most three cores "dirty"; a ring of a dirty
   core can add at most two marked copies (both copies of a clean recipient), or one unmatched newborn; rings of clean cores are synchronised and a clean core sees identical
   own neighbourhoods in both trees, so differences in the counts are confined to dirty cores, whose expected
   number grows at most like 3 e^{2u}. I checked the case analysis: (a) leaf of a dirty core rings: stays inside
   the dirty core; (b) dirty core picks a core neighbour: that neighbour becomes dirty; (c) dirty core picks a
   leaf: dirty newborn; (d) clean core picks a dirty neighbour: both copies of that neighbour gain a leaf, no new
   dirt; (e) clean core receives from a dirty neighbour in one tree only: counted under (b). Correct.
4. *Second moment.* The martingale increments are driven by vertex clocks at total rate n(T_v), whose conditional mean given T_0 is n(T_0)e^v, each event moving
   the conditional mean by at most 3 e^{2(b-v)}; the recursion for E (Z*_{j+1})^2 has additive terms of order
   e^{-(2 mu - 1) j} times a polynomial, summable because 2 mu > 1. This is the one place where the size of mu
   matters: the argument needs mu > 1/2 (certified: mu >= 0.5650). Positivity of the limit by the conditional
   second-moment bound, also using 2 mu > 1. Correct as far as I can check.
5. *Sections 4-5.* Borel–Cantelli on a mesh with the same variance bound; monotonicity of the counts between
   mesh points; transfer to fixed n through n e^{-tau_n} -> Y > 0; support from stars and paths. No gap found.

## Numerical look (own simulator; 12 single trajectories followed to t = 16, 20, 24)

    run   C/n^0.5657237 at t=16, 20, 24      fraction of cores of degree 2 / 3 / with no leaf at t = 24
    0     0.37470  0.38530  0.38273          0.38236  0.14484  0.14575
    1     0.44235  0.43984  0.44000          0.38345  0.14490  0.14495
    2     0.86170  0.86139  0.86120          0.38302  0.14517  0.14542
    3     0.90751  0.90294  0.90257          0.38339  0.14545  0.14547
    5     1.10101  1.09494  1.09330          0.38341  0.14562  0.14549
    8     0.38826  0.38291  0.38172          0.38421  0.14578  0.14582

Across the 12 runs, the last-window changes were below 1% and the t=24 core counts ranged from roughly 3e5 to 2e6; the displayed normalized counts differed across runs by up to a factor of about 3. The sampled degree-2, degree-3 and no-leaf fractions were near the pooled estimates 0.38345, 0.14523 and 0.14531 at t=24 (differences of order 1e-3), with larger scatter at t=16 (order 1e-2 in these runs). These are finite-time comparisons; the no-leaf fraction is an additional diagnostic.
Command from Fable folder `round2/`: `bin/hubmass 9700 12 16 20 24`; source `src/hubmass.cpp`, working output `outputs/hubmass.txt`, and build `build.sh` in that folder.

## Not done

The contour-shift estimate (2) and the uniformity of the error in (3)-(4) over all planted states were read,
not recomputed. No planted-error test (the claim is a proof).
