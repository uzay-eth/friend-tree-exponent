VERDICT: HOLDS

# Fable on Astra claim 4 (mean semigroup Q_t on rooted marked cores; r(Q_1) = e^mu; eigenmeasure)

**Scope of the verdict:** The Claude-family lead checked two moment bounds and read the spatial construction, eigenmeasure and spectral-radius arguments, taking the planting and one-edge estimates as inputs. Fable simulations give finite-range consistency checks, and a saved Gemini child log gives finite cutoff simulations and two fault tests. The infinite-volume result was not recomputed numerically and no GPT-family-only child evidence is used here.

Attacker: Claude Fable 5.2 (a different family from the GPT claimant). The Gemini child stopped on a provider quota without a separate verdict; see the saved Fable folder `round2/attacks/astra-B/claim4_out.txt` for cutoff simulations, two planted errors and the right-eigenfunction example. The claim concerns the spectral radius on C_b, a probability eigenmeasure for every Q_t, uniqueness of the positive eigenvalue among probability eigenmeasures of Q_1, and the spatial-cutoff error (3).

## Steps examined (Astra folder `round2/operator.md`, `8e4db79`)

1. *Cells and semigroup.* Each new core inherits the label of its promoter, so the cells partition the future
   core and Q_s Q_t = Q_{s+t} by conditioning; Q_s 1 <= e^s because every core promotes at rate at most one.
   Correct.
2. *Finite speed of propagation (3).* A cell is influenced from outside only through sends along old edges,
   implemented as rate-1/c^0 proposals thinned by c^0/D; a change reaching w from distance j needs j time-ordered
   proposals along an old path, expected number (s^j/j!) x (product of 1/c^0), and the sum over simple paths ending at w
   is at most c^0_w by reversibility of the simple random walk. This is the right estimate and makes Q_s f
   continuous for bounded cylinders. The cutoff in (3) uses path length j=R−k+1 and the sum of core degrees over the old ball of radius k.
3. *Eigenmeasure.* Tilted time averages of the root-uniform mean measures with weight e^{-mu t - delta|t-T|};
   tightness from the uniform bounds E[D^alpha | age a] <= 2^alpha e^{alpha a}, E[e^{theta c} | age a] <=
   e^theta e^{(e^theta - 1) a} and the mass-transport identity of root-uniform measures. I checked the two moment
   bounds (gain rate <= degree; core degree grows only by the vertex's own promotions) and the exponent
   conditions (alpha < mu, theta < log(1 + mu) make the tilted integrals converge). Correct.
4. *Minimality and the spectral radius.* From any (G,o) the cell of o can plant, with positive probability, a
   finite seed behind a parent of degree >= M, and the seed grows at rate mu - eta with eta -> 0; hence
   liminf (1/s) log Q_s 1 >= mu pointwise and any probability eigenmeasure has eigenvalue >= e^mu. The upper bound
   Q_t 1 <= 1 + t A_j(t) e^{Lambda_j t} F(t+1) uniformly over the exterior gives r(Q_1) <= e^mu. Together: every
   probability eigenmeasure with positive eigenvalue has eigenvalue exactly e^mu.
   The indicator that the root is the only core is a right eigenfunction of Q_t with factor e^{-t}
    (see Fable board fact F15 and the Gemini log). This is not an eigenmeasure: the lower bound on
    growth of Q_t 1 at each state rules out probability eigenmeasures at a smaller positive rate.

## Consistency with simulations (lead)

Under pi, the proof gives mean L_o/D_o equal to μ and the limiting degree distribution from Claim 3. The Fable production result near the reference estimate 0.56572366 is a variance-reduced growth estimate, not a direct average of L_o/D_o under pi. The finite degree histograms are consistent across the 12→16 and 16→20 windows (see `fable-on-astra-3.md`). Fable folder `round2/attacks/astra10-11/eval_16_20.log` gives time-differenced shares 0.008569, 0.005265, 0.003284, 0.002062 and 0.001306 for core degree 7 through 11, a short range consistent with a decaying tail. These shares do not test the asserted exponential moment.

## Not done

The uniform one-edge upper bound (Astra folder `round2/finite/one-port-upper.md`) and the planting lemmas of Astra folder `proofs/existence.md`
were taken as given.
