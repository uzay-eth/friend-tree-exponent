VERDICT: HOLDS

# Fable on Astra claim 3 (E[C_n] ~ c_2 n^mu; mean degree law p_k ~ mu A k^(-1-mu))

**Scope of the verdict:** The Claude-family lead read the sharp-mean and degree-law proofs, using the renewal construction, local operator and cited Tauberian theorem as inputs. The lead compared finite-time mean and degree statistics with the predictions; the simulations do not establish asymptotic rates or limits. No completed child calculation is relied on here, and no GPT-family-only check is used.

Attacker: Claude Fable 5.2 (a different family from the GPT claimant). The Gemini child stopped on a provider quota; no result from it is used here. The proof files are Astra folder `round2/fringe/asymptotics.md`, `round2/degree-law.md` and `round2/operator.md` (`8e4db79`). The claim concerns each fixed k≥2: E[N_{≥k}(n)] ~ c_k n^μ with c_k>c_{k+1}>0, and the ratios of mean point counts to mean core count tend to p_k, with p_k ~ μA k^{-1-μ} after the fixed-k limit.

## Steps examined

1. *Transform (3) and the pole.* For Re z > mu the generation series converges because r(K_{Re z}) < 1, giving
   the resolvent formula; eigenvalue 1 of K_mu is simple (strong positivity, compactness), with
   lambda'(mu) = l(K'_mu g) < 0 because the kernel is a positive measure in time. Correct.
2. *No other singularity on Re z = mu.* If K_{mu + i tau} v = v then |v| <= K_mu |v|, forcing equality and a
   constant phase of e^{-i tau t} v(1/2) against a component of the kernel that has a density bounded below on
   compact time intervals (the newborn, untouched up to time t, gives birth with mark 1/2 at rate 1/2). The
   density claim is right: from J_a the event "no ring of the root's leaf, no ring of the root, no killing up
   to t" has probability e^{-(2+a)t} and the direct-birth rate is then 1/2 with mark 1/2. So tau = 0. Correct.
   The positive-density argument here is analytic; the finite kernel simulation in `fable-on-astra-2.md` was a separate consistency check.
3. *Tauberian step.* M_t f_k is non-decreasing in t for the counts of degree >= k (degrees never decrease, new
   cores have degree 2), so Wiener–Ikehara applies to x -> M_{x/mu} f_k; local analytic continuation of the
   pole-subtracted transform is more than the theorem needs. Correct.
4. *Fixed n.* E_{T_2} N_{>=k}(t) = q^2 sum_n (n-1)(1-q)^{n-2} a_{n,k} with q = e^{-t} (size from an edge is the sum
   of two unit Yule processes, independent of the jump chain); Karamata-type transfer for the non-decreasing
   sequence a_{n,k}. Correct.
5. *Tail of p_k* (degree-law.md). pi(D >= k) = E_beta[e^{-mu T_k}] with T_k the hitting time of degree k by a
   newborn sampled at promotion, k e^{-T_k} -> W, giving A = E_beta W^mu; point masses from the exact flux identity
   mu pi(D >= k+1) = pi[(L + rho) 1{D = k}] and (L + rho)/k -> 1 except on an exponentially unlikely event. The flux
   identity is the stationarity equation for the indicator of degree >= k+1 (in my notation: the estimating
   equation of that feature) and holds exactly. Correct.

## Numerical tests (own simulators; none of these numbers is used by the claimant)

* **Mean.** For the continuous-time 3-path start, averaging seeds 301 and 302 (3 million runs each)
  gave E[S(t)] - 0.565723 E[C(t)] = 0.0937, 0.0963, 0.0955, 0.0945, 0.0940 at t = 5..9,
  with standard errors about 0.00051–0.00161. Here S = sum_u L_u/D_u. These finite offsets are
  consistent with the leading exponential and suggest, but do not establish, an additive correction
  of about -0.17. The reported amplitude estimate A_CT = 1.296 ± 0.004 (Fable folder `round2/RESULT.md`)
  gives c_2 = 2 A_CT/Gamma(3+μ) ≈ 0.725 ± 0.002 using the reference estimate 0.56572366; the underlying
  amplitude run is not linked here.
* **Degree law.** Time-differenced histograms over 40,000 runs between t = 16 and 20 (3.8e9 differenced core
  vertices): p_2 = 0.383456(8), p_3 = 0.145227(7), p_4 = 0.080368(5), p_5 = 0.052246(4), p_6 = 0.037181(4),
  p_8 = 0.022054(3), p_10 = 0.014849(2), p_20 = 0.004512(1). An independent sample at the earlier window 12 -> 16
  (100,000 runs) gives 0.383456(14), 0.145216(11), 0.080376(8), 0.052252(9), 0.037180(7) for k = 2..6: the same
  values within the reported sampling errors at these two finite windows.
* **One tail constant for both laws.** k^mu P(D >= k) = 0.8029, 0.7658, 0.7516, 0.7436, 0.7406 and
  p_k k^(1+mu)/mu = 0.872, 0.8015, 0.774, 0.757, 0.751 at k = 20, 50, 100, 200, 300 (p_k averaged over integer j with |j-k|≤floor(k/20)): both series decrease over the displayed range, consistent with a shared limiting constant; 0.735 is a rough extrapolation.
* **Single trajectories** (fable-on-astra-6.md): the reported finite-tree fractions for degrees 2 and 3 (Fable folder `round2/outputs/hubmass.txt`) are consistent with the pooled estimates.

From Fable folder `round2/`, the later-window command was `bin/degree_law 9500 0 40000 200 16 20 outputs/degree_law_1620.txt`. Evaluate either working histogram (`outputs/degree_law_1216.txt` or `outputs/degree_law_1620.txt`) with `.venv/bin/python src/degree_table.py FILE`. The earlier histogram’s simulator seed is not recorded here. The offset inputs are `outputs/beta_check_seed301.txt` and `outputs/beta_check_seed302.txt` in the same folder, for `bin/beta_check SEED 3000000 9 0.565723` through t=9. Fable folder `round2/build.sh` builds the simulators.

## Not done

Debruyne–Vindas’ Tauberian theorem was taken as cited. There was no planted-error test and no sensitivity study for alternative asymptotics.
