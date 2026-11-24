**Claim 10 (`6c8afec`): HOLDS**, conditional on Gold2, the repaired Claim 3 coboundary and positive mean growth, and Astra’s identification with the stationary core law.  
**Claim 11 (`93484ba`): HOLDS**, including the explicit finite-$`n`$ core-height tail and both almost-sure diameter bounds.

This was a cross-family audit: GPT-6.2 Sol Ultrafast checking Gemini-family Neon-2. In 1,000 independent trees at $`n=10^6`$, none of the eight marginal estimates was more than **1.284 run-level standard errors** from its nearest certified endpoint.

### Audit

| Item | Finding |
|---|---|
| **10 — signs, bracket, and stationary reduction** | Both snapshots bound the *same direction* of drift; only the observable forcing changes sign. For each fixed tree and witness, the residual is affine in $`\lambda`$, so taking the larger endpoint error covers the full interval even for a signed potential. At $`\lambda=\mu`$, adding $`\sigma R_{\psi, \mu}`$ cancels the empirical sum with the correct sign. Stopped Dynkin is dominated at fixed time by the Yule size, and the terminal contribution divided by the discounted occupation integral is $`O(1/t)`$ by the positive mean asymptotic. This uses more than the exponent bracket alone; the detailed proof supplies the dependencies that its opening paragraph understates. |
| **10 — analytic and bare-core tails** | Substitution gives the stated gain and promotion polynomials, including the mixed $`qx`$ term. The gain-rate extreme reduction is affine; the newborn-edge bound follows from the one-knot convex interpolation. For a negative Bernstein numerator bound, dividing by the **largest** denominator is the safe direction, and the script does that. The separate bare-core Lyapunov drift is bounded by $`40\,2^{-D_F}`$ per core at the $`D=5`$ seam and by $`-\mu`$ on the target tail. The resulting upper addend is strictly below $`10^{-479}`$. |
| **10 — identification and corollaries** | The nonfounder’s parent slot restores both degree counts, and Astra’s mean and empirical limits identify the degree scores and the bounded neighbor-reciprocal cylinder scores. The latter do require the spatial law; the one-sided fringe by itself does not record the exterior parent’s degree. Both degree-two balance equations check out, including the extra promotion contribution when the sole core has degree two. I also checked exact containment of the $`p_2`$ interval in the $`J+(2,2)`$ interval, rather than just overlap. |
| **11 — upper bound and tail** | The two physical jumps give an affine residual in the gain rate, so both extremes on all 5,148 finite rows suffice. At $`D=100`$ the gain uses the checked $`D=101`$ boundary. Beyond the cutoff the old core’s weight is identically one; the only positive contribution is the newborn term, bounded by $`zf(2,1)\le\lambda`$. The product-normalized potential is a supermartingale. Markov on dyadics gives a summable bound for every positive exponent; monotone core height and a countable sequence of thresholds give the claimed limsup. |
| **11 — embedded CMJ** | A real leaf has one neighbor, so only that neighbor can promote it. At each active core’s self-ring, the true leaf-choice probability dominates the tracked threshold. Removing the actual promoted leaf if designated, or an adapted designated leaf otherwise, preserves the leaf-stock inequality even on a physical promotion that produces no embedded child. A promoted leaf leaves its old driving set before becoming a driver. Disjointness **plus** the fresh marked clocks at stopping times gives the history-independent rates of the nonexplosive process of independent CMJ individuals. |
| **11 — speed, diameter, restarts, and clock** | The BRW has i.i.d. reproduction point processes, bounded offspring, mean $`>1`$, and finite positive tilts for displacements $`-\tau_i`$. Sibling ages need not be independent. The exact convex check and the first-emission lower bound give a strict speed inequality at a finite maximizer. The repaired binomial argument gives two surviving branches; the speed event can be taken simultaneously over all subtree roots. Extinction is detected in finite time on a failed trial, promotion hazard is at least one, and each fresh trial has the same conditional success probability. Finally $`ne^{-\tau_n}\to W>0`$, so $`\tau_n/\log n\to1`$. The Gamma envelope in the README should read $`\sigma_{v, J_{\max}}-s_v`$ because the ring times there are absolute; the calculation uses elapsed ages. |

### Replays and plants

All builds and runs used clean source archives under `round2/children/n2x/`; snapshot payloads were read in place and authenticated against both pinned digest lists before use.

- Strict build: `g++ -std=c++17 -O3 -fopenmp -Wall -Wextra -Werror -pedantic neon2-6c8afec/round2/neon2/stationary-marginals/verify_marginal.cpp -o verify_marginal` exited successfully. `python3 -OB replay_pair.py` and `python3 -OB replay_pair.py --also-bare-upper` each passed `--threads 1` to C++. All replayed fields matched `certificates_all.json`, including 35 exact tree-oracle checks per snapshot.

  | Snapshot | Error ticks at $`\mu_{\rm lo}, \mu_{\rm hi}`$ | Outward endpoint |
  |---|---:|---:|
  | `obs1_lower.bin` | `15388075, 11767580` | `0.38337348` |
  | `obs1_upper.bin` | `13356179, 14450017` | `0.38351583` |
  | `obs8_upper.bin` | `38490159, 15797774` | `0.14540145` |

  Each replay covered 49,204,104 finite motif rows per endpoint and 8,035,467 shared envelope checks. Separately, `python3 -OB bare_check.py` printed `bare-tail endpoint rows checked: 2621930` and `exact tail bound < 1e-479: True`.
- `python3 -OB neon2-93484ba/round2/neon2/height-diameter/cert_height.py` printed `HEIGHT_DIAMETER_EXACT_CHECKS_PASSED`: 5,148 upper rows, 7,500 CMJ rows, 23 trees / 59 rooted configurations, and 11/11 built-in negatives rejected. It reproduced the exact prefactor `861558416425/171115606632`, height endpoints `5000/5067` and `1.462248`, and diameter endpoints `10000/5067` and `2.924496`.
- `python3 -OB plants.py` rejected the tail-script mutation $`1/2\to1/3`$ with `tail oracle mismatch at D=1600, L=0`. Raising $`z`$ by $`10^{-5}`$ on the frozen height table was rejected: `far-tail condition failed: z*f[2][1]=6617573506752451/4000 > lam*S=1654390000000`. My first probe, a $`10^{-9}`$ increase, was accepted; it did not exhaust the certificate’s slack.

### Independent simulation

I built `sim.cpp` with `g++ -std=c++17 -O3 -Wall -Wextra -Werror sim.cpp -o sim`. The expanded run used `taskset -c 10 timeout 1200 ./sim 1 500` and `taskset -c 11 timeout 1200 ./sim 501 500`, under `nice -n 5` and a 2 GiB virtual-memory cap per worker. I observed about 80 MB RSS per worker. Both workers completed all 500 seeds. `python3 -OB sim_stats.py sim500_1.log sim500_501.log` reported **1,828,810 pooled cores** at $`n=10^6`$. The SE is for the pooled ratio, estimated from independent-run residuals $`X_i-\hat p C_i`$, where $`X_i`$ is the observable count in run $`i`$. The signed score uses the **nearest exact endpoint**, including when the estimate is inside.

| Observable | Estimate at $`10^6`$ | Run SE | Signed endpoint distance / SE | 12-run pilot distance / SE |
|---|---:|---:|---:|---:|
| $`p_2`$ | 0.38365604 | 0.00031562 | +0.444 | +1.031 |
| $`p_3`$ | 0.14514247 | 0.00026137 | −0.106 | +0.692 |
| $`p_4`$ | 0.08030304 | 0.00019886 | −0.140 | −0.233 |
| $`p_5`$ | 0.05201743 | 0.00015814 | −1.284 | −4.014 |
| $`p_6`$ | 0.03726631 | 0.00013953 | +0.444 | −0.519 |
| $`\pi(J)`$ | 0.25229740 | 0.00027145 | −0.012 | +0.233 |
| $`\pi(2,2)`$ | 0.13135864 | 0.00024070 | +0.548 | +1.509 |
| $`\pi(L_o=0)`$ | 0.14546071 | 0.00024840 | +0.239 | +1.799 |

The pilot is included in the 1,000 runs. Its $`p_5`$ discrepancy did not persist. The finite-size core-degree law at $`10^6`$ is not the stationary law, and its bias is unknown. A few-SE discrepancy would merit investigation, not refute the enclosure.

| $`n`$ | $`H/\log n`$, observed range; below / inside / above | $`\mathrm{diam}/\log n`$, observed range; below / inside / above |
|---:|---|---|
| $`10^5`$ | [0.43429448, 1.91089572]; 270 / 705 / 25 | [0.60801227, 3.04006137]; 471 / 528 / 1 |
| $`10^6`$ | [0.50667690, 1.80956034]; 138 / 853 / 9 | [0.86858896, 2.96767896]; 362 / 637 / 1 |

These comparisons are to **asymptotic** windows. Values on either side at a finite size do not violate them.

**Supplementary argument of my own.** `tilted-block.md` (`round2/children/n2x/tilted-block.md`) records a Biggins-free check of Claim 11’s lower endpoints. The certified convex inequalities give an exponentially tilted offspring-age law with mean strictly below $`5067/5000`$ and positive growth after undoing the tilt. The LLN then makes the expected number of generation-$`k`$ descendants born within $`k(5067/5000)`$ greater than one for some fixed $`k`$. Keeping these descendants every $`k`$ generations produces a supercritical block process; two surviving branches give the diameter lower bound.

**Scope:** I personally replayed three of the sixteen Claim 10 snapshots. I did not re-prove Gold2, Claim 3’s renewal and contraction inputs, or Astra’s spatial theorem. Neither verdict provides monotonicity past $`p_6`$, the degree-tail constant, or convergence of the normalized height and diameter to exact limits.
