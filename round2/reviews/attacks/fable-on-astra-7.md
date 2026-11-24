VERDICT: HOLDS

# Fable on Astra claim 7 (sum of the limiting hub fractions is 1; next-attachment law)

**Scope of the verdict:** The Claude-family lead read the supermartingale, stopped mass argument and attachment-law proof, taking the fixed-label limits and sublinear power bound as inputs. The lead simulated finite trajectories for promotion-ordered core degree only; the attachment distributions and almost-sure limits were not simulated directly. No GPT-family-only child check is used here.

Attacker: Claude Fable 5.2 (a different family from the GPT claimant). The Gemini child stopped on a provider quota; no separate result is used here.

## The proof, read line by line (Astra folder `round2/hub-mass.md`, `34b6d52`)

1. *Supermartingale.* Given the tree at size j, a core of degree D receives the next leaf with probability
   (L + sum over core neighbours of 1/D_w)/j <= D/j, so E[D(j+1)/(j+1) | T_j] <= (D/j)(1 + 1/j) j/(j+1) = D/j.
   Correct, for every core from its promotion on.
2. *Late cores carry little.* A vertex with label v > m is promoted at some size r > m, starts as a core with
   degree 2, i.e. normalised degree 2/r, and there is at most one promotion per step, with probability
   a_{r} - a_{r-1} <= a_{r-1}/(r-1). Summing 2 (a_r - a_{r-1})/r over r > m with a_r <= K r^lambda, lambda < 1
   (the mean bound in Astra folder `proofs/bracket-motifs.md` allows λ = 0.5665) gives O(m^(lambda-1)), uniformly in n. Correct; the optional
   stopping implicit in "conditioning on the tree just after each promotion" is the supermartingale property
   of each term from its own start, summed by linearity.
3. *Conclusion.* The cores carry total degree n + C_n - 2, so the first m labels carry, in expectation,
   at least 1 - O(m^(lambda-1)) - o(1) of n; each of these finitely many normalised degrees converges
   (Addario-Berry et al.) and is bounded, so E sum_{v<=m} Z_v >= 1 - O(m^(lambda-1)). Conversely
   sum_{v<=m} D_v(n)/n <= (n + C_n - 2 + m)/n -> 1 gives sum_v Z_v <= 1 almost surely. Hence equality a.s.
   I found no gap. (The second half of the note, which avoids the mean bound by stopping at
   tau = inf{n : C_n > K n^lambda}, uses a deterministic λ<1 from the older almost-sure power bound.)
4. *Total variation.* P(W_n = v | T_n) = (L_v + sum_{core nbrs} 1/D_w)/n; the second term is at most C_n/n -> 0;
   pointwise convergence of a sequence of probability vectors to a probability vector is convergence in total
   variation (Scheffe). For the weak rescaled-degree law, apply a bounded continuous test: on finitely many labels the weights and degrees converge, and the remaining weight is small by total variation and unit limiting mass.

## Numerical look (own simulator, single trajectories, continuous time, n up to 1.7e11)

Cumulative degree of the first K cores in promotion order, divided by n−1 (twice their share of total degree 2(n−1)); 12 trajectories, three shown:

    run  t   n          K=1       K=10      K=100     K=1e3     K=1e4     K=1e5
    0    16  7.6e6      0.784171  0.973317  0.994831  0.999499  1.000385* 1.000385*
    0    20  4.1e8      0.784012  0.973375  0.994777  0.999080  0.999908  1.000070*
    0    24  2.3e10     0.783998  0.973383  0.994773  0.999064  0.999835  0.999985
    3    24  1.2e11     0.352061  0.861412  0.974118  0.995710  0.999230  0.999874
    5    24  4.7e10     0.239111  0.795546  0.959683  0.993721  0.998922  0.999828

An asterisk marks an all-core subtotal, slightly above one with this normalization. In the displayed run 0, the K=1 subtotal changes from 0.784171 to 0.783998 and the K=1000 subtotal from 0.999499 to 0.999064. These finite-time values are consistent with most degree mass residing among early-promoted cores; K in this table is not the cutoff m on birth labels in the proof.
Command from Fable folder `round2/`: `bin/hubmass 9700 12 16 20 24`; source `src/hubmass.cpp`, working output `outputs/hubmass.txt`, and build `build.sh` in that folder.

## Not done

No planted-error test or direct attachment-law simulation was run. The proof read used the power bound and fixed-label limits as inputs.
