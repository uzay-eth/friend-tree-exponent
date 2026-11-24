# Astra rival attack on Neon Claim 2

**VERDICT: HOLDS.** Scope: all three posted subparts of Neon Claim 2, pinned to `ft2-neon@f7b6074ec93992c2375679e5057918d0cd37fb68`, `round2/open_questions/REPORT.md` §§1 and 3. This is a separate Astra-side attack, not the Neon-2 lead review. REPORT §2 (core-limit tails), simulations, and Claim 1 are outside scope.

**Local Astra attack commit:** `ft2-astra@989bfd4d569f3e4f64beebb76d7b8b9d6f67ec53` (`round2/attacks/neon2claim/` only). **Separate original Neon pin:** the Claim 1 production audit targets `ft2-neon@cec3855cc52ab4265ce3e36b8a7a38f1ae257f77`, with Astra attestation `de931f5` ([Claim 1 report](astra-on-neon-1.md)). This Neon-2 report targets the numbered Claim 2 at `f7b6074`.

| Subpart | Verdict | Finding |
|---|---|---|
| (i) Hub characterization | HOLDS | Independent first-degree-hit test: from a new hit at degree `d≥3`, `P(Z=0 | history at hit)≤1/d`, via a tagged birth–death family and a conditional large-stock urn bound. Checked the rejection-chain likelihood and simultaneous stock comparisons in the positive-probability freeze proof. Promotion does not suffice. |
| (ii) Exact means, bias, and separation | HOLDS | Independent integer donor-arrow enumeration reproduced both `n=8` endpoints and the gap. At `n=10`, the same argument improves the certified gap to `3337671095981771/82959593472000000` (~0.04023). Both directions of the reachable degree-tie rate reversal checked. Finite marginal ordering through `n=10` establishes no general limiting order. |
| (iii) Independent set and fixed-label forest | HOLDS | Directed conditional Borel–Cantelli gives hubhood of every neighbor of a finite-degree label. The singleton/deletion statement even holds for every diverging sublinear cutoff sequence on the same event. This does not describe a sampled moving-label component. |

The early mean separation is distinct from Astra's separated-label survival comparisons (`34db7ad`) and its subsequent eventual adjacent hub-probability result (`298c54c`). The fixed-label forest law is distinct from Astra's sampled-core and counted-component limits (Claims 19 and 20).

**Evidence:** `ft2-astra/round2/attacks/neon2claim/{README.md,enumerate_claim.py,exact.log}`. `taskset -c 3 nice -n 15 /usr/bin/time -f 'wall=%e s RSS=%M KiB' python3 -B enumerate_claim.py` passed, enumerating `362880` histories at size ten (`wall=3.90 s RSS=13696 KiB`). No packages, GMP run, or simulation used for this verdict.

**Portable replay:** `python3 -I -B replay.py` (Python 3.9+, standard library only), from the attack directory. A copy of just the four attack files, run outside both checkouts, printed `PASS: exact Neon Claim 2 transcript matched through n=10 (362880 size-ten histories).`
