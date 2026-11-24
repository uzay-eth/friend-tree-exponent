VERDICT: HOLDS

# Fable on Astra claim 1 (finite rational branching matrices G_{M,B}; s(G_{M,M^2}) -> mu with error bound (16))

**Scope of the verdict:** The Claude-family lead reconstructed one matrix by hand and read the seam budget, coupling and iteration. A Gemini child rebuilt six matrices exactly, tried six more, simulated the coupling on small examples and ran four planted controls; its saved output predates the quota stop. The limit assessment inherits the uniform newborn-branch and planting estimates, which the lead did not rederive. No GPT-family child evidence is used here.

Attacker: the Fable lead (Claude Fable 5.2, a different family from the GPT claimant), with the saved Gemini 4 Neon output in Fable folder `round2/checks/astra-claim1-recomputation-gemini-partial.txt`; the child did not deliver a separate verdict before the quota stop.

The claim in Astra folder `round2/finite/finite-operator.md` (`f7ab502`) concerns reachable rational matrices along `B=M²`, their explicit two-sided error (16), and the uniform comparison in Astra folder `round2/finite/one-port-upper.md`.

## What was recomputed

1. **(M,B) = (3,1) by hand (lead).** Types: founder; hub; singleton (D,L) = (2,1) with one inactive slot;
   singleton (2,0) with two slots. Founder: gain at rate 2 makes it a hub; promotion at rate 1 gives two (2,1)
   singletons after the cut: row (-3, 2, 2, 0). Hub: (0, 0, 1, 0). (2,1): gain at rate 1 to hub; promotion at rate
   1/2 leaves a (2,0) and a new (2,1): row (0, 1, -1, 1/2). (2,0): no transition. The class {hub, (2,1)} has matrix
   [[0,1],[1,-1]] with top eigenvalue (sqrt5 - 1)/2 = 0.6180339887. Agrees with the claim.
2. **All six matrices, entry by entry (Gemini child, own code from the written definition):** exact match for
   (3,1), (3,3), (3,9), (4,1), (4,2), (4,3): 4, 9, 27, 8, 23, 68 reachable types; spectral abscissae 0.618033988750,
   0.678337373185, 0.679280747970, 0.556123256049, 0.614850634815, 0.627558021428. Differences from the reference estimate 0.56572366: +0.0523, +0.1126, +0.1136, -0.0096, +0.0491, +0.0618. Six further matrices built by the child:
   (3,2) 0.670085, (3,4) 0.679188, (5,1) 0.532248, (6,1) 0.519650, (7,1) 0.511868, (10,1) 0.499961.
   With B = 1 the displayed values fall as M grows (0.618, 0.556, 0.532, 0.520, 0.512, 0.500),
   and the displayed fixed-M values rise in B. Only (3,9) is on the claimed sequence B = M², so these runs
   do not test its convergence.
3. **The error bound (16) evaluated (Gemini child).** For the six certified matrices the reported evaluations
   range from 84.00 to 107.65 in the free-parameter search and from 26,675.74 to 35,385.82 on the j=4
   schedule, compared with differences from the reference estimate of about 0.01 to 0.11. In the child's
   free-parameter search for total interval width below 0.0015, the estimated cutoff was M about
   10^(10^13795.5). This was a schedule-dependent calculation, not a global minimum. The six matrices
   do not supply an eight-digit enclosure.

## The proof, read by the lead

The seam budget s J + P <= C~ - 1 (section 3) is a correct potential argument (a size cut lowers P by exactly s; a hub
release cannot raise P). The window comparison (14) uses the same clean/dirty coupling as elsewhere in Astra's
work: exceptional rings of clean hubs have total rate at most 2C/M, sends across seams at most J <= C~/s, each
exceptional event dirties at most five cores and a dirty core dirties at most two per unit time, giving
(10/M + 2/s)(e^{2t} - e^t). The iteration over windows (15) needs the one-window bound conditionally on every
reachable state, which is how (14) is stated. I found no error. The uniform one-edge bounds U_j, L_j of
Astra folder `round2/finite/one-port-upper.md` were taken as inputs.

## Not done

The lead did not run planted errors. The saved Gemini output records detection of four selected faults: averaging over unique outcomes, dropping inactive slots, changing the hub newborn, and omitting a seam recipient's dirty mark. The implementation is in Fable folder `round2/attacks/astra-B/attack_claim1.py` (working files); for replay elsewhere, adjust its hard-coded input path to the local Astra folder `round2/outputs/finite-matrices.json`.
