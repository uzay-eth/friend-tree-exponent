VERDICT: HOLDS

# Fable on Astra claim 3, second attack (3b): reconstruction of the sharp mean and the degree law

**Scope of the verdict:** Claim 3 stands behind claims 6 and 7 in the owner's order of 10:05 UTC on 24 Nov, so it was reconstructed together with them. A GPT-6.2 Sol Ultrafast child (R3; the claimant's model family) received definitions and statements only and wrote a complete proof of the three parts of the claim; it was allowed to use the statements of claims 2 and 5 and the uniform mean lemmas of claim 6, which other children proved blind (see `fable-on-astra-6b.md`). A second GPT child compared the claimant's proofs (`fringe/asymptotics.md`, `degree-law.md`, `operator.md`) with that proof step by step and derived the nine groups of estimates of the claimant's texts that the blind proof does not write out. The Claude-family lead checked by hand the passage from continuous time to fixed n with its constants, verified the one slip reported below against the round-1 file it refers to, and had measured the constants of the claim earlier with its own simulator. The lead did not read the blind proof of the spatial eigenmeasure line by line. So: two independent GPT-family proofs, a GPT-family comparison, and Claude-family checks of the fixed-n passage and of the numbers.

Attacker: the Fable lead (Claude Fable 5.2). Files: Fable folder `round2/checks/astra6-reconstruction/blind-proof-claim3-S3.md`, `comparison-claim3.md`.

## 1. The route of the blind proof, and what the lead checked

1. Continuous-time means from the three-vertex path: Laplace transform of the mean through the renewal kernel on exact marks, a simple pole at mu with residue d_f, no other singularity on Re z = mu (non-lattice argument in the column of mark 1/2), Wiener-Ikehara for the non-decreasing mean. Gives E_{T_3} N_{>=k}(t) ~ d_{f_k} e^(mu t) with d_{f_k} > d_{f_{k+1}} > 0 and d_{f_k} -> 0.
2. Start from the edge: the first ring comes after an Exp(2) time H and always produces the path, so d*_k = d_{f_k} E e^(-mu H) = 2 d_{f_k}/(2 + mu). Checked.
3. Fixed n: the size from the edge has the law P(Y_t = n) = (n-1) x^2 (1-x)^(n-2), x = e^-t, independent of the shape; hence the exact identity sum_{m>=0} (m+1) a_{m+2} r^m = (1-r)^-2 E_{T_2} N_{>=k}(-log(1-r)) ~ d*_k (1-r)^(-mu-2) for a_n = E N_{>=k}(n). Karamata's theorem for the partial sums, then monotonicity of (m+1) a_{m+2} to pass from partial sums to terms: a_n ~ d*_k n^mu / Gamma(mu + 2). Checked by hand, including the exponent mu + 2 that comes from the initial size two. The claimant uses a tilted gamma measure instead; the two give the same constants c_k = d*_k/Gamma(mu+2).
4. p_k = (c_k - c_{k+1})/c_2 > 0, sum p_k = 1 because c_k -> 0. These are ratios of expectations, n first, then k, as the claim says.
5. Tail: representation of the limiting degree law through a newborn of exponential age, pi(D >= k) = E_beta[e^(-mu T_k); T_k < infinity] ~ A k^-mu with A = E[(W^deg)^mu], 0 < A < infinity, and the crossing-flux identity mu pi(D >= k+1) = pi[(k - c_o + rho_o) 1{D_o = k}] to turn the cumulative tail into the pointwise p_k ~ mu A k^(-1-mu). Not checked line by line by the lead.

Numbers (not part of the proof): the lead's own simulation gives c_2 = 0.725(2), the degree-tail constant A = 0.7335(5), and the ratio P(c = k+1)/P(c = k) -> 0.639(2) against 1/(1 + mu) = 0.6387 (FACTS F20, F23 on the board); the child's simulation gives c_2 = 0.7273(40).

## 2. Differences and remarks

* No estimate of the claimant's proof remained underived after the comparison (nine groups, each with a derivation in `comparison-claim3.md`).
* One slip in the claimant's text, found by the comparison child and verified by the lead against `proofs/upper-bound-banded.md` of the round-1 repository: `operator.md`, section 4, last paragraph, "the directed arc (w,u) carries r_w g(w,u) from row u to row w". In the convention of the round-1 certificate, g(D,c;d,k) has the focal vertex first, and a gain at w moves r_w g(D_u,c_u;D_w,c_w) from row u to row w; so the sentence should read r_w g(u,w). The table is not symmetric (the child exhibits the entry (2,1;2,2) = -12511645000 against 0 for its reverse). The sentence is an auxiliary remark on how the round-1 degree-band transfers integrate to zero under the eigenmeasure; no statement of claims 3 or 4 depends on the order of the arguments.
* One step stated without proof in the claimant's text (`operator.md` before section 5: "It applies with a locally finite exterior as well") was supplied by the comparison child (finite completion with uniform conditional bounds).
* The blind proof uses the capped lemma S6.2 for one spectral bound of the spatial semigroup; S6.2 does not depend on claim 3, so there is no circle between claims 3 and 6.

## 3. Not verified

The existence and identification of the eigenmeasure pi of the spatial semigroup on infinite rooted trees (tightness, the tilt, the mass-transport argument) and the uniform integrability behind 0 < A < infinity were read by GPT-family models only. The round-1 results are inputs.
