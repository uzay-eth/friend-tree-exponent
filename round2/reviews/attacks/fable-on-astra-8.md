VERDICT: HOLDS

# Fable on Astra claim 8 (a constructed family of finite local drift certificates is complete)

**Scope of the verdict:** The Claude-family lead rederived the two drift-error counts and read the finite grammar and completeness argument, taking Claim 6’s uniform renewal expansion as an input. No table or planted control was run; finite-time estimator measurements are context only. No GPT-family-only child check is used here.

Attacker: Claude Fable 5.2 (a different family from the GPT claimant). The Gemini child stopped on a provider quota; no separate result is used here. The claim is about the constructed rational tables obtained by averaging cut phases, not the original published descriptor family.

## Steps examined (Astra folder `round2/fringe/local-certificates.md`, `e01d4d3`)

1. *Projection and the drift error (1).* Phi_w averages, over the R phases of a periodic horizontal cut, the
   tabulated scores of the finite low components and hub particles. Against the true generator the auxiliary
   one differs only (i) at hubs: missing newborn at rate c/D and sends to core neighbours at rate 1/D each,
   total rate sum_high c/D <= 2C/M, each worth at most A(w) + B(w); (ii) on a suppressed low-low edge: two sends
   of rate at most 1/2, and each edge is suppressed in exactly one of the R phases, worth at most A(w) C / R on
   average. So |L Phi_w - Phi_{Gw}| <= (2[A+B]/M + A/R) C. I re-derived both counts; correct.
2. *Tests (2).* With w > 0 and Phi_w comparable to C, Phi_{Gw - lambda w} + D C <= 0 on every tree gives
   L Phi_w <= lambda Phi_w and hence mu <= lambda by the same product argument as in round 1; symmetrically for the
   lower test. Correct.
3. *"On every tree" as a finite linear programme.* Types record a vertex's own state, its parent's high flag
   and its low downward neighbourhood through R levels; the admissible typed trees form a finite grammar and
   the inequality over all derivations is equivalent to the existence of potentials v_i with
   v_i >= reward + sum of children's v. This is the standard reduction (finite suprema because any unbounded
   symbol could be inserted in a start context). Correct. It is the same kind of certificate as the round-1
   motif tables, with a different descriptor family.
4. *Approximation.* With w = e^{hG} c and M = R = ceil(e^{3h}), the coupling bound (4)-(5) and the uniform expansion
   E_T C(T_h) = e^{mu h} d_2 Psi(T) + O(C e^{sigma h}) (claim 6, which I attacked separately) give strict signs at
   distance O(e^{-(mu-sigma)h} + e^{-(1+mu)h}) from mu. The only input I could not verify independently is the
   uniformity in T of that expansion; given it, I found no gap in the steps I read, including the final rational approximation with strict drift margins.

## Remark on usefulness

With M = R = ceil(e^{3h}), the rate (6) contains e^{-(μ−σ)h}; the proof does not quantify σ and describes the state counts as growing rapidly. These are qualitative observations about table size; no size census or numerical gap estimate was obtained here.

## Not done

No table of the family was rebuilt and no planted-error test was run.
