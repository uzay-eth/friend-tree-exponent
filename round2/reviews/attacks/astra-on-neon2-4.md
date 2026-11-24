# Astra cold attack on Neon-2 Claim 4

**VERDICT: HOLDS — response transcendence.** Pin `ft2-neon2@dd711c233c59f86d9c7249a84b47fa1526ed7205`, `round2/neon2/obstruction/REPORT.md` §B (SHA-256 `49debc6252ac6e2a39b51b6cbb8621c592ddad17c481614c8c09515cd6c8d20a`). The affirmative verdict covers `U_s(J)`, the physical trigger `F(s)`, and all 361 permanent-mark entries `K_ij(s)` separately over `C(s)`, using their stated free/marked resolvent identities and positive semigroup bounds. The all-diameter eigenfunction obstruction in §A remains open as the submission states. No conclusion about the arithmetic of `μ` follows.

| Response | First `p`-singular order in `z=1/s` | Decisive term |
|---|---|---|
| `U_s(J)`, `p>=5` | `z^(p-1)`, valuation −1 | Terminal top-star promotion; feedback and lower-order Möbius terms cannot cancel it |
| `F(s)`, `p>=5` | `z^p`, valuation −1 | One more action: independent direct production from the top star; live reciprocal remains uncapped |
| `K_i,20`, `p>=23` | `z^(p-1)`, valuation −1 | Terminal top-star reward |
| `K_ij`, `3<=j<=19`, `p>=23` | `z^(p+j-2)`, valuation `−(j−1)` | Unique highest-pole route: promotion, `j−3` parent sends, shifted copy |
| `K_i,2`, `p>=23` | `z^(p+1)`, valuation −2 | Promotion then unshifted copy; the negative `−1/p` part survives in `1/20−1/p` |

The separate gain counts of the first degree-`p` vertex and the final small-degree top exclude earlier singular actions even on lineages where the degree-`p` vertex is not the top. The asymptotic expansion only needs finite path truncations and a bounded continuation at positive infinity; it need not converge as a power series. Any analytic algebraic relation would induce a formal one, descend to rational coefficients, and confine denominator primes to a fixed finite set. This contradicts the nonzero residues, without claiming meromorphic poles in the analytic variable.

**Evidence:** `ft2-astra/round2/attacks/neon2-transcendence/{README.md,replay.py}`. The independent standard-library replay passed eight exact coefficient checks, first poles at 5 and 7, and 741 large-prime leading entry residues (all 361 at both 23 and 29). Planted missing direct birth, wrong copy sign, missing shift, and capped live reciprocal were rejected. Large-prime checks are leading-residue replays, not exhaustive full coefficient censuses. No author files were edited; no contact or push.
