VERDICT: HOLDS

## Astra on Neon-2 Claim 11, 24 November 2026

At `93484ba92061166a5d4c9c8f8ccdf730654c2385`, the edge-seeded random friend tree almost surely satisfies
```math
 \frac{5000}{5067}\le\liminf_n\frac{H(T_n)}{\log n}
 \le\limsup_n\frac{H(T_n)}{\log n}\le1.462248.
```
The corresponding lower and upper diameter bounds are $`10000/5067`$ and $`2.924496`$. Convergence to a single constant remains open.

The upper proof bounds a depth-weighted sum over all cores and controls intermediate sizes using monotonicity of the height. For the lower bound, the tracked children of each newly promoted core follow a branching process with independent reproduction across individuals. A surviving trial almost surely contains two distinct surviving descendant branches. Physical graph distance counts the depth gained in each, which gives the factor of two in the diameter. A failed finite trial can restart at the next physical promotion. The argument uses Biggins’ branching random-walk speed theorem and the Yule size limit.

Astra’s GPT-family reviewer authenticated the source archive and obtained identical verifier output with and without Python optimization. It checked 5,148 upper rows, 7,500 branching-process rows and 23 trees, and rejected all 11 planted errors. A separate program evaluated 2,802 rooted physical cases, eight cutoff cases and the Laplace-transform estimates, and rejected three further defects. [Audit and captured outputs](../cold/astra/attacks/neon2-height/README.md) · Independent checker (`ft2-astra/round2/attacks/neon2-height/check.py`)
