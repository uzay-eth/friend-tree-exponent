# Neon D18 proof source

Neon, *Rigorous Certified Enclosure of μ via Unified Bernstein-3 Full-Tail Level-4 Motif Potentials*, `round2/proof-b3f-motif4.md`, source commit `cec3855cc52ab4265ce3e36b8a7a38f1ae257f77`, Git blob `3bbda54c3cbdbd235a6cf4af8ecf80ca75d16ae7`. The excerpts below preserve the transfer and comparison statements used here. To retrieve the full proof from a Git repository containing the pinned object:

```sh
git -C "$NEON_SOURCE" show cec3855cc52ab4265ce3e36b8a7a38f1ae257f77:round2/proof-b3f-motif4.md
```

## Motif definitions and transfers

Sections 2–4 define the eight-term potential `f,g,hV,hC,P4,Yt,Yb,S3`. The cubic Bernstein degree embedding reproduces the physical reciprocal degree exactly. A core gains a leaf at rate `L+ρ`, where `ρ` sums the reciprocal degrees of its core neighbours, or promotes one of its leaves at rate `L/D`. Section 3 transfers the centered receive terms from `f` to `g`, from `g` to `hV,hC`, and from those three-vertex motifs to the eligible four-vertex motifs. Each transfer adds and subtracts the same term, so the transfers cancel in the sum over a tree. Section 4 encloses the tail rows by rational Bernstein control bounds, including overflow and capped states.

> Because every transfer pair adds $`+X`$ to one motif and $`-X`$ to its neighbour motif, **the sum of all transferred terms over any tree $`T`$ is identically zero**.

## Core comparison bounds and code map

> Every motif in $`\Phi(T)`$ is charged to a unique core vertex $`v \in K(T)`$:
> - `f(v)` to $`v`$ ($`1`$ term);
> - `g(p, v)`, `hV(p, u, v)`, `qP4(a, b, c, v)` to bottom vertex $`v`$ (at most $`1`$ each, only if $`h_v = 1`$);
> - `qYb(u; w1, w2 -> v)` to bottom vertex $`v`$ (at most $`1`$ group of $`k_u - 1`$ siblings sharing $`(u, w_2)`$, only if $`h_v = 1`$);
> - `hC(v; w1, w2)` to apex $`v`$ ($`\binom{k_v}{2}`$ pairs when $`2 \le k_v \le K_T`$; $`0`$ otherwise);
> - `qYt(p, v; w1, w2)` to apex $`v`$ ($`\binom{k_v}{2}`$ pairs when $`v \in \mathcal{D}_{YT}`$);
> - `qS3(v; w1, w2, w3)` to apex $`v`$ ($`\binom{k_v}{3}`$ triples when $`v \in \mathcal{D}_{S3}`$).

Maximizing the grouped charge over the one-body states and overflow states gives the unrecharged `C0` ceiling. The separate lower comparison recharges the vertex by `(k−h)/4`, which sums to zero on every core tree. The tracked GMP verifier (`attacks/neon1/round2/cert/cert_b3f_verify_gmp.c`) checks the eight row families and the comparison constants. The tracked GMP receipt archive (`attacks/neon1/outputs/neon1-production-rehearsal.tar.gz`) contains `outputs/production-final/lo.log` and `outputs/production-final/hi.log` for the original D18 witnesses.
