# Smoothed signed scalar: physical and backward-law audit

Write `D=1+L+k`, `p=1/D`, `P=Σ_i 1/D_i`, and `b=L/D`. The summation
for `P` is over every actual core-child edge with multiplicity. Each
child uses its own full degree including the slot to this root. There is
no arrow from the host's inert parent. The `+1` in `D` and in `d=k+1`
counts this silent slot; it is not part of `P`.

For `g=V_(s+1)h`, intrinsic evolution of old children within `L⁰(Cg)`
cancels `C L⁰g`. Root gains occur at rate `L+P`; the children themselves
are unchanged, so the root-degree increment contributes `−(L+P)Cg/(D+1)`.
Each root-to-child send has rate `p` and changes just that old child,
giving `p C(g⁺−g)`. A promotion has rate `b` and adds the virgin child
with score `p g(J)`. Thus

```
B g = C(1−L⁰)g + L⁰(Cg)
    = η Cg + p Cg⁺ + bp g(J),
η = 1 − p − (L+P)/(D+1) = (k+1−P−p)/(D+1),
(s−L⁰)(V_s C h − Cg) = B g.
```

`P≤k/2` and `p≤1/2`, so `η≥0`; also `η≤1−p`. All three pieces of `B`
are positive and `B1=η kp+p(1−p)≤1−p≤1`. For bounded `h`, `g`, `Cg`,
and `Bg` are bounded; the displayed identity gives a bounded generator
source. Stopping discounted Dynkin and using nonexplosion gives the full
resolvent identity without assuming bounded total jump rates. Signed
starting measures follow by linearity. The exact checker expands host
and child arrows independently of the census formula and compares all
generator coefficients (not just one test vector).

## Backward-selected law and restart

This is an equality of discounted **endpoint measures**. For nonnegative
host source, weight its host trajectory segments by discounted occupation
and choose an arm with its actual factor in one of the three `B` channels.
The immediate `C` channel is a separate atom. After drawing the chosen
canonical shape, apply a fresh *intrinsic* trajectory with no incoming
arrow through its slot, for an independent `Exp(s+1)` duration. Its
expected endpoint reward times `1/(s+1)` is `g`. For the pre channel
start unchanged; for the post channel insert one leaf first; for
immigration start at `J`. In particular the `J` immediate channel is
zero, and the integrated `B` mass from the host `J` is exactly
`(s+1)V_s C1(J)`. For each side of a signed start the same construction
holds before subtraction.

This construction does not identify an interval of an already-conditioned
physical path with the independent exponential segment. In an actual
backward selection an outgoing arrow weights the previous histories by
their hazards and the no-earlier-change survival; a child attached at the
host observation time already has an age and a selected shape. For example,
two equiprobable histories with upward hazards `1/2`, `1/3`, competing at
total rate `3`, have conditional first-upward probabilities `10/17,7/17`
over the histories (before discount). Resetting these to `1/2,1/2` is wrong.
On an unconditioned forward strong-Markov restart retain the complete
physical state; for a backward gap carry the tilted endpoint law or prove
a new conditional kernel. In particular arbitrary restart with the prior
or an independent exponential **physical** gap is not justified.

The auxiliary construction is also discrete: a host jump discounted at
`s` has categorical mass `r/(s+λ)` and cemetery mass `s/(s+λ)`; integrate
the host state reward with residence factor `1/(s+λ)`. After the port
selection the new chain uses the independent `s+1` cemetery kernel.
If the host exits a finite domain, its remaining kernel is conditioned on
the actual exit shape. Time weights factor by the Markov and memoryless
properties, not by independence of the shapes at the two stages.

## Both finite seams and a matched dual

Let `v=αV_N` and `m` be the host's signed first-exit measure. The exact
factorization at a finite cut is

```
α V_s C h = α Cg + v Bg + m(V_s C h−Cg).
```

The subtraction in the seam is essential. For the scalar `U=V_s(b+CH)`
and `H=ΔU`, let `D(a)=α C+v B−m C` be the net smoothed coefficient,
where `α=ev_J+aΔ`. Then `U(J)=v b+m U+D(a)g−aH`. Use a freely chosen
finite dual `t` in smoothed coordinates: write `t g=a H+e g` with
`a=t V_(s+1),N` and the silent first-exit measure `e`. The identity is

```
U(J) = v b + m U + [D(a)+e−t] g.
```

Both seams are carried by the exact transpose. Iterate `t←D(a)` on the
silent start domain as a proposal; final arithmetic checks the identity
with a manufactured nonconstant signed forcing and nonzero boundary.
A small dual mismatch is observed, not assumed from a convergence theorem.
With a fixed **unmatched** `H` dual, separately pricing `g` and `H` can
instead worsen the interval by dropping their dependence.

## Short CPU-13 results

```
taskset -c 13 timeout 9s python3 -Bu round2/influence/zero/smoothed/check/check.py --cutoff 8 --extra 1 --scalar
taskset -c 13 timeout 9s python3 -Bu round2/influence/zero/smoothed/check/check.py --cutoff 6 --extra 1 --scalar --matched 5 --profile ft-zero-envelope-v2-s56572365-cap48-tight.txt
```

The first run passes 199 exact physical rows and rejects omission of the
`−p` in the numerator of `η` on 192 rows. The four unmatched signed `g`
coefficient masses (immediate, pre, post, birth) are approximately
`−.00297311,.04729972,.05661190,.19604722`. Re-expanding `g` in this
unmatched dual leaves interior mismatch TV `.036557910907` versus the
original killed mismatch `.000032086163`; no gain follows from the algebra
alone.

The second run uses the corrected v2 profile and its `g,U` callbacks. It
gives matched mismatch TV `≤.000000298685`, `g` width `≤.000000186908`,
host `U` width `≤.402669729605`, and full scalar
`[.792481469901,1.195151386414]`. The low queries after matching are:

| Seed | Exact outward `g` enclosure | Signed mismatch coefficient (approx.) |
|---|---|---:|
| `J` | `[-.071315953811,.444134559418]` | `−1.13654e−7` |
| `J⁺` | `[-.111481487458,.246665765915]` | `−4.3963e−8` |
| Bare root with `J` child | `[-.425627730072,1.119488346377]` | `−8.304e−9` |
| One leaf and `J` child | `[-.301260950085,.564459405148]` | `−8.107e−9` |
| Bare root with `J⁺` child | `[-.608325817734,1.399972208696]` | `−2.374e−9` |

This prices the matched smoothed uncertainty directly. It does not price
away the host seam, and it supplies no exponent sign.
