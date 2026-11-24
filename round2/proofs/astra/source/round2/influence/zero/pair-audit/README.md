# Independent audit of the six-count pair bound

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 14: Claude for transfer/scalar identities; GPT for full contraction and sign transfer; inputs include Claim 12 and the cap-48 continuation certificate. See the [Astra ledger](../../../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


`audit.py` starts with labelled planted trees. It enumerates rate-one source clocks, synchronizes each common directed edge at the minimum of its two marginals, then couples residual selections by the normalized product of the two residual laws. The root's inert slot is a selection with no growth. Projection is checked only after both selections from a clock have been applied. `corners.py` separately checks the resulting relaxed Bellman equations with exact rational arithmetic and upward integer ticks.

## Why the six-count projection gives an upper bound

Starting from a private A root leaf, neighbour sets of each common vertex stay nested: A contains B at even depth and B contains A at odd depth. At a common source the dominant side has the smaller per-edge rate. Its residual selection targets a private neighbour; the other selection is a common neighbour or the inert slot. The resulting unpaired append at a common vertex is therefore on the dominant side at the target's depth. A private source can hit a common target only along its attachment edge; it appends on the dominant side there. Synchronous appends at common targets add the same new label to both trees. This induction also allows private branches to grow without imposing a degree bound on them.

At the root, the residual for a particular A-private, B-common (including inert) selection pair has rate `1/(D_A D_B)`. The A and B projections may both change at this one clock. A common leaf sends at rate one on both sides. A common A-leaf/B-core sends in common at `p_B<=1/2` and appends privately in A at `1-p_B`. A clean common core sends only in common, at rate at most `1/2`. A dirty common core sends in common at `p_B` and privately in A at `p_A-p_B`; their sum is at most `1/2`. A private A core sends to the root at rate at most `1/2`, and a private leaf does so at rate one. All other sources either leave the six counts unchanged or grow inside an already dirty branch. The Bellman endpoint products cover the convex sets of these projected rates without assuming that all endpoints are physically attainable simultaneously.

Match clean common child shapes before taking the TV norm. The instantaneous bound is `c(p_B-p_A)+e(p_A+p_B)+v p_B+k p_A`, with reciprocals here taken at the root. Labelled dirt is permanent even when the two unlabelled rooted shapes subsequently coincide. The bound therefore can overcount both the instantaneous TV and cancellations across ages.

## Terminal inequalities, valid at arbitrary degree

For one tree with root degree `d` and `j` actual core children, set

```
F(d,j) = (j+1/s)/(1+s*d).
```

An increase of degree at fixed `j` can only lower `F`. Promotion raises `j` by one at rate `(d-1-j)/d` (the inert slot is excluded). Thus `(s-L)F-j/d >= 1/[d(1+s*d)]`. The bounded discounted generator comparison gives the one-side occupation bound `min(1/s,F)`. It strictly improves the old arity expression `(j+1/s)/(s*d)` before clipping. Two one-side occupation bounds can be added at a stopping boundary because the strong Markov property applies to the physical pair; neither side needs to be conditioned on its eventual future path.

Alternatively, let `m=v+e+l+k` and let `d0` be the starting B degree. The instantaneous bound is at most `min(2,2m/d0)`: the clean-core term costs at most `(l+k)/D_B`, and the remaining three classes are covered by the available `2m/D_B`. A residual root clock raises `m` by at most one, at total rate `(l+k)/D_A<=m/d0`. Each counted root child can otherwise generate at most one new counted child at rate at most one. Hence `E[m(t)]<=m(0) exp((1+1/d0)t)`, by localization and Grönwall. Integrate `min(2,2m(0) exp((1+1/d0)t)/d0)` against `exp(-st)`; the exact integral decreases in `d0` for fixed `m(0)`. Left discount and right toll endpoints give an upper Riemann sum, followed by the exact constant-toll tail. The one-side and counting bounds are independently valid, so their minimum is valid.

## Short exact checks

- The mixed fixture has all six classes and nine active projected successors. Splitting a simultaneous root selection into two separate clocks or omitting a private core's send fails the exact projected row check, even though splitting preserves the one-tree marginals.
- A held-out pair created by a residual root clock contains a deeper common subtree. At uniformization jump ages 0, 1, 2 the labelled support sizes are 1, 12, 157. Signed TV after aggregating by shape within each step is respectively `7/6`, `301/270`, `516853/486000`; the corresponding expected six-count bounds are `7/6`, `311/270`, `552803/486000`.
- Cap-six cross-replay checks all 14 reported seed states. An earlier producer initialized Bellman iteration at the fixed-only ratio, which need not be feasible when an LC source has compulsory outgoing rate. The replay found a maximum excess of `0.0071372136`, at seed B degree 2 with one common leaf. The revised producer starts at zero and all 14 seed bounds now match exactly.
- On three single-dirty held-out six-count states `(0,0,0,1,1,0)`, `(1,0,0,1,1,0)`, `(0,1,0,1,1,0)` in `(u,v,c,e,l,k)` order, the revised single-side boundary lowers the rational cap-six ceilings to `1.5791554986`, `1.2382110110`, `1.7726131257` from `1.6249025495`, `1.3794089541`, `1.8594863789`. These are statewise upper comparisons, not a uniform contraction on arbitrary post-dirty states or a sign for the full susceptibility.

Run from the repository root, with `PAIR_EXE` naming the compiled pair executable:

```
timeout -k 1s 5s taskset -c 8 python3 -Bu round2/influence/zero/pair-audit/audit.py
timeout -k 1s 8s taskset -c 8 python3 -Bu round2/influence/zero/pair-audit/corners.py --exe "$PAIR_EXE"
```
