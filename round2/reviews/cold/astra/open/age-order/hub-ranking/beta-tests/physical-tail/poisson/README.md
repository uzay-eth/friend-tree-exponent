# Poisson serial floor: one eighth and symbolic cutoff 9/2

Using the root's Poisson-stock lower lemma, the exact physical ledger gives strict conditional order on every old tree for `b≥20` and `1/8≤x<1`. The paired Fable splice from the [cold band review](../band/README.md) gives strict annealed order for all adjacent labels from `2→3`; founders tie. Independently of the finite splice, the symbolic cutoff is `(b+1)x≥9/2` for `b≥20`.

## Symbolic prices

At serial completion the tracked middle stock starts at one. The supplied lemma gives `e^(-1)(1−x)^b`; use `M=36787/100000`, strictly below the odd order-nine Taylor sum for `e^(-1)`. The pinned scaled ledger with unchanged locked prices already passes at `(b+1)x≥5`, hence at `x≥1/4` for `b≥20`.

For the sharper cutoff use `m=b+1`, `r=1/k`, `s=1−r`, `B=m−r`, `q=1−x` and `z=q/(xB)`. The normalized locked Beta tail is bounded by `z^s/Γ(r)` times `E[(1+(q/x)V)^(-s)]` for `V~Beta(1,B)`. The prefactor follows by Gamma recurrence and Jensen as in the scaled proof. If `mx≥Y=9/2`, then `z≤1/Y` by `mx≥r`. The degree-two Taylor polynomial is an upper bound on `(1+t)^(-s)` for every `t≥0`, by the negative third derivative. Using the first two moments of `V` gives the factor

```
1 − s z B/(B+1) · [1 − (1+s) z B/(B+2)]
  ≤ 1 − s z [20/21 − (1+s)/Y].
```

Here `B≥20`. The product with `z^s` increases on `[0,1/Y]`: with `c=20/21−(1+s)/Y`, its derivative has sign `1−(s+1)c z>0`. For `k<32`, evaluate at `1/Y` with a downward rational root and the pinned lower Gamma quadrature. For `k≥32`, retain the coarser original `Y0=4` envelope `2/[(k−1)Y0]`, which also controls the unbounded-degree tail of the skeleton ledger. The exact updated scaled ledger (including large-degree sparse/rich profiles) passes all 735 endpoint types with the inherited fees and reserve. Its checker is pinned at `85e9ce571c1b6a842e6d3aa8faa71dbda147b2df`; the frozen band and endpoint checkers keep their prior pins.

## Finite splice and limit of this sweep

For each label `20..34`, the replay checks closed boxes against the stock-sensitive physical continuation, clipping the last box at `9/[2(b+1)]`. The cuts before clipping are `1/8, 77/600, 2/15, 1/7, 1/6, 1/4`. In particular label 20 starts at `1/8` and meets the symbolic bound at `3/14`. For `b≥35`, already `36/8=9/2`. Early labels use the paired Fable receipts only.

At `b=20,x=1/10`, the unmodified ledger with the Poisson floor fails even with a singleton threshold box: the retained profile `(d,c,chain_ports,receive)=(2,1,0,0)` fails. At `b=20`, a narrow box starting at `1/9` failed a one-site chain check. Neither failure is evidence against the physical tail order. A different debit or stronger local continuation input is needed for this method to reach `1/10`.

## Replay

From the Astra root, let `P=round2/open/age-order/hub-ranking/beta-tests/physical-tail/poisson/check.py`. Use `nice -n 10 timeout 115s python3 -IB "$P" --stop 20` for the symbolic rows and `--finite-only --start B --stop E` for bounded finite replays (`--stop` is exclusive). Both symbolic rows passed. Finite label PASS lines were observed in separate invocations for `20..22`, `23..24`, `25..26`, and `27..34` (the first two invocations printed only those prefixes before stopping).
