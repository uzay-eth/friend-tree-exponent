# Erratum and signed transpose tail audit

The direct residence ledger here does not import signed surrogate profiles.
Earlier shifted terminal intervals that used those profiles are withdrawn
pending corrected full replay; see the adjoint README (`README.md`).

The previous chat subtotal counted the two- and three-core blocks twice.
The source intervals below are checked against the producer's
`checkpoint-ticks.tsv` (`physical-corrector/checkpoint-ticks.tsv`) and
`TAIL.md` (`physical-corrector/TAIL.md`). All numbers in the table are
integer ticks of `10^-12`.

| Disjoint source | Upper snapshot | Lower snapshot |
|---|---:|---:|
| Size at most 15 | [151379389,151379390] | [-132517199,-132517198] |
| At most two cores, size above 15 | [590972,590973] | [-465200,-465199] |
| Three cores, size above 15 | [-598504,-598503] | [1117746,1117747] |
| Correct subtotal | **[151371857,151371860]** | **[-131864653,-131864650]** |

The previously reported `[151364325,151364330]` and
`[-131212107,-131212102]` are withdrawn. The source owner transcribed the
checkpoint from the original tool output; there is no saved stdout for
those runs. This audit checks the recorded outward integers, not a fresh
run of the size-15 and large core-block jobs.

## Signed tail occupation by local motif

Let `α=ev_J+aΔ` for the chosen finite dual `a`. To price a disjoint mask
`B`, transport the signed discounted mass from `α`, merge on true
canonical destinations, and retain the signed residence row `v` only on
`B`. With `G=b−(s−L⁰−CΔ)Q`, form

```
v G = v b + κ Q,    κ = −s v + v L⁰ + v CΔ.
```

The mask is applied to residence before taking the transpose. In
particular the destination of an arrow or child derivative can lie
outside the mask. Pool the coefficients of each signed local motif
(`family`, oriented label tuple) before taking endpoints. Singleton,
smooth shift, interpolation seam and finite patch contributions are kept
in the same identity. The code checks exact agreement against direct
treewise `vG`, with positive and negative payments recorded separately
after forming the actual local residual. Taking absolute values before
forming `G` loses the principal cancellation: even after pooling equal
motifs the absolute ledger on this sample is about `.662`, whereas the
sum of the absolute net-tree payments is below `7.6e-6` on both snapshots.

```
taskset -c 13 timeout 8s python3 -Bu round2/influence/zero/adjoint/motif_tail.py --subtotal
taskset -c 13 timeout 9s python3 -Bu round2/influence/zero/adjoint/motif_tail.py --dual 6 --prefix 6 --last 7
taskset -c 13 timeout 9s python3 -Bu round2/influence/zero/adjoint/motif_tail.py --side lower --dual 6 --prefix 6 --last 7
```

The exact 48-tree size-7 band yields `[.000003084338,.000003084339]`
and `[-.000002809230,-.000002809229]`, respectively. These are validation
bands, not additions to the size-15 checkpoint. The same ordinary
frontier's uniform unit-source multiplier is bounded by `.708921094645`.
For scale, even ignoring the derivative correction, a blanket source
bound of `1e-3` at this multiplier would cost `7.09e-4`, several times
the checkpoint subtotal. A sharp ordinary-tail enclosure needs an
envelope for the *combined* residual or a joint occupation bound that
preserves the transfers between motifs; the individual motif norms are
not useful tail prices.

No sign claim follows from the finite bands or subtotal. The two
unresolved terms are discounted ordinary residual occupation beyond the
priced masks and the signed derivative mismatch paired with `ΔW`.
