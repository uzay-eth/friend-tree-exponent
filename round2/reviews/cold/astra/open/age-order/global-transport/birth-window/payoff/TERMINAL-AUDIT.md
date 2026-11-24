# Source audit: degree-two direct terminal

The tight retained row can absorb the payoff loss through `y=1/6400000000000000` (about `1.563e-16`) when charged its actual sources. This is **a local row certificate** for a degree-two old core with one old leaf and one direct core neighbour of degree `e≥31`, not an improvement of the all-tree cutoff. The other retained rows and short chains would need separate recalibration at this value.

Write `r=1/e`, `w=(1+r)/3`, `s=w/2`, and `t=1/4`. The two serial terms have weights `w,t` and post-placement old-host degrees three and two. The old-leaf edge has `s<t`, so its adverse direction is the new leaf at the old leaf and its **favourable** endpoint has parent degree three. Its extra charge is at most `(t−s)η_parallel(3)`. The direct core edge is either adverse at the other owner, or adverse here with favourable parent degree `e+1≥32`. The frozen favourable lower there was already zero, so **there is no new direct charge here**. The upper is still the frozen hub cap. The total loss paid by this row is

```
w η_serial(3) + t η_serial(2) + (t−s)η_parallel(3).
```

After subtracting this charge, the frozen reserve is concave in `r` on `[0,1/31]`: the serial terms and the old-leaf term are affine, and the direct cost is bounded with `max(B,s U(3,1))`. The exact endpoint replay gives loss `5.18536020e-6` against reserve `5.47868233e-6` at `r=0`, and loss `5.15374919e-6` against reserve `.00173635132` at `r=1/31`.

## Clock and amplitude audit

At the successful stock exit, apply the conditional Poisson–beta inequality on the actual finite tree. The fresh size-clock amplitude is independent of the *future embedded chain*, conditional on that exit history. The beta–gamma multiplication therefore compares remaining **amplitudes**, not shares rescaled by the random exit size. Chronological discount back to the original completion is `e^-τ`. The size amplitude at the original completion is used only through a marginal upper-tail union bound; it need not be independent of the exit time or its stock. Neither the foregoing local improvement nor the uniform proof transfers an inequality by identifying `S/m` with one.

The upper transform counts successful paths. Its return row retains a unit floor before the cutoff because an extra physical gain can turn a virtual failure into an early success. The new source accounting does not replace that transform with a conditional transform of the virtual clock given success.

Replay: `python3 -IB round2/open/age-order/global-transport/birth-window/payoff/terminal-audit.py`. This local result is separate from 1374’s large-size interior boxes; it supplies no coverage between the small-scale and interior certificates.
