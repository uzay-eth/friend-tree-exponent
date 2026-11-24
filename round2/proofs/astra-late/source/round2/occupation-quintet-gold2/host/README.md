> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical-floor and frontier audits and scoped Claude conditional reduction and arithmetic; original Neon-2 Gold2 witnesses and the four reward floors. See the [release index](../../../../README.md#claim-28) for exact scope and citations.

# Gold2 stationary host and census

**Both witness readers pass.** The D24 sided joint replay, with the Gold2 conditioned `f` floor and the disjoint edge floors listed below, gives the outward enclosure **`[0.565714160532366, 0.565736645246524]`**. This profile's lower trial at `99/175` has a certified margin lower bound of about **−170,037.73 ticks/core** and does not by itself certify the target. At fixed occupation model, an additional disjoint floor strictly above that deficit would suffice. The separate H3 pricing is not included in this host profile; the upper-side trial is a control.

## API and provenance

```python
from reader import open_host

with open_host('lower') as w:             # optionally work=another private Gold2 copy
    a, b = w.row_ticks(D, c)             # original inward singleton A+B*rho, no cap
    a, b = w.slack_ticks(D, c)           # same plus original cap_f
    finite, argmin, tail, floor = w.stock_floor(c)
    lo, hi = w.comparison_ticks(D, c)    # recharged lower; unrecharged upper
    higher = w.higher_ticks(D, c)        # pointwise sum of original higher cap slots
    cells = w.cells(24, 7)               # (D,k=children) -> (lo,hi,higher)
```

`c` counts all core neighbours, including the mother (`k=c−1` for a nonfounder). Local values are ticks (`S=2^40`). Lower/upper comparisons apply to different allocations because only the lower comparison includes `alpha*(c−2)`; do not interpret them as a pointwise interval of one recharged quantity. Summing the lower host comparison is conservative, since on a finite core tree `Σ(c−2)=−2`, and the solitary core is handled unshifted. `D=limit+1` is the LP label for **all** `D>limit`; at its final arity cell use **all** `k≥saturation`. Exact rows and local lookups are available at the actual physical degree, including arbitrary analytic degree. A tail floor on singleton slack uses the original analytic cap, not the minimum of sampled degrees.

`open_host` requires both completed hash-receipt lines, the detached Gold2 source commit `b435b4c824e2ba5835f58b85e67a27de7d6a0d49`, and the published receipt SHA-256 `5e11f9c0b6eb3bd43e880be608c37d6f0b1fdbcc632b91bb5e8037db52e13284`. It freshly authenticates the entire selected binary before mapping or constructing the layout; thereafter the layout, sparse reads and mapping use the authenticated open inode, checked again at exit. Wire parsing reuses pinned generic labelled-layout code only. No Gold1 value, receipt, discount or credit is loaded.

The local potential comparisons reduce the authenticated G2 triangular edge and higher **weight** blocks, not residual/envelope tables. H3 uses its middle, S3/S4 their hubs, P5 its center, B5 its hub, and each P4 contributes exactly half to each internal. Partner factors count actual allowed internal roles, not endpoint table indices or distinct census labels. Comparison extrema are reconstructed independently for all `D<4000`; both signs match the published finite extrema exactly. Higher cap cells use the published original family-wide maxima with the appropriate weighted multipliers and a valid global pointwise cap, so they can be looser than the verifier's unpublished per-block residual maxima. They are not a second claim on higher-family residuals. Analytic overflow uses separate singleton and convex G2 interpolation envelopes including infinity and all clamped arity classes; exact local higher support is absent once `D>=4000`.

## Durable composition

For a signed inward target shift `a=d*(target−lambda)`, the balance gives `d*(mu−lambda) M >= B−E+π[H−h]`; `B` is the sum of one original singleton floor and disjoint original `g` rebates, all in ticks. The joint objective prices `H−h+a*(M_max−u)` for `a>=0` (the `a<0` case uses `|a|*(l−M_min)`). The LP works in `π/μ`; the exact dual reprices every column and deducts a per-family repair before multiplying a nonnegative result by the original Gold2 lower `μ` bound. Separate/priced versions of the same inequality are compared, not summed.

| D24 quantities (ticks/core) | Lower | Upper |
|---|---:|---:|
| Denominator lower | 1,071,507,709,763.02 | 1,095,510,289,929.18 |
| Denominator upper | 1,359,401,542,039.25 | 1,360,177,372,416.95 |
| Conditioned `f` floor | 207,990.215716345 | 296,272.606830282 |
| Original higher cap savings (floor) | 16,626.091785895 | 15,892.212580932 |
| Disjoint edge floor included | 447,837.133169383 | 0 |

The lower edge total uses SHA-bound `edge18.log`, `absorb12.log` (supersedes `mixed20`), and [`payer/README.md`](../payer/README.md). D18 takes both paying-edge endpoints `≤18`; each FY takes paying-edge maximum `>18`; the absorber pays payer `≤5`, whereas the wider payer pays payer `>5` with donor `≤5`. The wider payer's exact normalized receipt is `1046834122639903380009/21990232555520000000000000000`. The replay multiplies by `2^40` once and asserts the exact tick value `1046834122639903380009/20000000000000000` (about `52,341.70613199517`); this corrects a factor-of-1,000 underestimate in the previous export. All edge receipt hashes are embedded in the conditioned export. The reader also permits `--edge-floor` for an independently supplied exact ticks/core floor. The conditioned receipt is SHA-1 Git-blob pinned; use its directed decimals as exact lower inputs. `--reduced-f` instead derives its own stopped newborn plus same-receive residual. It is weaker for this target and **replaces** the conditioned `f` path. No upper-side edge floor was supplied in these runs.

- [conditioned24.json](conditioned24.json): exact root inputs, signed trials, denominator and target controls, and all 307 cells per witness; SHA-256 `3f98385e0662b4a5aea71dcfef03cdbd66012de01f0083456e243592007de790`.
- [joint24.json](joint24.json): exact alternate stopped/reduced singleton controls, no edges; SHA-256 `104ad923d984e68a185ac5ca8678a2305bb9aa31fb510470a5fea27df514c99a`.

## Reproduce

From Astra root, use Python with NumPy and SciPy (the existing `.work/venv` here):

```sh
P=round2/occupation-quintet-gold2/host
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 5 nice -n 19 timeout 155s \
  /usr/bin/time -f 'wall=%e rss=%M exit=%x' round2/.work/venv/bin/python -Bu "$P/check.py" \
  --conditioned-f --completed-lower-edges --out "$P/conditioned24.json"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 5 nice -n 19 timeout 210s \
  /usr/bin/time -f 'wall=%e rss=%M exit=%x' round2/.work/venv/bin/python -Bu "$P/check.py" \
  --reduced-f --out "$P/joint24.json"
```

Both reported `PASS exact Gold2 host, overflow and joint inputs`; conditioned `wall=30.09 rss=227528 exit=0`; reduced `wall=54.16 rss=275876 exit=0` (its JSON digest was unchanged). Each sign checks 358 labelled host allocations and reduces 759,586,848 potential words. The sided LP has 601 rows and 5,542 columns. The runs did not invoke the full verifier or scan verifier envelope tables.
