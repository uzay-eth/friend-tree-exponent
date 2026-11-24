# GOLD2 upper FY: large payer, small donor

The widest completed replay gives a **floor of `10,114.29962382028` ticks/core** and a **matched ceiling of `50,074.33413881668` ticks/core**. They bound the same original GOLD2 upper-witness reward on the stated physical mask. At the fixed D24 comparison, even the ceiling covers only **0.43%** of the remaining sufficient credit for `μ<241/426`. This ceiling does not cover omitted payer degrees or arities.

## Mask and proof

Pay one original outward FY gap for each ordered pair of distinct physical neighbours at an exact G2 host outside M3. The wide mask uses `D_host≤251, c_host≤12`, payer degree `6..35` and arity at most `6`, donor degree `2..5` at all arities, and `max(D_host,D_payer)>18`. These gates apply at reward time. `D_payer>5` separates it from 2f74's upper small-payer credit, and the paying-edge maximum separates it from their degree-18 full edge credit. Payer and donor are roles, not equivalence classes of state labels. No H3 or GOLD1 term is used. Upper outward sign is `+1`; `--mode ceiling` bounds the same upper-witness reward above and does not switch witnesses.

The original 20-knot FY spline and the exact discrete donor reciprocal determine the gap. The replay uses singleton payer bins, so no monotonicity in payer degree is assumed. The rolling Bellman counts mother-payer, mother-donor, and both chronological two-child orders. The older neighbour evolves from its actual registration; births preserve the promoter's degree and increase its arity. Adjacent reciprocals enter its true gain rates, with other receives in their safe boxes. The bounds start at actual core birth. For the floor, use the original GOLD2 upper certified bound for discount and the lower bound for birth intensity; for the ceiling, use the reverse. Rounding is directed and every resulting continuation inequality is checked in integer arithmetic. Monotone degree and arity cuts justify stopping on exit. As in the frozen lower replay, founder-containing paths are excluded.

## Exact receipts

All runs SHA-authenticated the original `quintet_upper.bin` (`f1dcd7ad3e6df82a985eb3dd6fd793b5a3e38c78368d94e4b135ab174ab92876`) and the Gold2 receipt at source `b435b4c824e2ba5835f58b85e67a27de7d6a0d49`. Each pass checked four directed mask/rounding fixtures and 20 raw-knot interpolation and physical cap inequalities, including payer tail and the host's `251→252` seam. Input file identity was checked again at exit.

| Host degree/arity cap | Payer degree/arity cap | Direction | Exact ticks/core | Path states / reward labels | CPU user + sys; wall (s) |
|---|---|---|---:|---:|---:|
| `135/10` | `30/6` | floor | `185646528420795685827/25000000000000000` | `2686320 / 302400` | both passes: `79.98+4.11; 150.08` |
| `135/10` | `30/6` | ceiling | `366581807632754083097/10000000000000000` | `2686320 / 302400` | same replay |
| `251/12` | `35/6` | floor | `1011429962382027974409/100000000000000000` | `7319620 / 1456200` | `123.00+3.68; 186.66` |
| `251/12` | `35/6` | ceiling | `1251858353470416996307/25000000000000000` | `7319620 / 1456200` | `108.30+3.00; 141.46` |

From Astra root, run the commands below with prefix `taskset -c 11 nice -n 19 /usr/bin/time -f 'user=%U sys=%S wall=%e'`:

```sh
python3 -Bu round2/occupation-quintet-gold2/payer/upper/replay.py
python3 -Bu round2/occupation-quintet-gold2/payer/upper/replay.py --host 251 --stock 12 --receiver 35 --rstock 6 --mode floor
python3 -Bu round2/occupation-quintet-gold2/payer/upper/replay.py --host 251 --stock 12 --receiver 35 --rstock 6 --mode ceiling
```

## Scalar comparison at `241/426`

The fixed comparison uses the upper `f_floor`, `saved_plus_reduced`, `cap_error` and upper denominator `M_+` from [`host/conditioned24.json`](../../host/conditioned24.json), and the completed disjoint upper credits in [`edges/upper/replay.log`](../../edges/upper/replay.log) (D18) and [`edges/upper/absorb8.log`](../../edges/upper/absorb8.log) (small payer). These two edge floors sum to `422,385.43386298366` ticks/core. The additional sufficient credit in this separate comparison is

```
(28286843/50000000 − 241/426) M_+ + cap_error − f_floor − saved_plus_reduced − edge18 − small_payer
= 11,745,375.058991114 ticks/core (rounded).
```

The wide floor improves the fixed-comparison endpoint by `7.4360e−9`; replacing it by the matched mask ceiling would improve it by at most `3.6815e−8` at that same denominator. This measures the size of this FY mask against one certified comparison, not an upper bound on what a different occupation model, omitted FY masks, or other row families could prove. No scalar endpoint or release receipt was edited.
