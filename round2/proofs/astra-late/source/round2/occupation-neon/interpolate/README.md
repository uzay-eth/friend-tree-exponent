> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../README.md#claim-21) for exact scope and citations.

# New pointwise drift witnesses from Neon D18

**Repaired pointwise endpoints: `[0.565675148075, 0.565773278615]`.**
These are induced finite-tree witnesses, independent of stationarity.
The [verified composition](compose/README.md) with retained higher-order
occupation charges gives `[0.565676930574, 0.565771634792]`.

**Withdrawal.** The old upper endpoint `0.565770220167` and mixed upper
endpoint `0.565770220267` were invalid. The handwritten upper `k=1`
tail floor was 1,000 times too large. With parents `(-1,0,1,2,1,0,5)`
and leaves `(62,63,0,1,4,65,1)`, the old signed upper drift is strictly
negative. The composition replay reproduces it exactly, checks the same
tree at the repaired trial, and independently re-extracts both sets of
floors. The corrector now reads its row dump and both tail band tables
from hash-pinned sources.

## Pointwise proof

Use `ε=+1` on the lower side and `−1` on the upper side, and write
`δ=ε(s−λ)≥0`. Denote the checked nonnegative first-body row by `F(x,ρ)`.
All seven remaining signed physical motif charges stay nonnegative at
the **original** `λ`. For the `−δΦ` shift pay the pointwise denominator
charge `q(x)+t(k−h)` from the eight-motif envelope, using its exact
recharge cancellation and a raw all-degree cap `M`. Units here and in
the following equation are normalized potential units. They yield

```
ε(ℒ−s)Φ ≥ −δ M C + Σ_nonfounder [F(x,ρ)+δ(M−q(x)−t(k−1))].
```

On nonfounder census `x=(D,k)`, write `A` for the two physical census
changes with rates `L+ρ` and `L/D`. Backward solve for `z` with
`(s−A)z ≤ F+δ(M−q−t(k−1))`, checking the two exact physical-rate corners.
The charged row and the proposed transition use the **same** `ρ`.
The actual signed correction is `V(x)=z(x)−z(J)`. A promotion creating a
core adds no correction from the newborn since `V(J)=0`. The founder's
correction is zero, and the global signed constant `−δ M/s` pays its cap
cost without any founder residual assumption. Thus for the new potential

```
Φnew = Φ + ε (Σ_nonfounder V − δ M/s),
ε(ℒ−s)Φnew ≥ (s z(J)−δM)(C−1) ≥ 0.
```

The all-degree comparison to `C` follows from the existing strict `c0`
and the small uniform correction checked by the script (`<0.01` per
core, including the signed constant). `T3` remains positive. The usual
finite-tree drift product therefore proves the endpoints directly.

For the tail `D>64`, `z` depends only on `k`, so arbitrary gain rates
cancel **exactly**. The tail pays both the certified all-stock row floors
and the uniform charge deficit at each arity. For arity overflow it is
constant, and for the other arities check the two promotion corners
`max(0,1−(k+1)/65)` and `1`. The repaired upper margin is `+8.97e−13`
per core; replacing the tail by zero gives `−6.04e−8`. The falsely large
tail floor had inflated both the margin and the endpoint gain. The
central target remains uncertified by this construction.

## Interpolation obstruction and admissible mixing scale

For uncorrected `aΦ−+bΦ+` with `a,b≥0`, the residual at the common target
is the same combination of `R−=(ℒ−s)Φ−` and `R+=(ℒ−s)Φ+`. Exact physical
generator evaluation, independently assembled in two ways, gives:

| Physical fixture at `.56572365` | `R−/Φ−` | `R+/Φ+` |
|---|---:|---:|
| Six-cell thin train | −4.8461e−5 | −4.4088e−3 |
| Six-cell hub-16 train | +.166880 | −.00122657 |
| All-eight mixed fixture | −1.0933e−5 | −1.1703e−5 |
| Three-core chain, 14,001 leaves per core | −4.0344e−5 | −4.5806e−6 |
| Same chain, 1,000,000 leaves per core | +5.1296e−6 | +4.0901e−5 |

Thus **every nonzero positive mixture** has normalized residual below
`−4.5e−6` on one frozen-tail tree and above `+5.1e−6` on the other.
Tree recharges summing to zero cannot change these examples. The two
signs persist on the bulk homogeneous chains: inserting a core after the
endpoint neighborhoods separate changes the score and generator by a
fixed amount (motif plus one donor range); the replay checks the exact
cell differences at lengths 8, 9 and 10. Consequently changing only a
founder term or a correction supported on bounded leaf stocks away from
both frozen-tail stocks and their one-step images cannot fix either
sign at this target. This is specific to the stated frozen basis; it
does not exclude new tail-active motifs or a different census correction.

The optional mixed witness bounds `|(ℒ−s)Φother|≤K C` pointwise: assign
each motif to a core with bounded multiplicity, take absolute table
extrema, charge promotion changes and newborn copies, and use
`D·||ΔBernstein||₁≤6` on a pure cubic tail. At interior/activation rows
use the finite degree caps; the one-body exception uses `LF+KF+2`.
The conservative ceilings are `50661.783747` for the opposite upper
payload and `50499.006747` for the opposite lower payload. For margin
`m=s z(J)−δM>0`, use exact mixing weight `α=m/(2K)` and extra signed
constant `−αK/s`. This retains a pointwise margin on every tree.
This conservative weight does not explain the gain. Improving it requires local
signed residual cancellation (especially active in both frozen tails),
not a larger absolute all-eight bound.

## Replay

```sh
OPENBLAS_NUM_THREADS=1 taskset -c 8 nice -n 19 round2/.work/venv/bin/python -Bu round2/occupation-neon/interpolate/corrector.py b3f18_lo.bin b3f18_hi.bin
OPENBLAS_NUM_THREADS=1 taskset -c 8 nice -n 19 round2/.work/venv/bin/python -Bu round2/occupation-neon/interpolate/corrector.py b3f18_lo.bin b3f18_hi.bin --mix
OPENBLAS_NUM_THREADS=1 taskset -c 8 nice -n 19 round2/.work/venv/bin/python -Bu round2/occupation-neon/interpolate/extremes.py b3f18_lo.bin b3f18_hi.bin
```

The earlier production GMP sign verifications are dependencies. The
composition replay authenticates the payloads and re-extracts the
pointwise operands as well as checking the residual bounds.
