# Cold review: 34fa high-hub-limit / all arms at D≥101

**CORRECT in scope.** Add the two `all` floors below to `aa511af202d0d64cd1a08f900e0516965d9db952`. The outward enclosure is **`[0.565680279007146, 0.565768284181773]`**.

| Original D18 witness | Exact masked floor/core (book this) | Exact masked ceiling/core |
|---|---:|---:|
| Lower, mode 2 (`+`) | `52327439425129643019077301/100000000000000000000000000000000` | `6077184618479666155541/4000000000000000000000000000` |
| Upper, mode 1 (`−`) | `8091751687185104316481949937/50000000000000000000000000000000000` | `785167337696132124993/1000000000000000000000000000` |

Exact endpoints before outward rounding:

```text
14157243308580160981686611486065459691229739/25026934531690315946179968825000000000000000
44818646818747647275382808904477397643909921/79217319301603159409316399750000000000000000
```

The upper witness is already oriented: its nonnegative reward increases the upper witness's net numerator and **subtracts** from the upper exponent endpoint. In both lanes the increment uses the unchanged physical ceiling (`1357493815758092773/10^18`, `273436004214812983/(2·10^17)`). I checked both the exact endpoint increment and signed net-numerator formulas. Coefficients remain those at the *original* D18 trials (`.565673`, `.565775`); the frozen EDGE pivot trials and the discount enclosure serve different purposes.

## Physical mask and exactly-once registration

An owner is the original allocated `hV`, `P4`, or `Yt` row at a **nonfounder** apex with **current** `D≥101` and **at most** 12 actual core children (`k=12` included; core degree `k+1`). The entire first-child arm has at most seven ordinary vertices, including the child and its ordinary leaves, excluding the apex. Other arms need not be small. Each motif in scope has exactly one first child; `Yt` tips form an unordered pair. Identical arm shapes at distinct physical children count separately. All supported original coefficient terms are included, including centered transfers at their assigned owners, and unsupported occurrences contribute zero.

Let `H` be the continuation for not-yet-born arms and `A_j` the continuation of a live registered arm. On a promotion the host successor is `H(D,k+1)` **plus** `A_new(D,k+1,J)`; every already-registered arm advances its arity without obtaining a second birth reward. A hub degree gain changes `D` by one, preserves `k` and every arm's shape, and does not create a child. Thus the seam at `100→101` carries `H(101,k)+Σ_j A_j(101,k,x_j)` once. This is equivalent to registering every still-eligible pre-existing arm exactly once at first entrance in its actual entrance state, then registering each later arm at birth. `bulk.py` implements the equivalent earlier-registration formulation: zero reward below 101, each old registration already live, no extra seam injection, and `old` future-birth value zero above the seam. A birth while `D=101` is post-entrance; one while `D=100` is pre-entrance. No ages or internal slots are reset. Exits at size 8 or arity 13 are permanent; a promotion from `k=12` pays no newborn continuation and kills existing registrations. If there are already more than 12 children on threshold crossing, no arm can ever enter this mask. Founder-hosted charge is omitted; the stable-age/Palm identity is for a nonfounder newborn `(2,0)`. Direct founder births have at most linear expected count and their bounded masked contribution vanishes under mean-core normalization.

## Rates and arbitrary degree

On a retained arm `x`, write `p=1/D`, `L=D-k-1`, and `u=1/D_mother+Σ_other children 1/D_j`. The gain into the hub, including at degree 100, has **one** live rate `L+1/D_x+u`; the opposite arrow into `x` is `p` and leads to `x⁺`. Promotion is `Lp`. Intrinsic arm moves have their exact full-arm rates, with aggregate hazard `|x|−1/D_x` including size exits; the arm-to-hub arrow belongs to hub gain, not this hazard. Host-only gain includes the mother and **all** current children. The physical mother and sibling states may evolve invisibly in the compression, so the rates can be adapted, not independent. A rate in the holding denominator is the same one that multiplies its transition. Receive pressure in the full charge is the same physical `u` as in the gain; at `D≥101` the envelope bounds the charge for this actual `u∈[0,k/2]`, rather than modifying either arrow. At that seam the test function no longer changes with `D`, so the hub gain contributes zero to its continuation drift, regardless of its rate.

For fixed `k≤12` and full arm, potential weights have cubic Bernstein dependence on `18p`, with just one high slot. An apex degree increment replaces `p` by `p/(1+p)` and its weight difference is `O(p²)`, canceling the apparent pole in `L`. At promotion the only arity change is `k→k+1` at unchanged `D`. Below the apex all old slots and one-step intrinsic successors are inner (degree ≤8); the incoming tail-`g` centered transfer has the original cubic `p²`/`p³` form and the other relevant transfers have degree at most cubic in the high-slot reciprocal. Multiplying by `(1+p)^3` gives a polynomial of degree at most seven, including the affine exterior slope. Hence the reconstruction from eight rational points is an identity, not extrapolation from finite observations. Taking Bernstein extrema for the residual divided by `p` on `[0,1/101]` covers **every** true degree at and above 101 with no artificial upper cutoff. Exterior donor degree can itself be arbitrary. Extra direct checks below include a true degree above `10^71`.

At and above the seam the upper equation tests the full send interval `[0,1/101]`, with the same send in hold and gain. The lower equation kills on a dominating rate-`1/101` proposal, whether or not the physical send is accepted: on the actual path this amounts to replacing a nonnegative accepted successor by zero and adding nonnegative killing on rejected proposals. Both test the physical promotion interval `[1−(k+1)/101,1]`, including the exit. Each host promotion's newborn and future-host values use the same rate. Quotients reduce to corners by cross-multiplication; directed grid rounding preserves the inequalities. Nonnegative occupation permits discount `HI` for the floor and `LO` for the ceiling, then stable-age normalization with `LO` and `HI`, respectively. This uses the inherited physical stationary identity and original two-sign GMP certificates; neither is replaced by finite-case controls.

## Ledger ownership and the complement

| Earlier piece | Why high101 misses it |
|---|---|
| Original higher reserves (size ≤12; size 13–20 with 3/4 cores) | The **complete** planted apex fringe contains at least `D` ordinary vertices, even with a seven-vertex first arm. |
| `f/g` corrections, EDGE/memory, high-parent `g/hC`, paid `S3` and `Yb` | Different allocated charge owners. For example, an `hV` newborn potential term assigned to `g` is still paid in `g`, not as another live `hV`. |
| Frozen `hV` donor installment | Same row family, but apex is at most 48 for mother/apex components and at most 11 for middle/tip components. Whole-row pricing would overlap donor excess on the same occurrence. |
| Root P4 and wider low complementary P4 | Apex at most 12 and 16 respectively, on current descriptors. |
| Old `Yt` installment in `aa511af` | Same row family, but current apex at most 16 for the four booked donor roles (`--caps 16 5 5 5`). |
| Unpublished high-apex P4, `18≤D≤80` | Same charge family and possibly same directed component, but disjoint at reward time. Can coexist as a separately certified payment on this stated mask. The published degree-80 transport receipt is lower-witness only; the degree-48 signed contrast is not a second independent floor for each lane. Neither is added to the `aa511af+high101` bracket here. |

The `all` interval is for the sum over the **masked** families and arms. On that mask, the `old` and `young` objectives partition actual incidence; their separately optimized bounds can be summed as an alternative enclosure, not added to `all`. The outer sum over the **complement** is a distinct, still unpriced nonnegative remainder in each oriented witness: e.g. other degrees, larger arms and `k≥13`, net of already booked pieces. A masked ceiling does not control it (or the global signed action); dropping this remainder is valid only for the oriented floors. No extrapolated complement ceiling or unconditional all-arity price is asserted.

## Reproduction and pins

Companion: `high-hub-sol.py`. It checks 13 source pins, invokes the frozen `aa511af` verifier with full original payload hashing, checks labelled seam/exit controls, reconstructs both outer endpoints, and re-solves the Bellman inequalities on both signs with a separate recurrence and Bernstein implementation. It shares the pinned original row extraction and arm topology; the separate full-host check uses `Local.charge` and labelled projected arrows. Producer files were read-only.

```sh
P=round2/.venv/bin/python
R=round2/occupation-neon/high-hub-limit
L=b3f18_lo.bin U=b3f18_hi.bin
OUT=[local path]
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
nice -n 15 taskset -c 1 timeout 240s "$P" -EBu "$R/all-arms/exact_receipt.py" "$L" > "$OUT/lo.log"
nice -n 15 taskset -c 2 timeout 240s "$P" -EBu "$R/all-arms/exact_receipt.py" "$U" > "$OUT/hi.log"
nice -n 15 taskset -c 6 timeout 240s "$P" -EBu "$R/physical_check.py" "$L" "$U" > "$OUT/phys.log"
nice -n 15 taskset -c 1 timeout 240s "$P" -EBu round2/checks/high-hub-sol.py "$L" "$U" "$OUT" --solve lo
nice -n 15 taskset -c 2 timeout 240s "$P" -EBu round2/checks/high-hub-sol.py "$L" "$U" "$OUT" --solve hi
```

The fresh producer exact receipts contain the three family values and exact sum for both directions and all three roles, per witness. The lower and upper runs took 55.10 s and 36.60 s; the complete-host run took 3.35 s and matched 1,092 cases per sign (`1316/854/238` nonzero hV/P4/Yt occurrences per sign). Both companion solves passed (40.36/39.38 s): 356 labelled control events, 420 direct polynomial comparisons per witness, and all 36 exact family endpoints across both signs and Bellman directions. The frozen verifier authenticated both original 3,495,091,968-byte payloads and returned the exact `aa511af` endpoints.

| Scratch receipt | SHA-256 |
|---|---|
| `lo.log` | `2d50468e919e6316db6ad75d99b3dbb0af9b43131d942a7f0542cb5a4d481305` |
| `hi.log` | `b01626435491ef5abb8303e7c73dac76e9dfac837ec703571a82574cafe4233c` |
| `phys.log` | `2475139e2f3f736e51a39fc7c5b7c8cc53358803ff2fe58a28561eb3d4e66520` |
| `cold-lo.log` | `eeb679c12e75023143b543439091f535fc1781e84c5dcb8ffa43138471eabff1` |
| `cold-hi.log` | `b9e52f06a4ac123cd1bef86b0148eeb2945021978fef4161db61c243d04531d0` |
