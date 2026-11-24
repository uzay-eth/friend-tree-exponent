> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits and Claude replay of earlier freezes; original Neon-2 Gold1 witnesses and inherited stationary results. See the [release index](../../../../../../README.md#claim-25) for exact scope and citations.

# Lower-quintet FY: older paying neighbor

The widest completed cut certifies **an additional `4.04727235e−7/core`** on the lower witness, **83.6%** of the existing `175/5/8` FY floor (`4.84195503e−7/core`). Their disjoint sum is at least **`8.88922738e−7/core`**. The additional cut has a matching ceiling **`1.50212724e−6/core`**, ruling out `5e−6/core` **on this cut**. That ceiling is not a bound on all unpriced FY terms, on older donors, or on higher witness rows.

| Host / paying degree / donor degree; host / payer arity | New floor per core | Matching ceiling per core | CPU time, floor / ceiling | Path states visited per pass |
|---|---:|---:|---:|---:|
| `120/22/4; 8/5`, binned | `3.00175335e−7` | `1.18764734e−6` | `14.001 / 13.820 s` | 723,492 |
| Same, exact payer-degree reciprocal | `3.03570825e−7` | — | `13.247 s` | 723,492 |
| `175/36/5; 10/6`, binned | `4.04727235e−7` | `1.50212724e−6` | `67.142 / 74.378 s` | 4,284,000 |

Bounds are rounded outward. The narrow binned price already adds **62.0%** of the old floor and delivers **3.56×** the floor per CPU second of the widest cut (the bins cost about 1.12% of that narrow floor versus singleton bins). Widening raises the floor by roughly 35% at roughly 4.8× the floor CPU cost. This is an ROI comparison between measured extension runs and a credit comparison to the pinned baseline; the old receipt does not contain a directly comparable CPU timer. On the wide cut, the three contributions are `6.88372377e−8` (mother pays), `1.39177903e−8` (child pays, mother donates), and `3.21972207e−7` (child–child, sum of both directed orders), approximately. The matching ceiling shows that reaching `5e−6` requires a wider or different mask.

## Directed transport and disjointness

At an exact `G2` host outside `M3` with at least two neighbors, the receiver is the neighbor on the paying edge. The retained ordered occurrence requires **receiver degree `>5`**, donor degree `≤5` (or `≤4` in the small cut), and `max(host degree, receiver degree)>18`. The first gate excludes the old FY mask, which requires both other vertices to have degree at most 5. The last gate excludes the previous full edge slack. The complementary older-donor lane (donor `>5`, payer `≤5`) is also disjoint. No Jensen or `M3` payment enters this calculation. Footprint is unchanged by vertex age: `D_receiver>5` is tested **at reward time**, not at promotion or at a guessed age. Both degree and arity caps apply to physical states, and donor arity is unrestricted within its degree cut.

The birth transport starts with the physical newborn state `J`. For a mother–child pair, track the prospective mother from her birth and register her edge when she promotes the future host. For two children, track the host from birth and register an old edge at each promotion. At every later promotion by that host, the old-edge stages create the completed directed paths: one with the older child paying and one with the younger child paying. The older-receiver and older-donor stages use different caps and are propagated separately. A newborn paying neighbor stays at zero **reward** through degree 5, while its continuation follows every actual gain and promotion, including those before degree 6. Internal reciprocals in all transition rates are exact. Unseen neighbors only widen the receive hazard by `[0,(c−tracked)/2]`. No birth-rate reset is made on entry to a bin.

For each host, payer arity and physical payer degree in a reciprocal bin, take the exact signed knot maximum `max_k σ(ξ_k z−A(ξ_k))` at that **physical** slope `z=Δ^L g`, then subtract the physical donor evaluation. The extrema across payer degrees give reward lower and upper bounds within the bin. Exact cut gates are evaluated separately; bins cannot promote excluded states. The mean within-bin reward width on the wide cut is `122413.82` ticks over its 238,464 unweighted bands, not a stationary width.

For a floor, solve at `HI` with adverse low Bellman ratios and round down; for a ceiling solve at `LO` with adverse high ratios and round up. The same original three-stage balance cancels registrations against preceding stage outputs. All stopped states outside the degree and stock cuts have zero continuation, since neither coordinate decreases. The only remaining source is core birth, giving a lower factor `LO` or upper factor `HI`. For each state the exact adverse affine inequality is rechecked at the rounded continuation; there is no floating-point decision in the certificate. For both floor and ceiling the initial mother and host are populated from the physical core birth flux.

Implementation stores only the current and next host-degree layers, with current-layer promotions and donor/receiver coordinate jumps swept backwards. There is no stored dense three-dimensional history. This saves storage, not enumeration work: all supported lattice states swept per pass are listed above.

## Replay

From the unpacked project root, restore the private authenticated inputs as described in [the FY README](../README.md). Run serially on one niced core (change core 5 if unavailable):

```sh
taskset -c 5 nice -n 19 python3 -Bu round2/occupation-quintet/edges/fy/extended/replay.py \
  --host 175 --receiver 36 --donor 5 --stock 10 --rstock 6 --hash
taskset -c 5 nice -n 19 python3 -Bu round2/occupation-quintet/edges/fy/extended/replay.py \
  --host 175 --receiver 36 --donor 5 --stock 10 --rstock 6 --ceiling-only --hash
```

`--hash` reauthenticates the private binary without rerunning any production verifier. The measured runs reused the authenticated copy. Exact results are in [floor](cut36-floor.log), [ceiling](cut36-ceiling.log), [narrow binned](cut22.log), and [narrow exact](cut22-singletons.log). Each linked receipt ends with the exact corner-check and unchanged-copy `PASS`.
