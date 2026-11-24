> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits and Claude replay of earlier freezes; original Neon-2 Gold1 witnesses and inherited stationary results. See the [release index](../../../../README.md#claim-25) for exact scope and citations.

# Physical edge slack for the gold quintet

At maximum endpoint degree **18**, the original edge-cap residual pays at least **156.85265 ticks per core lower** and **337.42462 upper** (rounded down). These are **39.1%** and **35.7%** of the original `g` caps of 401 and 946. With the original `m_min` denominators alone these correspond to inward movements of approximately `1.5704e−9` and `2.1736e−9` in the respective endpoints. These floors are for an edge-only addition to the original singleton ledger; the higher per-core caps stay unchanged.

## Physical row and accounting

Use the verifier sign `σ=−1` lower, `+1` upper (opposite to the census inward sign), and keep the original nominal `λ`. For each **unordered physical edge** `{u,v}`, start with `−λ g(u,v)`. At endpoint `u` with state `(D,c)`, `L=D−c`, let `d=Δ_u^L g`, `z_w=1/D_w`, and `r_1` be the stored center when `u` is exact `G2`. Add:

* `(L+z_v)d` and, on promotion, `(L/D)[Δ_u^B g+H3(v;u^+,J)]`, with newborn `J=(2,1)` and the middle evaluated **after** promotion. The term is zero outside the supported middle. Both interpolation and jumps use exact physical degrees.
* `(z_v−r_1)Δ_u^L f` iff `u` is exact `G2`.
* If `u∈M3`, add `(c−1)r_1 d`; the remaining terms `(z_w−r_1)d` go to `H3`. Otherwise, for exact `G2` with `c≥2` add `(c−1)A_u(z_v)+Σ_{w≠v}[z_w d−A_u(z_w)]`, with the exact piecewise-linear FY function on the original knots. Else add `Σ_{w≠v}z_w d`.

Add the two endpoint contributions and multiply by `σ` to obtain `G`. Reversal does not change it; equal endpoint labels do not halve the physical edge. At each FY vertex the terms in `A_u` cancel when summed over its incident edges, by exchanging the ordered incident pairs. A promotion creates precisely one new edge: its *potential jump at creation* is in `f`, never in `G`. Every newborn `H3` is charged to the old edge left when deleting `J`. No additional `H3` cap residual is taken. The replay checks the exact identity `f_charged + g_charged = f_raw + g_raw − g_to_H3` in unsigned ticks, where the raw rows on the right already own their corresponding newborn motifs. This uses the same original centered singleton convention as the census and lifetime streams.

By the authenticated original `g` verification, `G≤cap_g` for all states (including overflow and both degree tails); set `R=cap_g−G≥0`. Pointwise the combined outward drift is bounded by `(cap_f+cap_g+higher_per_core) C − Σ_f slack_f − Σ_edges R`, since the positive `cap_g` is paid at most `C−1` times. The two slack sums charge different original row slots even though both use the same promotions as physical events. The edge calculation below only pays the masked part of `R` and does not reuse any singleton residual.

## Discounted physical certificate

Stop payoff when either edge endpoint exceeds the degree cap. Neither degree decreases. At a retained pair state let `V` be the nonnegative subsolution at the high discount `H=0.565748407438`. Its Bellman inequality is `H V ≤ R + L_edge V` for the *actual* endpoint hazard at every outside environment. For each endpoint a donor adds an affine expression in its reciprocal in each FY knot interval (or a single affine expression without FY). The worst sum over `c−1` independent box donors is attained by taking all at an extremizing knot. Checking the product of the two endpoint knot sets therefore suffices; the same donor affects its reward **and** the gain hazard. At each pair the next states have strictly greater degree or core arity, so the inequalities are solved backwards and rounded down. The replay checks every corner with exact rational arithmetic. Outside the mask use `V=0`.

To bound the birth input, solve another stopped Bellman subsolution `W` on the actual promoter’s state, with reward `b V((D,c+1),J)`, `b=L/D`. At a promotion the promoter moves to `(D,c+1)` regardless of the edge created; the new edge's value is booked once. Its gain uses actual receive in `[0,c/2]` and both endpoints of that rate interval are checked. From stable age, `π[bV_birth] ≥ μ W(J) ≥ LO W(J)`, since `μ≤H` and the focal typical core starts in `J`. To apply stationary balance to `V`, assign the **one edge value** to its younger endpoint, or equivalently sum `V` over unordered edges. The physical stationary identity is `μ π[V_edges] = π[L_edge V]+π[b V_birth]`; hence `π[R_edges] ≥ (H−μ)π[V_edges]+π[b V_birth] ≥ LO W(J)`. The birth mark is only for counting and does not enter `G`. Both test potentials are bounded per core, with bounded nonzero local rates on their stopped support. The founder carries no birth mass in the stable-age law.

## Replay

From the Astra root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 5 nice -n 19 timeout 100s \
  python3 -Bu round2/occupation-quintet/edges/replay.py
```

[Receipt](replay.log): both sides passed nine physical fixtures (49 edges each; includes FY, overflow and double-tail), plus every Bellman corner for degree caps 4, 8, 12, 15 and 18. The cap-18 calculation has 14,535 unordered pair classes and 170 singleton classes per side. Replay reads only the named private authenticated copies and golden receipt, checks headers and file stat stability, and does not rerun the original verifier or scan the giant row arrays. `--hash` optionally rehashes the copies. Provenance: private source `b04cc50c1e391e0420eb7b006ac781d484c68543`, quintet `PROOF.md` §3 and `verify_quintet.cpp` §§3–4. Tick scale is `2^40`.
