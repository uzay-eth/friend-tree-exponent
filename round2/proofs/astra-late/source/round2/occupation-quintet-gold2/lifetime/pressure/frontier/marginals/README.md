> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical-floor and frontier audits and scoped Claude conditional reduction and arithmetic; original Neon-2 Gold2 witnesses and the four reward floors. See the [release index](../../../../../../../README.md#claim-28) for exact scope and citations.

# Claim 10 conditional marginal overlay

**HOLD: the H3-dependent endpoint and gap numbers below are withdrawn.** The Claim 10 intervals, raw singleton roots, and marginal-only singleton/denominator comparisons are unaffected. See [`../repair/`](../repair/) for corrected lower-side accounting.

**Conditional on Neon-2 Claim 10 pending the independent Astra audit (863c and 2f74).** This replay does not verify the sixteen underlying marginal binaries or make an additional spectral assumption.

The overlay gives **`[0.565714780208554, 0.565736144184918]`**, versus matched Gold2-only **`[0.565714767317686, 0.565736152862047]`**. It moves the ends inward by about `1.2891e-8 / 8.6771e-9` in exponent units. The dominant gain is a tighter potential denominator (`1.343919880562e12 / 1.344347431578e12` ticks/core, down from `1.359401542040e12 / 1.360177372417e12`), not the roughly `964 / 793` ticks/core improvement in singleton floors. `241/426` still lies inside.

| Conditional result, ticks/core | Lower witness | Upper witness |
|---|---:|---:|
| Singleton floor | 333141.795910075 | 432682.450341019 |
| Singleton ceiling | 1509409.622162206 | 1607477.961544022 |
| Extra credit for 2e-6 from the cold-audited paired baseline: necessary / sufficient | 607566.676757462 / 2546650.824378763 | 760571.083535735 / 2541792.529890666 |
| Extra credit for 241/426: necessary / sufficient | 12616453.081887697 / 17362533.627379387 | 8045194.806163809 / 11352549.453463841 |

The necessary bounds are for the fixed non-singleton ledger even allowing the entire true original singleton and all higher-cap savings; the sufficient bounds use the existing floors and a new disjoint payment. Strict excess excludes a rational. The matched control gives sufficient budgets for `241/426` of 17579455.781670218 / 11497268.529874944 ticks/core. No overlap with the existing host, H3, FY or frontier credits: the overlay is a replacement, not an additional payment.

## Exact input and common normalization

The replay reads `round2/neon2/stationary-marginals/certificates_all.json` directly from Neon-2 commit `6c8afec5cd1b7f98bbcd7c6e7848c0d631a50f7d` (SHA-256 `cc8658e2e6354577716b44ff7f57bf6e3bb9831a1c7a7ac0ecf5c9bf717dbf9c`). From the JSON it uses the rational `certified_interval_exact` for five degree marginals, `J=(2,1)`, `(2,2)`, and all bare cores, and the two exact conditional-receive intervals. The exact Gold2 prior lies inside the Claim 10 input bracket. The board's displayed decimals are not used.

Claim 10 uses the physical stationary core-root marginal, with the existing identification of its occupation law with the physical core eigenmeasure (Neon-2 Claim 3 and the Astra operator/fringe identification). The lift here uses the same stationary measure: `pi B1=mu`, `beta=pi B/mu`, and `pi f=mu E_beta int exp(-mu t)f dt`. In the LP `eta=pi/mu`, so `eta(b)=1` and `eta(1)=1/mu`; `b=(D-c)/D` is the physical promotion intensity. Thus a certified marginal `[l,h]` becomes `l eta(1)<=eta(g)<=h eta(1)`, and a conditional receive `[l,h]` becomes `l eta(I)<=eta(I rho)<=h eta(I)`. At `J` only the mother port contributes. At `(2,2)` both actual core ports contribute. The receive identities follow from the exact physical two-state balances: `eta(J)*(mu+3/2)+eta(1_J rho)=1` and `mu eta(2,2)+eta(1_22 rho)=eta(J)/2`. These are already enforced by age balance with the lifted discount variable; the overlay adds the certified interval consequences.

For `D<=24` the bare indicator is exactly `c=D`. Every bare vertex with `D>24` is in the saturated overflow cell `(25,7)` in `(D,k)` labels; the converse is not used. Introduce an unknown nonnegative bare overflow mass bounded by the occupation of that cell. Use the sum of that mass and the exact bare cells in both inequalities. The Claim 10 JSON endpoint already includes its analytic bare-tail correction at 1600, so there is no separate credit or tail truncation here.

The dependencies run from the original Gold2 exponent prior and the established physical occupation identification to Claim 10, then from its exact intervals to this new overlay. Claim 10 does not depend on this overlay, the stronger Astra interval, or marked singleton prices. Both objective branches keep the exact original Gold2 prior. Tagged birth roots and old non-singleton payments are unchanged; the overlay constrains the feasible physical lift, it supplies no new additive reward. Each inequality gets its own slack class, with `1/LO` a valid bound for exact dual repair. The hidden bare mass is also bounded by `1/LO`.

## Replay

From `.`, on one niced core, run each `side` and `bound`:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 8 nice -n 17 timeout 240s \
  round2/.work/venv/bin/python -Bu round2/occupation-quintet-gold2/lifetime/pressure/frontier/marginals/replay.py \
  --side lower --bound floor --export
```

Use `--side upper` for the upper witness and `--bound ceiling` for the all-stock singleton upper comparison on either witness. Each output pairs the otherwise identical unconditioned-on-Claim-10 control with the constrained overlay. JSON exports include all eleven exact birth-mark roots keyed by `D,post-c`, original-Gold2 sources and discount, and exact objective values. `(0,0)` is the overflow mother mark. The roots and `results.control` do not rely on Claim 10; only `results["claim10 conditional"]` does.

All four runs completed with 2,800,083 physical states and 17,262 explicit corners per run. Claim 10 adds 21 rows and 22 columns (control 623/5564; overlay 644/5586). Receipt SHA-256 pins:

| File | SHA-256 |
|---|---|
| [`lower-floor.json`](lower-floor.json) | `95b05a81de95b5ef559414d006a23acc82f4790126db69d8509fe488277f979a` |
| [`upper-floor.json`](upper-floor.json) | `8b1b54a47a621e0f15006dfa5ec94aa26db53ce022e3502f793938ba271aae08` |
| [`lower-ceiling.json`](lower-ceiling.json) | `8810f46994605bc8839d18d13c5c5e79a40d3febd6d747ab94074dbf13473bfd` |
| [`upper-ceiling.json`](upper-ceiling.json) | `75395dd5ab62bd12373a5ab4b5bf4e404c4456c9365365cadfbab5b848293e4f` |
