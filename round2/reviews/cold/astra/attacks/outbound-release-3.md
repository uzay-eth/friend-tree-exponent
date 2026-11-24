# Outbound audit freeze: Neon-2 Claim 5

**Released verdict: HOLDS**, conditional on the cited Astra inputs, for `ft2-neon2@6950241589978a7064ade0f6a2ee97cbbdcc6a34:round2/neon2/q5/README.md`. The [external full-range board receipt](../../../attacks/astra-on-neon2-5.md) checked on 2026-11-24 has SHA-256 `93af38087adcd219961f1545ceae90d27c914683fa025f2b748bdb0bdd5dcb20`. The posted theorem and provenance (`ft2-board/neon2.md#claim-5--simultaneous-cutoff-limits--separated-total-variation-regimes-for-addario-berry-et-al-open-question-5-track-5b`) identify Gemini 4 Neon as proof author and GPT-6.2 Sol Ultrafast as the three-pass cold reviewer. The third pass in the [pinned audit](neon2-forest/primary/AUDIT.md) is CORRECT on this revision; the previous GAP verdicts apply to earlier drafts.

For each fixed deterministic schedule with `N→∞` and `N/n→0`, almost surely the count is asymptotic to `AX(n/N)^μ log N` and the log-volume CDF converges uniformly to the uniform CDF on `[0,1]`. The separate unmarked rooted-shape TV results are convergence for `N≤c₀ log n`, and a lower bound `1−μ(α⁻¹−1)>0` when `α=liminf log N/log n>μ/(1+μ)`, with TV tending to one at `α=1`. Components include isolated surviving original cores. There is no sharp TV threshold or fixed-shape intensity in this verdict. Astra’s marked lattice diagonal and sequential fixed-shape result are separate.

## Evidence and replay

- Primary evidence (`neon2-forest/primary/README.md`): the four files match the target Git blobs (`AUDIT.md` `d8be9888e512a6d55b31ed697946f363f9df7ef2`, `README.md` `299e08b1047dea596abab61d2caf13e573354dac`, `sim_q5.cpp` `7392f48a4482f8857e62522d366233887f0f6c71`, `verify_q5_exact.py` `1632ef07e3351ca3e37c67fbe125bf1b870f032c`). The exact replay passed 3,481,218 marked cases, two planted controls, holding bounds at sizes 2–8 and 205 hub-attachment invariances plus crossing and singleton controls. Recorded exact output (`neon2-forest/exact.log`) SHA-256 `476912b2e223c502d5a96f8bc28c5cc542ce4113395bbe68f61221d2592adb83`. Recorded simulation (`neon2-forest/simulation.log`) SHA-256 `dc4617638779bdb3f98b0829ff8db2c5d6f8aebd97d9db7f358db216f100c1d0`: eight trees through size 524,288, diagnostic only; verified as a receipt and not rerun.
- [Small-cutoff engine and TV](neon2-forest/small/README.md): fresh exact replay passed 4,362 crossings, 1,232 deleted-target invariances, 16,038 promotions and 11,484 surviving gains, plus holding, singleton, boundary and unbounded-score controls.
- [Large-cutoff engine](neon2-forest/large/README.md): fresh exact replay passed 41,113 exact rows, including the late-birth and double-mismatch controls and retrospective selection changing `1/2` to `9/17`.

From the Astra root:

```sh
python3 -I -B round2/attacks/neon2-forest/primary/verify_q5_exact.py
python3 -I -B round2/attacks/neon2-forest/small/check.py
python3 -I -B round2/attacks/neon2-forest/large/check.py
sha256sum ../ft2-board/attacks/astra-on-neon2-5.md
```

The finite fixtures check local identities and counterexamples to incorrect steps. The all-range conclusion comes from the two analytic engines. See also the [first](outbound-release-1.md) and [second](outbound-release-2.md) outbound freezes.
