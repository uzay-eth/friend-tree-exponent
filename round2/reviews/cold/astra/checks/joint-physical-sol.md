CORRECT (frozen composition and targeted clean-input preflight)

# Joint physical enclosure — independent Sol audit

**Frozen composition: CORRECT. Clean-input preflight: PASS after repair.** GPT‑6.2 Sol Ultrafast, independent of the calculation's authorship and within the GPT family after the owner’s policy change. The frozen certificates give `[0.565673570322, 0.565784239110]` outwards. The two previously absent imports are now in the manifest, and the cold entry checks them before production. The full joint clean production replay remains pending and is outside this verdict.

## Transfer and disjoint sources

Use the original prior `l=113/200`, `u=1133/2000`. Write `C` for core count, `B=LC` for its birth intensity, and normalize Laplace expectations by that of `C`. The signed witness inequality provides the raw residual `e` per core and the selected nonnegative physical slack. Summing the occupation inequalities at discount `u` gives slack at least `(u−L)T+vB`, where each promotion books a newborn once. Integration by parts trades `L` for the Laplace variable; the initial terms vanish after normalization, and `(u−μ)T` is nonnegative in the limit. The founder remainder vanishes and normalized birth intensity tends to `μ`. On a subsequence with limiting potential per core `m`, the positive raw comparison and the independently reviewed unsigned per-birth bound give `0<m≤2^40 μ y`. Hence

```
σ(λ−μ)m ≥ μv−e > 0
σ(λ−μ) ≥ (lv−e)/(l·2^40·y),  σ=−1 lower, +1 upper.
```

In the [five-pass construction](../five-next/split-audit.md), deleting an undirected leaf of an ordinary five-set recovers a unique receiving quartet, donor and contact; the donor is allowed to be oldest. Different deletions give different incidences. The nonnegative gain gap is `J(a,p)=max(a,0)/2−pa` for the signed quartet gain `a` and the actual donor reciprocal `p`. The whole-source candidates and summed per-deletion candidate are alternatives: at the current maximizing candidate the physical successor dominates that candidate's successor. Each following pass consumes only the actual-minus-worst gain-rate slack of that immediately preceding **combined auxiliary weight**. The intermediate slack terms telescope. Every enlarged connected set is booked using its unique youngest member, and each pass contributes its own singleton newborn value. This spends the original quartet-outside source once, so its archived credit must be removed.

The two original singleton runs target opposite directions of a high-degree singleton incidence. Mask 63 (`occupation-tail/continuation.md#additional-physical-sources`) uses the maximized singleton-row difference, maximized edge-row slack, edge-donor relaxation gap, maximized triple/quartet-row slack, and uncovered triple-donor relaxation gap. Maximized row slack still maximizes outside inputs and does not already contain the actual donor gap. The uncovered triple donor has a triple receiver and four-site support; the replacement has a quartet receiver and five-site support. The unsigned sheet bounds the nominal potential per birth on the original prior and consumes no signed slack. No curved, mask 64, six-site or signed stationary credit enters this sum. The inputs' separate reviews remain [mask 63](occupation-tail-hy4.md), [split](five-split-hy4.md) and [unsigned](unsigned-moments-hy4.md).

| | Lower | Upper |
|---|---:|---:|
| Nominal trial | `5656453/10000000` | `56581165/100000000` |
| Raw residual ticks/core | `115075` | `155608` |
| Archived Q4 removed | `7402100` | `7709975` |
| Five-pass replacement | `36290118` | `32703071` |
| Parent receiver | `11253928` | `8148713` |
| Child receiver | `301604` | `3226994` |
| Mask 63 | `15289028` | `17458365` |
| Total ticks/newborn | `63134678` | `61537143` |
| Unsigned `y` per birth | `2.024576014096` | `2.032667251118` |
| Exact endpoint | `1356995264017948174596364819/2398901654969101844480000000` | `2725370105716046104606536597/4816977775848776007680000000` |

## Executed checks and cold inventory

From the Astra root, `timeout -k 1s 14s taskset -c 8 python3 -B round2/checks/joint_physical_sol.py` printed both exact endpoints, `Frozen checks PASS: 25 split pins, 50 rank counts, 2 unsigned sheets and 10 planted failures`, and `Cold preflight: 74 static files, 25 local include edges, 4 isolated import probes; 7 named downstream outputs excluded from seed`. Exit 0, without `--overlay`. The ten rejected plants cover double mask 63, hidden archived Q4, wrong sign, inward rounding and wrong scale on each side. The checker imports no author arithmetic, hashes the frozen raw witnesses and dependent data, and binds the nominal profile headers and exact unsigned price-sheet bounds. It does not re-solve the dual or rescan production rows.

The original 72-file manifest failed with `ModuleNotFoundError: No module named 'negative'` because `mutations.py` also needs `corner_audit.py`. Both modules were already tracked at `210805b`. A temporary two-file overlay had passed the isolated probes. The runner has now added `round2/motif-high/negative.py` and `round2/motif-high/corner_audit.py` to the actual allowlist, and both preflight and cold entry call the isolated raw-import probes before production. The independent repeat above confirms the input repair without supplementing the manifest.

The independent inventory follows the handoffs: the physical stage generates both snapshots, both finite mbar tables and both analytic tail profiles; the split stage consumes the snapshots and emits its checked cold transcript; the six K7 generators compare to frozen gzip contents, and the unsigned checker consumes the same pinned compressed inputs plus rebuilt snapshot/mbar bytes. The two dual sheets and six gzip inputs are intentionally static. The three wrapper sources (`replay.py`, `arithmetic.py`, `check_tables_k7.py`) and its `inputs.txt` are entry-point files outside the copied allowlist. This targeted preflight is not a production replay, a cross-family verdict on the composition, or an eight-digit kernel check.
