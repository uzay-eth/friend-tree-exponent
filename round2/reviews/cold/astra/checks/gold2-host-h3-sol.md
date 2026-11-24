# Cold audit: GOLD2 D24 host plus matched 2417 H3

**HOLDS on both original GOLD2 witnesses.** The exact interval is **`[565714675151/10^12, 565736253459/10^12]`**, or **`[0.565714675151000, 0.565736253459000]`**. Both final one-sided restricted trials have strictly positive rational margins. The lower endpoint is strictly above `99/175`. The separately derived, unrestricted D24 denominator-only control is `[0.565714675049909, 0.565736253540236]` with outward rounding.

## Physical ledger and the sign of the denominator

Write `d=+1` lower or `−1` upper and `x=d(μ−λ)`. With the stationary original potential `M_π>0` in ticks/core, the bound is `x M_π ≥ B + π(H−h)`, where `B=F+G+C_H3−E`. Here `F` is the published GOLD2 conditioned singleton floor, `G` the eligible edge floors, `E` the original outward error (`19217/20062` ticks on lower/upper), `H` the original global higher cap (`18370/17240`), and `h` the capped host allocation. The matched H3 inputs are `1084423030952180046412413909/1550425090937500000000` and `3304204122287698944634328103/6201700363750000000000` ticks/core. They come from the two normalized threshold receipts multiplied by `2^40` once; the lifetime calculation already includes the single birth normalization. The separate `higher/` alternative and its mother-80–100 search are **not** in this ledger.

`B` is strictly positive on both signs, and the D24 host savings lower bound is nonnegative. Divide the positive net by the **upper** stationary potential bound, not the lower bound. The exact paired bounds are independently repriced from the 307 cut-24 cells with saturation seven. The lower cell allocation includes the original recharge whose tree sum is a fixed constant; the upper allocation is unrecharged. Both bound the same stationary potential. One does not require the lower and upper numbers of each individual cell to be ordered. The separate signed offset is `(B+saved)/M_hi`, added to the lower nominal and subtracted from the upper nominal. The replay prints the full rational endpoints before outward rounding.

## Cap preservation and disjointness

At a physical supported H3 middle `m`, the original GOLD2 signed row `P` is bounded by the original verifier's cap for that middle, using its physical endpoint interpolation. Adding a payment `0≤A≤[-P]_+` satisfies `P+A≤max(P,0)≤max(K_m,0)`. The same original per-middle count `C(c,2)` therefore still applies, with one payment per unordered path owned at its middle. The D24 reconstruction uses the published family maximum as a conservative ceiling on the unpublished per-middle maximum (408/233 ticks for H3); it does not assume all individual maxima equal it. The checker independently counts the physical host roles of all six higher families, including the P4 half at each internal and actual allowed partner arities, then compares every cap cell after taking the minimum with the original global higher cap. Consequently `H−h` is unchanged and paid once. The matched reward pays the original physical negative part on `19≤D_mother≤79,c_mother≤8`, middle `D≤5,c≤4`, child `D≤18,c≤3`, with the nonfounder mother gate; it is not an envelope-gap payment.

The lower side keeps D18 `g` and the two audited FY floors. D18 has both paying-edge endpoints `≤18`; each FY floor requires the paying-edge maximum `>18`. The absorbed small payer has degree `≤5`; the wide payer has degree `6..34`. The tests are disjoint at the actual charge time even if the same vertices can appear in another role. FY lives at non-M3 `g` hosts. The original `g→H3` transfer at an M3 middle is already assigned to H3 and carries no FY payment. The singleton floor is on its own original `f` row. No edge credit is booked on upper. The checker authenticates and reconstructs these exact inputs; the matched H3 and FY lifetime solvers were not rerun as part of this composition audit.

## Joint fixed trials

For `a=d(s−λ)>0`, price `H−h+a(M_global−u)` on each occurrence column with the original upper cell comparison `u`, and zero on auxiliary columns. The exact margin is `B−a M_global+rebate ≤ (x−a) M_π`. Use the original GOLD2 prior for the uniform trial at the denominator-only endpoint. A positive uniform margin permits an additional inward step of `margin/M_hi`; the saved uniform endpoints match this exactly.

For the stronger fixed trial, if the desired endpoint were false, the physical discount lies in `[LO,s]` on lower or `[s,HI]` on upper. Narrow only the shared physical LP's discount prior to that range, and apply its lower endpoint to convert the nonnegative dual objective from `η=π/μ` back to per-core units. The already certified birth/H3/FY inputs and upper cell comparisons stay valid. The strictly positive fixed-trial margin contradicts the adverse case, including equality. It is a certificate at that trial; its conditional margin is not subsequently divided to claim a stronger endpoint.

| Exact replay output (margins displayed in ticks/core) | Lower | Upper |
|---|---:|---:|
| Uniform margin | `136.02468951319864` | `107.97562535014244` |
| Restricted fixed-trial margin | `41.74286410186889` | `9.071193851267111` |

## Verification

Receipt SHA-256 `ae6f5cde2ef710b4196c3c445a3786580e5ccc9dfd9206b2962f60628866c6ed`; corrected host input SHA-256 `3f98385e0662b4a5aea71dcfef03cdbd66012de01f0083456e243592007de790`. The checker (`gold2-host-h3-sol.py`) reproduces the exact revised receipt without invoking the producer compositor. It freshly authenticates both GOLD2 payloads, reconstructs all host cells through the pinned host reader, independently counts cap roles, and for six objectives per sign proposes a dual then checks all 5,542 rational column residuals on the original 601 rows. The final uniform and restricted rebates, margins, endpoints, and directed decimals match exactly.

```sh
taskset -c 11 nice -n 19 env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  round2/.work/venv/bin/python -Bu round2/checks/gold2-host-h3-sol.py
```

Exit 0: `PASS independent GOLD2 H3+D24 0.565714675151000 0.565736253459000`.
