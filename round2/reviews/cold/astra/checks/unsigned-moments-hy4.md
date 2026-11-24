CORRECT

# Referee audit: V5 UNSIGNED per-newborn denominator (hy4, exact/astra)

Scope: the stronger immutable K7/144 pair
`round2/outputs/moments-v5-unsigned-fine-{lower,upper}-dual.json.gz`
(per newborn 2.024576014096 lower, 2.032667251118 upper) as first priority,
with the K5/80 pair `round2/outputs/moments-v5-unsigned-{lower,upper}-dual.json.gz`
as control. Supporting documents: `round2/moments/RESULTS.md`,
`round2/moments/README.md`. No signed donor, reserve, tip, star, floor or
lower-score file was opened at any point; no physical surplus was spent.

## Inputs disclosed
- Sheets (4): `round2/outputs/moments-v5-unsigned{,-fine}-{lower,upper}-dual.json.gz`.
- Per sheet the only data files opened are the sheet, `mc-{side}.bin`,
  `mc-mbar-{side}.i64`, `mc-mbar-{side}-tail.tsv`, and the three profiles
  `moments-v5-{side}-k{K}.tsv.gz`, `-score-k{K}.tsv.gz`, `-edge-score-k{K}.tsv.gz`.
  This is enforced by the audit hook in `round2/moments/check_unsigned.py`
  (`seen == allowed`) and by the `inputs` digest equality in
  `round2/moments/v5.py:546`. Both fine sheets carry `score_only=True`,
  `use_prior=False`, `edge_donor/lower_score/tip_edges/star_children/original_score/
  singleton_gap/tail_score` all off, and rows exactly
  `{M-per-core, M-per-birth}`; the LP interval is therefore the originally
  published `LOW=113/200, HIGH=1133/2000` (`v5.py:519-520,535`).
- Growth prior [0.565, 0.5665] confirmed: `round2/moments/occupation.py:8`.

## Commands run (repository root; individual pilots < 30 s except the two consolidated replays)
- `python3 -B round2/moments/check_unsigned.py round2/outputs/moments-v5-unsigned-fine-lower-dual.json.gz` -> PASS (22.9 s): M-per-core 1.146896986910 (margin 7.578e-11), M-per-birth 2.024576014096 (margin 4.15e-11), "Input-restricted replay passed".
- same for `-upper-` -> PASS (28.6 s): 1.151465636722 (9e-12), 2.032667251118 (3.46e-10).
- `bash round2/moments/replay.sh denominator unsigned-fine` -> PASS (47.0 s, one core, pinned): reproduces both certificates and prints outward intervals [0.565653724789, 0.565803233551] (three lifetime reserves) and [0.565660593037, 0.565795421997] (plus affine source).
- `python3 -B round2/checks/unsigned_moments_hy4.py unsigned-fine both` -> PASS (74 s, includes two model builds, duals, plants, census). `... unsigned lower` (K5 control) -> PASS (9.8 s).
- All commands are single-core; the two long ones are the exact documented replay and my one-shot consolidated script, not repeated pilots.

## What was re-derived and checked
1. **Palm/root and oriented-edge occupancy.** In `occupation.py`: `unit` sums the finite vertex colours plus the two tail colours to [1]=1 (`:89`); `child` gives incoming parent-child edges = vertex count (`:90-93`, root negligible per README §Occupation); `parent` gives outgoing edges = (c-1)v_s (`:94-98`); path/cherry/four-site marginals (`:99-123`). Each is an equality the physical measure satisfies, so the LP is a relaxation and an upper dual bound transfers to it.
2. **Price identification and signs.** `dual.py:verify` implements the Lagrangian lower bound on `c'x` for `Ax=b, Gx≥h, 0≤x≤U`: free equality prices, nonnegative inequality prices, plus the capped-variable correction. I reimplemented this from scratch in `unsigned_moments_hy4.py:myverify` and obtained bit-identical margins on all four sheets. Direction `-1` on `M-per-birth` proves `den ≤ y·nu`, i.e. an upper bound on potential per newborn, which is the correct direction for a denominator.
3. **Rounding/certification.** Prices are integer at scale 1e14; exact `Fraction` margins are strictly positive (1e-12..3e-10) on every row, and `IntegerProof` reconstructs the same inequality. The box correction covers rounded prices at the degenerate optimum.
4. **Parent/child orientation and all-degree tail.** `v5.py:_vertices` puts the parent reciprocal on the incoming edge and child reciprocals on outgoing edges; `_edges` carries `rin[b]=range(a)`, `rout[a]=range(b)`; `_tail` and the analytic caps cover every degree band through 1399 and the infinite band. `test_occupation.one_step` (direct physical one-step generator) passed 8/8 cyclic fixtures at K5 and K7, and my 7-tree census (degrees 2 … 2,000,003) satisfied `phi <= den <= mbar` pointwise once tail cells are included — confirming the reciprocal orientation and that the all-degree tail bound holds where finite-degree tables stop.
5. **Applicability to an independent positive numerator.** The certificate prices only the potential; the input guard and sheet flags confirm no reserve/donor/signed file is read. Combining with any separately proved positive per-newborn numerator `v` and per-core residual `e` gives `σ(λ-μ) ≥ (v - e/μ)/(S y) ≥ (v - e/0.565)/(S y)`; positivity requires `v > e/0.565`, checked in `denominator.py:37-38`. This spends none of the signed surplus.
6. **Exact substitution.** My independent `fractions` computation of `lam ± (credit/S - loss/(S·LOW))/y` reproduced the replay endpoints digit for digit (see Commands).

## Negative controls (plants)
On the fine upper M-per-birth sheet: baseline margin +3.46e-10; replacing the claimed ceiling by 2.0 gives -0.404; flipping one equality price's sign gives -2.2e-9; shifting one inequality row index by +1 gives -0.020; a negative inequality price and an out-of-range index are both rejected. The certificate is therefore correctly indexed, signed and tight — it proves the stated ceiling and no materially lower one.

## Limitations / open (not defects in the audited artifact)
- I did not regenerate the C++ row tables (`replay.sh profiles`) nor run `check_tables.py` (`unsigned-inputs`) in full; those are longer checks. Instead I relied on the embedded input-digest integrity and the physical-tree pointwise bound `phi ≤ den ≤ mbar`.
- Small finite trees need not satisfy the stationary growth/tail inequalities (`growth-low` etc.); those constraints are asymptotic, so the census validates the pointwise potential bound and the generator identities, not literal LP feasibility of each finite tree.
- `round2/moments/RESULTS.md` still quotes the K5/80 ceilings in its main table and does not yet cite the fine pair; that is a documentation lag only. The fine endpoints above are not claimed there.

## Verdict
The two rational unsigned-fine dual sheets are correct, stronger than the K5 controls, use the original [0.565,0.5665] prior, and are valid for any independent positive numerator without spending signed surplus. Standalone script: `round2/checks/unsigned_moments_hy4.py`; this report: `round2/checks/unsigned-moments-hy4.md`. No author files were modified; nothing published, committed or paid.
