# Neon-2 Claim 10: independent exact replay

This audit uses the pinned revision `6c8afec5cd1b7f98bbcd7c6e7848c0d631a50f7d`, an ancestor of `exact/neon2`, from `[local path]`. The source archive contains `round2/neon2/stationary-marginals` and its sibling `quintet-cert` dependencies. All ten marginal source and receipt files were byte-identical to the supplied `6c8afec.tar`.

**Result:** All sixteen authenticated snapshots passed. The full replayed JSON was byte-for-byte identical to the 57,761-byte Git-pinned receipt (SHA-256 `cc8658e2e6354577716b44ff7f57bf6e3bb9831a1c7a7ac0ecf5c9bf717dbf9c`). The replay took 360.16 s wall time and peaked at 1,648,856 KiB RSS. The independent endpoint and bare-core arithmetic passed on the replayed output. The manifest, endpoint errors and evidence digests are in `replay.json`.

## Reproduction

From the private scratch directory, with `source/` containing the extracted Git archive and `snapshots/` containing the sixteen read-only copies:

```sh
g++ -std=c++17 -O3 -fopenmp -Wall -Wextra -Werror -pedantic \
  source/round2/neon2/stationary-marginals/verify_marginal.cpp -o verify_marginal

# Authenticate all sixteen snapshots against the Git-pinned certificates_all.json
# before invoking the verifier. See replay.json for the input manifest.
nice -n 19 taskset -c 11,13 env LC_ALL=C python3 -B \
  source/round2/neon2/stationary-marginals/verify_all_marginals.py \
  --verify-bin ./verify_marginal --dir snapshots --threads 2 > replayed.json
cmp source/round2/neon2/stationary-marginals/certificates_all.json replayed.json
```

Each side reports 49,204,104 finite motif rows evaluated at each of the two discounts, and 8,035,467 discount-independent envelope checks. The batch includes the 16-panel rational analytic singleton tail and 35 concrete-tree oracle fixtures per side. The C++ verifier binds its in-memory input by SHA-256; the batch checks that digest against its own byte buffer.

## Signed coboundary and bare-core tail

For each endpoint, take the maximum of the two finite singleton caps and the analytic singleton cap. The core tree has at most `C` edges, so a positive edge cap costs at most `C` times that cap; a nonpositive edge cap can be dropped. The higher-motif cost uses the maximum of the sum of positive hosted charges per core vertex: `H3`, `S3` and `S4` at their centers, `P5` at its middle, `B5` at the degree-three interior host, and `P4` charged in halves at its two interior vertices. The one-vertex forcing has sign `−1` for lower and `+1` for upper. Thus the endpoint residual is bounded per core by `max(0, f + max(0,g) + higher)/2^40`. At fixed tree the signed residual is affine in the discount, so the larger endpoint cap suffices throughout the bracket. The independent receipt arithmetic checks both endpoint aggregates, the linear potential bounds and the resulting signed exact enclosures for every side.

At the true discount, Corollary 2.1 of the pinned Claim 3 (`ea6a093c247d6b9bfd38dc27346bfe3c54bcf012`) supplies `(L−μ)R = π(ψ)C−Aψ`. Adding `σR` cancels the physical reward with the correct sign. The combined potential is bounded by a constant times `C`; stopped Dynkin is justified on finite time intervals by the Yule bound, and the discounted expected core count approaches a positive constant. Its time integral grows linearly while the boundary stays bounded, giving `σ(π(ψ)−θ) ≤ ε`.

For bare cores the finite observable stops at `D = 1600`. Let `h=0` for `D<6`, and `h=2^(-max(1600−D,0))·4^(-L)` otherwise. Newborns contribute zero. For `6≤D<1600`, gain costs half the weight at rate at least `L`, and promotion adds three times the weight at rate `L/D≤L/6`; their sum is nonpositive. Above the cutoff gain costs three quarters of the weight and promotion obeys the same bound. At a bare vertex above or at the cutoff, discount contributes `−μ` and promotion vanishes. The only positive source is the degree-5 seam; checking `c=0,…,5` gives the maximum `40·2^-1600` at `c=5`. Summing and applying the same stationary argument gives the upper tail bound `(40/μ_lo)·2^-1600 < 10^-479`, added only to observable 8's upper side. The exact arithmetic check also exercises 38,616 endpoint-rate controls around the seams and at high degrees.
