> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits and Claude replay of earlier freezes; original Neon-2 Gold1 witnesses and inherited stationary results. See the [release index](../../../../README.md#claim-25) for exact scope and citations.

# Claim 6 quintet singleton census and sparse host bounds

Authenticated exports cover **4,750 exact singleton rows and 13 all-stock bands per witness**. The source is the detached Claim 6 checkout at `b04cc50c1e391e0420eb7b006ac781d484c68543`; only the two private golden snapshots are used. No D18/D20 reward or previous witness is an input.

## API and accounting

Import `census.py` from this directory. `Snapshot(path, 'lower'|'upper')` checks size, exact production header and full SHA-256, then uses read-only mmap for sparse reads. Use as a context manager if convenient. `verify=False` skips hashing for exploratory callers that already control authentication; the export CLI always authenticates. The source header, C++ verifier and Python tail script are checked against their Git blob pins by the CLI. `Snapshot.check_unchanged()` verifies file identity before export.

`c` means **core degree including the parent** at a nonfounder, not a child count. `L=D−c`. `f(D,c)` and `newborn_edge(D,c+1)` are unsigned potential ticks; `row_ticks(D,c)` returns exact inward `(A,B)` in ticks and `row(D,c)` divides both by `2^40`. With `rho=Σ_core_neighbor 1/D_neighbor`, the original inward singleton row is `A+B*rho`. Specifically,

```
Δf = f(D+1,c)−f(D,c)
birth = (L/D) [f(D,c+1)−f(D,c)+f(2,1)+g((D,c+1),(2,1))]
signed row = inward_sign * [rate*Δf + birth − λ f(D,c)] / 2^40
```

Centered rows use `rate=L+c*r1/S_r` and `B=0`. They transfer `(rho−c*r1/S_r)Δf` to existing `g`. Noncentered rows use `rate=L+rho`, with true `rho∈[0,c/2]`. This includes exact arity overflow below `d2`, the `d2` projection seam, and the `DF−1→DF` analytic seam. Promotion retains the old degree and uses the **post-promotion** census for the newborn edge exactly once. No `H3` newborn or FY term is charged to `f`.

Signs here are **+ lower / − upper**, the opposite of the verifier's signed residual caps. The frozen verifier caps `(f,g,higher)` are `(0,401,18785)` lower and `(-1053,946,19890)` upper in its orientation. Inward rows are at the original trials `56570237/10^8` and `56574828/10^8`. This is a row export, not an exponent update: a later stationary proof must debit disjoint residual ceilings and compose with its own comparison, or replace the original `f` price once. Do not add these rows to the frozen corrected endpoints as a second price for `f`.

## All-stock floor

`tail_floor(c, degree=None)` bounds both affine receive endpoints on the entire tail at fixed `c`. Set `x=1/D`, `q=L/D=1−cx`. Each exact normalized residual multiplied by `(1+x)^2` is a polynomial in `x`; use the exact `c+1` newborn-edge endpoint, rather than taking an extremum over all newborn arities. Convert to one-variable Bernstein on `[0,1/max(DF,c,degree)]`; the minimum coefficient divides by the largest denominator if nonnegative and by one otherwise. The replay checks both polynomial identities against direct rational rows at three integer degrees for every band. `stock_floor(c)` also scans every finite leaf stock below `DF=2000` at that fixed arity. All exported bands `c=0..12` are strictly positive. For example:

| Arity | Lower all-stock floor | Upper all-stock floor |
|---:|---:|---:|
| 1 | `1.5841787562e−9` | `3.5753279865e−9` |
| 2 | `4.2363564121e−10` | `1.7724070144e−9` |

This analytically covers unbounded **degree/leaf stock** at fixed arity. It does not certify every finite arity with an all-stock scan. The sparse reader can evaluate an individual row at any legal degree and arity.

## Host overlay

`HostCaps(snapshot, gold_receipt)` supplies `higher_ticks(D,c)` and `comparison_ticks(D,c)` for the parent's denominator lane. `higher_ticks` is a valid outward cap, not a fresh blockwise extraction: it uses exact state eligibility and combinatorial multiplicity under the degree-sum masks, times the already verified maxima for each family, clipped by the global host cap. The `P4/P5/B5` multiplicities were independently checked by explicit admitted-partner enumeration. Ineligible hosts have zero higher cap. No giant motif slabs are read to produce this overlay.

The comparison is the lower (recharged) and upper (unrecharged) bound on the verifier's local potential allocation. The lower and upper **need not be ordered at an individual host**, because only the lower includes `alpha*(c−2)`. In the one-core `c=0` case neither an edge nor recharge is present. Column extrema over `g` are exact sparse reads with zero included for the root; convexly combine column extrema in the tail. On higher hosts the optional comparison uses the conservative allowed weight bound times exact multiplicities, not a new reduction over the giant slabs. It is most useful off that support. The frozen global comparisons remain separate.

## Replay evidence

From Astra root:

```sh
R=round2/.work/replay-inputs
C=round2/occupation-quintet/census
for side in lower upper; do
  taskset -c 5 nice -n 19 timeout 65s /usr/bin/time -f 'wall=%e rss=%M exit=%x' \
    python3 -B "$C/replay.py" "$R/quintet_$side.bin" "$side" --source "$R/source" \
    --receipt "$R/gold.json" --compare "$C/$side-host.json"
done
```

The original exports used `--out` in place of `--compare`. For plain `{lower,upper}.json` omit `--receipt`. Plain export outputs: `rows 4750`, 13 stock bands, `wall=4.03 rss=88872 exit=0` lower and `wall=5.42 rss=88868 exit=0` upper. A later lower plain `--compare` reproduced the exact export (`wall=10.69 rss=90956 exit=0`). The final host exports reported cap ranges `0..18785` and `0..19890`, `wall=14.48 rss=117656 exit=0` lower and `wall=21.15 rss=117784 exit=0` upper. Each passed 11 exact physical cases and the fixed-arity polynomial checks. For wrong centering at the genuinely noncentered `(27,27)` overflow census, a planted frozen rate `c/4` instead of actual exterior endpoint `c/2` changes the inward row by `−4763483761203/2` ticks lower and `+540335198115/2` ticks upper. The source/audit checkout and golden binaries were left untouched. No full verifier scan was needed.
