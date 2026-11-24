> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT, Claude, and Gemini within their recorded scopes for hubhood, with GPT-only audits of the accepted positive-threshold proof; threshold addenda use Claim 24 or pinned Fable finite transports where stated. See the [release index](../../../../../../README.md#claim-26) for exact scope and citations.

# Claim 26 proof package: all-tree pre-placement hubhood

The [completed ledger](edge-dual/README.md#closing-the-remaining-low-endpoints), including the [low-end closure](interfaces/low-end.md), proves a strict hubhood gap conditional on each finite old tree of size `n≥2`. When the old tree has at least two cores the conditional gap is at least `371/[12500 n(n+1)]`. Stars and the size-two seed have separate positive bounds. Averaging at each size and using the separate founder calculation gives `p₁=p₂>p₃>p₄>⋯` in the original edge-seeded law. This package is for threshold zero.

## Cold reviews and rounding

- [42b5 full review](../../../../../../../../reviews/cold/astra/checks/preplacement-all-sol.md) and [fixture](../../../../checks/preplacement-all-sol.py).
- [89ae independent second review](../../../../../../../../reviews/cold/astra/checks/preplacement-all-second.md) and [fixture](../../../../checks/preplacement-all-second.py).
- [Inherited discharge review and input pins](../../../../../../../../reviews/cold/astra/checks/preplacement-discharge-sol.md), with the earlier checker retained.

**89ae’s +1 tick clarification:** `U` is rounded upward to a `10^-12` grid. An arbitrary rounded row need not satisfy `U ≤ 1/(j+1)` literally: analytic tail inequalities use the unrounded physical envelope, or allow one extra tick on the rounded comparison. The finite checks use the rounded rows directly. The high-exposure tail (`1/86 < .012`), the large-degree outgoing bound and the stocked toll comparison against `L(2)>1/4` retain their slack. This does not modify the pinned proof or executable arithmetic.

## Exact input inventory

Git blob IDs (SHA-1). The first six rows are included under `round2/open/age-order/global-transport/hostile/`. The final four rows are inherited in the same archive under `round2/open/age-order/`.

| Location | Note | Checker |
|---|---|---|
| `hostile/` | `4602cfd0688772210e90a562f01e05e7726f49c3` | `b5b63d953c0ae51a6d145fa3c00452c4493c1b95` |
| `hostile/discharge/` | `e76a611435ec1244e1601ba401589a7753cd7906` | `feccb888c61fbc58cb8fd93a65c5d7ed422e9575` |
| `hostile/mixed-stock/` | `dc88f8ca792a6def8cf189129b7d4edbfa8dd178` | `741caa3c965c324418b1d7f7b71ab00885901b1a` |
| `hostile/edge-dual/` | `3bfafb6a7067b9bc115fcaceb2d6c6817c25388d` | `439d34cb9d424c274baad5c6958a728c72fde216` |
| `hostile/interfaces/` | `5e23d428a78b98a2905276f9c44be6a74e546c3b` | `f2eb33b128991029c78aecd3e55568e1aba2ecd1` |
| `hostile/interfaces/low-end` | `5eff7df4586344f09985deadd6312b8128eb0823` | `f96ca178b1b94f7d0d0bb7f96cdf300f0ff68fa1` |
| `early/` | `22e27fc0c7e4a4754217b66be5a184ee95282c93` | `5ca7998a11fc514efcd1e90e859a42d8b1bf2b8c` |
| `preplacement-transfer/` | `2d6eae9e8f0fcb157ee270d94a3a8023cb3493bd` | `fa54080065cc4157f7043ab78f5a9448baccb202` |
| `adjacent/` | `b25578d819568143f158c9c3411e0ef846c94d4a` | — |
| `global-transport/` | `cf321f9e2abe32cf992cc1e9c3f79c8632552f4a` | — |

## Archive replay

From the root of the unpacked commit, use Python 3 with assertions enabled. All scripts use standard-library inputs only.

```sh
set -e
h=round2/open/age-order/global-transport/hostile
python3 -IB "$h/check.py"
python3 -IB "$h/discharge/check.py" 9
python3 -IB "$h/discharge/check.py" --stocked 7
python3 -IB "$h/discharge/check.py" --pairs 6
python3 -IB "$h/mixed-stock/check.py"
python3 -IB "$h/interfaces/check.py"
python3 -IB "$h/interfaces/check.py" 9
python3 -IB "$h/interfaces/low-end.py"
python3 -IB "$h/edge-dual/check.py"
python3 -IB round2/checks/preplacement-discharge-sol.py
python3 -IB round2/checks/preplacement-all-sol.py
python3 -IB round2/checks/preplacement-all-second.py
```

The producer’s `FAIL local ledger only` lines are expected controls for the superseded ledger before transfers and suppression. The final exact checks must pass.
