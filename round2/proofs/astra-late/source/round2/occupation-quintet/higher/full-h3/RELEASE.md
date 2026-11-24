> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits and Claude replay of earlier freezes; original Neon-2 Gold1 witnesses and inherited stationary results. See the [release index](../../../../../README.md#claim-25) for exact scope and citations.

# Original quintet H3 cohort addendum

The [cohort proof](README.md), [sibling-port proof](../../diagnostic/h3-port/README.md), and [cold audit](../../../checks/quintet-h3-cohort-sol.md) certify these exact stationary floors in normalized potential units per physical core:

| Original witness | Three-component H3 floor |
|---|---|
| Lower | `677678669854000165705787542011/1151564514912178502778811842560000000` |
| Upper | `236418900054226499063534772363/359863910910055782118378700800000000` |

The sum pays the sibling-box chain, the disjoint direct slice, and the complementary sibling increment. The original higher host cap still applies after clipping. This addendum gives no composed exponent endpoint.

## Source authentication

Paths are relative to `round2/occupation-quintet/` unless otherwise stated. The audit checks the pinned producer inputs, the two original private payload digests, and the original quintet receipt. `screen.py` loads the inherited `round2/attacks/neon2-quintet/transfer/check.py` (Git blob `57328f8e69ab053699d40d011560114aedd5a57d`). The original source checkout is `b04cc50c1e391e0420eb7b006ac781d484c68543`.

| Source | SHA-256 |
|---|---|
| `higher/full-h3/README.md` | `fbc0663d4cb4fc03546b00b13da65cd661d3dfc68477b61757556c08e018d604` |
| `higher/full-h3/cohort.py` | `4fe6d9fe0b9faeae5aae85dffe29f2f3f8fa8ea8cb1dbdb64e17d4561ea514ab` |
| `higher/README.md` | `b6ce7c49ddd7083d76c7497baf8e7dcd9443d868b16def31fd3350f874e7e996` |
| `higher/screen.py` | `4c9ff5dc3366d33d7f349fd14b7c5db7e5632b61b1880db91f9e13c196a33504` |
| `higher/price.py` | `66b09ba51b34b82371814fea48fd66a85a74370981b587f23e3b815b9d01bb54` |
| `diagnostic/h3-port/README.md` | `a2759d0ef436c3893d6250ad6c9ed5db65c095b859ad2919a2241fc7b39e56bc` |
| `diagnostic/h3-port/price.py` | `913088ec708fcb6fc6dcd107822bfc8dfedc86b0b57898bb364febbaa509e06b` |
| `diagnostic/h3-port/measure.cpp` | `aa8abfa683047d979a2b844e4c675eefd359a1815d40cfaaa655fd1fdabe0c03` |
| `diagnostic/h3-port/measure.py` | `50b6bbbbd48150b6934183fd8b9a574ea2e23ec382f4790e601d5f75097978b5` |
| `diagnostic/h3-port/measure.sh` | `0e07d12df2d568ae2452bd3ab903682789667e238f49d9c9bff2cb9a61981bf2` |
| `diagnostic/census.cpp` | `7f2042ee18b4b7fdbf374e6cc19bc03eab40a6e7ec3fd8bab5a6b2265a85b6b1` |

## Replay

From a clean archive, mount the original read-only private inputs at `round2/.work/replay-inputs` and run:

```sh
python3 -IBu round2/checks/quintet-h3-cohort-sol.py --rerun --side lower
python3 -IBu round2/checks/quintet-h3-cohort-sol.py --rerun --side upper
```

These commands check the three exact floor components and their sum for each witness. The empirical companion `diagnostic/h3-port/measure.sh` additionally uses the existing out-of-tree `round2/outputs/numeric4-census-ct-48125000{0,1,2,3,4,5,6,7}.tsv` and `round2/outputs/numeric4-root-tree.tsv`. Its sample outputs are diagnostic and do not enter the stationary certificate.
