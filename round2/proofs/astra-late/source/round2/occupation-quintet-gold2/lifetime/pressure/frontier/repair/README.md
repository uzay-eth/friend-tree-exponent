> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical-floor and frontier audits and scoped Claude conditional reduction and arithmetic; original Neon-2 Gold2 witnesses and the four reward floors. See the [release index](../../../../../../../README.md#claim-28) for exact scope and citations.

# Corrected Gold2 lower trial at `99/175`

The old H3-dependent frontier and Claim 10 endpoint compositions are **on HOLD**. This replay uses only the corrected and disjoint original lower-witness credits at `99/175`:

- D18 and both FY sectors: `4478371331693831604075238101/10000000000000000000000` ticks/core (pinned `host/conditioned24.json`).
- Corrected H3 sum: `919076249761939897477242027/65029541606195200000000` ticks/core, from the two disjoint corrected masks in `round2/checks/gold2-lower-h3-alt-sol.md`. This is the only H3 payment in the objective.

For `delta=99/175−lambda`, original `M0=m_max`, higher cap `H`, local host comparisons `U,h`, and marked singleton continuation `V`, the exact repaired dual lower-prices `pi[b V(D,c+1) + H−h + delta*(M0−U)]`. The trial margin adds the corrected edge and H3 payments once, and subtracts original `E+delta*M0` once. The mark value replaces the old singleton; it does not add to any prior singleton floor or its residual. All trials use the original exact Gold2 prior. The Claim 10 column restrictions are conditional pending independent Astra audit by 863c and 2f74.

| Replay | Gold2-only remaining credit (ticks/core) | Conditional Claim 10 remaining credit (ticks/core) | Claim 10 joint price increment (ticks/core) |
|---|---:|---:|---:|
| `standard` | 30884.736359813 | 20737.619337905 | 10147.117021908 |
| `D1100-M5-T3-C12` | 18493.560217849 | 8302.968462437 | 10190.591755411 |
| `D1750-M5-T3-C12` | **14384.069081350** | **4179.179585182** | 10204.889496168 |

A positive number in either gap column is a sufficient size of a new disjoint credit (take a strict excess), not a lower bound on the actual unknown physical reward. The higher-degree solves use the same actual marked degree/stock at maternal birth and the two real newest children; the eleven root values are combined pointwise by maximum with the original valid marks and then substituted once. D1750 certifies an improvement over standard of at least 16500.667278463 / 16558.439752723 ticks/core in the joint trial price. No degree-2/3 reward is added independently of these roots. The corrected ledger alone still does not exclude `99/175`.

Each JSON receipt records exact roots, trial margins, grid integer duals, and one rational repair per column class. Dual reconstruction uses zero for absent rows and classes. If `n` is the integer dual, each exact column satisfies `dot(n,A_i) <= GRID*score_i/SCALE + defect[kind_i]`. Every kind has total mass at most `1/LO` under the original physical lift, so the repaired stationary joint floor is `SCALE*(LO*n_birth−sum(defects))/GRID`. The reference Gold2 cells and root input are pinned in `replay.py`, and the optional marginal overlay uses the original Claim 10 pin. The underlying SciPy solution is only a proposal; every listed column is checked exactly.

Replay from `.`:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 8 nice -n 17 timeout 60s round2/.work/venv/bin/python -Bu round2/occupation-quintet-gold2/lifetime/pressure/frontier/repair/replay.py --export
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 8 nice -n 17 timeout 165s round2/.work/venv/bin/python -Bu round2/occupation-quintet-gold2/lifetime/pressure/frontier/repair/replay.py --probe --live-cut 1100 --export
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 8 nice -n 17 timeout 480s round2/.work/venv/bin/python -Bu round2/occupation-quintet-gold2/lifetime/pressure/frontier/repair/replay.py --probe --live-cut 1750 --export
```

All three exited zero. Receipt SHA-256: [`standard.json`](standard.json) `14485fe7b373042a0ef5b11fddf33652c7e0edf58c97b78026d53ac81519514c`; [`D1100-M5-T3-C12.json`](D1100-M5-T3-C12.json) `8bbb894e9f02ea7b2ded8f5443f7f040be7573d0b3a62d4a3fab4b177c048476`; [`D1750-M5-T3-C12.json`](D1750-M5-T3-C12.json) `5ffd94a84339edcaf48984aa3de19e5bbd79decda38bfacfb93db0a7871c0a78`.
