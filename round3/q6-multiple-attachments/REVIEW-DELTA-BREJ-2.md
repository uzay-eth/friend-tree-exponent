> **Public checkpoint export note:** Edited export of `multiple-astra/brej-delta2-audit/REVIEW.md` (frozen original SHA-256 `4d8ef8c5770c47e4a81dd39b7b802936e433963fde393cb06abaa3c793da20df`). Changes are limited to local path/session/link redaction. All mathematical content, formulas, labels, errata, and scoped audit verdicts are preserved. **Line numbers cited in errata and audits refer to the frozen original source bytes listed in [`../PROVENANCE.md`](../PROVENANCE.md), not line numbers in this reformatted export.**

HOLDS — express Astra audit of the pinned Q6 rejection correction. All six maximum-degree rows, their stated seed domains, and the incorporated repairs to the earlier Astra review pass.

# Scope and evidence

Neon incorporated repairs derived in Astra's earlier review. This is Astra's subsequent, express audit of that incorporation, performed after the physical three-slot work at `ee20cc05899d685f2064f54a18bd005c3098a704`. The accepted correction is local `multiple-neon/DELTA-BREJ-2.md`, SHA-256 `4f5d8e2c39dc2a3d5179a8340653a76a0ad9e275d398a97198811526637cb64d`. Its author verifier is `q6/verify_delta2_brej.py` at `d3e7a73bad1568e94930c5fdd7c689928bf21328`, Git blob `565a1c30e455486e178ce15dd47e043843565bdc`, SHA-256 `12707502be5923c51b846752e3ca74ce89b94d9f361982b2194a4315aadb5a9d`. I verified both input hashes before reading the correction.

The earlier Astra report at `585e159c62dd8e349a4201d00491887de76cbd93` has SHA-256 `d17dc108858dcf17bea66babf8218501e7cab8cad95c38afd19bc95c9c3e8980`. I checked the relevant cross-references in the first repair, SHA-256 `b6d1b7e6e3ec046ccac7ec099c618f0262a667924c2f190e9a6494498b3487b5`. The original rejection statements serve as previously reviewed dependencies. An earlier run reported a passing replay of Neon's verifier. The separate check.py and REPLAY.txt here are my own exact calculations. The quantitative rate addendum remains outside this verdict. Human review is pending.

## Six-row verdict

At m=2, p is the single-walk candidate law, M=max p, S=∑p². For each of the four non-rejection rows, z denotes the largest coordinate of **that row's** limiting vector Z (called h for collapse). Their reviewed mass results give M→z>0 and S→∑Z²>0 a.s. by total-variation convergence. The strict upper endpoint z<1 uses the previously established dispersion of hub mass or, for A-simple, p_v≤1/2. Each degree column below is the a.s. limit of its maximum divided by n.

| Rule | Finite connected seed covered | Native degree | Projected simple degree | Verdict |
|---|---|---:|---:|---|
| A-simple | Simple, minimum degree ≥2 | 2z | 2z | HOLDS |
| A-multi | Loopless multigraph, minimum multidegree ≥2 | 2z | z | HOLDS on the stated pure-leaf theorem domain |
| B-multi | Loopless multigraph with at least one edge | 2z | 2z−z² | HOLDS, including the ordinary edge |
| B-simple-col | Simple with at least one edge | 2z−z² | 2z−z² | HOLDS, including the ordinary edge |
| B-rej-pair | Simple with at least two vertices | 0 | 0 | HOLDS, with M,S→0 |
| B-rej-seq | Simple with at least two vertices | 0 | 0 | HOLDS, with M,S→0 |

A-simple requires two distinct neighbors at every mediator. A-multi draws two incident edges with replacement, so the process itself starts from an ordinary edge. The table's minimum-degree restriction specifies the seed domain used to prove its projected coefficient by Theorem 3(1). B-multi's extended seed result supplies total-variation convergence from the ordinary edge. Under either rejection rule an ordinary-edge seed forces the first two selected labels to be the two seed vertices. Under both rules degree-two density tends to one by the earlier reviewed conclusion, and simple edge count grows by exactly two per birth.

For A-simple the conditional inclusion probability is 2p_v. For both multiplicity-retaining models the expected number of copies gained is 2p_v. With independent slots the chance of gaining a distinct neighbor is g(p_v)=2p_v−p_v², both when the copies are stored and when a duplicate is collapsed. With a common A-multi mediator x having row entry r_x(v), the projected increment has mean avg_x g(r_x(v)). Its excess above p_v is avg_x r_x(v)(1−r_x(v)), bounded uniformly in v by one quarter of the non-pure mediator fraction. That fraction tends to zero on the stated domain. Both rejection processes stay simple and their uniform degree limit follows directly from reviewed Theorem 6.

The moving-maximum argument also HOLDS. Separate off each bounded seed or birth degree, pad old-label increments by zero before birth, and center by their conditional means. Azuma and a union over the n labels bound the probability of an error exceeding ηn by 2n exp(−η²n/8). Borel–Cantelli on rational thresholds gives uniform o(n) error. The conditional mean vectors converge in sup norm to 2Z, g(Z), or Z by the preceding formulas; their Cesàro averages have the same limit. Maximizing over the active labels converges to the supremum of the full limit vector. A nonnegative unit-mass vector attains its positive maximum at a finite label. This argument covers changing maximizers.

## Repair-by-repair disposition

| Earlier issue | Express verdict on the accepted correction |
|---|---|
| Smallest separating seed | HOLDS. The P3 center has inclusion probabilities 8/9 and 14/15. K2 forces the same pair. With minimum degree ≥2, K3 has uniform masses and the four-vertex diamond is a smallest witness. |
| Sequential relative expansion | HOLDS in its corrected scope. The exact arbitrary-label additive remainder is at most (M²+SM)/[2(1−M)]. At p_i=S a general relative expansion fails; the recorded nonuniform six-vertex example and its blow-ups retain this equality and strict positive sequential excess. The factor is an inclusion statistic. The endpoint-law drift also has a neighbor-renormalization loss. |
| Strict comparison at every maximizer | HOLDS for positive nonuniform p of length ≥3, including tied maxima. Conditional on the first label J, the inclusion chance is h(i)=1, h(j)=M/(1−p_j) elsewhere. Tuple rejection reweights J by (1−p_J)/(1−S). The difference of the inclusion probabilities is Cov(h(J),p_J)/(1−S)>0. Unequal masses give strictly ordered h values; ties contribute zero to the covariance. At a maximum, the nonnegative terms p_j(M−p_j)/[2(1−M)(1−p_j)] give the displayed lower and upper bounds. |
| Theorem 6 stopping initialization | HOLDS. On Y_N>b the stop is immediate. Multiplying the stopped variable and its sum of increments by the initial measurable guard 1{Y_N≤b} justifies the bound b+Δ_max while retaining every infinite-survival outcome. The unconditional bound max(N,b+Δ_max) works as well, since Y_N≤N. |
| Theorem 7 floor, signed debit, caps | HOLDS. At degree ≥339 and zero leaf neighbors the weight is 10/49, below the finite-table floor. Using 10/49 for the global comparison covers the K₃₄₀ seed. Monotonicity and nonnegative weights give 0≤R⁻≤2w, with R⁻=0 at ℓ=0. The tight λ at δ=1/20000 is exactly 39223999/39998000, below the printed λ and 99/100. All 336 finite caps have zero increment. The earlier exact outward seam bound exceeds the stated strict threshold 0.00002619601 for the 57,288 off-diagonal entries. |
| Theorem 7 dyadic step | HOLDS. The first repair selects λ≤49/50+ε/4 and tests at 49/50+ε/2, giving summable probabilities. Monotone interpolation leaves a vanishing dyadic factor at the target 49/50+ε. Astra's ε/2, 3ε/4 pair also has both required strict gaps. Use the global floor, take the countable union over deterministic stopping starts, then a countable sequence of ε decreasing to zero. |
| Theorem 7 core degree | HOLDS. Subtract degree two for each degree-two vertex and degree one for each remaining degree-one seed from 2e_n. Writing s₁ for the latter count gives 2e_n−2(n−C_n−s₁)−s₁=2n+2C_n+O_seed(1), so the core fraction of all degree tends to 1/2. |

The prior WRONG verdict remains accurate for the original printed common coefficient. This HOLDS applies to the accepted corrected table and specified repairs. Replay the delta-only finite checks with `python3 -B round3/q6-multiple-attachments/check_delta2.py` (portable export with local-host pin stub removed). They enumerate 47 old-label cases across the six rules and compare the repaired constants using exact arithmetic. The almost-sure limits use the argument above and the previously reviewed prerequisites.
