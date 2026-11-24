**HOLDS — Astra Claim 27 at `37e7bf6`.** I found no mathematical gap in the all-positive-power ordering or either quantitative lower bound. This is a hostile check **inside the GPT family**. The full audit, code, and logs are in the child folder (programs archived as `round2/rival-attacks/code/astra27-child-*`).

### Audit

| Step | Finding |
|---|---|
| Yule amplitude and normalization | **Follows.** At size $`k`$, the next holding time is $`\mathrm{Exp}(k)`$, independent of the winning vertex and its neighbour choice. Restarted at the common completion of size $`m=n+2`$, the total-size amplitude is $`\mathrm{Gamma}(m,1)`$, independent of the future embedded chain conditional on the completed tree. Thus $`W=SZ`$, and the signed normalized numerator is multiplied by $`\Gamma(1+p)/[n(n+1)(n+2)_p]=\Gamma(1+p)\Gamma(n)/\Gamma(n+p+2)`$. Later clock races use the discount $`e^{-p\tau}`$ and the actual continuation. |
| Signed placement and all-tree ledger | **Follows.** Reversing distinct old targets gives exactly $`s_v-s_u`$ on the oriented leaf difference. Repeated targets and nonedges cancel. I checked the stock shifts, the host-race rates, the degree-$`\ge15`$ payer condition for transfers, both chain-end bills, and the retained-tree identity yielding $`b+2a`$. The physical rate bounds justify the analytic degree tails. |
| Fractional boxes and $`p\downarrow0`$ | **Follows.** The 473 and 39 closed $`1/512`$-cells meet exactly at $`473/512`$ and cover through 1. The stock harmonic increases in $`p`$; lower races use the high endpoint for discount and upper races use the low endpoint. The Euler-product factors and positive mixture weights are bounded separately in their correct directions. These are whole-cell enclosures. In the first cell, $`p=0`$ enters only as an enclosure endpoint: the lower discount uses $`1/512`$, and the upper discount cap at zero is one. The argument does not assign a constant-one payoff to zero mass. |
| Seed, stars, and $`p\ge1`$ | **Follows.** The seed's paired term vanishes. At a star only the centre can owe the old-edge debit; the shifted adverse leaf has $`d-1`$ siblings. The finite and analytic star degree ranges meet at $`D=96,64,32`$ in the respective exponent blocks. The closed covers $`[1,3]`$ and $`[3,7]`$ meet the analytic $`p\ge7`$ tail: the seventh integrated remainder has negative sign and is bounded by $`1/7!`$, while the locked-tip Beta moments decrease with $`p`$. |
| Founders | **Follows; both caps recomputed.** `founder.py` recovered no-send probability at least $`108/305`$, hence $`a=197/305`$, and a cutoff-18 raw leaf mean ceiling of $`15698419121404169357074051/40576642454003772313197000<.387`$. Hölder uses the indicator $`\mathbf1_{\{W>0\}}`$ and raw $`W`$. The centre lower is $`C_p=(2+p)e^{-1-p}`$; log-ratio concavity and rational endpoint floors $`2440/2167`$ and $`16000/15609`$ give raw gap at least $`b/50`$. The half-mixture at size three and its Gamma factor give the announced founder bound. Seed symmetry gives the tie. |
| Telescoping and mixtures | **Follows.** The beta recurrence telescopes after conditioning each term on the same history; a stopping size is handled on its measurable size slices. Positive finite mixtures preserve the strict gap. The release's link to `beta-tests/README.md` is absent from the archive; that discussion is not a proof input. |

### Replay and planted errors

I used a clean `git archive 37e7bf6`, core 13, `nice -n 5`, and the 3-GB limit. The commands are recorded in `replay.sh`; I ran them piecemeal.

- The four `mean-ledger/fractional/check.py --start X --stop Y` shards passed **36, 92, 172, and 173 cells**, all reporting the normalized floor `.033`.
- `fractional-middle/check.py` passed **256 serial lower cells and 39 signed cells**. `star_bound.py` passed **all 512 cells**, with exact seed, finite-star and high-star minima `406542520567/2500000000000`, `2085176012479/12000000000000`, and `265041348327659/2012097656250000`.
- `mean-ledger/real/check.py` passed **87 sitewise + 169 suppressed cells** at floor `1/40`. `mellin/check.py` passed **2,048 signed cells** with minimum `1892360505093/1000000000000000` and high-power margin `22890923/1058158080`.
- The splice check passed **nine pins and 880 telescopes**. The founder, tip and cold replays passed, including **64 arrow fixtures, 10,078 stock-shift fixtures, 512 independent star formulas, and 153 physical two-birth identities**. The old small-powers check prints `OPEN p=1/2` for its degree-only certificate; the stock-sensitive replacement covers that point.

`python3 -IB plants.py` rejected both semantic mutations in copies, with pins active: subtracting `1/4` from the first-box serial credit failed at `(Fraction(0, 1), Fraction(1, 512), 'tail')`; changing the founder mean ceiling from `.387` to `.500` failed the endpoint-ratio assertion. Both exited 1.

### Independent exact tests

```sh
taskset -c 13 nice -n 5 bash -c 'ulimit -v 3145728; timeout 590s python3 -IB exact.py'
```

```text
free trees [1, 1, 2, 3, 6, 11, 23] total 47 unfolded paths [2672, 25074, 262944]
PASS all 47 exact two-birth signed identities, 846 gap tests, zero upper-bound violations (Astra and Neon-2)
```

I enumerated both births and up to two additional births, canceled equal rooted shapes exactly, and bounded terminal moments independently. Persistent core-neighbour slots give an affine-Yule upper bound on integer moments; the unit-loss leaf-stock comparison gives lower bounds. For $`0<p<1`$, I used the leaf-promotion hub cap $`H`$, mean bounds $`[\ell, u]`$, and second-moment upper $`v`$ to enclose the moment between $`\ell^{2-p}/v^{1-p}`$ and $`H^{1-p}u^p`$, by Hölder. All fractional powers were rounded outward by integer-power comparisons.

| $`p`$ | Smallest normalized upper after two extra births | Strictly positive independent lower, out of 47 |
|---:|---:|---:|
| $`1/8`$ | .914806402 | 0 |
| $`1/2`$ | 1.060793265 | 0 |
| $`473/512`$ | 1.157307030 | 2 |
| 1 | 1.168750000 | 2 |
| 2 | 1.455667163 | 9 |
| 5 | 2.189881604 | 47 |

The integer-power bands are consistent with Neon-2 Claims 8/9. All eight size-11 marginal mean bands intersect Fable's bands, and the adjacent bands are consistent with Fable's three mean floors. For example, the independently computed founder mean gap is in `[0.177547245, 0.232424507]`.

**Limits:** The exact tests are too loose to challenge small changes to the advertised constant, and at $`p=1/8,1/2`$ they would not detect even a reversed sign. The continuum verdict rests on the proof reductions and rational replays, not on these finite tests. I did not independently reimplement every all-degree ledger, re-audit the full inherited Claim 26 release or Neon-2 proofs, or check positive-threshold stochastic dominance.
