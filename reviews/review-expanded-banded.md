# Cold review of the expanded banded upper bound

Reviewer: Claude Opus 5.5 (Claude family). The proof was written by GPT-6.2 Astra (GPT family).

Object: [`proofs/upper-bound-banded.md`](../proofs/upper-bound-banded.md), the expanded text with the sections "The generator on a labelled history", the displayed local inequality and the derivation of $`Gq`$. I read it against its sources: [`bounds-transfer.md`](../proofs/bounds-transfer.md), [`bounds-hub.md`](../proofs/bounds-hub.md), and the integer verifier [`verify/banded-upper/check.cpp`](../verify/banded-upper/check.cpp). I did not edit the proof.

**Overall verdict: CORRECT.** Every step follows from its inputs. There are four minor changes, listed at the end, and none of them affects the bound.

Verdicts: **SOUND** means the step follows as written. **MINOR** means the step is right but its wording or a reference needs a small fix.

## Step verdicts

1. Target $`\lambda = 5657/9983`$ and $`17/30 - \lambda = 1/299490`$: since $`17 \cdot 9983 - 30 \cdot 5657 = 1`$ and $`30 \cdot 9983 = 299490`$. SOUND.
2. Reproduction command and the SHA-256 of `recipe.txt.gz`: the digest `584bb93c…0f60a` matches the shipped file. SOUND.
3. Band rows, with missing cutoffs equal to 2, the first band starting at $`D = 2`$ and the last ending at $`H = 2400`$: matches `Band` and the `geometry.txt` rows, `check.cpp` 151–158. SOUND.
4. Cutoffs nonincreasing in $`c`$ and across bands, the cutoff 620 at $`c = 1`$, and core-degree caps 20, 7, 5: matches `geometry.txt` and is enforced at `check.cpp` 155, 226 (with the caps at 89). SOUND.
5. Exact labels for $`d < h_c`$ and the witness condition at $`c = 1`$ ($`d \le 14`$, $`2 \le k \le 4`$, exact $`(e, j)`$ for $`e < 6`$): matches `witness` and `label`, 107–121. SOUND.
6. Witness kept on a near gain and replaced by the fresh $`(2,1)`$ on a near conversion: matches 250–251. SOUND.
7. Witness at the focal vertex's birth taken as the minimum label from the creator's old core neighbours: the creator has 1 to 3 old core neighbours and degree at most 14. Band 0 then has cutoffs $`620, 240, 45 \ge 6`$, so $`e < 6`$ is exact, including $`(5,5)`$ under cap 5. Matches `newborn`, 262–270. SOUND.
8. Labels defined along the history, ties by earliest birth, drift in the history filtration: well defined, since no tie changes the recorded data. SOUND.
9. Pool marks ($`h_c`$, powers of two above $`h_c`$ from 4, and $`H`$), only the largest mark kept, and each mark a lower bound on the true degree: matches 102–104 and `index` at 137–141. Marks stay lower bounds because degrees never decrease and both resets and projections move to marks no larger than the old bound. SOUND.
10. Projection on a central conversion or band crossing, the witness kept on a crossing at $`c = 1`$, and the new $`(2,1)`$ inserted: matches 229–239. SOUND.
11. Central rates $`L + \sum 1/d_w`$ and $`L/D`$. Near rates $`(d-k)/d`$ and gain in $`[d-k+1/D,\ d-k+1/D+(k-1)/2]`$. Saturated near rates $`[1/D,\ d-(k+1)/2+1/D]`$. Witness rates $`[e-j+1/d,\ e-j+1/d+(j-1)/2]`$ and $`(e-j)/e`$. Checked against the integer endpoints at 304–319. SOUND.
12. Exhaustion of local transitions (own $`(L, c)`$, one near type, the far witness; pooled gains change nothing; the witness differs from the centre in a tree): SOUND.
13. Generator on a labelled history: $`a_x`$, $`b_x`$; ties between converted leaves do not matter; at $`T_3`$, $`a = 2`$ and $`b = 1`$; total rate $`n`$. SOUND.
14. Identity for $`n\,\mathbb E[\Delta\Phi \mid \mathcal F_n] - \lambda\Phi`$ with $`I_u = \{u\} \cup S_u \cup J_u`$ (and $`\{u\}`$ above the cutoff): every pair of an event and a changed weight is counted once, and newborns go to their creator. SOUND.
15. Signed cancellation: $`g(d_w,k_w;D,c)`$ is added at the focal gain and $`g(D,c;d_w,k_w)`$ is subtracted at the near gain, both with the rate of the same vertex. Witness variants share a key, and keys with an endpoint of degree at least $`H`$ are zero. Matches 196–201, 280–281, 288, 311. MINOR (change 2).
16. Arrow cancellation and its bounds: exact $`(U - a\mathbf 1_{d_w \ge T})/d_w`$; pooled $`\max(0, U - a\mathbf 1_{b \ge T})/b`$, which is valid for $`b \ge T`$ and for $`b < T`$ because $`a \ge 0`$; sender charge at the maximum over consistent types. Matches 212–221, 243, 293–294, 323. SOUND.
17. Displayed amended row, and the equality of the $`g_{wu}`$ coefficients $`L_u + \rho_u = a_u`$ at both ends: SOUND.
18. Replacing the newborn by its trial reward, sign-selected rate endpoints, and multiplicities: valid termwise even for dependent rates. SOUND.
19. Rate intervals with capped core degree (the upper expression decreases in the true $`k'`$), with a witness ($`k`$ exact, $`k-2`$ unknown terms), and for the far witness: SOUND.
20. Boundary: $`q = L/D`$ rises by $`(1-q)/(D+1)`$ on a gain and falls by $`1/D`$ on a conversion. So $`Gq = (L+\rho)(1-q)/(D+1) - q/D \le 1-q`$, since $`L + \rho \le D`$. SOUND.
21. Tail drift is at most $`(B + a_*)(1-q) + sq`$, with arrow charge $`a_* c/D = a_*(1-q)`$ and no signed transfer. The affine bound in $`q`$ reduces to $`B + a_* \le \lambda A`$ at $`q = 0`$ and $`s \le \lambda(A+B)`$ at $`q = 1`$. Checked in integers at 202–205. The rounded-up seam value $`A + \lceil B(L+1)/H \rceil`$ is safe because the rate is nonnegative. SOUND.
22. Order of construction (decreasing $`D`$, then $`c`$, then reverse lexicographic bag order) puts every non-newborn successor first: SOUND, and the checker enforces it in `get`, 142–146.
23. The row test $`\sum \lceil\cdot\rceil + \text{fee} \le \lfloor \lambda z \rfloor`$ and the monotonicity in $`z`$: matches 338–350. The floating proposal only picks the first trial. SOUND.
24. Actual newborn weights must be at most the trial rewards: checked at 369–371, with minimum reserve 4,340,525. SOUND.
25. Explicit form $`\lceil \max\{r_-(v-z), r_+(v-z)\} \rceil`$ and the targets of each kind of term: matches 286–323. The function `up` rounds negative values toward $`+\infty`$ correctly. SOUND.
26. Near and far transitions strictly advance the sorted bag: type indices increase with $`d`$, then $`k`$, then $`(e, j)`$, and pool entry is at the first mark. Reprojection happens only on central events. SOUND.
27. Conclusion: $`mC_n \le \Phi_n`$ and $`\mathbb E C_n \le (\Phi_3/m)\prod_{j=3}^{n-1}(1+\lambda/j)`$. SOUND.
28. Reported numbers: root $`1539738114735`$, $`m = 3908512529`$, 226,545,617 rows, 63,309 transfers, slack 3, `PASS upper 5657/9983`. These match my run (below). SOUND.
29. Optional strict improvement with reserve $`e`$: $`\Phi + eC`$ has drift at most $`\lambda\Phi/n`$, and $`\Phi \le MC`$ gives $`\lambda M/(M+e)`$. This uses $`\mathbb E[\Delta C \mid \mathcal F_n] = S_n/n`$ and the reserve in the tail inequality too. SOUND (it is stated as optional and not used numerically).
30. Scale $`10^{-12}`$ and homogeneity: SOUND.
31. Correspondence table: the content of every row is right. Some line ranges are off by a few lines. MINOR (change 1).

## Required changes

1. Correspondence table. `up` is at lines 28–31, not 26–28; the anchor `#L26` lands on `narrow`. The newborn reward lookup is `solve` 296–299, not 298–301. `Recipe` spans 32–69. The other ranges are correct for the current `check.cpp`.
2. After "so both lookups exist", add the reason the exact labels are current. A neighbour that is exact and unsaturated now has never been pooled or saturated. Degrees and core degrees never decrease, and cutoffs and caps never increase along a vertex's history. So both rows read the true current types.
3. Optional: give the arrow threshold from the recipe, $`T = 8`$ (18 arrow entries).
4. Optional: in the correspondence table, note that the compile-time hooks `FRIEND_TREE_OBSERVE_ROW` and `FRIEND_TREE_NO_MAIN` (lines 352–354 and 379–385) are inactive in the production build. Only `verify/tree-drift/` defines them.

## Evidence

```text
$ sha256sum verify/banded-upper/recipe.txt.gz
584bb93c2797138ae3de54cb85c7d8e487bb7d7cc5320c021b5a64ec0cf0f60a  verify/banded-upper/recipe.txt.gz
$ g++ -std=c++17 -O2 -fsanitize=undefined -fno-sanitize-recover=all verify/banded-upper/check.cpp -o check_banded
$ ./check_banded verify/banded-upper/geometry.txt recipe.txt 5657 9983
verified rows 226545617, signed transfers 63309, minimum integer drift slack 3
scale 1000000000000, minimum newborn reserve 4340525, weight range [3908512529, 1764716280715], centre 1539738114735
PASS upper 5657/9983
exit status 0, real 85.37 s
```

The recipe header has `scale 1000000000000` and `threshold 8`, and 18 `arrow` records.
