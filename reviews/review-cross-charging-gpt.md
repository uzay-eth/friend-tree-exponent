VERDICT: CORRECT

The stated bounds on `E[C_n]` follow for every `n >= 3`. The two rational tables satisfy the drift inequalities at the stated exponents, and the proposed constants round in the safe directions. Two zero-rate conventions in the row formula should be made explicit. The general formulation of Lemma 3 also needs domain qualifications, which the supplied lower certificate satisfies.

## 1. Proof audit

Line numbers refer to `xcharge.md`. The audit uses the written rates, state space, and inequalities.

| Lines | Status | Finding |
| --- | --- | --- |
| 1–9 | VALID | The subsequent cancellation is an identity on each tree for arbitrary multipliers. The method that found them plays no part in the proof. |
| 11–23 | VALID | Every chosen recipient W is core or a leaf of a core vertex. For a core recipient w, a chosen leaf of w contributes 1 to the scaled rate, and a chosen core neighbour z contributes `1/D_z`. For a leaf recipient, its sole possible source is its core neighbour w. Each such leaf contributes `1/D_w`, giving total `L_w/D_w` at w. A gain adds one leaf to w. A conversion replaces a leaf neighbour of w by a core neighbour at unchanged degree and creates a degree-2 vertex with one leaf and one core neighbour. The newborn has non-TOP class `(1,1)`, which meets the degree and core-neighbour bounds. Existing degrees and core-neighbour counts are nondecreasing. |
| 25–47 | VALID | The detailed, coarse, and beyond states partition the core. Existing TOP vertices remain TOP because either their degree or their core-neighbour count crossed a nondecreasing threshold. Each TOP vertex has degree at least `d_T = min(Dcap,c1+1) >= 2`. For `L >= Lmax`, the degree alone meets the TOP threshold whenever a class is needed. A vertex with `c=0` has no core edge and needs no class. Capping L in f leaves a finite range. |
| 49–70 | NEEDS REPAIR | Read the detailed row as omitting the self-conversion term when `L=0`, the TOP-credit term when S contains no TOP, and both class-dependent payments and receipts when `c=0`. Otherwise the displayed formula can evaluate an undefined `f(-1,...)`, an undefined `S - TOP`, or an undefined `b=cls(u)` even at zero rate. The event discussion in §§1 and 3 gives these omission rules. The neighbour conversion term already has an explicit zero-rate convention. When a conversion leaves the detailed range, the f lookup discards S in the coarse case and uses G in the beyond case, as defined in §2. |
| 72–79 | NEEDS REPAIR | The coarse and beyond row is correct with the same omission of the self-conversion term at `L=0`. A coarse or beyond vertex records no neighbour classes, so there is no state change in its row from a neighbour's gain or conversion. The exact credit decomposition uses total core credit `rho`. A beyond conversion keeps its own value at `f(G)` and its newborn sees TOP. |
| 81–95 | VALID | With zero multipliers, a gain changes the recipient's own L and only the recorded classes of that recipient at its detailed core neighbours. A conversion changes the creator, contributes the newborn once, and can update the creator's class only at detailed existing core neighbours. Existing neighbours retain their own pairs `(L,c)`. The gain rate of a non-TOP neighbour w of u splits into `l + 1/D_u + rho'(w ∖ u)`, and the conversion rate is `l/(l+k)`. The creator's own conversion row already includes its newborn. |
| 96–100 | VALID | Each ordered core edge `w → u` out of a TOP w contributes `+kappa(cls(u))/D_w` at w and the same term with a minus sign at u in `rho_top(u)`. A detailed TOP giver uses its multiset S and a coarse or beyond giver uses its actual neighbour classes. A detailed TOP, coarse, or beyond receiver subtracts using `kappa(TOP)`. Both orientations of a TOP–TOP edge are accounted for individually. |
| 102–106 | VALID | Fix w of class `a ∈ A2`. For every ordered pair of distinct core neighbours `z,u` of w, the charge is `nu(a,cls(u))/D_z`. The row of u subtracts it as part of the rate for w gaining from neighbours other than u. The row of w adds it within the credit from z after removing one occurrence of `cls(z)` from `S_w`. The argument counts multiplicities even when different neighbours share a class. If z is TOP, the addition is in the TOP-credit term and removal still excludes the source occurrence. If u is TOP, the charge is zero, so a coarse or beyond u needs no subtraction and a detailed TOP u subtracts zero. A w outside A2 creates no nu charge by definition. |
| 108–109 | VALID | The identity follows from the preceding calculation. The author's stated software check was outside my evidence. |
| 111–119 | VALID | A TOP neighbour contributes a reciprocal degree in `(0,1/d_T]`, giving the closed interval `[0,n_top/d_T]`. A specific neighbour of class `(l,k)` has exactly `k-1` other core neighbours, each with degree at least 2. This gives `[0,m_a(k-1)/2]` for each class sum and zero at `k=1`. Identical classes have identical coefficients in the detailed row. For coarse and beyond states, the same TOP bound and ordinary core bound give `rho-rho_top <= (c-n_top)/2`. Each of the c payment terms lies between the class minimum and maximum. These bounds need no independence or simultaneous attainment. |
| 123–131, 143–147 | VALID | For a scalar in `[0,hi]` times e, the maximum is `max(0,hi*e)`. In a coarse row the two credit bounds are weighted by `c-n_top` and `n_top`. On dividing those weights by c, the result is a convex combination of `(c/2)*max(0,Df)` and `(c/d_T)*max(0,Df-kappa(TOP))`, bounded by the stated three-way maximum. The coarse payment bound follows from every multiplier being at most Kp. Detailed payments are exact signed terms. |
| 133–136, 147, 152 | VALID | For every actual `L>Lmax`, both L and `L-1` index the terminal f block, including the first tail point `L=Lmax+1`. The creator's own class before and after conversion is TOP if used. Thus every f difference and uncertainty coefficient is fixed within a family, the self-gain vanishes, the conversion rate is `1-c*y`, and the remaining rate factors are constants or multiples of y. In a coarse family `Df=0` and the upper credit maximum vanishes because `kappa(TOP)>=0`. The lower coarse minimum vanishes because `Df=0`. The right-hand side is constant. Affinity therefore carries the two endpoint inequalities to every interior y, with y=0 interpreted as the limiting endpoint. |
| 137, 147–149 | VALID | In a beyond row, the term `-kappa(TOP)*rho_top` is nonpositive. The remaining upper bound is `(L/D)*f(newborn TOP) + (c/D)*Kp`, a convex combination. The two beyond inequalities bound it by `lambda*G`. |
| 139–149 | VALID | Sum the upper row bounds and apply the drift identity, then take expectations to get `E[Phi_(n+1)] <= (1+lambda/n) E[Phi_n]`. The assumptions imply `lambda>0` by the newborn beyond inequality. Start the iteration from the deterministic three-vertex path with `Phi_3=f(2,0,{})`. The pointwise inequality `C_n<=Phi_n` uses `f>=1` also at G. |
| 151–156, 161–163 | VALID | The minimum of an interval term is `min(0,hi*e)`. The lower coarse credit is `rho*Df` because `kappa(TOP)=0`, and `rho<=c/2` gives the stated minimum. Coarse payments are nonnegative since every kappa is nonnegative. A beyond vertex has no own change or TOP debit, and its newborn and payments are nonnegative, giving a row of at least `delta*f(G)=0`. |
| 158–159 | NEEDS REPAIR | The general lower statement should specify `delta>=0` for the intended growth regime and `sup f>0` for the displayed quotient. The all-zero lower certificate otherwise gives `0/0`. If arbitrary real deltas are read into the definition, the product conclusion can be false. Take f equal to 1 on the box and tail and zero at G, with all multipliers zero and `delta=-12`. Every lower row bound is nonnegative, but the proposed product gives 6 at n=5, where every five-vertex tree has at most 3 core vertices. At `delta=11/20` with the supplied positive supremum, the step is valid. The expectation recursion can be iterated with positive factors and `Phi_n<=C_n*sup f`. |
| 165–171 | VALID | LP provenance is immaterial to exact inequalities. The raw text tables satisfy the sign constraints without clipping. The checks below use independent code. |
| 173–182 | VALID at the supplied delta | The cited almost-sure theorem asks for `0<delta<=1`, as holds here. Finiteness and the lower drift give its first two conditions. At most `Dcap-1` existing neighbours can see a changed non-TOP class, plus the recipient, with at most one newborn. Positivity on the finite list of newborn states gives the required positive minimum. For a vertex with positive f, one has `c<=cmax`. Its own gain rate when `L<Lmax` is bounded by `Lmax+cmax/2`, its conversion rate by 1, and the relevant per-neighbour gain and conversion rates by `Dcap` and 1. Adding gives the stated R. If this section is used as a general assertion, state the cited theorem's exponent range. The expectation conclusion is independent of this section. |

## 2. Independent exact certificate check

`certificate.py` decodes each hexadecimal token by applying a power of 2 to its integer significand. `verify_rows.py` independently generates every allowed multiset and every box pair from the degree and core constraints, and checks that this universe is exactly the set of table keys. The c=0 states start at L=2, c=1 at L=1, and other states at L=0 because a core vertex has degree at least 2. There are 93 tail families, each tested at both endpoints. Class indices in the evidence use the text-file ordering. A final field of `None` means coarse and an empty tuple means no core neighbours.

At the claimed exponents each table passed 5,665 box and 186 tail inequalities. The upper table passed its two explicit beyond conditions. For the lower table I also checked the two nonnegative coefficients of its beyond row separately, reusing the minimum kappa already checked by the sign conditions. The smallest upper drift slack was `11258999068417/112589990684262400`. The smallest lower drift slack was exactly zero. The combined sign and drift minimum is zero for each table. Kp was recomputed, and the lower equality conditions were checked exactly in addition to the displayed inequality counts.

Reducing the upper exponent to `579/1000` produced 5,668 failed drift checks. Increasing the lower exponent to `552/1000` produced 5,747. Changing only the upper `f(2,0,{})` to 1 in memory preserved the sign constraints and broke its box inequality.

Command and output from `rows.log`

```text
$ taskset -c 0 python3 -B checks/fable-xcharge-gpt/verify_rows.py | tee checks/fable-xcharge-gpt/rows.log
Class indices: {0: (1, 1), 1: (0, 2), 2: (2, 1), 3: (1, 2), 4: (0, 3)} TOP: 5
certificate: upper, exponent=29/50
  box    count=5665 failures=   0 min=11258999068417/112589990684262400 (9.99999999999e-05) at (0, 2, (5, 5))
  tail   count= 186 failures=   0 min=11622735139569/14073748835532800 (0.000825845002308) at (61, 3, (0, 0, 4), '0')
  beyond count=   2 failures=   0 min=929818811201/1125899906842624 (0.00082584500234) at newborn(TOP)
  sign   count=5667 failures=   0 min=0 (0) at ('f', (0, 2, (5, 5)))
  drift tests=5853, all inequalities=11520, PASS
certificate: lower, exponent=11/20
  box    count=5665 failures=   0 min=0 (0) at (0, 2, (5, 5))
  tail   count= 186 failures=   0 min=15996238495629/18014398509481984 (0.000887969614262) at (61, 0, (), '0')
  beyond count=   2 failures=   0 min=0 (0) at minimum kappa
  sign   count=5672 failures=   0 min=0 (0) at ('f', (0, 2, (5, 5)))
  drift tests=5853, all inequalities=11525, PASS
exponent control: upper, exponent=579/1000
  box    count=5665 failures=5573 min=-836833706109201/112589990684262400 (-0.00743257638644) at (59, 0, ())
  tail   count= 186 failures=  94 min=-1046046162604753/140737488355328000 (-0.00743260502109) at (61, 3, (0, 0, 4), '0')
  beyond count=   2 failures=   1 min=-167367386016051/22517998136852480 (-0.00743260502105) at newborn(TOP)
  sign   count=5667 failures=   0 min=0 (0) at ('f', (0, 2, (5, 5)))
  drift tests=5853, all inequalities=11520, FAIL
exponent control: lower, exponent=69/125
  box    count=5665 failures=5654 min=-1868170698331/703687441776640 (-0.00265483023772) at (60, 0, ())
  tail   count= 186 failures=  93 min=-59713318342587/22517998136852480 (-0.00265180403603) at (61, 0, (), '0')
  beyond count=   2 failures=   0 min=0 (0) at minimum kappa
  sign   count=5672 failures=   0 min=0 (0) at ('f', (0, 2, (5, 5)))
  drift tests=5853, all inequalities=11525, FAIL
In-memory perturbation: upper f(2, 0, ()) 4064395905426029/562949953421312 -> 1
table control: upper, exponent=29/50
  box    count=5665 failures=   1 min=-19585578232166161/879609302220800 (-22.2662245416) at (2, 0, ())
  tail   count= 186 failures=   0 min=11622735139569/14073748835532800 (0.000825845002308) at (61, 3, (0, 0, 4), '0')
  beyond count=   2 failures=   0 min=929818811201/1125899906842624 (0.00082584500234) at newborn(TOP)
  sign   count=5667 failures=   0 min=0 (0) at ('f', (2, 0, ()))
  drift tests=5853, all inequalities=11520, FAIL
```

## 3. Direct drift by ordered moves

`direct_drift.py` rebuilds the core, classes, states, and potentials on each successor tree and then enumerates each ordered `(V,W)` with weight `1/(n*degree(V))`. Successor potentials are reused for moves that have the same W. The sums of potential values use exact integers in units of a common dyadic denominator, and all probability and slack operations use fractions. The direct program uses the text parser but imports no row equation or multiplier calculation.

The unlabelled generator adds a leaf to every vertex of each previous tree and deduplicates by the least rooted encoding over all roots. Every tree has a leaf whose removal gives a predecessor. Sorting the rooted child encodings recursively characterizes the rooted isomorphism class, so taking the least encoding over roots gives the claimed unrooted deduplication. The one-vertex tree has no move. At n=2 the upper inequality has zero right-hand side and positive drift, as expected outside the stated n≥3 domain.

In each applicable family the relative slack divides the signed difference by `exponent*Phi(T)`. The numerator is target minus drift in upper mode and drift minus target in lower mode. The exact minimum for every family, with its witness, is in the output. The boundary tests exercise detailed, coarse, and beyond vertices with more than 60 leaves, neighbouring beyond vertices, and conversions that can cross from 12 to 13 core neighbours.

Command and output from `direct.log`

```text
$ taskset -c 0 python3 -B -u checks/fable-xcharge-gpt/direct_drift.py | tee checks/fable-xcharge-gpt/direct.log
Domain n>=3. n=1 has no friend move. Outside domain at n=2: Phi=('0', '0'), nE[Delta]=('4064395905426029/281474976710656', '7031833084929469/2251799813685248')
unlabelled counts n=1..11: [1, 1, 1, 2, 3, 6, 11, 23, 47, 106, 235]
all unlabelled: trees=434, n=3..11, ordered moves=7914, max L=10, max c=5, violations=0
  upper: min relative=100204372259656/581185360279379281 (0.000172413792755) at n=8/tree=13; zero targets=0
  lower: min relative=49944344079853/54938778488190300 (0.000909090909085) at n=8/tree=13; zero targets=0
Random seed=992623
random friend: trees=320, n=3..400, ordered moves=122922, max L=376, max c=12, violations=0
  upper: min relative=46671163043097/270692745944453047 (0.000172413792916) at #91/n=7; zero targets=0
  lower: min relative=79632879644297/87596167609367923 (0.000909090909084) at #135/n=25; zero targets=0
stars: trees=11, n=3..400, ordered moves=1694, max L=399, max c=0, violations=0
  upper: min relative=11622621764571/67411206236733896 (0.000172413793098) at n=59; zero targets=0
  lower: min relative=19927152947299/21919868242265246 (0.000909090909081) at n=60; zero targets=0
paths: trees=11, n=3..400, ordered moves=1694, max L=2, max c=2, violations=0
  upper: min relative=4216691603497/24456811300492647 (0.000172413793102) at n=4; zero targets=0
  lower: min relative=50076244053936/55083868459376999 (0.00090909090909) at n=4; zero targets=0
brooms: trees=6, n=14..351, ordered moves=1502, max L=200, max c=2, violations=0
  upper: min relative=3581184932119255126/5780778078293629350201 (0.000619498773974) at spine=150/leaves=200; zero targets=0
  lower: min relative=74911696679340099/6838277021696144926 (0.0109547619147) at spine=3/leaves=10; zero targets=0
caterpillars with bare interior: trees=4, n=11..390, ordered moves=1156, max L=180, max c=2, violations=0
  upper: min relative=11953087564861156868/6415500680221818231593 (0.00186315739966) at spine=80/left=130/right=180; zero targets=0
  lower: min relative=353684670591674423/11263785050927787807 (0.0314001615791) at spine=9/left=10/right=11; zero targets=0
high core degree: trees=6, n=25..242, ordered moves=1004, max L=61, max c=120, violations=0
  upper: min relative=2221251791761735669/230309882309742374359 (0.00964462214771) at arms=120/leaves=1/hubs=1; zero targets=0
  lower: min relative=111294036513341315/17982873312090390528 (0.00618889065067) at arms=120/leaves=1/hubs=1; zero targets=0
more than 60 leaves: trees=7, n=62..201, ordered moves=1498, max L=75, max c=14, violations=0
  upper: min relative=11622735140014/67411863812025839 (0.000172413793104) at isolated core, 61 leaves; zero targets=0
  lower: min relative=761725642649/835042648650815 (0.000912199686902) at isolated core, 61 leaves; zero targets=0
Total in-domain violations=0
```

## 4. Initial values, suprema, and roundings

At n=3, the sole core vertex is the centre and has state `(2,0,{})`. The exact f values and maxima are in the output. The upper maximum is attained at a core endpoint with 60 leaf neighbours and sole core neighbour of class `(0,2)`. This is possible on a three-vertex core path with a leaf at the other endpoint. The lower maximum is attained at the centre of a 61-vertex star.

For `s = lambda` or `delta`, Gamma recurrence gives
`prod_(m=3)^(n-1) (1+s/m) = [2/Gamma(3+s)] * Gamma(n+s)/Gamma(n)`.
Use the upper prefactor `f(x_0)` and lower prefactor `f(x_0)/sup f`, respectively.

`check_constants.py` encloses the constants using rational arithmetic. For a unit-rate Gamma random variable of shape x, Jensen with `t^s` gives `Gamma(x+s)/Gamma(x) <= x^s`. Applying Jensen at shape `x+s` with power `1-s` gives the complementary bound `x/(x+s)^(1-s) <= Gamma(x+s)/Gamma(x)`. The script uses both at `x=20000`, encloses the required rational powers by integer comparisons, and outward-rounds the recurrence product after each multiplication on a rational grid. This certifies that the lower coefficient exceeds 0.5022 and the upper coefficient is below 3.974.

Command and output from `constants.log`

```text
$ taskset -c 0 python3 -B checks/fable-xcharge-gpt/check_constants.py | tee checks/fable-xcharge-gpt/constants.log
Jensen at N=20000, rational grid=1/100000000000000000000000000000000000
upper: f(2,0,())=4064395905426029/562949953421312; sup f=4649094056001841/562949953421312; maximizers=[(60, 1, (1,))]
  certified coefficient interval: [3.973695817631798, 3.973744216839823]
  safe decimal bound 3.9740: PASS
lower: f(2,0,())=7031833084929469/4503599627370496; sup f=7970861646212325/4503599627370496; maximizers=[(60, 0, ())]
  certified coefficient interval: [0.502205275797162, 0.502211490540453]
  safe decimal bound 0.5022: PASS
lower newborn minimum over all classes=1
```

## 5. Conclusion

The exact certificates and Lemmas 1–3 at the supplied exponents prove

```text
0.5022 * Gamma(n + 0.55) / Gamma(n) <= E[C_n] <= 3.974 * Gamma(n + 0.58) / Gamma(n)
```

for every integer n≥3. The direct move calculations gave zero in-domain violations. The omissions in the written row notation and the qualifications for a general lower lemma have the explicit repairs above and do not change any of the checked inequalities.
