# Cold challenge: first-defect near-one constant

**Holds** for each finite restart and tag, hence for every fixed label. Cross-read against 89ae's proof (`open/age-order/near-one-constant/README.md`), reviewed draft blob `739aa16a9eebe747e0e3306ed33a8e0ecca3af6a`; the writer corrected a display typo and added the ratio-CDF consequence in `3fb152f173f862b5122c972952e92832640ed80d`. Credit Fable Claim 5(e) for the original fixed-label two-sided orders and the physical deficit comparison. The additional result here is a coefficient and a conditional entrance law. The coefficient is a negative moment of a specified branching-system limit; replacing that limit by a Gamma law is not justified. Gamma factors enter through the exact prefix survival and a separate comparison urn.

## 1. Exact first-defect Palm law

Fix a size-`n0` restart, mark `o`, and put `a=Σ_{w~o}1/D_w`, `β=n0−a`. Let `M` be the size **before** the first attachment outside `o`. Along the all-tag prefix at size `m`, `r_o=m−β`, since each tag birth contributes a new reciprocal one and preserves all old neighbour degrees. Consequently

```
s_m = P(M≥m) = Γ(n0) Γ(m−β) / [Γ(a) Γ(m)],
w_m = P(M=m) = (β/m) s_m,
c = Γ(n0)/Γ(a),          w_m ~ c β m^(−β−1).
```

Use independent rate-one vertex clocks, each choosing uniformly among genuine neighbours. At the exit, a directed arrow `u→z`, `z≠o`, has conditional mass `1/(β D_u)`. All arrows internal to the initial nontrivial components of the tree with `o` removed keep their weights: root degrees count their slot towards `o`. Their total rate is `β−1`. In addition, the **whole** tag clock has rate one. Its core-target arrows have combined conditional probability `C_0/(β D_o(m))`; each creates a leaf in the chosen existing branch. Its other arrows each promote a singleton neighbour into a virgin two-vertex branch. The latter have limiting mass `1/β`. All stocked singleton leaves are equivalent after forgetting chronological labels. For a birth restart, `n0=v`, `a=1/d`, and the send `o→parent` has exact probability `1/(β D_o(m))`, rather than zero or a term removed from the first-defect hazard. For the initial edge there are no nontrivial components, `a=β=1`, and the virgin exit mass is exactly one.

The port convention matters: a branch root's send towards `o` attaches a new leaf **at `o`** in the physical tree. It is a null transition of that branch, not a branch birth, and does not promote an existing singleton. After the exit at physical size `m+1`, the future clocks can be restarted by the strong Markov property; conditioning on `M` and the exit arrow does not constrain subsequent rings.

## 2. Complete autonomous limit and tail control

Retain every nontrivial component of the tree minus `o`, with a single port at its root. In the limiting initial seed, either modify the indicated old component by the chosen intrinsic exit arrow, or add a virgin component to the unmodified old ones. Inside each branch, independent rate-one clocks choose their actual neighbours plus the one port at the root. Port choices are ignored within the component. Independently, a rate-one Poisson process starts independent virgin branches at its arrival times. This omits neither old components untouched by the arrow nor immigrants arriving arbitrarily late. Conditional on the limiting exit type, the branches and immigration clocks are independent. Mixing over types need not preserve independence.

Couple on `[0,T]` using identical intrinsic branch clocks, choices (including port), and tag rings. Their evolutions agree until the physical tag first selects a core. On agreement, the number of cores is at most `C_0+1+P_t` for post-exit tag ring count `P_t` and the degree is at least its order-`m` initial value. Hence the failure probability is `O((1+T)^2/m)` from the Poisson mean count and the sum of the successive core-target probabilities. The pre-exit tag-to-core type adds only `O(1/m)` under a coupling of initial types. Many real leaves can send to `o`, but these sends only increment the physical star reservoir, never alter a branch's intrinsic clocks, and can only reduce its core-selection hazard. Both processes are nonexplosive by the linear size bound.

Here is a uniform infinite-time argument, in unit-rate clock time started at the exit. Let `N(t)=size−D_o`, and let `k_i` be the sizes of the nontrivial components. Exactly `N=1+Σ_i(k_i−1)`. In either process every outside action increments `N` by one. Intrinsic branch `i` acts at rate `k_i−1/D_i`, with `D_i` counting the port; the physical tag always acts outside, and its ideal counterpart is the Poisson immigration clock. Thus the intensity is `N+A`, `A=Σ_i(1−1/D_i)`. At all times `0≤A≤C_0+1+P_t` even on the failed-coupling paths. The initial deficit after exit is fixed. The drift identity gives `E[N(t)+A(t)]≤K e^t` with `K` uniform over `m` and the exit arrow. For the integrating factor,

```
e^(−t) N(t) = N(0)+∫_0^t e^(−s)A(s) ds + martingale(t).
```

The drift tail in `L¹` is at most `K(1+T)e^(−T)`, and the martingale tail in `L²` is at most `sqrt(K e^(−T))` by the compensated jump isometry. This proves convergence and gives uniform `L¹` control of future increments. In the ideal system the same limit is `V=Σ_i W_i+Σ_j e^(−τ_j)W'_j`, where `W` is a branch's `lim e^(−t)(k(t)−1)` and `τ_j` are immigration times. Each fixed branch limit exists by the same discount argument; the sum agrees in `L¹` with the whole-system limit by the expected tail of arrivals `O(e^(−T))` and bounded expected discounted virgin size. There is at least one nontrivial component at exit. Its size minus one has growth rate at least its state, so Yule comparison makes its limit positive a.s.; the first branch transition time is an independent exponential factor in its limit, proving atomlessness even after adding independent branches and immigrants. Positivity also follows from the inverse-moment bound below.

For the physical system let `G_m=lim e^(−t) size(t)`. Total intensity is physical size, so `G_m~Gamma(m+1,1)` (shape, rate), the sum of `m+1` independent single-founder Yule limits. In particular `G_m/m→1` in probability; independence from the branch system is neither true by construction nor needed. The joint finite-horizon coupling and the uniform tail bound therefore give, with the exit type retained,

```
Y_m := m(1−Z_o) = [lim e^(−t)N(t)] / (G_m/m)  ⇒  V.
```

## 3. Inverse UI and the full endpoint mixture

After the first defect, `D_o` is unchanged, size increased by one, and its physical deficit `δ=D_o−r_o` did not decrease. Its fractional upper urn has red mass `r_o` and blue mass at least `β+1`. More simply retain the **predefect** shift `δ0`: the red mass can start at `m−β` and blue at `β+1` at size `m+1`, and the physical shifted degree stays below that urn under common uniforms. Thus `1−Z_o` is stochastically above `B_m~Beta(β+1,m−β)` uniformly over the exit arrow. For `β<p<β+1`,

```
E[Y_m^(−p) | exit arrow] ≤ m^(−p) Γ(β+1−p) Γ(m+1)
                                 / [Γ(β+1) Γ(m+1−p)] = O(1).
```

This is a bound on the true physical share, not just on an intrinsic branch or an unnormalized clock limit. It proves inverse UI at power `β`, `0<E[V^(−β)]<∞`, and `P(Y_m<x)≤K x^p`. In particular, at a birth restart the correct blue gap is `N_post+1−1/d=β+1` with `N_post=(m+1)−D_o`; confusing it with the pre-exit value loses the margin. Sizes with `εm<η` contribute at most `O(η^(p−β))` after dividing the endpoint mixture by `ε^β`; those with `εm>R` contribute at most `O(R^(−β))`. The middle is a Riemann sum using `Y_m⇒V` jointly with its type and atomlessness. Integrating `β x^(−β−1) 1{V<x}` over `x` gives exactly `V^(−β)`. Hence

```
P(Z_o>1−ε | restart) ~ [Γ(n0)/Γ(a)] E[V^(−β)] ε^β.
```

At fixed nonfounder label `v`, only the positive-probability birth sector `d=2` contributes to exponent `v−1/2`; sectors with `d>2` have strictly larger exponents. Sum their finite restart laws before conditioning: the leading coefficient is `Γ(v)/Γ(1/2)` times the birth-law average of `1{d=2}E[V^(−(v−1/2))|restart]`. Founders are exchangeable, with exponent one and restart factor `Γ(2)/Γ(1)=1`.

## 4. Conditional entrance law

For a fixed restart, retain exit type `J`, scaled pre-exit size `X=εM` and final amplitude `Y=M(1−Z_o)`. Conditional on the strict endpoint event their limit has joint measure

```
 β x^(−β−1) 1{0<y<x} dx P((V,J)∈d(y,j)) / E[V^(−β)].
```

Equivalently `(Y,J)` has the `V^(−β)` tilt of its limiting Palm law, and `X/Y` is independent Pareto of exponent `β` on `[1,∞)`. In particular `(1−Z_o)/ε` has limiting CDF `u^β`, `0≤u≤1`, independent of the tilted pair; entrance type does **not** retain its original Palm weights in general. The same weighting applies to restart states in the `d=2` sector for unconditional fixed-label conditioning. Replacing `M` by the post-exit size `M+1` has no effect on these limits. No claim is made about convergence of a microscopic chronological label on the promoted stock leaf. These fixed-label limits do not supply a uniform-in-label near-one window or resolve intermediate-threshold age order.

## Replay

`near-one-constant-sol.py` enumerates real directed arrows, including tag-to-core choices, and checks branch accounting, the postdefect inverse gap, and prefix telescoping on every tagged history through size six with seven prefix lengths. A solvable ordinary Pólya urn checks the exact first-defect mixture with censored prefix residual and the Gamma normalization independently of the physical branch construction. The code is a finite control; the infinite-horizon and inverse UI steps above establish the limiting statement.

```text
nice -n 10 timeout 12s python3 -IB round2/checks/near-one-constant-sol.py
PASS 6104 prefix cases: directed Palm weights, all exit types, branch deficit and urn gap
PASS 40 exact urn-mixture truncations and four Gamma survival-ratio controls
```
