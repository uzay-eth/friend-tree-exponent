# Independent root-ring Mellin / PGFL check

## Transfer for an arbitrary prescribed interface

Fix a core with initial leaf stock `L`, `d` nonleaf slots (including the planted parent), and a list of labelled existing child shapes. Keep **every** tick of its independent rate-one clock. In a gap of duration `t` there is no local promotion, send to a child, or parent selection. Leaf births form rate-one Yule families. Condition on imposed leaf impulses at times `v`, and on the upward send point processes `N_i` of all existing children. An individual entering at time `u` contributes a Yule family of age `t-u`. Thus the conditional generating factor at stock test `x^L` is

```
f_t(x)^L * product_v f_(t-v)(x) * product_i product_{u in N_i} f_(t-u)(x),
f_h(x) = x exp(-h) / [1-(1-exp(-h))x].
```

An independent externally prescribed Poisson input of rate `gamma(u)` replaces the imposed product by `exp(integral gamma(u)(f_(t-u)(x)-1)du)`; a deterministic input list uses the product literally. If exact input marks are correlated with the children, condition jointly rather than multiplying marginal factors.

Between monitored root ticks each child evolves with no input from this root except the specified impulses in its own interface. An upward selection is nonmutating for the child and records its send. For an arm with prescribed incoming schedule `I`, define the **operator**

```
B_I[w](t) g(Y) = E_Y^I [product_{u in N_up} w(u) * g(Y_t)].
```

Use child mask `w(u)=f_(t-u)(x)` in the root gap. Recursively, a parent-slot selection at an arm's tick supplies its externally specified mask at that tick, with no arm state change. On a silent planted root use mask 1. In a fixed gap masks can include additional nonnegative weights (e.g. a second observable); multiply masks at the same send. A prescribed downward send is an **add-one-leaf impulse at that time**, applied to the selected arm's endpoint configuration, not to all arms or to the monitored root. Interior marked times partition a gap by the semigroup composition, not by fresh independent copies of the gap duration. The operator carries child endpoint laws: a scalar PGFL with `g=1` cannot generally evaluate the next ring. Initial branches need not be virgin or identical. Given fixed shared times and a fixed action history, the arm integrations use the same symbolic variables, and factor on tensor endpoint tests (or are integrated against their tensor product kernel). The realized root choices reweight those variables and do not imply posterior arm independence after averaging.

Let `G(y)` denote the resulting pre-ring stock series, with child configurations carried as coefficients. For each particular old child or the parent, its choice is the operator

```
P_d G(x) = integral_0^1 z^(d-1) G(zx) dz.
```

For promotion it is

```
Q_d G(x) = integral_0^1 z^d G'(zx) dz.
```

On a term `y^L` these are `x^L/(L+d)` and `L x^(L-1)/(L+d)`, respectively. For promotion decrement the stock, increase `d` by one, and append a fresh virgin arm; **the total degree is unchanged by the promotion**. For a particular old-child choice apply its impulse; for the parent multiply by its outgoing mask. With masks one and endpoint test 1 the sum of choices is mass-preserving, as `d P_d G(x)+x Q_d G(x)=G(x)` before child-state updates. For higher derivatives carry divided derivatives `G^(j)/j!`. Their product rule is convolution of nonnegative jets. Chain rules at scaled arguments and Yule maps also have nonnegative coefficients on `0<=x<1` (and at `x=1` whenever moments are finite). Keep jets conditional on the shared time and integrate only after taking products; a derivative of a product of separately time-averaged arm responses is generally wrong. The same operators work for nested order windows, including any number of child send positions before the root ring.

An independent observation clock of rate `s` gives `V_s b = E[b(X_Exp(s))]/s`. Equivalently integrate each root gap against `exp(-(1+s)t) dt`; before the root tick the instantaneous toll can be evaluated by the same pre-ring stock pressure with no promotion/state change. A complete root step is therefore a gap kernel, followed by one of the ring operators; weighted exterior output contributes exactly at the parent branch. Independent killing adds its rate to `s` when it is constant.

## Truncation certificates

All bounds here apply to a **single intrinsic** discounted toll in `[0,1]`, or the absolute comparison under a coupling up to first disagreement for such tolls. For a `[0,1]` terminal observation rather than a resolvent, remove the factor `1/s`. Nonnegative masks at most one preserve the bounds; arbitrary unbounded masks or derivative tests need separate moment charges.

* **Root tick horizon.** After the `n`th root tick, the universal continuation charge is `(1+s)^(-n)/s`, from the actual rate-one clock, including misses. This is optimal on the class of bounded tolls (constant one attains it). At `.565728`, 27 ticks suffice for `1e-5`; at 28 the upper bound is `0.000006243`. A tighter state-specific alternative is to propagate a continuation interval at the stop instead of `[0,1/s]`.
* **Repeated contacts to a fixed old edge.** From initial source degree `D0`, the predictable contact rate to that slot is at most `1/D0`. Successive discounted first-arrival bounds give `(1+s D0)^(-k)/s` at the `k`th contact, whether or not earlier contacts changed the child or returned to the root. For a specified set of `m` old slots the factor becomes `m/(m+s D0)` per contact. A sum over different stop policies is valid by the union charge, with a separate bound for each tracked core. At initial `D0=4`, 12 repetitions of one edge cost at most `0.000001214`; a broad set of old slots may cost nearly the full root-clock tail.
* **Compress misses as operators.** A run of choices from a specified set of old slots and the parent is a positive subkernel `M` on full retained histories. In a stationary interface sum `I+M+M^2+...` before the next complementary event. In nonstationary time use time-ordered iteration on the enlarged state including absolute time and mask schedule. The root-clock discount makes it contractive on bounded tests, and the preceding slot bound sharpens it. A parent-slot miss is a scalar geometric factor only when its gap/stock/interface update truly repeats; choosing an old core is an impulse even though it does not change the monitored root's degree. The checker includes a stationary weighted example with exact remainder. A continuation-defect bound can be multiplied by the discounted exit mass, rather than paying the entire `1/s` on every boundary.
* **One old depth path.** Remove one deep incoming old edge, retaining its degree slot. With coupled common minimum-rate arrows, no root discrepancy can occur until a first omitted send crosses that edge and a first differential send crosses each successive old rootward edge. Growth on the tree cannot bypass the path. With starting source degrees `D_i` the bound is `product_i (1+s D_i)^(-1)/s`, by conditional Laplace bounds at each first crossing. Repeated returns after the first root disagreement are already charged by the residual toll range. For a single chain of minimum degree two, 16 edges suffice for `1e-5`. For an antichain, telescope and sum the actual path charges; there is no general claim that branch multiplicity is absorbed by this single-path factor. Signed old-edge continuation (`influence/intrinsic/README.md`) offers sharper path factors and a signed terminal when useful.
* **Depth of new promotions anywhere.** With `C` old cores, let depth count the new edges in a descendant chain after the initial skeleton. Every core promotes at rate at most one, so the expected number of depth-`h` births by `t` is at most `C t^h/h!` by induction. For a first depth-`h` exit split at `t=alpha h`, giving the bound `min(1,C(e alpha)^h + exp(-s alpha h))/s` on remaining toll. This is geometric in newly created depth, even with all branches retained, but slower in the pilot: the rational comparison with `alpha=.3` and one initial core gives `0.000007710` at depth 74. An output that observes only the root can additionally delay discrepancy to the first passage back along an old skeleton; the pilot does not use that refinement for this bound.

These local and depth bounds do not give geometric convergence in **total recursive core-ring count**: there can be many simultaneously active arms and missed rootward influence depends on their timing. In particular, a hard all-core-ring cap is at least as restrictive as the hard core count considered in the core-cut obstruction (`sol-hub-tail/core-cut/README.md`) on the common physical histories. An implementation intended to be practical should preserve passive arm output kernels, aggregate repeated contacts by subkernels, and charge a joint stopping policy by its discounted mass and the accuracy of its continuation. No full `V` or `Delta` at `1e-5` is inferred from the local tail constants.

## CPU 7 independent check

```sh
timeout 14s taskset -c 7 python3 -uB round2/sol-bipartite/root-clock-check/check.py
```

`oracle` races individually labelled leaf and core arrows. `advance` independently sums the coefficients of the Yule PGFs, integrating a shared gap exactly via a beta recurrence and using the exact monomial Mellin integrals at its next core or exogenous input tick. Refining a root gap at each nonroot tick partitions the same histories; it does not give separate durations to siblings. Both checks use the same explicit **physical-event** budget and payoff zero on premature budget exhaustion, to make the comparisons rational and finite. They compare a discounted bounded endpoint test after the last tracked root ring, without claiming that the physical-event-budget tail is small. A separate rational test compares the root-gap impulse PGFL with independent stock enumeration and a rigorous remainder, and a polynomial test checks derivative product jets and the missing chain-rule `z` fault.

| Fixture | Root ticks | Physical event budget | Common value, ceiling in `1e-9` units |
|---|---:|---:|---:|
| Virgin | 2 | 6 | 61923739 |
| One leaf, two distinct old arms (one contains a grandchild) | 2 | 6 | 7013015 |
| Bare fork with a nested old arm | 3 | 5 | 385949 |

The old-pair and bare-fork fixtures include external root immigration. All include ring-dependent exterior masks `(2/5,3/7,4/9)`, a nonconstant test of final child stocks and a planted parent slot. The replay rejects wrong physical denominator, promotion without stock decrement, idle not consuming a tick, omitted upward gain, omitted downward impulse, same-age treatment of scheduled immigrants, and missing Mellin chain-rule factor. The directly coordinating thread was sent the fixture conventions and output; its backward-jet implementation uses a different cutoff axis.
