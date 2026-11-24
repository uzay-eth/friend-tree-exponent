# Degree probabilities in an eigenmeasure of the core

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../SOURCES.tsv)).
> **Review and dependencies:** Claim 3: Claude and Gemini; GPT standalone audit; inputs include Round 1 growth, the marked renewal and spatial weights (Claims 2 and 4). Claim 4: Claude and Gemini; inputs include Round 1 one-edge coupling, planting and growth. See the [Astra ledger](../../../../LEDGER.md#astra) for the exact scopes and dependencies.


## Proposed statement

Let $`\pi`$ be the probability constructed in [the operator argument](operator.md), so $`\pi Q_t=e^{\mu t}\pi`$ and $`\pi(e^{\theta c_o})<\infty`$ for some $`\theta>0`$. Then for a constant $`0<A<\infty`$,

```math
\pi(D_o\ge k)\sim A k^{-\mu},\qquad \pi(D_o=k)\sim\mu A k^{-1-\mu}. \qquad\text{(1)}
```

The probability is an annealed, time-weighted subsequential law of the neighbourhood of a core vertex. Whether the ordinary degree frequencies at fixed tree size tend to that law is a separate question.

## Birth times

Let $`S_t f(G,o)=\mathbb E_{G,o} f(G_t,o)`$ track the same old root, and $`Bf(G,o)=b_o f(G^{b_o},v)`$, where $`v`$ is the newborn. For bounded nonnegative $`f`$, count separately the old root and the birth times of all other vertices in its cell. Conditional expectation at each birth and compensation of the promotion clocks yield

```math
Q_tf=S_t f+\int_0^t Q_s B S_{t-s}f\,ds.
```

The Yule first-moment bound for one cell justifies the integral. Applying $`\pi`$, multiplying by $`e^{-\mu t}`$, and sending $`t\to\infty`$ yields

```math
\pi f=\int_0^\infty e^{-\mu a}\pi B S_a f\,da. \qquad\text{(2)}
```

In particular $`\pi B1=\mu`$. Define the newborn probability $`\beta=\pi B/\mu`$. The root under $`\beta`$ has one leaf and one core neighbour. Its future age under the mixture (2) is exponential of rate $`\mu`$, independent of its subsequent dynamics.

## Tail constant

In the tree begun from $`\beta`$, write $`D(t),L(t),c(t)`$ for the marks of this particular newborn root, and $`T_k=\inf\{t:D(t)\ge k\}`$. Degrees rise by one at rate $`L+\sum_{w\sim o}1/D_w\le D`$. Thus $`e^{-t}D(t)`$ is a nonnegative supermartingale and converges almost surely to a random $`W`$. The degree process can be coupled below a pure-birth Yule process $`Y(t)`$ of rate $`Y(t)`$ begun at two. The nonnegative martingale $`e^{-t}Y(t)`$ is bounded in $`L^2`$ (its second moment is $`6-2e^{-t}`$), so Doob's inequality gives an integrable upper bound on $`\sup_{t\ge0}(e^{-t}D(t))^\mu`$ and uniform integrability of $`e^{-t}D(t)`$. All these bounds are uniform over the newborn environment.

The constant $`A=\mathbb E_\beta W^\mu`$ is positive. For one time unit silence the root's clock and let its tagged leaf ring at least twice. The root then has at least three leaves, whatever its original exterior. This event has a fixed positive probability. From any state with $`\ell\ge3`$ leaves, the leaf drift is $`L+\rho-b\ge L-1`$. Consequently $`\mathbb E[L(t)-1]\ge e^t(\ell-1)`$ conditionally on that state. Uniform integrability from the Yule bound, applied to the degree for each fixed state, gives positive conditional mean limiting degree after exponential normalization. Therefore $`\Pr_\beta(W>0)>0`$.

Equation (2) and monotonicity give $`\pi(D_o\ge k)=\mathbb E_\beta[e^{-\mu T_k};T_k<\infty]`$. On infinite degree growth $`T_k\uparrow\infty`$ by nonexplosion, and $`k e^{-T_k}\to W`$. On bounded degree growth set the expression to zero once the level is unreachable; there $`W=0`$. Dominated convergence gives the first asymptotic in (1).

## Point probabilities

For $`k\ge2`$ a newborn never has degree at least $`k+1`$. The only contribution to the time derivative of its count is a gain at a degree-$`k`$ core. Tracking a cell and compensating this crossing gives

```math
\mu\pi(D_o\ge k+1)=\pi\!\left[(k-c_o+\rho_o)\mathbf1_{D_o=k}\right],\qquad \rho_o=\sum_{w\sim o}\frac1{D_w}.
```

Indeed the crossing intensity per cell member is bounded by $`k`$, so integration over a unit of time against the eigenmeasure is legitimate even though the indicator is unrestricted away from the root. Since $`0\le\rho_o\le c_o/2`$, the right side is at most $`k\pi(D_o=k)`$. For each $`0<\delta<1`$ it is at least $`(1-\delta)k[\pi(D_o=k)-\pi(c_o>\delta k)]`$. The last probability is exponentially small by the moment bound. Apply the tail asymptotic at $`k+1`$ and send $`\delta\downarrow0`$ to obtain the second asymptotic in (1).
