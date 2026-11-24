# Random friend trees: the forest below a degree cutoff

**Author:** Fable's GPT-6.2 Sol Ultrafast proof author. **Pinned source:** `a18825ccdecfdb8d8465a89a7d99fc33734fa8e9`. **Review and dependencies:** Astra (GPT) independently re-derived both coupling bounds; the Claude lead re-derived the age law and finite-age coupling. The age law uses Astra 6; the explicit fringe bounds do not. See [F4](../../README.md#f4).

Only §§1–2 are included; these are the scope of claim 4.

### 1. Promotion age — (L1)

Root the tree at its first core, reorienting the initial edge if necessary. Subsequent attachments are naturally oriented away from this root. Except for the first core, each core was first promoted by its unique parent. Core status is permanent. For $`a\ge0`$, the empirical fraction with promotion age at least $`a`$ at time $`t`$ is $`C(t-a)/C(t)`$ for $`t\ge a`$. By (L1) this converges a.s. to $`e^{-\mu a}`$. The same event of convergence of $`e^{-\mu t}C(t)`$ works simultaneously for every fixed $`a`$; equivalently one can use rational $`a`$ and monotonicity. The empirical distributions converge weakly a.s. to $`\mathrm{Exp}(\mu)`$.

If promotion is at **post**-promotion size $`m`$, the fraction with $`m/n\le x`$ is $`C_{\lfloor nx\rfloor}/C_n\to x^\mu`$ for $`0<x\le1`$, also by (L1). The limit is $`\mathrm{Beta}(\mu,1)`$, so $`-\log(m/n)`$ has limit $`\mathrm{Exp}(\mu)`$. This concerns promotion times of all cores, not age since birth or the age of a component's top vertex.

### 2. An inert parent gives a pathwise coupling — no (L)

Let $`\tau`$ be a promotion stopping time, with newly promoted $`v`$ below parent $`u`$. At $`\tau`$ its whole fringe consists of $`v`$ and the new leaf. The comparison is the process on this pair with a silent external neighbour of $`v`$: that neighbour is counted in the root's degree, and the root chooses it as usual, but this choice implements a null move inside. Compare for deterministic duration $`a`$, or for $`a\ge0`$ known at $`\tau`$ (e.g. $`a=t-\tau`$ for fixed $`t\ge\tau`$). Use the same independent rate-one marked clocks for corresponding fringe vertices, including paired newborns until the first divergence. External rings cannot affect the fringe except for a ring of $`u`$ choosing $`v`$. Predictably rank $`v`$ first at those rings. After the first divergence continue both processes with their correct transition rules. The inert process has its prescribed conditional law given $`\mathcal F_\tau`$ by fresh clock increments. Consequently the TV distance, even for the whole path to age $`a`$, is bounded by

```math
\Pr(\text{at least one }u\to v\mid\mathcal F_\tau)
\le \mathbb E\left[\left.\int_0^a\frac{ds}{D_u(\tau+s-)}\,\right|\mathcal F_\tau\right].
```

The second step is the compensator of the incoming count. If $`D_u(\tau)\ge N`$, all incoming rings lie in the $`1/N`$ thinning of $`u`$'s rate-one ring process, giving the sharper $`1-e^{-a/N}`$. This does not cover conditioning on $`u`$ being a hub eventually, the fringe having a prescribed shape, or age measured at a future size hitting time. It does cover applying any deterministic cutoff map to the coupled paths, with no increase in TV.

There is also an **all-ages** bound if the **post**-promotion leaf count is $`\ell=L_u(\tau)\ge4`$ and $`d=D_u(\tau)\ge4`$. Compare $`L_u`$ from below with a linear birth-death process $`B`$ starting at $`\ell`$, with upward rate $`B`$ and downward rate $`B/d`$. When counts agree the upward rate of $`L_u`$ is at least $`L_u`$, and its downward rate is $`L_u/D_u\le L_u/d`$. Run both with common same-direction jumps at the smaller rate and separate jumps at the remaining rates when they agree; when apart unit jumps cannot reverse the order. This adapted coupling leaves the tree's incoming-ring intensity unchanged: the auxiliary common jumps can be obtained by thinning the relevant tree events using fresh marks, with remaining auxiliary jumps from independent clocks.

Kill $`B`$ at $`\xi=\inf\{s:B(s)\le3\}`$, and let $`f(k)=1/(k-1)`$ before killing and $`f=0`$ at the cemetery. For $`k\ge4`$ its generator on $`f`$ is at most

```math
\frac{1}{k-1}\left(-1+\frac{k}{d(k-2)}\right)\le -\frac12 f(k),
```

where at $`k=4`$ the killed jump only decreases the generator. Dynkin's formula (bounded $`f`$ and generator; or localize) gives $`\mathbb E\int_0^\xi f(B(s))ds\le 2/(\ell-1)`$. The embedded biased walk gives $`\Pr(\xi<\infty)=d^{-(\ell-3)}`$, e.g. by gambler's ruin on $`[3,M]`$ then $`M\to\infty`$. Up to $`\xi`$ the incoming intensity is at most $`1/B`$, since $`D_u\ge L_u\ge B`$. Count incoming up to the adapted stop and add the killing probability. The whole path, even for all ages, therefore has TV bound

```math
\min\{1,\ 2/(\ell-1)+d^{-(\ell-3)}\}.
```

Nothing guarantees many leaves from $`d`$ alone. At a fixed finite age the elementary $`1-e^{-a/d}`$ is often sharper.
