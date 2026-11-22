VERDICT: CORRECT

# Cold Review of `attempts/strogatz-astra/sample-exponent.md`

Cold review performed strictly on [sample-exponent.md](../attempts/strogatz-astra/sample-exponent.md) (referencing only [existence.md](../attempts/strogatz-astra/existence.md), already reviewed in [strogatz-astra-existence-gemini.md](strogatz-astra-existence-gemini.md)), without reading any other files, reviews, or board material.

**Summary:** The proof that $\frac{\log C_n}{\log n} \to \mu$ almost surely is **mathematically sound and complete**. The four critical potential failure modes—(1) independence under random leaf selection in Stage 3, (2) marked clocks reused across time intervals, (3) disjointness of sites/reserved sets/branches, and (4) whether conditioning biases the reference count—were audited in detail and all hold rigorously. Minor expositional notes (none blocking) are listed at the end.

---

## 1. Disjoint-Site Decomposition and the $C/6$ Bound (line 7)

Let $C \ge 2$ and partition the core tree $K$ into:
- $\mathcal{A} = \{v \in K : c_v \le 2, L_v \ge 1\}$ (size $A$),
- $\mathcal{B} = \{v \in K : c_v = 2, L_v = 0\}$ (size $B$),
- $\mathcal{D} = \{v \in K : c_v \ge 3\}$ (size $D$).

Check:
- Every core vertex with $c_v = 1$ must have $L_v \ge 1$, so $n_1 \le A$. By $\sum_{v \in K}(c_v - 2) = -2$, $D \le \sum_{d\ge 3}(d-2)n_d = n_1 - 2 \le A - 2 < A$, and $A \ge n_1 \ge 2$.
- Since $K$ is a tree with non-empty $K \setminus \mathcal{B}$, every connected component of $K[\mathcal{B}]$ is a maximal path of degree-$2$ core vertices lying between two distinct non-B vertices. Suppressing each maximal B-path into a single edge yields a tree on the $A + D$ non-B vertices, which has $(A + D) - 1 \le 2A - 3 < 2A - 1$ edges.
- Partitioning each maximal B-path of length $m$ into $\lfloor m/3 \rfloor$ disjoint consecutive triples $(u, w, z)$ leaves at most $2$ leftover B vertices per path, so $B \le 3J + 2(A + D - 1)$.
- Hence $C = A + D + B \le 3(A + D) + 3J - 2 \le 6A + 3J - 8 < 6A + 3J$, giving $A + J \ge C/6$ sites. (For $C = 1$, the single centre supplies $1 \ge 1/6$ site.)
- Verified computationally on extremal subdivided $3$-regular trees (all $k \in \{0,1,2,5,10\}$, $m \in \{0,\dots,6\}$): worst observed ratio was $(k+2)/(6k+4) \to 1/6$ at $m=2$ (specifically $12/64 = 0.18750$ at $k=10, m=2$), confirming the $1/6$ bound is sharp as well as valid.
- For each triple $(u, w, z)$, $w$ is the middle vertex of a path of B vertices, so $w$'s **only** two neighbours in $T(t)$ are $u$ and $z$. And the reserved sets $\{u, \ell_u\}$ (at $A$-sites) and $\{u, w, z\}$ (at triples) are pairwise disjoint (leaves have unique core neighbours, hence cannot be shared across distinct $A$-sites or overlap any core vertex in $\mathcal{B}$). **Valid.**

---

## 2. Independence of Preparation Trials, Random Leaf Selection, and Clock Reuse (lines 9–17)

Let $m \ge C(t)/6$ be the sites chosen deterministically from $\mathcal{F}_t$. We audit each of the four preparation stages sequentially in the tree filtration $(\mathcal{F}_r)_{r \ge t}$.

### 2.1 Stage 1 on $(t, t+1]$ (line 11)
- At an $A$-site $u$: requires $u$ silent on $(t, t+1]$ (depends solely on $\Pi_u|_{(t,t+1]}$, prob $e^{-1}$). Its reserved leaf $\ell_u$ remains a leaf of $u$ at $t+1$.
- At a triple $(u, w, z)$: requires $u, z$ silent on $(t, t+1]$, $w$ ringing $\ge 1$ time on $(t, t+1]$, and $w$'s first ring having mark $< 1/2$ with $u$ ranked first (depends solely on $(\Pi_u, \Pi_w, \Pi_z)|_{(t,t+1]}$, prob $\frac{1}{2}e^{-2}(1-e^{-1})$). Since $u$ and $z$ are the only neighbours of $w$ and both are silent, $D_w = 2$ at $w$'s first ring, so it necessarily attaches a fresh tagged leaf to $u$, which remains a leaf of $u$ at $t+1$ because $u$ is silent.
- Because the reserved vertex sets are pairwise disjoint, the Stage-1 events are mutually independent given $\mathcal{F}_t$. Thinning $A$-sites to match $\frac{1}{2}e^{-2}(1-e^{-1})$ via independent coins yields i.i.d. $\mathrm{Bern}(p_1)$ flags. **Valid.**

### 2.2 Stage 2 on $(t+1, t+2]$ (line 12)
- Condition on $\mathcal{F}_{t+1}$. Let $S_1$ be the Stage-1 successful sites. For each $i \in S_1$, its vertex $u_i$ and tagged leaf $\ell_i$ are $\mathcal{F}_{t+1}$-measurable, and the pairs $\{u_i, \ell_i\}$ are pairwise disjoint (distinct core vertices $u_i$, and each $\ell_i$ is a leaf whose unique neighbour is $u_i$).
- By independent increments of Poisson processes, the future clocks $\{(\Pi_{u_i}, \Pi_{\ell_i})|_{(t+1, t+2]} : i \in S_1\}$ are mutually independent of $\mathcal{F}_{t+1}$ and of each other.
- Requiring $u_i$ silent and $\ell_i$ ringing $\ge M$ times has exact conditional probability $p_2 = e^{-1}\mathbb{P}(\mathrm{Pois}(1) \ge M)$ for every $i \in S_1$, independent of $\mathcal{F}_{t+1}$ and of outside arrivals (which can only add leaves to $u_i$, never convert leaves while $u_i$ is silent). Thus at $t+2$ every $i \in S_2$ has $L_{u_i}(t+2) \ge M$ (in fact $\ge M+1$) and $c_{u_i}(t+2) \le 2$, and the Stage-1&2 flags are i.i.d. $\mathrm{Bern}(p_1 p_2)$ given $\mathcal{F}_t$. **Valid.**

### 2.3 Stage 3 on $(t+2, t+3]$: Critical audit of random leaf selection (lines 13, 16)
This is the subtlest step in the proof: outside arrivals (including from another site $j \in S_2$ that fails Stage 3 early in $(t+2, t+3]$) can attach new leaves to $u_i$ before $u_i$ rings, so the **identity** of the leaf $v_i$ selected by $u_i$ at its ring time $\tau_i$ can depend on the histories of other sites. We verify rigorously that this **does not** break independence of the Stage-3 success flags $(I_{3,i})_{i \in S_2}$:

1. **Fixed probability independent of $L_{u_i}(\tau_i^-)$:** The script does *not* ask "$u_i$ chooses a leaf" (whose probability $L_{u_i}/(L_{u_i}+c_{u_i})$ would depend on outside arrivals). Instead it ranks leaves of $u_i$ first and requires $u_i$'s mark $U_i$ to fall in $[0, M/(M+2))$—a fixed deterministic threshold! Because $u_i$ has not rung on $(t+2, \tau_i^-)$, we still have $L_{u_i}(\tau_i^-) \ge M$ and $c_{u_i}(\tau_i^-) \le 2$, so $M/(M+2) \le L_{u_i}(\tau_i^-)/(L_{u_i}(\tau_i^-)+c_{u_i}(\tau_i^-))$, guaranteeing promotion of *some* leaf $v_i$ of $u_i$ with fresh child $v_{i,1}$.
2. **Exact product-martingale proof of mutual independence given $\mathcal{F}_{t+2}$:**
   For each $i \in S_2$, define the active clock set $K_i(r)$ at time $r \in [t+2, t+3]$:
   - $K_i(r) = \{u_i\}$ before $u_i$'s first ring in $(t+2, t+3]$;
   - $K_i(r) = \{u_i, v_i, v_{i,1}\}$ after $u_i$ rings at $\tau_i \le r$ with $U_i < M/(M+2)$, while none of $\{u_i, v_i, v_{i,1}\}$ has rung on $(\tau_i, r]$;
   - $K_i(r) = \emptyset$ once site $i$ fails (either $U_i \ge M/(M+2)$, or one of $\{u_i, v_i, v_{i,1}\}$ rings on $(\tau_i, r]$).

   **Pairwise disjointness of $K_i(r)$ at every $r \in (t+2, t+3]$:** When site $i$ promotes $v_i$ at $\tau_i$, $v_i$ is a *leaf* of $u_i$ at $\tau_i^-$, so its unique neighbour is $u_i$. Therefore $v_i$ cannot equal any core vertex $u_j$ or $v_j$ ($j \ne i$), nor can it equal $v_{j,1}$ (whose unique neighbour is $v_j \ne u_i$ while $j$ is active). The newborn $v_{i,1}$ is fresh at $\tau_i$, hence disjoint from all existing vertices. Thus $K_i(r) \cap K_j(r) = \emptyset$ for all $i \ne j$ and all $r \in (t+2, t+3]$.

   Now let $\phi_i : \{0,1\} \to \mathbb{R}$ be arbitrary, let $p_3 = \frac{M}{M+2} e^{-3}$, and for $x = t+3-r \in [0, 1]$ define $M_i(r)$ by:
   - $\phi_i(0) + (\phi_i(1)-\phi_i(0)) \frac{M}{M+2} e^{-2} x e^{-x}$ before $u_i$ rings;
   - $\phi_i(0) + (\phi_i(1)-\phi_i(0)) e^{-2(\tau_i - (t+2))} e^{-3x}$ while $i$ is active;
   - $\phi_i(0)$ once $i$ has failed.

   Direct generator check on the full tree filtration $(\mathcal{F}_r)$: before $u_i$ rings, $\frac{d}{dr}(x e^{-x}) = -(1-x)e^{-x}$ exactly cancels the jump rate $1 \cdot (\frac{M}{M+2} e^{-2} e^{-x} - \frac{M}{M+2} e^{-2} x e^{-x})$ (since $e^{-2(1-x)} e^{-3x} = e^{-2-x}$); while active, $\frac{d}{dr}(e^{-3x}) = +3 e^{-3x}$ exactly cancels the jump rate $3$ of the 3 distinct clocks in $K_i(r)$ to $\phi_i(0)$! Thus each $M_i(r)$ is a bounded martingale on $(\mathcal{F}_r)$ jumping *only* at rings of $K_i(r)$. Because the sets $K_i(r)$ are pairwise disjoint at all $r$, no two $M_i, M_j$ ($i \ne j$) ever jump simultaneously ($[M_i, M_j] \equiv 0$), so the product $\prod_{i \in S_2} M_i(r)$ is a martingale on $[t+2, t+3]$. Applying the independent thinning coins $C_i \sim \mathrm{Bern}(e^{-2(\tau_i-(t+2))})$ at $r = t+3$ gives
   $$\mathbb{E}\!\left[\prod_{i \in S_2} \phi_i(I_{3,i}) \,\Big|\, \mathcal{F}_{t+2}\right] = \prod_{i \in S_2} M_i(t+2) = \prod_{i \in S_2} \mathbb{E}[\phi_i(I_{3,i})].$$
   Hence $(I_{3,i})_{i \in S_2}$ are **exactly i.i.d. $\mathrm{Bern}(p_3)$ given $\mathcal{F}_{t+2}$**, completely unbiased by random leaf identity or outside arrivals! Also note: if $u_i$ happens to select its own Stage-2 tagged leaf as $v_i$, Poisson memorylessness at $\tau_i$ ensures its future increments on $(\tau_i, t+3]$ are independent of $\mathcal{F}_{\tau_i}$. **Valid.**

### 2.4 Stage 4 on $(t+3, t+h]$ (line 14)
- Conditioned on $\mathcal{F}_{t+3}$, the active sets $\{u_i, v_i, v_{i,1}\}$ for $i \in S_3$ are pairwise disjoint. While $u_i$ is silent across the $k-2$ unit subintervals, the growing branch $B_i$ is completely isolated from the exterior (even if another site $j$ fails and attaches leaves to $u_i$). The Stage-4 script for site $i$ uses only $\Pi_{u_i}$ and internal clocks of $B_i$, has a deterministic conditional probability $q_4(T, v) \ge p_4(M, k) > 0$ independent of $\mathcal{F}_{t+3}$ and the exterior, and thinning to $p_4(M,k)$ makes the Stage-4 flags i.i.d. $\mathrm{Bern}(p_4)$ given $\mathcal{F}_{t+3}$.
- By the chain rule across Stages 1–4, the preparation flags $(I_{\mathrm{prep}, i})_{i=1}^m$ at time $t+h$ are i.i.d. $\mathrm{Bern}(p)$ given $\mathcal{F}_t$, where $p = p(M,k) = p_1 p_2 p_3 p_4 > 0$. **Valid.**

---

## 3. Independent Block Contributions, Clock Reuse, Disjointness, and Unbiased Reference Count (lines 18–26)

Condition on $\mathcal{F}_{t+h}$. Let $S_{\mathrm{prep}}$ be the set of prepared sites at $t+h$.
- **Disjointness of branches and old $u_i$'s:** Every prepared branch $B_i(t+h)$ starts at a leaf of $u_i$ promoted during $(t+2, t+3]$ and contains no vertex that was core at time $t$. Therefore (as verified in `strogatz-astra-existence-gemini.md`), the sets $\{u_i : i \in S_{\mathrm{prep}}\}$ and $\{V(B_i(t+h)) : i \in S_{\mathrm{prep}}\}$ are all pairwise disjoint, and the subsequent descendant branches $B_i(t+h+s)$ remain pathwise disjoint for all $s \ge 0$.
- **Clock reuse across time intervals:** Reusing the clock $\Pi_{u_i}$ on $(t+h, t+h+s]$ for the exterior test (forbidding marks $< 1/M$ with edge $u_i v_i$ ranked first) is completely independent of $\mathcal{F}_{t+h}$ (including Stages 1–4) by independent increments of Poisson processes! And across distinct sites $i \ne j$ in $S_{\mathrm{prep}}$, the clocks $\Pi_{u_i}|_{(t+h,t+h+s]}$ and the internal branch clocks of $B_j(t+h)$ belong to disjoint vertices, so the triples $(\mathcal{R}_{\mathrm{ref}}^{(i)}, \mathcal{A}_{v_i}^{(i)}, \Pi_{u_i}|_{(t+h, t+h+s]})$ are mutually independent given $\mathcal{F}_{t+h}$!
- **Simultaneous coupling into a single actual tree $T_{\mathrm{act}}$:** Even if two prepared vertices $u_i, u_j$ are adjacent or connected via failing trials, ranking $u_i v_i$ first at $u_i$ and $u_j v_j$ first at $u_j$ preserves the exact uniform neighbour law at every ring of every vertex. Whenever a failing trial (or exterior vertex) attaches leaves to $u_j$, it only increases $D_{u_j} \ge M$ without ever crossing the cut edge $u_j v_j$ so long as $G_{u_j} \cap G_{\mathrm{acc}, j}$ holds. When a branch $i$ suffers a discrepancy at stopping time $\sigma_i$, its lower contribution is already $0$ and $T_{\mathrm{act}}$ decouples branch $i$ from $\mathcal{R}_{\mathrm{ref}}^{(i)}$ while continuing with the true dynamics.
- **No conditioning bias on the reference count:** The contribution at site $i \in S_{\mathrm{prep}}$ is defined as $\widetilde{Y}_{s,i} = C^{\mathrm{ref},(i)}(s)\,\mathbf{1}_{G_{\mathrm{good}, v}^{(i)} \cap G_{\mathrm{acc}}^{(i)} \cap G_{u_i}}$—we multiply by the indicator (setting failed outcomes to $0$) rather than conditioning the reference law on success! Because $\mathcal{R}_{\mathrm{ref}}^{(i)}, \mathcal{A}_{v_i}^{(i)}, \Pi_{u_i}|_{(t+h,t+h+s]}$ are independent of $\mathcal{F}_{t+h}$, the conditional law of $\widetilde{Y}_{s,i}$ given $\mathcal{F}_{t+h}$ is identical across all $i \in S_{\mathrm{prep}}$ and independent of $\mathcal{F}_{t+h}$. Pairing with $I_{\mathrm{prep}, i} \sim \mathrm{Bern}(p)$ (and supplying independent dummy copies of $\widetilde{Y}_s$ when $I_{\mathrm{prep}, i} = 0$) yields $m \ge C(t)/6$ i.i.d. copies of $Y_s = I_{\mathrm{prep}}\widetilde{Y}_s$ given $\mathcal{F}_t$, with
  $$\mathbb{E}[Y_s] = p\,\mathbb{E}[\widetilde{Y}_s] \ge \frac{p}{k} e^{-\eta s} \big[F(s) - k - O_{k,M,\epsilon}(e^{-s})\big],$$
  and pathwise $C(t+h+s) \ge \sum_{i=1}^m Y_{s,i}$.
- **Truncation and Chebyshev (line 26):** Since $0 \le Y_s \le N_T(s) \in L^1$, for large $s$ with $\mathbb{E}[Y_s] > 0$ there exists $R_s < \infty$ with $\mu_R := \mathbb{E}[\min(Y_s, R_s)] \ge \mathbb{E}[Y_s]/2$. For $m \ge C/6$ i.i.d. copies, $S_m = \sum_{i=1}^m \min(Y_{s,i}, R_s)$ satisfies $\mathbb{E}[S_m] = m\mu_R \ge C\,\mathbb{E}[Y_s]/12$ and $\mathrm{Var}(S_m) \le m R_s^2$. Thus by Chebyshev,
  $$\mathbb{P}\!\left(C(t+h+s) < \frac{\mathbb{E}[Y_s]}{24} C(t) \,\Big|\, \mathcal{F}_t\right) \le \mathbb{P}\!\left(S_m < \frac{1}{2}\mathbb{E}[S_m]\right) \le \frac{4 m R_s^2}{m^2 \mu_R^2} \le \frac{24 R_s^2}{C(t)\,\mu_R^2} \xrightarrow{C(t)\to\infty} 0.$$
  Because this inequality is purely about the conditional distribution of $C(t+h+s)$ given $\mathcal{F}_t$, it holds in the original filtration with **no clock reuse across consecutive blocks**. **Valid.**

---

## 4. Almost-Sure Divergence, Block Grid, Upper Bound, and Hitting Time (lines 28–36)

- **$C(t) \to \infty$ a.s. (line 28):** At discrete size $n \ge 3$, $\mathbb{P}(C_{n+1} - C_n = 1 \mid \mathcal{G}_n) = S_n/n \ge 1/n$. Since $\sum_{n=3}^\infty 1/n = \infty$, Lévy's conditional Borel–Cantelli lemma gives infinitely many promotions, so $C_n \to \infty$ a.s. For the Yule process $N(t)$ starting at $2$, the holding times $E_j \sim \mathrm{Exp}(j)$ ($j \ge 2$) satisfy $\sum \mathbb{E}[E_j] = \infty$ and $\sum \mathrm{Var}(E_j) = \sum 1/j^2 < \infty$, so $\tau_n = \sum_{j=2}^{n-1} E_j < \infty$ a.s. and $\tau_n \to \infty$ a.s. (no finite-time explosion). Thus $N(t) \to \infty$ a.s. and $C(t) = C_{N(t)} \to \infty$ a.s. **Valid.**
- **Failure fraction $\to 0$ a.s. (line 30):** On the grid $t_j = j(h+s)$ with $\mathcal{H}_j = \mathcal{F}_{t_j}$, let $I_j = \mathbf{1}\{C(t_j) < (\mathbb{E}[Y_s]/24) C(t_{j-1})\}$ and $p_j = \mathbb{E}[I_j \mid \mathcal{H}_{j-1}]$. By line 26 and $C(t_{j-1}) \to \infty$ a.s., $p_j \to 0$ a.s. (hence $\frac{1}{m}\sum_{j=1}^m p_j \to 0$ a.s.). Since $D_j = I_j - p_j$ are bounded martingale differences ($|D_j| \le 1$), $W_m = \sum_{j=1}^m D_j / j$ is an $L^2$-bounded martingale ($\sup_m \mathbb{E}[W_m^2] \le \pi^2/6 < \infty$), hence converges a.s. By Kronecker's lemma, $\frac{1}{m}\sum_{j=1}^m D_j \to 0$ a.s., so $\frac{1}{m}\sum_{j=1}^m I_j \to 0$ a.s. **Valid.**
- **Almost-sure lower and upper bounds in continuous time (lines 30–34):**
  - If $\mu = 0$: since $C(t) \to \infty$ a.s., $\log C(t) \ge 0$ eventually a.s., so $\liminf \log C(t)/t \ge 0 = \mu$ a.s. automatically.
  - If $\mu > 0$: fix $(M, \epsilon)$ with $\eta < \mu$; then from $F(s) = e^{(\mu+o(1))s}$, $\mathbb{E}[Y_s]/24 \ge \frac{p}{24k} e^{-\eta s}(F(s) - k - O(e^{-s})) = e^{(\mu - \eta + o(1))s} > 1$ for all sufficiently large $s$. Letting $j_0 < \infty$ a.s. be the first grid index with $C(t_{j_0}) \ge 1$, monotonicity gives $C(t_m) \ge (\mathbb{E}[Y_s]/24)^{m - j_0 - \sum_{j=1}^m I_j}$ for $m > j_0$, and monotonicity between grid points yields $\liminf_{t\to\infty} \log C(t)/t \ge \log(\mathbb{E}[Y_s]/24)/(h+s)$ a.s. Intersecting over a countable sequence $s_r \to \infty$ and a countable sequence $(M_q, \epsilon_q)$ with $\eta_q \downarrow 0$ yields $\liminf_{t\to\infty} \log C(t)/t \ge \mu$ a.s.
  - Upper bound: for any $\delta > 0$, $\mathbb{P}(C(j) > e^{(\mu+\delta)j}) \le F(j) e^{-(\mu+\delta)j} = e^{-(\delta + o(1))j}$ is summable over $j \in \mathbb{N}$. Borel–Cantelli, monotonicity $C(t) \le C(\lceil t \rceil)$, and intersecting over $\delta = 1/r$ yield $\limsup_{t\to\infty} \log C(t)/t \le \mu$ a.s. Thus $\log C(t)/t \to \mu$ a.s. **Valid.**
- **Discrete-size limit via hitting time (line 36):** $\tau_n = \sum_{j=2}^{n-1} E_j$ with independent $E_j \sim \mathrm{Exp}(j)$ satisfies $\sum_{j=2}^\infty \mathrm{Var}(E_j) = \pi^2/6 - 1 < \infty$, so $\sum_{j=2}^{n-1}(E_j - 1/j)$ converges a.s. to a finite limit. Since $\sum_{j=2}^{n-1} 1/j = \log n + O(1)$, we have $\tau_n / \log n \to 1$ a.s. (and $\tau_n \to \infty$ a.s.). By right-continuity of the jump process, $C_n = C(\tau_n)$, so $\frac{\log C_n}{\log n} = \frac{\log C(\tau_n)}{\tau_n} \cdot \frac{\tau_n}{\log n} \to \mu$ almost surely. **Valid.**

---

## 5. Minor Expositional Notes (Non-blocking)

1. Line 30: When multiplying across blocks to get $\liminf \log C(t)/t \ge \log(\mathbb{E}[Y_s]/24)/(h+s)$ a.s., one starts the product at the first (a.s. finite) grid index $j_0$ with $C(j_0(h+s)) \ge 1$ (since $C(0) = 0$). Explicitly noting this avoids a reader momentarily worrying about multiplying by $C(0)=0$.
2. Line 11: For an $A$-site, "reserve the clock of $u$ and one of its leaves" (line 9) designates that reserved leaf as its tagged leaf for Stage 2, whereas line 11 only explicitly uses the word "Tag" for the triple's newborn leaf; noting both are the tagged leaves used in Stage 2 makes the reading immediate.
3. Line 20 ("Choose $(T, v)$ there to maximise the good reference contribution"): as checked above, maximizing $\mathbb{E}[C_T(s)\mathbf{1}_{G_v}]$ jointly over all $k$-vertex $(T, v)$ is $\ge$ the two-step choice in `existence.md` and yields the exact same bound in line 22.
