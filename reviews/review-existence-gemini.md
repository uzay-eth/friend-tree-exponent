VERDICT: CORRECT

# Cold Review of `attempts/strogatz-astra/existence.md`

Cold review performed strictly on the statement and proof in `existence.md`, without reading any other files in `attempts/strogatz-astra/`, other reviews, or board material.

**Summary:** The proof that $\lim_{n \to \infty} \frac{\log \mathbb{E}[C_n]}{\log n}$ exists is **mathematically sound and complete**. All probability estimates, planting uniformity and disjointness claims, clock couplings, quantifier orders, and de-Poissonization arguments check out rigorously. A few minor expositional notes (none affecting validity) are collected at the end.

---

## 1. Continuous-Time Embedding and Yule Second Moment (lines 5–6)

- **Jump-chain equivalence:** Assigning each vertex an independent rate-$1$ Poisson clock that selects a uniformly random neighbour and attaches a fresh leaf to it gives total birth rate $j$ at tree size $j$, with the ringing vertex uniform over the $j$ vertices and independent of the $\mathrm{Exp}(j)$ holding time. Thus the embedded skeleton is identically the discrete random friend tree $(T_n)_{n \ge 2}$, independent of $(N(t))_{t \ge 0}$. **Valid.**
- **Second moment of $N_T(t)$:** Starting from size $k$, the pure-birth generator $\mathcal{G}_{\mathrm{Yule}} g(j) = j(g(j+1)-g(j))$ yields $\frac{d}{dt}\mathbb{E}[N_T(t)] = \mathbb{E}[N_T(t)]$ ($\implies k e^t$) and $\frac{d}{dt}\mathbb{E}[N_T(t)^2] = 2\mathbb{E}[N_T(t)^2] + k e^t$. Differentiating $m_2(t) = k(k+1)e^{2t} - k e^t$ gives $m_2(0) = k^2$ and $m_2'(t) = 2m_2(t) + k e^t$. **Exact.**

---

## 2. Eligibility Bound, Planting Uniformity, and Disjointness (lines 8–24)

### 2.1 At least $C/4$ eligible vertices (lines 10–12)
- When $C \ge 1$, every path between two non-leaves consists solely of non-leaves, so the core vertices induce a tree.
- For $C = 1$, the unique core vertex has $c = 0 \le 2$ and $L \ge 2 \ge 1$, so it is eligible ($1 \ge 1/4$).
- For $C \ge 2$, let $E$ be the number of eligible core vertices, $I$ the number of type-B vertices ($c=2, L=0$) with no B neighbour, and $n_d$ the number of core vertices of core degree $d$. Every core vertex has $c \ge 1$; vertices with $c = 1$ must have $L \ge 1$ (since $c + L \ge 2$), hence are eligible ($n_1 \le E$). Thus the ineligible core vertices are partitioned into $I$ and $\{c \ge 3\}$.
  - Because $I$ is an independent set in the core tree and each vertex in $I$ has core degree $2$, their $2I$ incident core edges are pairwise distinct among the $C-1$ core edges, so $I \le (C-1)/2$.
  - On the core tree, $\sum_{v \in K} (c_v - 2) = -2$ implies $n_{\ge 3} \le \sum_{d \ge 3}(d-2)n_d = n_1 - 2 < n_1 \le E$.
  - Combining $C - E = I + n_{\ge 3} \le \frac{C-1}{2} + (E - 2)$ yields $E \ge \frac{C+5}{4} > \frac{C}{4}$ (and even the weaker $n_{\ge 3} \le n_1$ gives $E \ge \frac{C+1}{4} > \frac{C}{4}$).
- Note that the constant $1/4$ is asymptotically tight (subdividing every edge of a 3-regular core tree once with a type-B vertex gives $E/C = (k+2)/(4k+3) \to 1/4$). **Valid.**

### 2.2 Explicit 4-stage planting construction and uniformity of $(h, p)$ (lines 14–23)
Fix $M \ge 2$, $k \ge 2$, and any rooted unlabelled tree $(T, v)$ on $k$ vertices. Take an eligible vertex $u$ at time $t$ (so $c_u(t) \le 2$).
- **Stage 1 (unit $[t, t+1]$):**
  - If $L_u(t) \ge 1$, keeping $u$ silent on $(t, t+1]$ has probability $e^{-1}$. Since only rings of $u$ can convert leaves of $u$, $c_u$ cannot increase and no leaf of $u$ can be lost.
  - Otherwise $u$ is type B with a type-B neighbour $w$ (degree $2$, neighbours $u$ and $z \ne u$, since the graph is a simple tree). Keeping $u$ and $z$ silent on $(t, t+1]$ prevents *any* arrival at $w$ throughout $(t, t+1]$ (only neighbours of $w$ can redirect to $w$). Thus before $w$'s first ring, $D_w = 2$ deterministically, and requiring $w$ to ring at least once with its first ring choosing $u$ has probability $e^{-2}(1-e^{-1})/2 > 0$. Any further rings of $w$ or arrivals at $u$ from outside only add leaves to $u$ or $z$; since $u$ is silent, $u$ retains the created leaf and keeps $c_u(t+1) = 2$.
  - In both cases, at $t+1$ we have $L_u(t+1) \ge 1$ and $c_u(t+1) \le 2$, with conditional probability $\ge \frac{1}{2}e^{-2}(1-e^{-1})$. **Valid.**
- **Stage 2 (unit $[t+1, t+2]$):**
  - Choose one leaf $\ell$ of $u$ deterministically from $\mathcal{F}_{t+1}$ (fixed vertex ordering). Keeping $u$ silent and requiring $\ell$ to ring $\ge M$ times on $(t+1, t+2]$ has conditional probability $e^{-1}\mathbb{P}(\mathrm{Pois}(1) \ge M) > 0$. While $u$ is silent, $\ell$ stays a leaf (its only neighbour is $u$), so every ring of $\ell$ adds a leaf to $u$. Hence $L_u(t+2) \ge M+1 \ge M$ and $c_u(t+2) \le 2$. Exterior arrivals at $u$ only increase $L_u$. **Valid.**
- **Stage 3 (unit $[t+2, t+3]$):**
  - Let $\tau_1$ be the first ring time of $u$ after $t+2$ (an $\mathcal{F}$-stopping time). Before $\tau_1$, $u$ has not rung since $t+2$, so $L_u(\tau_1^-) \ge M$ and $c_u(\tau_1^-) \le 2$. On $\{\tau_1 \le t+3\}$, $u$ chooses a leaf $v$ with conditional probability $L_u(\tau_1^-)/(L_u(\tau_1^-)+c_u(\tau_1^-)) \ge M/(M+2)$, creating a new leaf $v_1$ attached to $v$, while preserving $D_u(\tau_1) = D_u(\tau_1^-) \ge M$.
  - Conditioned on $\mathcal{F}_{\tau_1}$ on this event, the clocks of $u, v, v_1$ on $(\tau_1, t+3]$ are independent rate-$1$ Poisson processes, so the probability that all three remain silent on $(\tau_1, t+3]$ is $e^{-3(t+3-\tau_1)}$. Integrating against the $\mathrm{Exp}(1)$ density $e^{-(\tau_1-(t+2))}$ on $(t+2, t+3]$ gives conditional probability $\ge \frac{M}{M+2} e^{-1} e^{-2} > 0$. Because $u, v, v_1$ are all silent on $(\tau_1, t+3]$, no exterior or internal birth touches $\{v, v_1\}$ before $t+3$. Also, $v$ was a leaf at $\tau_1^- > t$ and $v_1$ was born after $t$, so neither was core at $t$. **Valid.**
- **Stage 4 ($k-2$ unit subintervals $[t+3, t+k+1]$):**
  - Since $T$ is a tree on $k \ge 2$ vertices, the root $v$ has at least one neighbour $w_1$ in $T$, matched to $\{v, v_1\}$. Order the remaining $k-2$ vertices of $T \setminus \{v, w_1\}$ in parent-before-child order $w_2, \dots, w_{k-1}$.
  - At step $j \in \{2, \dots, k-1\}$ (branch size $m = j \le k-1$), the desired attachment vertex $y_j$ already lies in the connected branch of size $m \ge 2$, hence has an internal neighbour $z_j$ in the branch. Keeping $u$ silent throughout all $k-2$ units completely isolates the branch from exterior arrivals (the only boundary vertex is $u$). In the $(j-1)$-th unit subinterval, requiring $z_j$ to ring *exactly once*, choose $y_j$, all other $m-1$ internal vertices to stay silent for the unit, and the newborn $w_j$ to stay silent for the remainder of the unit gives conditional probability at least $\frac{1}{m} e^{-(m+2)} \ge \frac{1}{k+1} e^{-(k+1)} > 0$ (indeed $D_{z_j} \le m \le k-1 < k+1$ even when $z_j = v$).
  - Total time $h = k + 1$ and lower bound $p = p(M, k) > 0$ depend **only** on $M$ and $k$, uniformly over all $(T, v)$ of size $k$, all eligible $u$, all $t \ge 0$, and all exterior configurations. **Valid.**

### 2.3 Pathwise disjointness and no independence needed (line 23, lines 48–52)
- Suppose $u_1 \ne u_2$ are distinct core vertices at time $t$ whose plantings succeed at $t+h$, creating cut edges $u_1 v_1$ and $u_2 v_2$ with branches $B_1(t+h)$ and $B_2(t+h)$ containing no vertices that were core at $t$. For any $s \ge 0$, since $T(t+h)$ is a subtree of $T(t+h+s)$, $u_2$ lies on the $u_1$-side of $u_1 v_1$ (i.e., $u_2 \notin B_1(t+h+s)$) and $u_1 \notin B_2(t+h+s)$. If some vertex $x$ belonged to $B_1(t+h+s) \cap B_2(t+h+s)$, the unique simple path $u_2 \leadsto u_1$ in $T(t+h+s)$ avoids $B_1 \cup B_2$, and concatenating it with $u_1 \to v_1 \leadsto x$ (which lies in $B_1$ after $u_1$ and therefore avoids $u_2$) yields a simple path from $u_2$ to $x$ not traversing $u_2 v_2$, contradicting $x \in B_2(t+h+s)$.
- Therefore, pathwise on every realization, $C(t+h+s) \ge \sum_{u \in \mathcal{E}(t)} \mathbf{1}_{E_u} C^{\mathrm{int}}(B_u(t+h+s))$, since every vertex of internal degree $\ge 2$ in $B_u$ has degree $\ge 2$ in $T(t+h+s)$. By linearity of conditional expectation given $\mathcal{F}_t$, no independence among the events $\{E_u\}$ (even when two type-B vertices share a type-B neighbour $w$ in Stage 1, making their Stage 1 events mutually exclusive) is ever needed. **Valid.**

---

## 3. Clock Coupling Between Planted Branch and Freestanding Copy (lines 25–34)

Condition on $\mathcal{F}_{t+h}$ on the event $E_u$ (which is $\mathcal{F}_{t+h}$-measurable), so at time $t+h$ the edge $uv$ separates a branch $B_0 \cong (T, v)$ from an exterior containing $u$ with $D_u(t+h) \ge M$.
- Construct the freestanding reference process $(T_{\mathrm{ref}}(r))_{r \in [0, s]}$ from $T$ using independent rate-$1$ clocks and i.i.d. $\mathrm{Unif}[0,1)$ internal choice marks $\mathcal{R}_{\mathrm{ref}}$; attach an independent family $\mathcal{A}_v$ of i.i.d. $\mathrm{Unif}[0,1)$ acceptance marks to the rings of $v$; and generate all exterior clocks and marks $\mathcal{R}_{\mathrm{ext}}$ independently of $(\mathcal{R}_{\mathrm{ref}}, \mathcal{A}_v)$.
- Up to the stopping time $\sigma$ of the first discrepancy ($u$ chooses $v$ or $v$ chooses $u$), the internal branch $B_{\mathrm{act}}(r)$ equals $T_{\mathrm{ref}}(r)$:
  - Every internal vertex $x \ne v$ has the exact same neighbours in $T_{\mathrm{act}}$ and $T_{\mathrm{ref}}$.
  - At a ring of $v$ with current internal degree $D_v = D_v^{\mathrm{ref}}(r^-)$, $v$ has total degree $D_v + 1$ in $T_{\mathrm{act}}$. Accepting the uniform internal choice with probability $D_v/(D_v+1)$ (via $A_i < D_v/(D_v+1)$) and otherwise selecting $u$ produces the exact uniform choice over all $D_v + 1$ neighbours of $v$.
  - At a ring of $u$, ranking $uv$ first with independent mark $W \sim \mathrm{Unif}[0,1)$, since $D_u^{\mathrm{act}}(r^-) \ge D_u(t+h) \ge M$ forever, $W \ge 1/M \implies W \ge 1/D_u^{\mathrm{act}}(r^-)$, so $u$ cannot choose $v$. The event $G_u = \{\text{no ring of } u \text{ on } (t+h, t+h+s] \text{ has } W \in [0, 1/M]\}$ has exact probability $e^{-s/M}$ and is independent of $(\mathcal{R}_{\mathrm{ref}}, \mathcal{A}_v)$.
  - After $\sigma$ (if $\sigma \le s$), decoupling $T_{\mathrm{act}}$ preserves the exact true law of $T_{\mathrm{act}}$.
- On $G_u \cap G_{\mathrm{acc}}$, $\sigma > s$ and $B_{\mathrm{act}}(s) = T_{\mathrm{ref}}(s)$, so every reference core vertex is an actual core vertex in $B_{\mathrm{act}}(s)$. Taking expectations and conditioning on $\mathcal{R}_{\mathrm{ref}}$ gives the exact lower bound
  $$\mathbb{E}[C^{\mathrm{int}}(B_{\mathrm{act}}(s)) \mid \mathcal{F}_{t+h}] \ge e^{-s/M} \, \mathbb{E}\!\left[ C_T(s) \prod_{\text{rings of } v \text{ in } [0,s]} \frac{D_v}{D_v + 1} \right].$$
**Valid.**

---

## 4. High-Degree Vertex Estimate and Quantifiers (lines 35–43)

- **Chernoff constants (checked numerically in Python):**
  - $A = 8$: $A \log A - A + 1 = 8 \ln 8 - 7 \approx 9.6355 > 5$.
  - $b(\epsilon) = 7/\log(1/\epsilon)$: $b(\epsilon) \to 0$ as $\epsilon \downarrow 0$, $b(\epsilon) > \epsilon$ for small $\epsilon$, and $b \log(b/\epsilon) - b + \epsilon = 7 - \frac{7(\ln\ln(1/\epsilon) + 1 - \ln 7)}{\ln(1/\epsilon)} + \epsilon \to 7 > 5$ (numerically $6.126, 6.001, 6.032, 6.149$ at $\epsilon = 10^{-2}, 10^{-3}, 10^{-4}, 10^{-6}$).
  - For $X \sim \mathrm{Pois}(as)$ and $x > a$, $\mathbb{P}(X \ge xs) \le e^{-s(x\log(x/a)-x+a)}$ gives $\le e^{-5s}$ per original vertex for both $[0, \epsilon s]$ (threshold $bs$) and $[0, s]$ (threshold $As$).
- **Degree lower tail (line 39):**
  - In any $k$-vertex tree ($k \ge 2$), each original vertex $w$ has at least one original neighbour $n(w)$, which can be ranked first in $w$'s internal choice ordering using the internal mark $U_{w,i}$ (separate from the acceptance mark $A_{v,i}$ at $v$, so both orderings coexist without altering any marginal law).
  - If all $k$ original degrees remain $< M$ throughout $[0, \epsilon s]$, then at every ring of an original vertex $w$ with $U_{w,i} \in [0, 1/M]$, we have $1/M < 1/D_w^{\mathrm{ref}}$, so $w$ necessarily selects $n(w)$ and increases $\sum_{x \in V(T)} D_x^{\mathrm{ref}}$ by $+1$. Since the sum starts at $2(k-1) > 0$ and remains $< kM$, there can be at most $kM$ such marks on $[0, \epsilon s]$.
  - The total count $Y$ of such marks across the $k$ independent original clocks on $[0, \epsilon s]$ is $\mathrm{Pois}(k \epsilon s / M)$ with mean $\ge 6s$ (since $k \ge 6M/\epsilon$). By stochastic monotonicity of Poisson distributions, $\mathbb{P}(Y \le kM) \le \mathbb{P}(\mathrm{Pois}(6s) \le kM) = e^{-6s}\sum_{j=0}^{kM} (6s)^j/j! \le e^{-4s}$ for all sufficiently large $s$ (with $k, M$ fixed).
  - Combining the $2k$ upper tails and the lower tail gives $\mathbb{P}(\text{bad}) \le 2k e^{-5s} + e^{-4s} = O_{k,M,\epsilon}(e^{-4s})$, **uniform over all $k$-vertex trees $T$**. **Valid.**
- **Cauchy–Schwarz & acceptance weight (lines 41–42):**
  - $\mathbb{E}[C_T(s) \mathbf{1}_{\mathrm{bad}}] \le \sqrt{\mathbb{E}[N_T(s)^2]\,\mathbb{P}(\mathrm{bad})} \le \sqrt{k(k+1)}\, e^s \cdot O_{k,M,\epsilon}(e^{-2s}) = O_{k,M,\epsilon}(e^{-s})$, uniform in $T$.
  - On any good history assigned to original vertex $v$: before $\epsilon s$, $v$ rings $\le bs$ times and $D_v \ge 1 \implies D_v/(D_v+1) \ge 1/2 = e^{-\log 2}$; after $\epsilon s$, $v$ rings $\le As$ times and $D_v \ge M \implies D_v/(D_v+1) \ge M/(M+1) = (1+1/M)^{-1} \ge e^{-1/M}$ (since $e^{1/M} \ge 1 + 1/M$). Multiplying by $e^{-s/M}$ yields weight $\ge e^{-(b\log 2 + (A+1)/M)s} = e^{-\eta s}$.
  - Partitioning the good event by the assigned $v \in V(T)$ (fixed tie-breaking) and taking the best $v = v(T, s)$ gives $\mathbb{E}[C_T(s) \mathbf{1}_{G_v}] \ge (\mathbb{E}[C_T(s)] - O_{k,M,\epsilon}(e^{-s}))/k$. Because the comparison lower bound is an expectation of a nonnegative random variable on the reference space, restricting to $G_v \in \sigma(\mathcal{R}_{\mathrm{ref}})$ is valid and requires choosing $(T, v)$ only from the laws before planting. **Valid.**

---

## 5. Iteration, Quantifier Order, and De-Poissonization (lines 44–58)

- **Stopping at size $k$ (lines 46–47):** Let $\tau_k = \inf\{r \ge 0 : N(r) = k\}$. On $\{\tau_k > s\}$, $C(s) < k$. On $\{\tau_k \le s\} \in \mathcal{F}_{\tau_k}$, $s \le \tau_k + s$, so monotonicity of $C$ and the strong Markov property at $\tau_k$ give $\mathbb{E}[C(s)\mathbf{1}_{\{\tau_k \le s\}}] \le \mathbb{E}[\mathbf{1}_{\{\tau_k \le s\}}\mathbb{E}[C(\tau_k+s)\mid \mathcal{F}_{\tau_k}]] \le \max_{|V(T)|=k}\mathbb{E}[C_T(s)]$. Hence $F(s) \le k + \max_T \mathbb{E}[C_T(s)]$. **Valid.**
- **Supermultiplicative inequality (lines 48–52):** Choosing $T$ maximizing $\mathbb{E}[C_T(s)]$ and its pigeonhole root $v$, and summing over eligible vertices at time $t$ (note $|\mathcal{E}(t)| \ge C(t)/4$ holds trivially when $C(t)=0$ too) gives, for all $t \ge 0$ and all $s \ge s_1(k, M, \epsilon)$ (so that $F(s) - k - O_{k,M,\epsilon}(e^{-s}) > 0$):
  $$F(t + h + s) \ge F(t)\, B_s, \qquad B_s := \frac{p}{4k} e^{-\eta s} \big(F(s) - k - O_{k,M,\epsilon}(e^{-s})\big).$$
  **Valid.**
- **Divergence $F(s) \to \infty$ (line 54):** For $C = 1$, $\sum_u L_u/D_u = 1$; for $C \ge 2$, the core tree has at least two leaves ($c_u = 1, L_u \ge 1$), each contributing $L_u/(L_u+1) \ge 1/2$, so $\sum_u L_u/D_u \ge 1$. (At $N = 2$, the jump rate from $C = 0$ to $C = 1$ is $2 > 1$, so in fact $F'(r) \ge 1$ for all $r \ge 0$ and $F(s) \ge s \to \infty$.) **Valid.**
- **Block iteration and quantifier order (lines 54–56):**
  - Fix $M \ge 2$ and small $\epsilon > 0$; this fixes $A=8$, $b(\epsilon)$, $k \ge 6M/\epsilon$, $h = k+1$, $p(M,k)$, and $\eta = b(\epsilon)\log 2 + 9/M$.
  - Fix any large $s \ge s_1(k,M,\epsilon)$ and any $t_0 > 0$ (where $F(t_0) \ge 1 - e^{-2t_0} > 0$). Iterating gives $F(t_0 + m(h+s)) \ge F(t_0) B_s^m$ for all $m \ge 0$. For $t \ge t_0$, setting $m = \lfloor (t-t_0)/(h+s) \rfloor$ and using monotonicity of $F$ gives $t^{-1}\log F(t) \ge t^{-1}\log F(t_0) + (m/t)\log B_s$; since $m/t \to (h+s)^{-1}$ regardless of the sign of $\log B_s$, $\underline{\mu} := \liminf_{t\to\infty} t^{-1}\log F(t) \ge (h+s)^{-1}\log B_s$.
  - Because $s \le F(s) \le 2e^s$, both $\underline{\mu}$ and $\bar{\mu} := \limsup_{s\to\infty} s^{-1}\log F(s)$ lie in $[0, 1]$. Taking $s_j \to \infty$ with $s_j^{-1}\log F(s_j) \to \bar{\mu}$, we have $F(s_j) - k - O(e^{-s_j}) = (1-o(1))F(s_j)$ and $h, p, k$ fixed, so $(h+s_j)^{-1}\log B_{s_j} \to \bar{\mu} - \eta$, giving $\underline{\mu} \ge \bar{\mu} - \eta$.
  - Crucially, $k$ does not appear in $\eta = b(\epsilon)\log 2 + (A+1)/M$. Sending $M \to \infty$ and $\epsilon \downarrow 0$ yields $\underline{\mu} = \bar{\mu} =: \mu \in [0, 1]$. **Valid; quantifier order is completely non-circular.**
- **De-Poissonization (line 58):**
  - Independence of $N(t)$ and $(T_n)_{n\ge 2}$ gives $F(t) = \mathbb{E}[a_{N(t)}]$, where $N(t) = G_1 + G_2$ for i.i.d. $G_1, G_2 \sim \mathrm{Geom}(e^{-t})$ on $\{1,2,\dots\}$ (verified by Yule forward ODE).
  - At $t = \log n$ ($n \ge 2$), $\mathbb{P}(N(\log n) \ge n) \ge \mathbb{P}(G_1 \ge n) = (1 - 1/n)^{n-1} \ge e^{-1}$, so monotonicity of $a_n$ gives $a_n \le e F(\log n)$, hence $\limsup_n (\log a_n)/\log n \le \mu$.
  - At $t = (1-\alpha)\log n$ ($0 < \alpha < 1$, $q = n^{-(1-\alpha)}$), using $a_m \le a_n$ for $m \le n$ and $a_m \le m$ for $m > n$ gives $F(t) \le a_n + \mathbb{E}[N(t)\mathbf{1}_{\{N(t)>n\}}]$. The exact tail equals $(1-q)^{n-1}(n(n+1)q + 2(n+1)(1-q) + 2(1-q)^2/q) \le O(n^2)\exp(-(1-o(1))n^\alpha)$ (checked numerically in Python: $3.565 \times 10^{-2}$ at $\alpha=0.5, n=100$ and $3.851 \times 10^{-10}$ at $\alpha=0.5, n=1000$), which decays faster than any power of $n$ and is $o(1)$. Since $F((1-\alpha)\log n) \to \infty$, $a_n \ge (1-o(1))F((1-\alpha)\log n)$ for large $n$, giving $\liminf_n (\log a_n)/\log n \ge (1-\alpha)\mu$. Letting $\alpha \downarrow 0$ yields $\lim_n (\log a_n)/\log n = \mu$. **Valid.**

---

## 6. Minor Expositional Notes (Non-blocking)

None of the following are mathematical gaps; they are tiny wording clarifications that could improve readability:
1. Line 14 ("and a set of all rooted unlabelled trees $(T,v)$ on $k \ge 2$ vertices"): clearer as "Fix integers $M \ge 2$ and $k \ge 2$, and consider the finite set of all rooted unlabelled trees $(T,v)$ on $k$ vertices."
2. Line 21 ("prescribe one ring of an internal neighbour"): implicitly means *exactly one* ring in that unit subinterval (as explicitly stated in Stage 3). Also, the chosen source has degree at most $k-1$ (so $\le k+1$ holds with room to spare).
3. Line 41 ("bad histories contribute $O_k(e^{-s})$"): noting explicitly that the implied constant is uniform over all $k$-vertex trees $T$ (and depends only on $k, M, \epsilon$) makes the subsequent maximization over $T$ immediately transparent.
4. Line 58: noting explicitly that $F((1-\alpha)\log n) \to \infty$ (from line 54) is what lets one drop the $+o(1)$ inside $\log$ even without assuming $\mu > 0$.
