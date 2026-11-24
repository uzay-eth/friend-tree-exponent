# Generator and finite-window ratio identities

**Author:** Fable lead (Claude Fable 5.2; retained in the final Claude Opus 5.5 report). **Pinned source:** `a18825ccdecfdb8d8465a89a7d99fc33734fa8e9`. **Review and dependencies:** Astra (GPT) independently derived the generator and window identities. See [F1](README.md#f1) for the finite-window assumptions and review scope.

### 2.1 Continuous-time form and the generator identity

Give every vertex an independent rate-1 clock; when the clock of x rings, x picks a uniform neighbour y and y
receives a new leaf. The jump chain is the random friend tree and the size is a Yule process, independent of
the jump chain. Only the core (vertices of degree >= 2) and the number of leaves on each core vertex matter:
for a core vertex v write L_v for its number of leaf neighbours, c_v for its number of core neighbours,
D_v = L_v + c_v. A ring changes this structure in exactly one of two ways:

* **G_v, "v gains a leaf":** a leaf of v rings (L_v leaves, each picks v with probability 1), or a core
  neighbour j of v rings and picks v (probability 1/D_j). Total rate g_v = L_v + sum_{j ~ v, core} 1/D_j.
* **K_v, "v converts a leaf":** v rings and picks one of its leaves (probability L_v/D_v). That leaf becomes a
  core vertex with one leaf and the single core neighbour v; v now has L_v - 1 leaves and c_v + 1 core
  neighbours. Rate k_v = L_v/D_v.

(For n >= 3 every leaf hangs on a core vertex, so there is no other case.) Let F(u) be a vector-valued function
of the neighbourhood of radius 2 of a core vertex u, of polynomial growth in the numbers of core neighbours, and
Phi(T) = sum over core u of F(u). Then by Dynkin's formula (the tree is a non-explosive Markov chain whose size
has all moments)

    d/dt E[Phi(T(t))] = E[(G Phi)(T(t))],   (G Phi)(T) = sum_v g_v [Phi(T after G_v) - Phi(T)] + sum_v k_v [Phi(T after K_v) - Phi(T)].   (1)

The right-hand side is a finite explicit sum over the tree. The code computes it as sum_u a_u with
a_u = sum over events of rate x [F_after(u) - F(u)] + k_u F(newborn child of u) (the same terms, grouped by the
vertex whose feature changes), and b_u = F(u). With F = 1: Phi = C and G Phi = S := sum_u L_u/D_u, the identity
d/dt E[C] = E[S] used in round 1.

### 2.2 What a ratio of such sums estimates

Fix a weight vector theta and put den(t) = sum_u theta.b_u, num(t) = sum_u theta.a_u, so E[num(t)] = d/dt E[den(t)]
by (1). Write beta(t) = E[num(t)] - mu E[den(t)]. Then for t1 < t2 with E den(t2) != E den(t1), identically,

    (E num(t2) - E num(t1)) / (E den(t2) - E den(t1)) = mu + (beta(t2) - beta(t1)) / (E den(t2) - E den(t1)).   (2)

So the time-differenced ratio equals mu up to the change of beta across the window divided by the growth of
E[den]. If E[den(t)] = A e^{mu t} + B + r(t) with A != 0 and r(t), r'(t) -> 0, then beta tends to the constant -mu B
and, for windows of fixed length, the correction in (2) tends to zero faster than 1/E[den(t2)]. Two caveats, both
stated by the cold check: (i) A != 0 is needed (an observable supported on configurations of vanishing frequency,
such as the indicator that the initial vertex is still the only core vertex, has ratio -1 at all times); for the
observables sane and full, den/C is 0.783 and 0.865 at each of the checkpoints 12, 16, 20, 24 (`data/offsets_main.txt`).
(ii) The size of the correction in (2) is measured (2.5, 2.7), not proved.

The estimator is the sample version of (2) over independent runs:
mu_hat = sum_runs [num(t2) - num(t1)] / sum_runs [den(t2) - den(t1)], with a delete-one-run jackknife error.
For theta = (plain equation only) this is the differenced S/C estimator. Its standard error is about 320 times
larger than for the fitted theta below (2.5e-6 against 7e-9 on the primary window, section 3.1).
