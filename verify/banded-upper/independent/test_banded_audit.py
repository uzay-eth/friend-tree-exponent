"""Independent exact-rational (fractions.Fraction) audit and negative-control suite
for the multi-band upper certificate specification (attempts/strogatz-astra/banded/BOUND.md)
and our standalone C++ verifier (verify_banded.cpp).

Tests:
1. Exact-tree brute-force verification of the multi-band state machine & BOTH cancellations:
   - Band-dependent cutoffs h_c(D) with central gains crossing band boundaries
   - Witness preservation on c=1 cross-band gains, witness init at birth, witness swap on w-conversion
   - Signed edge transfer cancellation g(D,c; d,k) across all mutual-observable pairs
   - Nonnegative arrow transfer cancellation a(D,c) with threshold T, exact vs pooled receiver
     (b < T vs b >= T), and consistent-type maxima at finite senders (exact, saturated k=K_c, pooled)
   - Analytic tail f(L, D) = A + B*(L/D) for D >= H with tail arrow bound a_star*(1-q)
   - Checked on ALL 436 histories n=3..7 plus 55 random & adversarial multi-band trees!
2. Wrong-rate / broken-cancellation negative controls in the exact-tree model (all must fail).
3. Black-box CLI negative controls on ./verify_banded using altered copies of smoke-recipe.txt
   and bands-smoke.txt (all must fail/reject!).
"""
import os, random, subprocess, sys, tempfile
from fractions import Fraction as Fr

WITNESS_CAP = 6
BRIDGE_DEG = 14
WITNESS_CORE = 4

def k_limit(c):
    if c == 1: return 20
    if c == 2: return 7
    return 5

class BandedSpec:
    def __init__(self, bands):
        # bands: list of (d_end, [h1, h2, ...])
        self.bands = []
        d_start = 2
        for d_end, hs in bands:
            assert d_end > d_start
            self.bands.append((d_start, d_end, [2] + list(hs)))
            d_start = d_end
        self.H = d_start
        self.h1 = self.bands[0][2][1]

    def band_of(self, D):
        for bi, (ds, de, _) in enumerate(self.bands):
            if ds <= D < de:
                return bi
        raise ValueError(f"D={D} outside [2, {self.H})")

    def cap(self, D, c):
        if c <= 0: return 2
        hs = self.bands[self.band_of(D)][2]
        return hs[c] if c < len(hs) else 2

    def cap_b(self, b, c):
        if c <= 0: return 2
        hs = self.bands[b][2]
        return hs[c] if c < len(hs) else 2

    def pool_thresholds(self, b, c):
        hc = self.cap_b(b, c)
        th = [hc]
        m = 4
        while m < self.H:
            if m > hc:
                th.append(m)
            m *= 2
        if self.H > hc:
            th.append(self.H)
        return th

    def pool_mark(self, b, c, deg_lb):
        th = self.pool_thresholds(b, c)
        best = th[0]
        for x in th:
            if x <= deg_lb:
                best = x
        return best

    def label_plain(self, b, c, d, k):
        hc = self.cap_b(b, c)
        if d >= hc:
            return ('P', self.pool_mark(b, c, d), 0, 0, 0)
        return ('E', d, min(k, k_limit(c)), 0, 0)

    def label_bridge(self, b, c, d, k, wd, wk):
        if c != 1 or d > BRIDGE_DEG or d >= self.cap_b(b, 1) or k < 2 or k > WITNESS_CORE:
            return self.label_plain(b, c, d, k)
        if wd >= WITNESS_CAP:
            return ('E', d, k, WITNESS_CAP, 0)
        return ('E', d, k, wd, wk)

    def sort_key(self, lbl):
        kind, d, k, wd, wk = lbl
        if kind == 'E':
            if wd == 0: return (0, d, k, 0, 0)
            if wd >= WITNESS_CAP: return (0, d, k, 2, 0)
            return (0, d, k, 1, wd * 100 + wk)
        return (1, d, 0, 0, 0)

    def canonical_bag(self, b, c, items):
        hc = self.cap_b(b, c)
        lst = sorted(items, key=self.sort_key)
        out = list(lst)
        for j in range(len(out) - 1):
            if out[j][0] == 'P' and out[j][1] > hc:
                out[j] = ('P', hc, 0, 0, 0)
        return tuple(sorted(out, key=self.sort_key))

    def edge_observable(self, D, c, d, k):
        if not (2 <= D < self.H and 2 <= d < self.H and 1 <= c <= D and 1 <= k <= d):
            return False
        if d >= self.cap(D, c) or D >= self.cap(d, k):
            return False
        Kc, Kk = k_limit(c), k_limit(k)
        return (k < Kc or (k == Kc and d == Kc)) and (c < Kk or (c == Kk and D == Kk))


class BandedHistoryTree:
    def __init__(self, spec, adj, pool_lb=None, witness=None, witness_pooled=None):
        self.spec = spec
        self.adj = [list(a) for a in adj]
        self.pool_lb = pool_lb if pool_lb is not None else {}
        self.witness = witness if witness is not None else {}
        self.witness_pooled = witness_pooled if witness_pooled is not None else {}

    def core_vertices(self): return [u for u, a in enumerate(self.adj) if len(a) >= 2]
    def core_nbrs(self, u): return [w for w in self.adj[u] if len(self.adj[w]) >= 2]
    def deg(self, u): return len(self.adj[u])
    def cdeg(self, u): return sum(1 for w in self.adj[u] if len(self.adj[w]) >= 2)

    def nbr_label(self, u, w):
        D = self.deg(u)
        if D >= self.spec.H:
            return None
        b = self.spec.band_of(D)
        c = self.cdeg(u)
        plb = self.pool_lb.get(u, {}).get(w, None)
        if plb is not None:
            return ('P', self.spec.pool_mark(b, c, plb), 0, 0, 0)
        dw, kw = self.deg(w), self.cdeg(w)
        if dw >= self.spec.cap_b(b, c):
            return ('P', self.spec.pool_mark(b, c, dw), 0, 0, 0)
        if c == 1 and dw <= BRIDGE_DEG and 2 <= kw <= WITNESS_CORE:
            x = self.witness[u]
            assert x != u and x in self.adj[w] and self.deg(x) >= 2
            if self.witness_pooled.get(u, False) or self.deg(x) >= WITNESS_CAP:
                return self.spec.label_bridge(b, 1, dw, kw, WITNESS_CAP, 0)
            return self.spec.label_bridge(b, 1, dw, kw, self.deg(x), self.cdeg(x))
        return self.spec.label_plain(b, c, dw, kw)

    def state(self, u):
        D = self.deg(u); c = self.cdeg(u)
        if D >= self.spec.H:
            return ('TAIL', D, c)
        b = self.spec.band_of(D)
        raw = [self.nbr_label(u, w) for w in self.core_nbrs(u)]
        return ('FIN', D, c, self.spec.canonical_bag(b, c, raw))

    @staticmethod
    def init_from_tree(spec, adj):
        ht = BandedHistoryTree(spec, adj)
        for u in ht.core_vertices():
            D, c = ht.deg(u), ht.cdeg(u)
            if D >= spec.H: continue
            b = spec.band_of(D)
            hc = spec.cap_b(b, c)
            cn = ht.core_nbrs(u)
            pooled = [(w, spec.pool_mark(b, c, ht.deg(w))) for w in cn if ht.deg(w) >= hc]
            if pooled:
                ht.pool_lb[u] = {}
                max_w = max(pooled, key=lambda pair: (pair[1], pair[0]))[0]
                for w, mk in pooled:
                    ht.pool_lb[u][w] = mk if w == max_w else hc
            if c == 1:
                w = cn[0]; dw, kw = ht.deg(w), ht.cdeg(w)
                if dw < hc and dw <= BRIDGE_DEG and 2 <= kw <= WITNESS_CORE:
                    others = [x for x in ht.core_nbrs(w) if x != u]
                    x = min(others, key=lambda x: (ht.deg(x), ht.cdeg(x), x))
                    ht.witness[u] = x
                    ht.witness_pooled[u] = (ht.deg(x) >= WITNESS_CAP)
        return ht

    def _project_pool_dict(self, b_new, c_new, old_dict, nbrs):
        hc_new = self.spec.cap_b(b_new, c_new)
        nd = {}
        for x in nbrs:
            if x in old_dict:
                nd[x] = self.spec.pool_mark(b_new, c_new, old_dict[x])
            elif self.deg(x) >= hc_new:
                nd[x] = self.spec.pool_mark(b_new, c_new, self.deg(x))
        above = [(mk, x) for x, mk in nd.items() if mk > hc_new]
        if len(above) > 1:
            above.sort()
            keep_x = above[-1][1]
            for x in list(nd.keys()):
                if x != keep_x and nd[x] > hc_new:
                    nd[x] = hc_new
        return nd

    def step(self, W):
        spec = self.spec
        n = len(self.adj)
        new_adj = [list(a) for a in self.adj] + [[W]]
        new_adj[W].append(n)
        new_plb = {u: dict(d) for u, d in self.pool_lb.items()}
        new_wit = dict(self.witness)
        new_wpool = dict(self.witness_pooled)

        if len(self.adj[W]) >= 2:
            # Core vertex W gains a leaf: D_W -> D_W + 1, c_W unchanged.
            Dw_old = self.deg(W)
            Dw_new = Dw_old + 1
            cW = self.cdeg(W)
            # 1. Update W's own state if crossing a band boundary or entering tail!
            if Dw_new >= spec.H:
                new_plb.pop(W, None); new_wit.pop(W, None); new_wpool.pop(W, None)
            elif Dw_old < spec.H:
                b_old = spec.band_of(Dw_old)
                b_new = spec.band_of(Dw_new)
                if b_new != b_old:
                    if cW >= 2:
                        nd = self._project_pool_dict(b_new, cW, self.pool_lb.get(W, {}), self.core_nbrs(W))
                        if nd: new_plb[W] = nd
                        else: new_plb.pop(W, None)
                    # Note: for cW == 1, h_1 is constant across bands so pooled marks and witnesses are preserved!
            # 2. Update every core neighbour u of W
            for u in self.core_nbrs(W):
                Du = self.deg(u)
                if Du >= spec.H: continue
                cu = self.cdeg(u)
                hc = spec.cap(Du, cu)
                if u not in new_plb or W not in new_plb[u]:
                    if Dw_new >= hc:
                        new_plb.setdefault(u, {})[W] = hc
                if cu == 1 and (Dw_new > BRIDGE_DEG or Dw_new >= hc):
                    new_wit.pop(u, None); new_wpool.pop(u, None)
            # 3. Distance-2 witness update
            for u, x in list(new_wit.items()):
                if x == W and Dw_new >= WITNESS_CAP:
                    new_wpool[u] = True
        else:
            # Leaf W of core vertex p is promoted to newborn core (L=1, c=1, D=2)!
            p = self.adj[W][0]
            Dp = self.deg(p)
            cp_old = self.cdeg(p)
            cp_new = cp_old + 1
            # 1. Newborn W init (in band 0, c=1)
            if Dp >= spec.h1:
                new_plb[W] = {p: spec.pool_mark(0, 1, Dp)}
            elif Dp <= BRIDGE_DEG and 1 <= cp_old < WITNESS_CORE:
                p_nbrs = self.core_nbrs(p)
                lbls = [(spec.sort_key(self.nbr_label(p, x)), x) for x in p_nbrs]
                lbls.sort()
                chosen_x = lbls[0][1]
                new_wit[W] = chosen_x
                was_pooled = chosen_x in self.pool_lb.get(p, {})
                new_wpool[W] = was_pooled or (self.deg(chosen_x) >= WITNESS_CAP)
            # 2. Update p itself (c: cp_old -> cp_new, same Dp -> same band)
            new_wit.pop(p, None); new_wpool.pop(p, None)
            if Dp < spec.H:
                bp = spec.band_of(Dp)
                # Include newborn W (degree 2) in projection in case h_{cp_new}(bp) == 2!
                nd = self._project_pool_dict(bp, cp_new, self.pool_lb.get(p, {}), self.core_nbrs(p))
                if spec.cap_b(bp, cp_new) <= 2:
                    nd[W] = 2 # newborn (2,1) is pooled at mark 2 when h == 2!
                if nd: new_plb[p] = nd
                else: new_plb.pop(p, None)
            # 3. Update every core neighbour u of p
            for u in self.core_nbrs(p):
                Du = self.deg(u)
                if Du >= spec.H: continue
                if self.cdeg(u) == 1 and Dp < spec.h1 and Dp <= BRIDGE_DEG:
                    if 2 <= cp_new <= WITNESS_CORE:
                        new_wit[u] = W; new_wpool[u] = False
                    else:
                        new_wit.pop(u, None); new_wpool.pop(u, None)

        return BandedHistoryTree(spec, new_adj, new_plb, new_wit, new_wpool)


def verify_banded_drift(ht, f_fn, g_edge, a_arrow, T, tail_A, tail_B,
                        bug_asym_edge=False, bug_drop_sat_arrow=False, bug_bad_pool_arrow=False):
    """Exact rational verification of brute-force drift == summed exact amended drifts <= upper bound."""
    spec = ht.spec
    n = len(ht.adj)
    core = ht.core_vertices()
    phi0 = sum(f_fn(ht.state(u)) for u in core)

    bf_drift = Fr(0)
    for V in range(n):
        dV = len(ht.adj[V])
        for W in ht.adj[V]:
            ht2 = ht.step(W)
            phi1 = sum(f_fn(ht2.state(u)) for u in ht2.core_vertices())
            bf_drift += Fr(1, dV) * (phi1 - phi0)

    a_star = max([Fr(0)] + list(a_arrow.values()))
    edge_net = Fr(0)
    arrow_net = Fr(0)
    sum_exact = Fr(0)
    sum_upper = Fr(0)

    for u in core:
        st_u = ht.state(u)
        fu = f_fn(st_u)
        Du = ht.deg(u); cu = ht.cdeg(u); Lu = Du - cu
        cn = ht.core_nbrs(u)
        rho_u = sum(Fr(1, ht.deg(w)) for w in cn)
        r_u = Lu + rho_u

        if Du >= spec.H:
            # Analytic tail check
            q = Fr(Lu, Du)
            # True outgoing arrow charge at tail sender (if Du >= T)
            out_a_true = sum(a_arrow.get((ht.deg(w), ht.cdeg(w)), Fr(0)) for w in cn) * Fr(1, Du) if Du >= T else Fr(0)
            arrow_net += out_a_true
            # True drift at tail vertex u
            fu_gain = f_fn(ht.step(u).state(u))
            row_ex = r_u * (fu_gain - fu) + out_a_true
            if Lu > 0:
                lv = next(v for v in ht.adj[u] if len(ht.adj[v]) == 1)
                htc = ht.step(lv)
                row_ex += Fr(Lu, Du) * (f_fn(htc.state(u)) + f_fn(htc.state(lv)) - fu)
            # Analytic upper bound: (B + a_star)*(1 - q) + s * q
            s_H = f_fn(('FIN', 2, 1, (('P', spec.H, 0, 0, 0),)))
            row_up = (tail_B + a_star) * (1 - q) + s_H * q
            assert row_ex <= row_up, f"Tail upper bound violated at u={u}: {row_ex} > {row_up}"
            sum_exact += row_ex
            sum_upper += row_up
            continue

        # Finite vertex u (Du < H)
        fu_gain = f_fn(ht.step(u).state(u))
        out_edge = Fr(0)
        for w in cn:
            lbl = ht.nbr_label(u, w)
            if lbl[0] == 'E' and (lbl[2] < k_limit(cu) or (lbl[2] == k_limit(cu) and lbl[1] == k_limit(cu))):
                ev = g_edge.get((lbl[1], lbl[2], Du, cu), Fr(0))
                if bug_asym_edge and ev != 0:
                    ev += Fr(1)
                out_edge += ev
        U = (fu_gain - fu) + out_edge
        edge_net += r_u * out_edge

        a_here = a_arrow.get((Du, cu), Fr(0))
        # Incoming arrow subtractions (exact) vs upper bound
        row_ex = Lu * U
        row_up = Lu * U
        for w in cn:
            dw_t = ht.deg(w)
            sub_true = a_here if dw_t >= T else Fr(0)
            arrow_net -= sub_true * Fr(1, dw_t)
            row_ex += (U - sub_true) * Fr(1, dw_t)

            lbl = ht.nbr_label(u, w)
            if lbl[0] == 'E':
                assert lbl[1] == dw_t
                row_up += (U - sub_true) * Fr(1, dw_t)
            else:
                b_mk = lbl[1]
                assert b_mk <= dw_t
                if bug_bad_pool_arrow:
                    # BUG: subtract a_here even when b_mk < T!
                    adj = U - a_here
                else:
                    adj = U - (a_here if b_mk >= T else Fr(0))
                row_up += max(Fr(0), adj) * Fr(1, b_mk)

        # Outgoing arrow additions at sender Du >= T
        if Du >= T:
            for w in cn:
                a_true = a_arrow.get((ht.deg(w), ht.cdeg(w)), Fr(0))
                arrow_net += a_true * Fr(1, Du)
                row_ex += a_true * Fr(1, Du)

                lbl = ht.nbr_label(u, w)
                if lbl[0] == 'E':
                    _, dw, kw, _, _ = lbl
                    if kw < k_limit(cu) or bug_drop_sat_arrow:
                        a_mx = a_arrow.get((dw, kw), Fr(0))
                    else:
                        a_mx = max((a_arrow.get((dw, k2), Fr(0)) for k2 in range(k_limit(cu), dw + 1)), default=Fr(0))
                else:
                    b_mk = lbl[1]
                    a_mx = max((v for (dd, kk), v in a_arrow.items() if dd >= b_mk), default=Fr(0))
                assert a_true <= a_mx, f"Arrow max failed: true={a_true} > max={a_mx} for lbl={lbl}"
                row_up += a_mx * Fr(1, Du)

        # Central conversion
        if Lu > 0:
            lv = next(v for v in ht.adj[u] if len(ht.adj[v]) == 1)
            htc = ht.step(lv)
            c_term = Fr(Lu, Du) * (f_fn(htc.state(u)) + f_fn(htc.state(lv)) - fu)
            row_ex += c_term
            row_up += c_term

        # Neighbour & witness moves
        for w in cn:
            lbl = ht.nbr_label(u, w)
            if lbl[0] == 'P':
                assert ht.step(w).state(u) == st_u
                continue
            _, dw, kw, wd, wk = lbl
            dw_t, kw_t = ht.deg(w), ht.cdeg(w)
            Lw_t = dw_t - kw_t
            rw_t = Lw_t + sum(Fr(1, ht.deg(y)) for y in ht.core_nbrs(w))

            in_e = Fr(0)
            if kw < k_limit(cu) or (kw == k_limit(cu) and dw == k_limit(cu)):
                in_e = g_edge.get((Du, cu, dw, kw), Fr(0))
            edge_net -= rw_t * in_e

            v_wg = (f_fn(ht.step(w).state(u)) - fu) - in_e
            row_ex += rw_t * v_wg
            if wd > 0:
                rw_lo = Fr(dw - kw) + Fr(1, Du) + (Fr(1, wd) if wk > 0 else Fr(0))
                rw_hi = Fr(dw - kw) + Fr(1, Du) + Fr(1, wd) + Fr(kw - 2, 2)
            else:
                rw_hi = Fr(1, Du) + Fr(2*dw - kw - 1, 2)
                rw_lo = Fr(dw - kw) + Fr(1, Du) if kw < k_limit(cu) else Fr(1, Du)
            assert rw_lo <= rw_t <= rw_hi
            row_up += (rw_hi if v_wg > 0 else rw_lo) * v_wg

            if Lw_t > 0:
                lw = next(v for v in ht.adj[w] if len(ht.adj[v]) == 1)
                if kw == k_limit(cu):
                    assert ht.step(lw).state(u) == st_u
                else:
                    cw = Fr(dw - kw, dw) * (f_fn(ht.step(lw).state(u)) - fu)
                    row_ex += cw
                    row_up += cw

            if wk > 0:
                x = ht.witness[u]
                Lx = wd - wk
                rx_t = Lx + sum(Fr(1, ht.deg(y)) for y in ht.core_nbrs(x))
                v_xg = f_fn(ht.step(x).state(u)) - fu
                row_ex += rx_t * v_xg
                rx_lo = Fr(wd - wk) + Fr(1, dw)
                rx_hi = Fr(wd - wk) + Fr(1, dw) + Fr(wk - 1, 2)
                assert rx_lo <= rx_t <= rx_hi
                row_up += (rx_hi if v_xg > 0 else rx_lo) * v_xg
                if Lx > 0:
                    lx = next(v for v in ht.adj[x] if len(ht.adj[v]) == 1)
                    cx = Fr(wd - wk, wd) * (f_fn(ht.step(lx).state(u)) - fu)
                    row_ex += cx
                    row_up += cx

        assert row_ex <= row_up, f"Finite row upper bound failed at u={u}: {row_ex} > {row_up}"
        sum_exact += row_ex
        sum_upper += row_up

    assert edge_net == 0, f"Edge cancellation failed: net={edge_net}"
    assert arrow_net == 0, f"Arrow cancellation failed: net={arrow_net}"
    assert sum_exact == bf_drift, f"Exact drift mismatch: {sum_exact} != {bf_drift}"
    assert bf_drift <= sum_upper, f"Upper bound failed: {bf_drift} > {sum_upper}"


def run_exact_tree_tests():
    # 3-band geometry with H=22, T=6 so cross-band gains (D=9->10, D=15->16), witness (d<=14),
    # pooled arrows (b<T and b>=T), saturated arrows (k=K_c), and tail (D>=22) ALL get hit!
    spec = BandedSpec([
        (10, [18, 12, 8, 5, 4]),
        (16, [18, 9,  6, 4, 3]),
        (22, [18, 6,  4, 3, 2]),
    ])
    T = 6
    rng = random.Random(777)
    tail_A, tail_B = Fr(12, 10), Fr(6, 10)
    f_cache = {}
    def f_fn(st):
        if st[0] == 'TAIL':
            _, D, c = st
            return tail_A + tail_B * Fr(D - c, D)
        if st not in f_cache:
            f_cache[st] = Fr(rng.randint(80, 250), 100)
        return f_cache[st]
    g_edge = {}
    for D in range(2, spec.H):
        for c in range(1, min(D, 5) + 1):
            for d in range(2, spec.cap(D, c)):
                for k in range(1, min(d, 5) + 1):
                    if spec.edge_observable(D, c, d, k):
                        g_edge[(D, c, d, k)] = Fr(rng.randint(-40, 40), 100)
    a_arrow = {}
    for D in range(2, spec.H):
        for c in range(1, min(D, 8) + 1):
            a_arrow[(D, c)] = Fr(rng.randint(0, 60), 100)

    # 1. All 436 histories n=3..7
    ht3 = BandedHistoryTree(spec, [[1], [0, 2], [1]])
    layer = [ht3]
    n_hist = 0
    for n in range(3, 8):
        for ht in layer:
            verify_banded_drift(ht, f_fn, g_edge, a_arrow, T, tail_A, tail_B)
            n_hist += 1
        if n < 7:
            layer = [ht.step(W) for ht in layer for W in range(n)]
    print(f"TREE PART 1 PASS: {n_hist} exhaustive histories n=3..7 verified (exact drift + edge & arrow cancellation)!")

    # 2. Random friend trees + adversarial multi-band trees
    adv = []
    for _ in range(45):
        ht = BandedHistoryTree(spec, [[1], [0, 2], [1]])
        for _ in range(rng.randint(5, 70)):
            V = rng.randrange(len(ht.adj))
            ht = ht.step(rng.choice(ht.adj[V]))
        adv.append(("rft", ht))

    # Adversarial tree crossing band 0->1 (D=9->10), band 1->2 (D=15->16), tail (D=21->22),
    # with c=1 cross-band witness, saturated k=K_3=5 sender arrow, and pooled b<T vs b>=T arrows!
    adj = []
    def av(): adj.append([]); return len(adj) - 1
    def lk(a, b): adj[a].append(b); adj[b].append(a)
    def al(v, k):
        for _ in range(k): lv = av(); lk(v, lv)

    hub = av() # will have D = 21 (H-1)
    # Vertex b01 at D=9, c=2 (band 0 top, crosses into band 1 on gain!)
    b01 = av(); lk(hub, b01)
    nb_p = av(); lk(b01, nb_p); al(nb_p, 10) # deg 11 (exact in band 0 h_2=12, pooled in band 1 h_2=9!)
    al(b01, 7) # D=9, c=2
    # Vertex b12 at D=15, c=1 with witness (crosses band 1->2 preserving witness!)
    w_br = av(); lk(hub, w_br)
    b12 = av(); lk(w_br, b12); al(b12, 14) # D=15, c=1!
    x_br = av(); lk(w_br, x_br); al(x_br, 3) # x_br has d=4, k=1
    al(w_br, 2) # w_br has d=5, k=3
    # Vertex u_sat at D=10 (>= T=6), c=3 with saturated core nbr w_sat (d=7, k_true=7 > K_3=5!)
    u_sat = av(); lk(hub, u_sat)
    w_sat = av(); lk(u_sat, w_sat)
    for _ in range(6):
        z = av(); lk(w_sat, z); al(z, 2)
    z2 = av(); lk(u_sat, z2); al(z2, 2)
    al(u_sat, 7) # D=10, c=3; h_3(band 1)=6 so make w_sat deg 5 or check band 0 (D=8, c=3, h_3=8)!
    # Also put a band-0 sender at D=8 (>=T=6), c=3 (h_3=8) with saturated w_sat2 (d=7, k_true=7 > K_3=5)!
    u_sat0 = av(); lk(hub, u_sat0)
    w_sat0 = av(); lk(u_sat0, w_sat0)
    for _ in range(6):
        z = av(); lk(w_sat0, z); al(z, 2)
    z3 = av(); lk(u_sat0, z3); al(z3, 2)
    al(u_sat0, 5) # D=8, c=3, h_3(band 0)=8 > 7!
    # Vertex u_pool_ltT at D=16, c=4 (h_4(band 2)=3 < T=6) with pooled nbr of degree 4 (< T) and degree 8 (>= T)!
    u_plt = av(); lk(hub, u_plt)
    p4 = av(); lk(u_plt, p4); al(p4, 3) # deg 4 (< T=6, pooled at mark 3 or 4 < T)
    p8 = av(); lk(u_plt, p8); al(p8, 7) # deg 8 (>= T=6, pooled at mark 8 >= T)
    p2 = av(); lk(u_plt, p2); al(p2, 1)
    al(u_plt, 12) # D=16, c=4!

    al(hub, max(0, spec.H - 1 - len(adj[hub])))
    ht_adv = BandedHistoryTree.init_from_tree(spec, adj)
    adv.append(("adv_multi_band_H-1", ht_adv))
    adv.append(("adv_cross_band_0_1", ht_adv.step(b01)))
    adv.append(("adv_cross_band_1_2_c1_witness", ht_adv.step(b12)))
    adv.append(("adv_enter_tail_H", ht_adv.step(hub)))
    adv.append(("adv_deep_tail_H+2", ht_adv.step(hub).step(hub)))

    for name, ht in adv:
        verify_banded_drift(ht, f_fn, g_edge, a_arrow, T, tail_A, tail_B)
    print(f"TREE PART 2 PASS: {len(adv)} random & adversarial multi-band trees verified identically!")

    # 3. Exact-tree negative controls (deliberate bugs MUST fail!)
    nc_tree = []
    # NC-T1: asymmetric edge transfer
    caught = False
    try:
        verify_banded_drift(ht_adv, f_fn, g_edge, a_arrow, T, tail_A, tail_B, bug_asym_edge=True)
    except AssertionError as e:
        caught = True
        nc_tree.append(("NC-T1_asym_edge", str(e)))
    assert caught, "NC-T1 MISSED"

    # NC-T2: omit max over k' >= K_c at saturated sender arrow (make a_arrow[(7,7)] huge!)
    a_big_sat = dict(a_arrow)
    a_big_sat[(7, 7)] = Fr(500, 1)
    a_big_sat[(7, 5)] = Fr(0, 1)
    caught = False
    try:
        verify_banded_drift(ht_adv, f_fn, g_edge, a_big_sat, T, tail_A, tail_B, bug_drop_sat_arrow=True)
    except AssertionError as e:
        caught = True
        nc_tree.append(("NC-T2_sat_sender_arrow", str(e)))
    assert caught, "NC-T2 MISSED"

    # NC-T3: subtract arrow on pooled neighbour with b < T (unsafe when true deg < T and U > 0!)
    # Minimal clean tree: u0 at D=16, c=4 (h_4(band 2)=3 < T=6) with 4 core nbrs each of degree 4 (< T),
    # so all 4 have pooled mark b=3 < T. True credit = 4*(U/4) = U, buggy upper = 4*max(0, U-a)/3 << U!
    adj_t3 = [[]]
    for _ in range(4):
        w = len(adj_t3); adj_t3.append([0]); adj_t3[0].append(w)
        for _ in range(3):
            lv = len(adj_t3); adj_t3.append([w]); adj_t3[w].append(lv)
    for _ in range(12):
        lv = len(adj_t3); adj_t3.append([0]); adj_t3[0].append(lv)
    ht_t3 = BandedHistoryTree.init_from_tree(spec, adj_t3)
    f_cache[ht_t3.state(0)] = Fr(1, 1)
    f_cache[ht_t3.step(0).state(0)] = Fr(500, 1)
    a_big_lt = dict(a_arrow)
    a_big_lt[(16, 4)] = Fr(450, 1)
    caught = False
    try:
        verify_banded_drift(ht_t3, f_fn, g_edge, a_big_lt, T, tail_A, tail_B, bug_bad_pool_arrow=True)
    except AssertionError as e:
        caught = True
        nc_tree.append(("NC-T3_unsafe_pool_arrow_b_lt_T", str(e)))
    assert caught, "NC-T3 MISSED"

    print("TREE PART 3 Negative controls (all caught):")
    for k, msg in nc_tree:
        print(f"  {k}: CAUGHT ({msg})")


def run_cli_negative_controls(exe, bands_src, rec_src):
    with open(bands_src) as f: bands_txt = f.read()
    with open(rec_src) as f: rec_lines = f.read().splitlines()

    with tempfile.TemporaryDirectory() as td:
        def run_case(name, num, den, b_txt, r_lines, expect_fail=True):
            bp = os.path.join(td, "b.txt")
            rp = os.path.join(td, "r.txt")
            with open(bp, "w") as f: f.write(b_txt)
            with open(rp, "w") as f: f.write("\n".join(r_lines) + "\n")
            res = subprocess.run([exe, str(num), str(den), bp, rp], capture_output=True, text=True)
            if expect_fail:
                assert res.returncode != 0, f"{name} UNEXPECTED PASS!\n{res.stdout}"
                msg = (res.stderr.strip() or res.stdout.strip().splitlines()[-1])
                print(f"  {name}: REJECTED (rc={res.returncode}: {msg})")
            else:
                assert res.returncode == 0, f"{name} UNEXPECTED FAIL!\n{res.stdout}\n{res.stderr}"

        print("\nCLI Negative controls on verify_banded (smoke geometry):")
        # 1. Stricter exponent 56/100 (smoke recipe is calibrated at 3/5 = 0.60, should fail at 0.56!)
        run_case("NC-CLI-1_stricter_lambda_56/100", 56, 100, bands_txt, rec_lines)
        # 2. Reduce reward 0 dramatically -> fails newborn consistency!
        bad_rew = [ln if not ln.startswith("reward 0 ") else "reward 0 1000" for ln in rec_lines]
        run_case("NC-CLI-2_understated_reward", 3, 5, bands_txt, bad_rew)
        # 3. Unobservable nonzero edge transfer (sender deg >= h_c: e.g. D=20, c=2 has h_2=8, d=10 >= 8!)
        bad_edge = rec_lines + ["edge 20 2 10 2 12345"]
        run_case("NC-CLI-3_unobservable_edge", 3, 5, bands_txt, bad_edge)
        # 4. Negative arrow entry
        bad_arr_neg = rec_lines + ["arrow 11 2 -500"]
        run_case("NC-CLI-4_negative_arrow", 3, 5, bands_txt, bad_arr_neg)
        # 5. Arrow entry at D >= H (=48)
        bad_arr_tail = rec_lines + ["arrow 48 2 500"]
        run_case("NC-CLI-5_tail_arrow_D>=H", 3, 5, bands_txt, bad_arr_tail)
        # 6. Missing reward line
        bad_miss = [ln for ln in rec_lines if not ln.startswith("reward 5 ")]
        run_case("NC-CLI-6_missing_reward", 3, 5, bands_txt, bad_miss)
        # 7. Varying h_1 across bands
        bad_b1 = bands_txt.replace("28 30 ", "28 28 ")
        run_case("NC-CLI-7_varying_h1", 3, 5, bad_b1, rec_lines)
        # 8. Witness cutoff h_3 < 6 at D < 16 (keep monotone across bands so specifically witness rule triggers!)
        bad_b2 = "16 30 10 5 5 4\n28 30 8 5 4 3\n48 30 6 5 4 3\n"
        run_case("NC-CLI-8_witness_cutoff_h3<6", 3, 5, bad_b2, rec_lines)
        # 9. Massive corrupted arrow entry (arrow 2 1 inflated 100x breaks consistency at 3/5!)
        bad_arr_huge = [ln if not ln.startswith("arrow 2 1 ") else "arrow 2 1 50000000000000" for ln in rec_lines]
        run_case("NC-CLI-9_inflated_arrow_breaks_bound", 3, 5, bands_txt, bad_arr_huge)
        # 10. Non-monotone cutoffs across bands
        bad_b3 = "16 30 10 5 5 4\n28 30 8 6 4 3\n48 30 6 6 4 3\n"
        run_case("NC-CLI-10_nonmonotone_bands", 3, 5, bad_b3, rec_lines)
        # 11. Corrupted edge transfer (make edge 2 1 2 2 huge positive -> breaks consistency at 3/5!)
        bad_e_huge = [ln if not ln.startswith("edge 2 1 2 2 ") else "edge 2 1 2 2 50000000000000" for ln in rec_lines]
        if bad_e_huge == rec_lines:
            bad_e_huge = rec_lines + ["edge 2 1 2 2 50000000000000"]
        run_case("NC-CLI-11_corrupted_edge_breaks_bound", 3, 5, bands_txt, bad_e_huge)


if __name__ == '__main__':
    import argparse
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description="Independent exact-tree and CLI negative-control audit.")
    ap.add_argument("--exe", default=os.path.join(here, "verify_banded"),
                    help="Path to compiled verify_banded executable")
    ap.add_argument("--bands", default=os.path.join(here, "fixtures", "bands-smoke.txt"),
                    help="Path to smoke bands geometry file")
    ap.add_argument("--recipe", default=os.path.join(here, "fixtures", "smoke-recipe.txt"),
                    help="Path to smoke recipe file")
    args = ap.parse_args()
    run_exact_tree_tests()
    run_cli_negative_controls(args.exe, args.bands, args.recipe)
    print("\nALL BANDED VERIFIER SELF-TESTS & NEGATIVE CONTROLS PASSED!")
