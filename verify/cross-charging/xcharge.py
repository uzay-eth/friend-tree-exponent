#!/usr/bin/env python3
"""Cross-charged local drift certificates for the random friend tree (statement: proofs/xcharge.md).

State of a core vertex u.
  class of a core vertex w with c_w >= 1 (as seen by its neighbours):
      cls(w) = (L_w, c_w) if L_w + c_w < Dcap and c_w <= c1, else TOP.
  A function of (L_w, c_w) only; TOP is absorbing (D and c never decrease).
  detailed state (c_u <= c1):       x_u = (L_u, c_u, S_u), S_u = multiset of the classes of u's core neighbours;
  coarse state (c1 < c_u <= cmax):  x_u = (L_u, c_u);
  beyond (c_u > cmax):              one state G  (f = G in upper mode, f = 0 in lower mode).
  f(L, .) = f(Lmax, .) for L > Lmax.
Potential Phi = sum_u f(x_u).  Multipliers (real numbers; they enter through identities valid on every tree,
so they never change the total drift):
  kappa(b), b a class:
      sum_u kappa(cls u) * (credit u gets from its TOP neighbours) = sum_{w TOP} (1/D_w) sum_{v ~ w} kappa(cls v)
  nu(a, b), a a non-TOP class with k >= 2, b a non-TOP class (nu(a, TOP) = 0):
      sum_u sum_{w ~ u, cls w = a} nu(a, cls u) * (credit w gets from its neighbours other than u)
        = sum_{w: cls w = a} sum_{z ~ w} (1/D_z) sum_{u ~ w, u != z} nu(a, cls u)

usage:
  xcharge.py MODE Dcap c1 Lmax cmax LO HI STEPS [OUT.pkl]     MODE = upper | lower   (float LP bisection)
  xcharge.py test [NTREES] [SEED]                             brute-force check of the drift identity
"""
import sys, time, itertools, pickle, random
from fractions import Fraction as Fr
import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog


def _add(d, key, coef):
    if key is None or coef == 0:
        return
    d[key] = d.get(key, 0) + coef


class Scheme:
    def __init__(self, Dcap, c1, Lmax, cmax):
        """Lmax: an integer (f(L, .) = f(Lmax, .) for L >= Lmax), or a list of block starts
        [0, 1, ..., s, b_1 < b_2 < ... < b_m]: f is constant in L on each block [start_i, start_{i+1} - 1]
        and on the last block [b_m, infinity).  Blocks of more than two elements must start at L >= Dcap + 1."""
        if isinstance(Lmax, int):
            starts = list(range(0, Lmax + 1))
        else:
            starts = list(Lmax)
        assert starts[0] == 0 and all(x < y for x, y in zip(starts, starts[1:]))
        self.starts = starts
        self.nL = len(starts)
        Lmax = starts[-1]
        assert Lmax >= Dcap >= 3 and cmax >= c1 >= 1
        for i in range(self.nL - 1):
            assert starts[i + 1] - starts[i] <= 2 or starts[i] >= Dcap + 1
        import bisect as _b
        self._bis = _b.bisect_right
        self.Dcap, self.c1, self.Lmax, self.cmax = Dcap, c1, Lmax, cmax
        self.classes = [(D - k, k) for D in range(2, Dcap) for k in range(1, min(D, c1) + 1)]
        self.TOP = len(self.classes)
        self.ncls = self.TOP + 1
        self.cidx = {a: i for i, a in enumerate(self.classes)}
        # smallest possible degree of a TOP vertex (c >= 1, D >= 2, and D >= Dcap or c > c1)
        self.top_min_deg = min(Dcap, c1 + 1)
        self.A2 = [i for i, (l, k) in enumerate(self.classes) if k >= 2]
        self.shapes = {c: [tuple(s) for s in itertools.combinations_with_replacement(range(self.ncls), c)]
                       for c in range(0, c1 + 1)}
        self.c11 = self.cidx[(1, 1)]

    def cls(self, L, c):
        return self.cidx.get((L, c), self.TOP) if (L + c < self.Dcap and c <= self.c1) else self.TOP

    def blk(self, L):
        """start of the block containing L"""
        return self.starts[self._bis(self.starts, L) - 1]

    def fkey(self, L, c, S):
        """variable key of f at (L, c, S); S is ignored for coarse states."""
        L = self.blk(L)
        if c > self.cmax:
            return 'G'
        if c > self.c1:
            return ('f', L, c)
        return ('f', L, c, S)

    @staticmethod
    def repl(S, a, b):
        S = list(S); S.remove(a); S.append(b)
        return tuple(sorted(S))

    def state_keys(self):
        for c in range(0, self.c1 + 1):
            for S in self.shapes[c]:
                for i, L in enumerate(self.starts):
                    last = self.starts[i + 1] - 1 if i + 1 < self.nL else L + 2
                    if last + c >= 2:
                        yield ('f', L, c, S)
        for c in range(self.c1 + 1, self.cmax + 1):
            for L in self.starts:
                yield ('f', L, c)

    # ------------------------------------------------------------------ one row
    def row(self, L, c, S, R, y=None):
        """Row of the vertex state with L leaves, c core neighbours and neighbour classes S (None if c > c1).
        y = None: actual L (any L >= 0).  y given: 'all L > Lmax' row evaluated at 1/D = y (affine in y).
        Returns dict(fx, main, top, nbr, coarse, pay):
          main   : exact part, dict key -> coefficient
          top    : None or (hi, e):  adversarial factor in [0, hi] times the expression e   (credit from TOP nbrs)
          nbr    : list of (a, hi, e): adversarial factor in [0, hi] times e   (outer credit of class-a neighbours)
          coarse : None or (c, e1, e2): coarse receiver, n_nt * [0,1/2] * e1 + n_top * [0,tau] * e2, n_nt + n_top = c
          pay    : coefficient p >= 0 such that the vertex pays p * (sum over its core neighbours v of kappa(cls v))
                   and this sum is not known from the state (coarse and beyond states).
        """
        TOP, Lmax, c1, cmax = self.TOP, self.Lmax, self.c1, self.cmax
        if y is None:
            Lr = L; D = L + c; LD = R(L, D); iD = R(1, D); tail = L > Lmax
        else:
            Lr = Lmax + 1; LD = 1 - c * y; iD = y; tail = True
        main, top, nbr, coarse, pay = {}, None, [], None, 0
        fx = self.fkey(Lr, c, S)
        has_up = self.blk(Lr + 1) != self.blk(Lr)        # f changes when L -> L + 1
        fup = self.fkey(Lr + 1, c, S)
        if c > cmax:
            # beyond the box: f = G; conversion keeps G, newborn sees a TOP creator
            if Lr >= 1:
                _add(main, ('f', 1, 1, (TOP,)), LD)
            e2 = {}; _add(e2, ('k', TOP), -1)
            coarse = (c, {}, e2)
            pay = iD
            return dict(fx=fx, main=main, top=top, nbr=nbr, coarse=coarse, pay=pay)
        if c > c1:
            if has_up:
                _add(main, fup, R(Lr, 1)); _add(main, fx, -R(Lr, 1))
            e1 = {}
            if has_up:
                _add(e1, fup, 1); _add(e1, fx, -1)
            e2 = dict(e1); _add(e2, ('k', TOP), -1)
            coarse = (c, e1, e2)
            if Lr >= 1:
                _add(main, self.fkey(Lr - 1, c + 1, None), LD); _add(main, fx, -LD)
                _add(main, ('f', 1, 1, (TOP,)), LD)
            pay = iD
            return dict(fx=fx, main=main, top=top, nbr=nbr, coarse=coarse, pay=pay)
        # ---------------- detailed state
        bx = self.cls(Lr, c) if c >= 1 else None
        ntop = S.count(TOP)
        # own leaf gains: exact rate = L + sum over non-TOP neighbours of 1/D_a
        if has_up:
            exact = R(Lr, 1)
            for a in S:
                if a != TOP:
                    l, k = self.classes[a]
                    exact += R(1, l + k)
            _add(main, fup, exact); _add(main, fx, -exact)
        # x (class bx in A2) receives credit from neighbour a: its other neighbours b see x gain a leaf
        if bx is not None and bx != TOP and bx in self.A2:
            for a in set(S):
                if a == TOP:
                    continue
                l, k = self.classes[a]
                rest = list(S); rest.remove(a)
                for b in rest:
                    if b != TOP:
                        _add(main, ('nu', bx, b), S.count(a) * R(1, l + k))
        # credit from TOP neighbours: total in [0, ntop * tau]
        if ntop:
            e = {}
            if has_up:
                _add(e, fup, 1); _add(e, fx, -1)
            _add(e, ('k', bx), -1)
            if bx != TOP and bx in self.A2:
                rest = list(S); rest.remove(TOP)
                for b in rest:
                    if b != TOP:
                        _add(e, ('nu', bx, b), 1)
            top = (ntop * R(1, self.top_min_deg), e)
        # own conversion
        if Lr >= 1:
            Sc = tuple(sorted(S + (self.c11,)))
            _add(main, self.fkey(Lr - 1, c + 1, Sc), LD); _add(main, fx, -LD)
            _add(main, ('f', 1, 1, (self.cls(Lr - 1, c + 1),)), LD)
        # neighbour events
        for a in sorted(set(S)):
            if a == TOP:
                continue
            m = S.count(a)
            l, k = self.classes[a]
            fg = self.fkey(Lr, c, self.repl(S, a, self.cls(l + 1, k)))
            rate = m * (l + iD)
            _add(main, fg, rate); _add(main, fx, -rate)
            if k >= 2:
                e = {}
                _add(e, fg, 1); _add(e, fx, -1)
                if bx != TOP:
                    _add(e, ('nu', a, bx), -1)
                nbr.append((a, R(m * (k - 1), 2), e))
            if l >= 1:
                rate = R(m * l, l + k)
                _add(main, self.fkey(Lr, c, self.repl(S, a, self.cls(l - 1, k + 1))), rate); _add(main, fx, -rate)
        # payments of a TOP giver (exact: it knows the classes of its neighbours)
        if bx == TOP:
            for b in S:
                _add(main, ('k', b), iD)
        return dict(fx=fx, main=main, top=top, nbr=nbr, coarse=coarse, pay=pay)

    def row_Ls(self):
        """values of L whose rows have to be checked: in a block [a, b] the rows of a, a+1, b-1, b (the rows of
        the interior a < L < b are affine in 1/D: f does not change, the vertex is TOP, the newborn sees a TOP
        creator, the conversion stays in the block); the last block [Lmax, infinity) needs Lmax, Lmax + 1 and
        the limit 1/D -> 0."""
        out = []
        for i, a in enumerate(self.starts[:-1]):
            b = self.starts[i + 1] - 1
            out.extend(sorted({a, min(a + 1, b), max(b - 1, a), b}))
        out.extend([self.Lmax, self.Lmax + 1])
        return out

    def all_rows(self, R):
        """all rows of the certificate."""
        Ls = self.row_Ls()
        for c in range(0, self.cmax + 1):
            shapes = self.shapes[c] if c <= self.c1 else [None]
            for S in shapes:
                for L in Ls:
                    if L + c >= 2:
                        yield ('box', L, c, S), self.row(L, c, S, R)
                yield ('tail0', c, S), self.row(None, c, S, R, y=R(0, 1))

    # ------------------------------------------------------------------ LP
    def build_lp(self, mode):
        upper = (mode == 'upper')
        R = lambda p, q: p / q
        keys = list(self.state_keys())
        if upper:
            keys.append('G')
        nf = len(keys)
        kap = [('k', b) for b in range(self.ncls)]
        nus = [('nu', a, b) for a in self.A2 for b in range(self.TOP)]
        allkeys = keys + kap + nus + (['Kp'] if upper else [])
        idx = {k: i for i, k in enumerate(allkeys)}
        nv = len(allkeys)
        ri, rj, rv = [], [], []
        lam_rows, lam_cols = [], []
        st = dict(nrow=0, naux=0)
        sgn = 1.0 if upper else -1.0          # every row is stored as  expr <= 0

        def put(d, scale):
            for k, v in d.items():
                if k == 'G' and not upper:
                    continue
                ri.append(st['nrow']); rj.append(idx[k]); rv.append(scale * v)

        for name, r in self.all_rows(R):
            auxs = []
            if r['top'] is not None:
                auxs.append([r['top']])
            for (a, hi, e) in r['nbr']:
                auxs.append([(hi, e)])
            if r['coarse'] is not None:
                c, e1, e2 = r['coarse']
                if upper:
                    auxs.append([(c / 2, e1), (c / self.top_min_deg, e2)])
                else:
                    auxs.append([(c / 2, e1)])        # kappa(TOP) = 0 in lower mode
            aux_ids = []
            for _ in auxs:
                aux_ids.append(nv + st['naux']); st['naux'] += 1
            put(r['main'], sgn)
            if upper and r['pay']:
                # the vertex has c core neighbours; sum of kappa over them <= c * Kp
                c = name[2] if name[0] == 'box' else name[1]
                ri.append(st['nrow']); rj.append(idx['Kp']); rv.append(sgn * r['pay'] * c)
            for a in aux_ids:
                ri.append(st['nrow']); rj.append(a); rv.append(sgn)
            fx = r['fx']
            if not (fx == 'G' and not upper):
                lam_rows.append(st['nrow']); lam_cols.append(idx[fx])
            st['nrow'] += 1
            for a, items in zip(aux_ids, auxs):
                for hi, e in items:
                    if not e:
                        continue
                    put(e, sgn * hi); ri.append(st['nrow']); rj.append(a); rv.append(-sgn); st['nrow'] += 1
        if upper:
            # beyond the box, all (L, c) with c > cmax:  (L/D) f(nb_TOP) + (c/D) Kp <= lambda G
            for k in [('f', 1, 1, (self.TOP,)), 'Kp']:
                ri.append(st['nrow']); rj.append(idx[k]); rv.append(1.0)
                lam_rows.append(st['nrow']); lam_cols.append(idx['G']); st['nrow'] += 1
            for b in range(self.ncls):
                ri.append(st['nrow']); rj.append(idx[('k', b)]); rv.append(1.0)
                ri.append(st['nrow']); rj.append(idx['Kp']); rv.append(-1.0); st['nrow'] += 1
        ntot = nv + st['naux']
        A0 = sp.csr_matrix((rv, (ri, rj)), shape=(st['nrow'], ntot))
        Alam = sp.csr_matrix((np.ones(len(lam_rows)), (lam_rows, lam_cols)), shape=(st['nrow'], ntot))
        self.lp = dict(A0=A0, Alam=Alam, idx=idx, nv=nv, nf=nf, naux=st['naux'], keys=allkeys, mode=mode)
        return self.lp

    def solve(self, x, Fcap=1e3, Mcap=1e3, verbose=False):
        """feasibility at lambda = x (upper) or delta = x (lower); returns the solution vector or None."""
        lp = self.lp; upper = (lp['mode'] == 'upper')
        idx, nv, nf, naux = lp['idx'], lp['nv'], lp['nf'], lp['naux']
        A = (lp['A0'] - x * lp['Alam']) if upper else (lp['A0'] + x * lp['Alam'])
        n = nv + naux
        lo = np.zeros(n); hi = np.full(n, np.inf)
        cost = np.zeros(n)
        if upper:
            lo[:nf] = 1.0; hi[:nf] = Fcap
            for b in range(self.ncls):
                lo[idx[('k', b)]] = -Mcap; hi[idx[('k', b)]] = Mcap
            lo[idx[('k', self.TOP)]] = 0.0
            lo[idx['Kp']] = 0.0; hi[idx['Kp']] = Mcap
            cost[:nf] = 1.0 / nf
        else:
            lo[:nf] = 0.0; hi[:nf] = Fcap
            for b in range(self.ncls):
                lo[idx[('k', b)]] = 0.0; hi[idx[('k', b)]] = Mcap
            hi[idx[('k', self.TOP)]] = 0.0
            lo[nv:] = -np.inf; hi[nv:] = 0.0                   # aux <= 0 (min with 0)
            cost[:nf] = -1.0 / nf
        for k in lp['keys'][nf:nv]:
            if k[0] == 'nu':
                lo[idx[k]] = -Mcap; hi[idx[k]] = Mcap
        Aeq = beq = None
        if not upper:
            Aeq = sp.csr_matrix(([1.0], ([0], [idx[('f', 1, 1, (self.TOP,))]])), shape=(1, n)); beq = [1.0]
        res = linprog(cost, A_ub=A, b_ub=np.zeros(A.shape[0]), A_eq=Aeq, b_eq=beq,
                      bounds=np.column_stack([lo, hi]), method='highs',
                      options=dict(presolve=True, disp=verbose))
        return res.x if res.status == 0 else None


def bisect(S, mode, lo, hi, steps, verbose=True):
    """upper: lo infeasible, hi feasible -> smallest feasible lambda; lower: lo feasible, hi infeasible."""
    best = None
    for _ in range(steps):
        mid = (lo + hi) / 2
        t0 = time.time()
        x = S.solve(mid)
        ok = x is not None
        if verbose:
            print(f"  {mode} {mid:.6f} {'feasible' if ok else 'infeasible'}  ({time.time()-t0:.1f}s)", flush=True)
        if (mode == 'upper') == ok:
            hi = mid
        else:
            lo = mid
        if ok:
            best = (mid, x)
    return lo, hi, best


# ---------------------------------------------------------------------- brute-force test of the drift identity
def random_tree(n, rng):
    adj = [[1, 2], [0], [0]]
    while len(adj) < n:
        v = rng.randrange(len(adj)); w = rng.choice(adj[v])
        adj.append([w]); adj[w].append(len(adj) - 1)
    return adj


def tree_states(S, adj):
    """state key and row arguments of every core vertex."""
    deg = [len(a) for a in adj]
    core = [u for u in range(len(adj)) if deg[u] >= 2]
    Lc = {}
    for u in core:
        c = sum(1 for w in adj[u] if deg[w] >= 2)
        Lc[u] = (deg[u] - c, c)
    out = {}
    for u in core:
        L, c = Lc[u]
        Sm = tuple(sorted(S.cls(*Lc[w]) for w in adj[u] if deg[w] >= 2)) if c <= S.c1 else None
        out[u] = (L, c, Sm)
    return deg, core, Lc, out


def test_identity(ntrees=40, seed=1, verbose=True):
    """For random trees and random values of f, kappa, nu:  n E[Phi(T') - Phi(T) | T]  ==  sum over core
    vertices of (row evaluated with the TRUE values of the uncertain rates and payments)."""
    rng = random.Random(seed)
    R = lambda p, q: Fr(p, q)
    worst = 0
    for it in range(ntrees):
        Dcap = rng.choice([3, 4, 5]); c1 = rng.choice([1, 2, 3]); Lmax = rng.choice([Dcap, Dcap + 1, 7])
        if rng.random() < 0.5:
            # block structure: singletons 0..Dcap+1, then blocks
            st = list(range(0, Dcap + 2)); x = Dcap + 1
            for _ in range(rng.randrange(1, 4)):
                x += rng.randrange(1, 5); st.append(x)
            Lmax = st
        cmax = c1 + rng.choice([0, 1, 3])
        S = Scheme(Dcap, c1, Lmax, cmax)
        mode = rng.choice(['upper', 'lower'])
        val = {}
        def V(key):
            if key == 'G':
                return val.setdefault('G', Fr(rng.randrange(1, 50), 7) if mode == 'upper' else Fr(0))
            if key == ('k', S.TOP) and mode == 'lower':
                return Fr(0)
            if key not in val:
                val[key] = Fr(rng.randrange(-30, 60), rng.randrange(1, 9))
            return val[key]
        ev = lambda d: sum(co * V(k) for k, co in d.items())
        n = rng.randrange(3, 90)
        adj = random_tree(n, rng)
        deg, core, Lc, st = tree_states(S, adj)
        phi = lambda stt: sum(V(S.fkey(L, c, Sm)) for (L, c, Sm) in stt.values())
        base = phi(st)
        # exact drift by enumeration of all moves (V, W)
        drift = Fr(0)
        for w in range(n):
            pw = sum(Fr(1, deg[v]) for v in adj[w])          # n * P(W = w)
            adj2 = [list(a) for a in adj]; adj2[w].append(n); adj2.append([w])
            drift += pw * (phi(tree_states(S, adj2)[3]) - base)
        # sum of rows with the true rates
        tot = Fr(0)
        for u in core:
            L, c, Sm = st[u]
            r = S.row(L, c, Sm, R)
            t = ev(r['main'])
            cn = [w for w in adj[u] if deg[w] >= 2]
            rho_top = sum(Fr(1, deg[w]) for w in cn if S.cls(*Lc[w]) == S.TOP)
            rho_nt = sum(Fr(1, deg[w]) for w in cn if S.cls(*Lc[w]) != S.TOP)
            if r['top'] is not None:
                assert rho_top <= r['top'][0]
                t += rho_top * ev(r['top'][1])
            for (a, hi, e) in r['nbr']:
                x = Fr(0)
                for w in cn:
                    if S.cls(*Lc[w]) == a:
                        x += sum(Fr(1, deg[z]) for z in adj[w] if deg[z] >= 2 and z != u)
                assert 0 <= x <= hi
                t += x * ev(e)
            if r['coarse'] is not None:
                cc, e1, e2 = r['coarse']
                t += rho_nt * ev(e1) + rho_top * ev(e2)
            if r['pay']:
                t += r['pay'] * sum(V(('k', S.cls(*Lc[w]))) for w in cn)
            tot += t
        err = abs(tot - drift)
        worst = max(worst, err)
        if verbose:
            print(f"  tree {it}: n={n} Dcap={Dcap} c1={c1} Lmax={Lmax} cmax={cmax} {mode}: drift={float(drift):.4f} "
                  f"rows={float(tot):.4f} diff={err}")
        assert err == 0, "DRIFT IDENTITY FAILS"
    print(f"drift identity holds exactly on {ntrees} random trees (random f, kappa, nu)")


def parse_L(arg):
    """'30' -> 30 (singletons 0..30, constant beyond);  '12:14,17,21,28' -> singletons 0..12, then blocks
    starting at 14, 17, 21, 28 (so [13] is the block 13..13, [14,16], [17,20], [21,27], [28, infinity))."""
    if ':' not in arg:
        return int(arg)
    s0, rest = arg.split(':')
    return list(range(0, int(s0) + 1)) + [int(x) for x in rest.split(',')]


if __name__ == '__main__':
    if sys.argv[1] == 'test':
        test_identity(int(sys.argv[2]) if len(sys.argv) > 2 else 40, int(sys.argv[3]) if len(sys.argv) > 3 else 1)
        sys.exit(0)
    mode = sys.argv[1]
    Dcap, c1, cmax = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[5])
    Lmax = parse_L(sys.argv[4])
    lo, hi, steps = float(sys.argv[6]), float(sys.argv[7]), int(sys.argv[8])
    t0 = time.time()
    S = Scheme(Dcap, c1, Lmax, cmax)
    lp = S.build_lp(mode)
    print(f"{mode}: Dcap={Dcap} c1={c1} Lmax={Lmax} cmax={cmax}: {S.ncls} classes, {lp['nf']} f-variables, "
          f"{lp['nv']-lp['nf']} multipliers, {lp['naux']} aux, {lp['A0'].shape[0]} rows; built in {time.time()-t0:.1f}s",
          flush=True)
    if steps == 0:
        # no bisection: solve once at HI (upper) / LO (lower) and dump that solution
        v = hi if mode == 'upper' else lo
        x = S.solve(v)
        best = (v, x) if x is not None else None
        print(f"{mode} at {v}: {'feasible' if x is not None else 'INFEASIBLE'}  ({time.time()-t0:.1f}s)", flush=True)
    else:
        lo, hi, best = bisect(S, mode, lo, hi, steps)
        print(f"{mode} threshold in ({lo:.6f}, {hi:.6f})")
    if len(sys.argv) > 9 and best is not None:
        val, x = best
        sol = {k: float(x[i]) for k, i in lp['idx'].items()}
        pickle.dump(dict(mode=mode, params=(Dcap, c1, Lmax, cmax), value=val, sol=sol), open(sys.argv[9], 'wb'))
        print("solution at", val, "written to", sys.argv[9])
