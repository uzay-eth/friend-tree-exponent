#include <bits/stdc++.h>
using namespace std;

// Fast C++ test of Planted-Tree Quotient Bellman Enclosures:
// Tests BOTH Depth-1 AND Depth-1+Bare-Chain quotient on (s - L0 - C Delta) w = b!
// Goal: see how tight [s_lo, s_hi] gets as a function of depth & bare-chain length!

struct Child {
    int d, c;
    bool operator<(const Child& o) const { return d != o.d ? d < o.d : c < o.c; }
    bool operator==(const Child& o) const { return d == o.d && c == o.c; }
};

struct State {
    bool fine;
    int d, c;
    vector<Child> chs;
    bool operator<(const State& o) const {
        if (fine != o.fine) return fine < o.fine;
        if (d != o.d) return d < o.d;
        if (c != o.c) return c < o.c;
        return chs < o.chs;
    }
};

void gen_combs(const vector<Child>& pool, int k, int start, vector<Child>& cur, vector<vector<Child>>& out) {
    if ((int)cur.size() == k) { out.push_back(cur); return; }
    for (int i = start; i < (int)pool.size(); ++i) {
        cur.push_back(pool[i]);
        gen_combs(pool, k, i, cur, out);
        cur.pop_back();
    }
}

int main() {
    const int D0 = 24, C0 = 6;
    const int D1 = 8, C1 = 4;
    const int D_ch = 8, C_ch = 3;

    vector<Child> T_ch;
    map<Child, int> ch_id;
    for (int d = 2; d <= D_ch; ++d)
        for (int c = 1; c <= min(d, C_ch); ++c) {
            ch_id[{d, c}] = (int)T_ch.size();
            T_ch.push_back({d, c});
        }
    auto clamp_ch = [&](int d, int c) -> int {
        return ch_id[{min(d, D_ch), min(c, C_ch)}];
    };

    vector<State> states;
    map<State, int> sidx;
    for (int d = 2; d <= D1; ++d) {
        for (int c = 1; c <= min(d, C1); ++c) {
            vector<vector<Child>> combs;
            vector<Child> cur;
            gen_combs(T_ch, c - 1, 0, cur, combs);
            for (auto& chs : combs) {
                State st{true, d, c, chs};
                sidx[st] = (int)states.size();
                states.push_back(st);
            }
        }
    }
    for (int d = 2; d <= D0; ++d) {
        for (int c = 1; c <= min(d, C0); ++c) {
            if (d <= D1 && c <= C1) continue;
            State st{false, d, c, {}};
            sidx[st] = (int)states.size();
            states.push_back(st);
        }
    }
    int Nst = (int)states.size();
    int J_idx = sidx[{true, 2, 1, {}}];
    cout << "Depth-1 C++ prototype: Nst=" << Nst << " states\n";

    // Precompute transition indices
    vector<int> gain_idx(Nst), prom_idx(Nst, -1);
    vector<vector<int>> fiber(T_ch.size());
    struct ChTrans {
        int ch_cid; // child coarse id
        int dj, cj;
        int gain_to; // -1 if clamped at D_ch
        int prom_to; // -1 if clamped at C_ch or dj==cj
    };
    vector<vector<ChTrans>> ctrans(Nst);

    for (int i = 0; i < Nst; ++i) {
        const auto& st = states[i];
        fiber[clamp_ch(st.d, st.c)].push_back(i);
        // gain_idx
        if (st.fine) {
            if (st.d + 1 <= D1) gain_idx[i] = sidx[{true, st.d + 1, st.c, st.chs}];
            else gain_idx[i] = sidx[{false, min(st.d + 1, D0), st.c, {}}];
            if (st.d > st.c) {
                if (st.c + 1 <= C1) {
                    auto nchs = st.chs;
                    nchs.push_back({2, 1});
                    sort(nchs.begin(), nchs.end());
                    prom_idx[i] = sidx[{true, st.d, st.c + 1, nchs}];
                } else {
                    prom_idx[i] = sidx[{false, st.d, min(st.c + 1, C0), {}}];
                }
            }
            for (int j = 0; j < (int)st.chs.size(); ++j) {
                int dj = st.chs[j].d, cj = st.chs[j].c;
                int gto = -1, pto = -1;
                if (dj < D_ch) {
                    auto nchs = st.chs; nchs[j].d = dj + 1; sort(nchs.begin(), nchs.end());
                    gto = sidx[{true, st.d, st.c, nchs}];
                }
                if (cj < C_ch && dj > cj) {
                    auto nchs = st.chs; nchs[j].c = cj + 1; sort(nchs.begin(), nchs.end());
                    pto = sidx[{true, st.d, st.c, nchs}];
                }
                ctrans[i].push_back({ch_id[{dj, cj}], dj, cj, gto, pto});
            }
        } else {
            gain_idx[i] = sidx[{false, min(st.d + 1, D0), st.c, {}}];
            if (st.d > st.c) prom_idx[i] = sidx[{false, st.d, min(st.c + 1, C0), {}}];
        }
    }

    auto solve = [&](double s, bool upper_mu) -> double {
        vector<double> v(Nst, 1.0), vn(Nst, 1.0);
        vector<double> minD(T_ch.size()), maxD(T_ch.size());
        for (int it = 0; it < 120; ++it) {
            double gMinD = 1e9, gMaxD = -1e9;
            for (size_t k = 0; k < T_ch.size(); ++k) {
                double mn = 1e9, mx = -1e9;
                for (int idx : fiber[k]) {
                    double dv = v[gain_idx[idx]] - v[idx];
                    if (dv < mn) mn = dv;
                    if (dv > mx) mx = dv;
                }
                minD[k] = mn; maxD[k] = mx;
                if (mn < gMinD) gMinD = mn;
                if (mx > gMaxD) gMaxD = mx;
            }
            double md = 0.0;
            for (int i = 0; i < Nst; ++i) {
                const auto& st = states[i];
                int d = st.d, c = st.c, L = d - c;
                if (d == D0) {
                    // Use generous diagnostic tail at D0=24: upper_mu => 1/s, lower_mu => (1 - c/(d*(s+1)))/s?
                    // Actually for exact diagnostic of INTERIOR Depth-1 gap, use (1 - 0.55*c/d)/s vs sound endpoints:
                    vn[i] = (1.0 - (upper_mu ? 0.40 : 0.70) * c / d) / s;
                    continue;
                }
                double b = (double)L / d;
                double gv = v[gain_idx[i]];
                double pv = (L > 0) ? v[prom_idx[i]] : 0.0;
                if (!st.fine) {
                    double xi = (upper_mu == (gv >= v[i])) ? (c - 1) * 0.5 : 0.0;
                    double cd = ((double)(c - 1) / d) * (upper_mu ? gMaxD : gMinD);
                    double tot_r = L + xi + b;
                    double rhs = b + cd + (L + xi) * gv + b * pv;
                    vn[i] = rhs / (s + tot_r);
                } else {
                    double rg = L, cd = 0.0, erhs = 0.0, er = 0.0;
                    for (const auto& ct : ctrans[i]) {
                        if (ct.dj < D_ch) rg += 1.0 / ct.dj;
                        else if (upper_mu == (gv >= v[i])) rg += 1.0 / D_ch;
                        cd += (upper_mu ? maxD[ct.ch_cid] : minD[ct.ch_cid]) / d;
                        if (ct.gain_to >= 0) {
                            double nv = v[ct.gain_to];
                            int Lj = ct.dj - ct.cj;
                            double r_lo = Lj + 1.0 / d;
                            double r_hi = r_lo + ((ct.cj < C_ch) ? (ct.cj - 1) * 0.5 : (ct.dj - 1) * 0.5);
                            double rj = (upper_mu == (nv >= v[i])) ? r_hi : r_lo;
                            erhs += rj * nv; er += rj;
                        }
                        if (ct.prom_to >= 0) {
                            double nv = v[ct.prom_to];
                            double rj;
                            if (ct.dj < D_ch) rj = (double)(ct.dj - ct.cj) / ct.dj;
                            else {
                                double r_lo = max(0.0, (double)(D_ch - ct.cj) / D_ch);
                                double r_hi = 1.0;
                                rj = (upper_mu == (nv >= v[i])) ? r_hi : r_lo;
                            }
                            erhs += rj * nv; er += rj;
                        }
                    }
                    double tot_r = rg + b + er;
                    double rhs = b + cd + rg * gv + b * pv + erhs;
                    vn[i] = rhs / (s + tot_r);
                }
                double diff = fabs(vn[i] - v[i]);
                if (diff > md) md = diff;
            }
            v.swap(vn);
            if (md < 1e-12) break;
        }
        return v[J_idx];
    };

    double lo = 0.50, hi = 0.60;
    for (int k = 0; k < 35; ++k) {
        double m = 0.5 * (lo + hi);
        if (solve(m, false) > 1.0) lo = m; else hi = m;
    }
    double s_lo = lo;
    lo = 0.54; hi = 0.62;
    for (int k = 0; k < 35; ++k) {
        double m = 0.5 * (lo + hi);
        if (solve(m, true) < 1.0) hi = m; else lo = m;
    }
    double s_hi = hi;
    cout << fixed << setprecision(8) << "Depth-1 Bellman quotient: mu in [" << s_lo << ", " << s_hi << "] (width = " << (s_hi - s_lo) << ")\n";
    return 0;
}
