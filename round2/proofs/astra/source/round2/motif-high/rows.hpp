// One generator formula on embedded connected sets. The only table-specific choices
// are the source-array bindings and the frozen layout in snapshot.hpp.
#pragma once
#include "snapshot.hpp"
#include <gmpxx.h>

namespace mh {
inline std::string decimal(Wide n) {
    if (n == 0) return "0";
    bool negative = n < 0;
    if (negative) n = -n;  // preflight bounds keep n strictly above INT128_MIN
    std::string s;
    while (n) { s.push_back(char('0' + int(n % 10))); n /= 10; }
    if (negative) s.push_back('-');
    std::reverse(s.begin(),s.end());
    return s;
}

struct Grid {
    // Inward ambiguity is resolved by exact GMP arithmetic in the scan.
    static constexpr I B = 4096;
    Wide up = 0;
    int noninteger = 0;
    void add(Wide n, I d = 1) {
        if (n == 0) return;
        n *= B;
        if (d == 1) { up += n; return; }
        Wide quot = n/d;
        bool rem = n != quot*d;
        up += quot + (n > 0 && rem);
        noninteger += rem;
    }
    static I ceil_grid(Wide x) {
        Wide q = x/B + (x > 0 && x%B != 0);
        if (q < INT64_MIN || q > INT64_MAX) throw std::runtime_error("row result exceeds int64");
        return I(q);
    }
};
struct Rational {
    mpq_class sum = 0;
    void add(Wide n, I d = 1) {
        if (!n) return;
        mpq_class t(mpz_class(decimal(n)), mpz_class(d));
        t.canonicalize();
        sum += t;
    }
    I ceil() const {
        mpz_class q;
        mpz_cdiv_q(q.get_mpz_t(),sum.get_num_mpz_t(),sum.get_den_mpz_t());
        require(q.fits_slong_p(),"exact row result exceeds int64");
        return q.get_si();
    }
};

struct Cut { int donor, contact; };
struct Extension { Family f{}; Labels from{}; };

inline Extension extend_shape(Family f, int vertex, bool parent) {
    int n = shapes[f].n;
    Labels parents = shapes[f].parent;
    if (parent) { parents[0] = n; parents[n] = -1; }
    else parents[n] = vertex;
    for (int larger = V; larger <= T; ++larger) {
        const auto &shape = shapes[larger];
        if (shape.n != n+1) continue;
        Labels map{0,1,2,3};
        do {
            bool same = true;
            for (int j = 0; j <= n; ++j) {
                int pp = shape.parent[j];
                if (parents[map[j]] != (pp < 0 ? -1 : map[pp])) same = false;
            }
            if (same) return {Family(larger), map};
        } while (std::next_permutation(map.begin(),map.begin()+n+1));
    }
    throw std::logic_error("new rooted shape missing");
}
struct Plan {
    Family f;
    int n;
    Labels children{}, degree{};
    std::array<Labels,4> neighbors{};
    Cut cuts[3]{};
    int ncut = 0;
    std::array<Extension,4> child{};
    Extension parent{};
    explicit Plan(Family family) : f(family), n(shapes[f].n) {
        auto &p = shapes[f].parent;
        for (int j = 1; j < n; ++j) {
            ++children[p[j]];
            neighbors[p[j]][degree[p[j]]++] = j;
            neighbors[j][degree[j]++] = p[j];
        }
        for (int j = 0; j < n; ++j) {
            if (degree[j] == 1) cuts[ncut++] = {j,neighbors[j][0]};
            if (n < 4) child[j] = extend_shape(f,j,false);
        }
        if (n < 4) parent = extend_shape(f,0,true);
    }
};

class Checker {
    const Snapshot &s;
    std::array<Plan,6> plans{Plan(V),Plan(C),Plan(P),Plan(F),Plan(Y),Plan(T)};
    struct Fraction { Wide n; I d; };

    static I changed(const Snapshot &s, Family f, Labels &x, int i, int to) {
        int old = x[i];
        x[i] = to;
        I answer = s.value(f,x);
        x[i] = old;
        return answer;
    }
    static Labels projected(const Snapshot &s, const Extension &ext, const Labels &x) {
        Labels out{};
        for (int j = 0; j < 4; ++j)
            out[j] = ext.from[j] == 3 ? c4_lvl_inner(&s.c.L4,1,0) : s.s3[x[ext.from[j]]].coarse.a;
        return out;
    }
    bool covered(const Extension &ext, const Labels &x) const {
        return s.active(ext.f,projected(s,ext,x)) >= 0;
    }
    Fraction birth(const Extension &ext, const Labels &post) const {
        // Delete the NEW child to recover its unique old embedding. Internal vertices
        // of an active extension all project to a single fine physical state.
        Labels x = projected(s,ext,post);
        int act = s.active(ext.f,x);
        if (act < 0) return {0,1};
        int var[2]{}, m = 0;
        Mass mass[2]{};
        for (int j = 0; j < 4; ++j) if (ext.from[j] != 3) {
            Mass a = s.s3[post[ext.from[j]]].coarse;
            if (a.nb != 0) {
                if (m == 2) throw std::logic_error("three interpolated old endpoints at birth");
                var[m] = j; mass[m++] = a;
            }
        }
        if (m == 0) return {s.value(ext.f,act,x),1};
        if (m == 1) {
            Wide w = Wide(mass[0].na)*s.value(ext.f,act,x);
            x[var[0]] = mass[0].b;
            return {w + Wide(mass[0].nb)*s.value(ext.f,act,x), mass[0].na+mass[0].nb};
        }
        Wide w = 0;
        for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) {
            x[var[0]] = i ? mass[0].b : mass[0].a;
            x[var[1]] = j ? mass[1].b : mass[1].a;
            w += Wide(i ? mass[0].nb : mass[0].na)*(j ? mass[1].nb : mass[1].na)*s.value(ext.f,act,x);
        }
        return {w, I(mass[0].na+mass[0].nb)*(mass[1].na+mass[1].nb)};
    }

    struct Envelope { I m, r; int den; };
    Envelope envelope(Family f, const Labels &x, int act, int donor) const {
        const I *const *e = nullptr;
        size_t idx = 0, n3 = s.c.L3.NG, n4 = s.c.L4.NG;
        const size_t a = act;
        if (f == V) {
            e = s.e23 + (donor == 0 ? 0 : 2);
            idx = a*n3 + x[donor == 0 ? 2 : 0];
        } else if (f == C) {
            e = s.e23;
            idx = size_t(s.hu3[x[0]])*n3 + x[donor == 1 ? 2 : 1];
        } else if (f == P || f == F) {
            bool second = donor == 3;
            e = (f == P ? s.ep : s.ef) + (second ? 2 : 0);
            idx = a*n4 + x[!second ? 3 : f == P ? 0 : 1];
        } else if (f == Y && donor != 0) {
            e = s.eyv;
            idx = (a*s.c.L4.npars+s.c.L4.p_to_pi[x[0]])*n4 + x[donor == 2 ? 3 : 2];
        } else {
            // hC is symmetric. Either ordering of the two physical arguments gives a
            // validated envelope. Their average is valid and symmetric at the corners.
            int yact = f == Y ? act : s.yi4[x[0]];
            int i = f == Y || donor == 1 ? 2 : 1;
            int j = f == Y || donor != 3 ? 3 : 2;
            size_t base = size_t(yact)*n4*n4;
            size_t ij = base + size_t(x[i])*n4 + x[j], ji = base + size_t(x[j])*n4 + x[i];
            return {s.sign*(s.eyc[0][ij]+s.eyc[0][ji]), s.eyc[1][ij]+s.eyc[1][ji], 2};
        }
        return {s.sign*e[0][idx],e[1][idx],1};
    }
    template<class Sum> void sources(Sum &out, const Plan &p, const Labels &x, int act) const {
        const auto &v = s.states(p.f);
        for (int i = 0; i < p.ncut; ++i) {
            auto [j,a] = p.cuts[i];
            bool from_parent = shapes[p.f].parent[a] == j;
            int cent = p.f == C ? 4 : p.f == V ? (from_parent ? 2 : 3) : (from_parent ? 0 : 1);
            I z = s.cent[cent][v[x[a]].center];
            const State &donor = v[x[j]];
            I d = donor.type == 1 ? 1 : donor.d;
            I num = (donor.type == 1 ? 0 : s.h.S_r) - z*d;
            Envelope e = envelope(p.f,x,act,j);
            out.add(Wide(num)*(e.m + (num >= 0 ? e.r : -e.r)),d*s.h.S_r*e.den);
        }
    }

public:
    explicit Checker(const Snapshot &snapshot) : s(snapshot) {}
    template<class Sum> void evaluate(Family f, Labels x, Sum &out) const {
        const Plan &p = plans[f-V];
        const Cert4vLevel &lev = s.level(f);
        const auto &v = s.states(f);
        int act = s.active(f,x);
        if (act < 0) throw std::runtime_error("row outside active support");
        I w = s.value(f,act,x);
        out.add(-Wide(s.h.P)*w,s.h.Q);
        sources(out,p,x,act);
        for (int i = 0; i < p.n; ++i) {
            int old = x[i];
            const State &u = v[old];
            if (u.type == 1) continue;
            I delta = changed(s,f,x,i,u.up)-w;
            if (u.type == 2) {
                // A tail is an undirected endpoint, and its sole internal neighbour is inner.
                const State &adj = v[x[p.neighbors[i][0]]];
                if (p.degree[i] != 1 || adj.type != 0) throw std::logic_error("tail at branching vertex");
                int base = lev.D0-u.k-1;
                if (delta >= 0) out.add(delta);
                else if (u.k <= lev.KT)
                    out.add(Wide(std::max(0,base*adj.d+1))*delta,(lev.D0+1)*adj.d);
                if (u.k <= lev.KT) {
                    I bdelta = changed(s,f,x,i,u.bud)-w;
                    if (bdelta >= 0) out.add(bdelta);
                    else out.add(Wide(std::max(0,base))*bdelta,lev.D0);
                }
                continue;
            }
            out.add(Wide(u.l)*delta);
            for (int j = 0; j < p.degree[i]; ++j) {
                int a = p.neighbors[i][j], dlabel = x[a];
                const State &donor = v[dlabel];
                if (donor.type == 0) out.add(delta,donor.d);
                if (donor.type == 2) {
                    // Same coordinate appears as a reciprocal and in the gain jump.
                    x[a] = lev.INF_ID;
                    I d0 = changed(s,f,x,i,u.up)-s.value(f,x);
                    x[a] = dlabel;
                    out.add(std::max(delta,d0),lev.D0);
                }
            }
            int uncentered = 0;
            if (i == 0) {
                if (p.n == 3 && covered(p.parent,x)) out.add(Wide(s.cent[0][u.center])*delta,s.h.S_r);
                else ++uncentered;
            }
            int free_children = u.k-p.children[i];
            if (free_children < 0) throw std::runtime_error("infeasible child count");
            if (free_children > 0) {
                if (p.n == 3 && covered(p.child[i],x))
                    out.add(Wide(free_children)*s.cent[1][u.center]*delta,s.h.S_r);
                else uncentered += free_children;
            }
            if (delta > 0) out.add(Wide(uncentered)*delta,2);
            if (u.l > 0) {
                x[i] = u.bud;
                I bdelta = s.value(f,x)-w;
                Fraction newborn{0,1};
                if (p.n < 4) newborn = birth(p.child[i],x);
                x[i] = old;
                out.add(Wide(u.l)*(Wide(bdelta)*newborn.d+newborn.n),I(u.d)*newborn.d);
            }
        }
    }
    bool valid(Family f, const Labels &x) const {
        auto &v = s.states(f);
        auto &p = plans[f-V];
        for (int i = 0; i < p.n; ++i) {
            if (x[i] < 0 || size_t(x[i]) >= v.size()) return false;
            const State &u = v[x[i]];
            if (u.type != 1 && u.k < p.children[i]) return false;
        }
        return s.active(f,x) >= 0;
    }
};
} // namespace mh
