#include "rows.hpp"
#include <fstream>
#include <iostream>
#include <sstream>

using namespace mh;

namespace {
struct Stats {
    I max = INT64_MIN;
    Labels witness{};
    size_t rows{}, positive{}, fallback{};
    std::vector<I> by_anchor;
    I fork = 0;
};
Family parse_family(const std::string &t) {
    for (int f = V; f <= T; ++f) if (t == names[f]) return Family(f);
    throw std::runtime_error("unknown family " + t);
}

template<class Visit>
void rows(const Snapshot &s, Family f, Visit visit) {
    const auto &c = s.c;
    const auto &g = s.level(f);
    if (f == V) {
        for (int p = 0; p < g.npars; ++p) for (int u : s.hu_state)
            for (int w = 0; w < g.NG; ++w) visit(Labels{g.pars[p],u,w});
    } else if (f == C) {
        for (int u : s.cu_state) for (int a = 0; a < g.NG; ++a)
            for (int b = a; b < g.NG; ++b) visit(Labels{u,a,b});
    } else if (f == P) {
        for (int t = 0; t < c.NP4_ACT; ++t) for (int a = 0; a < g.npars; ++a)
            for (int b = 0; b < g.NG; ++b) visit(Labels{g.pars[a],c.p4_u1[t],c.p4_u2[t],b});
    } else if (f == F) {
        for (int t = 0; t < c.NF4_ACT; ++t) for (int a = 0; a < g.NG; ++a)
            for (int b = 0; b < g.NG; ++b) visit(Labels{c.f4_u1[t],a,c.f4_u2[t],b});
    } else if (f == Y) {
        for (int t = 0; t < c.NY4_ACT; ++t) for (int a = 0; a < g.npars; ++a)
            for (int b = 0; b < g.NG; ++b) for (int d = b; d < g.NG; ++d)
                visit(Labels{g.pars[a],c.y4_u[t],b,d});
    } else {
        for (int t = 0; t < c.NT4_ACT; ++t) for (int a = 0; a < g.NG; ++a)
            for (int b = a; b < g.NG; ++b) for (int d = b; d < g.NG; ++d)
                visit(Labels{c.t4_u[t],a,b,d});
    }
}
int anchor(Family f) {
    static constexpr int idx[] = {-1,2,0,3,3,1,0};
    return idx[f];
}
size_t row_count(const Snapshot &s, Family f) {
    const auto &z = s.z;
    size_t lens[] = {0,z.n_hV,z.n_hC,z.n_hP4,z.n_hF4,z.n_hY4,z.n_hT4};
    return lens[f];
}
Stats scan(const Snapshot &s, const Checker &check, Family f, bool gmp_all) {
    Stats st;
    st.by_anchor.assign(s.level(f).NG,INT64_MIN);
    rows(s,f,[&](const Labels &x) {
        Grid g;
        check.evaluate(f,x,g);
        I hi = Grid::ceil_grid(g.up), lo = Grid::ceil_grid(g.up-g.noninteger);
        if (hi != lo || gmp_all) {
            Rational exact;
            check.evaluate(f,x,exact);
            I e = exact.ceil();
            require(e >= lo && e <= hi,"integer filter disagrees with exact rational");
            hi = e;
            ++st.fallback;
        }
        if (hi > st.max) { st.max = hi; st.witness = x; }
        ++st.rows;
        if (hi > 0) ++st.positive;
        auto &a = st.by_anchor[x[anchor(f)]];
        a = std::max(a,hi);
        if (f == F && hi > 0) {
            Wide charged = Wide(s.s4[x[0]].k-1)*hi;
            require(charged <= INT64_MAX,"fork charge exceeds int64");
            st.fork = std::max(st.fork,I(charged));
        }
    });
    require(st.rows == row_count(s,f),"incomplete family");
    std::cout << "family=" << names[f] << " rows=" << st.rows << " positive=" << st.positive
              << " exact_fallback=" << st.fallback;
    if (st.rows) {
        std::cout << " max_ceil=" << st.max << " corner=";
        for (int j = 0; j < shapes[f].n; ++j) std::cout << (j ? "," : "") << st.witness[j];
    } else std::cout << " empty";
    std::cout << std::endl;
    return st;
}
I nonnegative(const Stats &st, int anchor) {
    return std::max<I>(0,st.by_anchor[anchor]);
}
Wide aggregate(const Snapshot &s, const std::array<Stats,NFamily> &a) {
    // A single deepest endpoint determines either chain. In a fork the fixed
    // grandparent of its deepest endpoint has k-1 options for the side child.
    Wide linear = std::max<I>(0,a[V].max) + Wide(std::max<I>(0,a[P].max)) + a[F].fork;
    Wide branch = 0;
    for (int i = 0; i < s.c.L3.INNER_N; ++i) {
        auto &u = s.s3[i];
        if (u.k < 2) continue;
        I c = nonnegative(a[C],i), y = 0, t = 0;
        if (u.d < s.c.L4.D0) {
            int j = c4_lvl_inner(&s.c.L4,u.l,u.k);
            y = nonnegative(a[Y],j); t = nonnegative(a[T],j);
        }
        Wide pairs = u.k*(u.k-1)/2;
        Wide charge = pairs*(Wide(c)+y) + pairs*(u.k-2)/3*t;
        branch = std::max(branch,(charge+u.k-1)/u.k);
    }
    std::cout << "chain_and_fork_units=" << decimal(linear)
              << " per_child_branch_units=" << decimal(branch) << '\n';
    return linear+branch;
}
void query(const Snapshot &s, const Checker &check) {
    std::string line;
    while (std::getline(std::cin,line)) {
        if (line.empty()) continue;
        std::istringstream in(line);
        std::string tag, extra;
        require(bool(in >> tag),"missing family");
        Family f = parse_family(tag);
        Labels x{};
        for (int i = 0; i < shapes[f].n; ++i) require(bool(in >> x[i]),"missing corner coordinate");
        require(!(in >> extra) && check.valid(f,x),"invalid requested corner");
        Rational r;
        check.evaluate(f,x,r);
        std::cout << r.sum << '\n';
    }
    require(!std::cin.bad(),"input error");
    (void)s;
}
}

int main(int argc, char **argv) try {
    if (argc < 2) throw std::runtime_error("usage: mh-check FILE [--query | --family V|C|P|F|Y|T] [--pq P Q] [--gmp] [--profile PATH] [--budget UNITS]");
    Snapshot s(argv[1]);
    bool gmp = false, bounded = false, queries = false;
    I budget = 0; Family only = Edge;
    std::string profile;
    for (int i = 2; i < argc; ++i) {
        std::string flag = argv[i];
        if (flag == "--query") queries = true;
        else if (flag == "--gmp") gmp = true;
        else if (flag == "--pq" && i+2 < argc) {
            auto positive = [](const std::string &arg) {
                size_t used;
                long n = std::stol(arg,&used);
                require(used == arg.size() && n > 0 && n <= INT_MAX,"invalid rational exponent");
                return int(n);
            };
            s.h.P = positive(argv[++i]);
            s.h.Q = positive(argv[++i]);
        }
        else if (flag == "--profile" && i+1 < argc) profile = argv[++i];
        else if (flag == "--family" && i+1 < argc) only = parse_family(argv[++i]);
        else if (flag == "--budget" && i+1 < argc) {
            std::string v = argv[++i]; size_t len;
            budget = std::stoll(v,&len); require(len == v.size() && budget >= 0,"invalid budget"); bounded = true;
        } else throw std::runtime_error("invalid argument " + flag);
    }
    Checker check(s);
    if (queries) {
        require(!bounded && !gmp && profile.empty() && only == Edge,"invalid combination of query flags");
        query(s,check);
        return 0;
    }
    std::cout << "mode=" << (s.sign == 1 ? "UPPER" : "LOWER") << " P=" << s.h.P
              << " Q=" << s.h.Q << " S=" << s.h.S << " Sr=" << s.h.S_r << std::endl;
    std::array<Stats,NFamily> st{};
    std::ofstream prof;
    if (!profile.empty()) {
        prof.open(profile); require(prof.good(),"cannot write profile");
        prof << "family\tanchor_label\tdegree\tchildren\tkind\tmax_ceil\n";
    }
    for (int fi = V; fi <= T; ++fi) {
        Family f = Family(fi);
        if (only != Edge && only != f) continue;
        st[f] = scan(s,check,f,gmp);
        if (prof.is_open()) {
            auto &v = s.states(f);
            for (int j = 0; j < s.level(f).NG; ++j) if (st[f].by_anchor[j] != INT64_MIN)
                prof << names[f] << '\t' << j << '\t' << v[j].d << '\t' << v[j].k << '\t'
                     << v[j].type << '\t' << st[f].by_anchor[j] << '\n';
        }
    }
    if (prof.is_open()) { prof.close(); require(prof.good(),"incomplete profile"); }
    if (only != Edge && bounded && st[only].max > budget) {
        std::cerr << "family " << names[only] << " exceeds budget " << budget << '\n';
        return 2;
    }
    if (only == Edge) {
        Wide r = aggregate(s,st);
        std::cout << "higher_residual_per_core=" << decimal(r) << "/" << s.h.S << std::endl;
        if (bounded && r > budget) {
            std::cerr << "certified residual exceeds budget " << budget << '\n'; return 2;
        }
    }
    return 0;
} catch (const std::exception &e) {
    std::cerr << "mh-check: " << e.what() << '\n'; return 1;
}
