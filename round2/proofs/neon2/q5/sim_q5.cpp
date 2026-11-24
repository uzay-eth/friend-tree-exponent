// g++ -std=c++17 -O3 sim_q5.cpp -o sim_q5 && ./sim_q5
// Diagnostic Monte Carlo, not an asymptotic test. AX is calibrated using a
// distinct intermediate degree cutoff; shape identities use 64-bit hashes.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <stdexcept>
#include <unordered_map>
#include <vector>

using U = uint64_t;
using Hist = std::unordered_map<U, int>;
constexpr double mu = 0.5658;

U mix(U z) {
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

struct Random {
    std::mt19937_64 gen{0x5A51957ULL};
    int pick(int n) {
        U bound = static_cast<U>(n), value;
        // The accepted interval has length a multiple of n.
        do { value = gen(); } while (value < -bound % bound);
        return static_cast<int>(value % bound);
    }
    double unit() { return (1.0 + static_cast<double>(gen() >> 11)) * 0x1.0p-53; }
};

struct Node {
    int parent;
    std::vector<int> children;
};

struct Tree {
    std::vector<Node> x;
    bool planted;
    explicit Tree(bool inert) : x{{-1, {}}}, planted(inert) {
        add(0);
        if (!inert) add(0);
    }
    void add(int p) {
        x[p].children.push_back(static_cast<int>(x.size()));
        x.push_back({p, {}});
    }
    int size() const { return static_cast<int>(x.size()); }
    int degree(int v) const { return x[v].children.size() + (v != 0 || planted); }
    void ring(Random& rng) {
        int v = rng.pick(size());
        int j = rng.pick(degree(v));
        int p = j < static_cast<int>(x[v].children.size()) ? x[v].children[j] : x[v].parent;
        if (p >= 0) add(p);
    }
};

struct Forest {
    Hist volumes, shapes;
    int count = 0, hubs = 0;
    std::vector<int> sizes;
};

U shape(int size, U child_sum) {
    return mix(0x9e3779b97f4a7c15ULL * size + child_sum);
}

Forest forest(const Tree& t, int N) {
    Forest f;
    std::vector<int> sz(t.size());
    std::vector<U> sum(t.size());
    for (int v = t.size() - 1; v >= 0; --v) {
        int d = t.degree(v), p = t.x[v].parent;
        if (d > N) { ++f.hubs; continue; }
        ++sz[v];
        U hash = shape(sz[v], sum[v]);
        if (p >= 0 && t.degree(p) <= N) {
            sz[p] += sz[v];
            sum[p] += hash;
        } else if (d >= 2) {
            ++f.count;
            ++f.volumes[sz[v]];
            ++f.shapes[hash];
            f.sizes.push_back(sz[v]);
        }
    }
    return f;
}

void sample_inert(Random& rng, int N, int samples, Hist& volumes, Hist& shapes) {
    Tree t(true);
    std::vector<int> sz;
    std::vector<U> sum;
    for (int k = 0; k < samples; ++k) {
        for (auto& v : t.x) v.children.clear();
        t.x.resize(2);
        t.x[0].children.push_back(1);
        double remaining = std::log(N) * rng.unit();
        for (;;) {
            remaining += std::log(rng.unit()) / t.size();
            if (remaining <= 0) break;
            t.ring(rng);
        }
        sz.assign(t.size(), 1);
        sum.assign(t.size(), 0);
        for (int v = t.size() - 1; v > 0; --v) {
            int p = t.x[v].parent;
            sz[p] += sz[v];
            sum[p] += shape(sz[v], sum[v]);
        }
        ++volumes[t.size()];
        ++shapes[shape(t.size(), sum[0])];
    }
}

double mass(const Hist& h, U k, int n) {
    auto p = h.find(k);
    return p == h.end() ? 0.0 : double(p->second) / n;
}

double tv(const Hist& p, int pn, const Hist& q, int qn) {
    double overlap = 0;
    for (auto [key, cnt] : p) overlap += std::min(double(cnt) / pn, mass(q, key, qn));
    return 1 - overlap;
}

double uniform_ks(std::vector<int>& sizes, int N) {
    std::sort(sizes.begin(), sizes.end());
    double result = 0;
    for (size_t i = 0; i < sizes.size(); ++i) {
        double x = std::min(1.0, std::log(sizes[i]) / std::log(N));
        result = std::max({result, x - double(i) / sizes.size(), double(i + 1) / sizes.size() - x});
    }
    return result;
}

struct Row {
    double ratio = 0, same_cutoff_ratio = 0, ks = 0, volume_tv = 0;
    double shape_tv = 0, shape_witness = 0, size_witness = 0, count = 0, predicted = 0;
    int runs = 0;
};

struct Ref {
    Hist v, s, v_holdout, s_holdout;
};

int main() {
    const auto start = std::chrono::steady_clock::now();
    constexpr int reps = 8, samples = 1800;
    const std::vector<int> sizes{1 << 14, 1 << 17, 1 << 19};
    const std::vector<double> alphas{0.25, 0.36, 0.50, 0.70, 0.90, 1.0};
    std::map<std::pair<int, double>, Row> rows;
    std::map<int, Ref> refs;
    Random physical, inert;
    inert.gen.seed(0x10F07EULL);
    for (int rep = 0; rep < reps; ++rep) {
        Tree t(false);
        t.x.reserve(sizes.back());
        for (int n : sizes) {
            while (t.size() < n) t.ring(physical);
            // Estimate AX from the cumulative degree count at k ~ n^.32.
            int k = std::max(3, static_cast<int>(std::pow(n, .32)));
            int calibration_hubs = 0;
            for (int v = 0; v < n; ++v) calibration_hubs += t.degree(v) >= k;
            double ax = calibration_hubs / std::pow(double(n) / k, mu);
            if (ax <= 0) throw std::runtime_error("calibration has no hubs");
            for (double a : alphas) {
                int N = a == 1.0 ? static_cast<int>(n / std::log(n)) : static_cast<int>(std::pow(n, a));
                auto f = forest(t, N);
                if (!f.count) continue;
                auto [it, inserted] = refs.try_emplace(N);
                Ref& ref = it->second;
                if (inserted) {
                    sample_inert(inert, N, samples, ref.v, ref.s);
                    sample_inert(inert, N, samples, ref.v_holdout, ref.s_holdout);
                }
                Row& r = rows[{n, a}];
                double L = std::log(N);
                r.ratio += f.count / (ax * std::pow(double(n) / N, mu) * L);
                r.same_cutoff_ratio += f.hubs ? f.count / (f.hubs * L) : 0;
                r.count += f.count;
                r.ks += uniform_ks(f.sizes, N);
                r.volume_tv += tv(f.volumes, f.count, ref.v, samples);
                r.shape_tv += tv(f.shapes, f.count, ref.s, samples);
                for (auto [key, cnt] : f.shapes) {
                    const double p = double(cnt) / f.count;
                    if (p > mass(ref.s, key, samples)) r.shape_witness += p - mass(ref.s_holdout, key, samples);
                }
                int min_v = static_cast<int>(std::ceil(f.count * std::log(L) / L)) + 1;
                for (auto [key, cnt] : f.volumes) {
                    if (key > static_cast<U>(min_v)) r.size_witness += double(cnt) / f.count - mass(ref.v_holdout, key, samples);
                }
                r.predicted += std::max(0.0, 1 - mu * (std::log(n) / L - 1));
                ++r.runs;
            }
        }
    }
    std::cout << "Means over " << reps << " physical trees, " << samples << "+" << samples
              << " reference draws per cutoff; 1.00 means N=floor(n/log n).\n"
              << "Plugin TV has sampling bias; holdout event gaps are lower-bound diagnostics.\n"
              << "n alpha M M/(AX*(n/N)^mu*logN) M/(hub*logN) KS volTV shapeTV shapeGap volGap volAsymLB runs\n";
    std::cout << std::fixed << std::setprecision(3);
    for (auto [key, r] : rows) {
        double d = r.runs;
        std::cout << key.first << ' ' << key.second << ' ' << r.count / d << ' ' << r.ratio / d << ' '
                  << r.same_cutoff_ratio / d << ' ' << r.ks / d << ' ' << r.volume_tv / d << ' '
                  << r.shape_tv / d << ' ' << r.shape_witness / d << ' ' << r.size_witness / d << ' '
                  << r.predicted / d << ' ' << r.runs << '\n';
    }
    std::cout << "Seconds: " << std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count() << '\n';
}
