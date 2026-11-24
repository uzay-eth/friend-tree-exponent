#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

struct Vertex {
    uint64_t leaves;
    double time;
    std::vector<size_t> core;
};

struct Simulation {
    std::mt19937_64 rng;
    std::vector<Vertex> v;
    double now;

    double uniform() { return (double(rng() >> 12) + 0.5) * 0x1p-52; }
    uint64_t index(uint64_t n) {
        // Multiply-high with rejection of the short residue classes.
        __uint128_t x = __uint128_t(rng()) * n;
        while (uint64_t(x) < n && uint64_t(x) < (uint64_t(0) - n) % n)
            x = __uint128_t(rng()) * n;
        return uint64_t(x >> 64);
    }
    double exp1() { return -std::log(uniform()); }

    void advance(size_t i, double t) {
        auto &w = v[i];
        const double dt = t - w.time;
        w.time = t;
        const auto k = w.leaves;
        if (dt <= 0 || !k) return;
        const double odds = std::expm1(dt);
        if (k < 8) {
            const double q = -std::expm1(-dt);
            for (uint64_t j = 0; j < k; ++j) {
                const double u = uniform();
                if (u < q) w.leaves += uint64_t(std::log(u) / std::log(q));
            }
        } else {
            const double intensity = std::gamma_distribution<double>(double(k), odds)(rng);
            w.leaves += std::poisson_distribution<uint64_t>(intensity)(rng);
        }
    }

    uint64_t observe(double t) {
        while (now <= t) {
            size_t i = index(v.size());
            advance(i, now);
            const auto j = index(v[i].leaves + v[i].core.size());
            if (j < v[i].core.size()) {
                const size_t target = v[i].core[j];
                advance(target, now);
                ++v[target].leaves;
            } else {
                --v[i].leaves;
                v[i].core.push_back(v.size());
                v.push_back({1, now, {i}});
            }
            now += exp1() / v.size();
        }
        uint64_t n = v.size();
        for (size_t i = 0; i < v.size(); ++i) {
            advance(i,t);
            n += v[i].leaves;
        }
        return n;
    }

    void run(double t1, double t2) {
        v.clear();
        const double first = exp1() / 2;
        if (first > t2) {
            double y = std::gamma_distribution<double>(2, std::exp(-t2))(rng);
            std::cout << "2,0,2,0," << y << '\n';
            return;
        }
        v.push_back({2, first, {}});
        now = first + exp1();
        const auto n1 = first <= t1 ? observe(t1) : 2;
        const auto c1 = first <= t1 ? v.size() : 0;
        const auto n2 = observe(t2);
        const auto c2 = v.size();
        // Future total population, conditionally on the entire observed past, is Yule(n2).
        double y = std::gamma_distribution<double>(double(n2), std::exp(-t2))(rng);
        std::cout << n1 << ',' << c1 << ',' << n2 << ',' << c2 << ',' << y << '\n';
    }
};

int main(int argc, char **argv) {
    if (argc != 5) return 2;
    Simulation sim;
    sim.rng.seed(std::strtoull(argv[1], nullptr, 10));
    int reps = std::atoi(argv[2]);
    double t1 = std::atof(argv[3]), t2 = std::atof(argv[4]);
    if (reps < 1 || t1 < 0 || t2 < t1) return 2;
    std::cout << std::setprecision(12) << "n1,c1,n2,c2,yule_limit\n";
    for (int i = 0; i < reps; ++i) sim.run(t1, t2);
}
