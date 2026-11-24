#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

using U = uint64_t;
struct Core { U leaves; double stamp; std::vector<size_t> nbr; };

int main(int argc, char** argv) {
    if (argc != 3) return 2;
    const U reps = std::strtoull(argv[1], nullptr, 10);
    std::mt19937_64 rng(std::strtoull(argv[2], nullptr, 10));
    auto pick = [&](U n) { return std::uniform_int_distribution<U>(0, n-1)(rng); };
    auto wait = [&]() { return std::exponential_distribution<double>(1)(rng); };
    auto refresh = [&](Core& c, double t) {
        const double dt = t-c.stamp;
        c.stamp = t;
        if (dt <= 0 || c.leaves == 0) return;
        const auto init = c.leaves;
        if (init <= 4) {
            std::geometric_distribution<U> g(std::exp(-dt));
            for (U i = 0; i < init; ++i) c.leaves += g(rng);
        } else {
            const double rate = std::gamma_distribution<double>(init, std::expm1(dt))(rng);
            c.leaves += std::poisson_distribution<U>(rate)(rng);
        }
    };
    const std::array<double, 7> edges{0, .01, .02, .04, .06, .08, .10};
    std::array<std::array<U, 6>, 2> counts{};
    std::array<double, 2> sum{}, sum_n{};
    double paired2 = 0;
    std::vector<Core> c;
    for (U run = 0; run < reps; ++run) {
        c.clear();
        c.push_back({2, 0, {}});
        double next = wait(), prev = 0;
        for (int k = 0; k < 2; ++k) {
            const double t = 13+k;
            while (next <= t) {
                const size_t i = pick(c.size());
                refresh(c[i], next);
                const U target = pick(c[i].nbr.size()+c[i].leaves);
                if (target < c[i].nbr.size()) {
                    auto& v = c[c[i].nbr[target]];
                    refresh(v, next);
                    ++v.leaves;
                } else {
                    --c[i].leaves;
                    c[i].nbr.push_back(c.size());
                    c.push_back({1, next, {i}});
                }
                next += wait()/c.size();
            }
            U n = c.size();
            for (auto& v : c) { refresh(v, t); n += v.leaves; }
            const double x = c.size()/std::pow(n, .56572366);
            sum[k] += x;
            sum_n[k] += n;
            if (k == 1) paired2 += (x-prev)*(x-prev);
            prev = x;
            for (size_t b = 0; b+1 < edges.size(); ++b)
                counts[k][b] += (edges[b] < x && x <= edges[b+1]);
        }
    }
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "reps=" << reps << ", paired RMS=" << std::sqrt(paired2/reps) << '\n';
    for (int k = 0; k < 2; ++k) {
        std::cout << "t=" << 13+k << ", mean size=" << sum_n[k]/reps << ", mean Xproxy=" << sum[k]/reps << '\n';
        for (size_t b = 0; b+1 < edges.size(); ++b) {
            const double a = edges[b], z = edges[b+1], p = double(counts[k][b])/reps;
            std::cout << a << " " << z << " count=" << counts[k][b] << ", density=" << p/(z-a)
                << " +/- " << std::sqrt(p*(1-p)/reps)/(z-a)
                << ", kappa=1 " << (std::pow(z, 1/.56572366)-std::pow(a, 1/.56572366))/(z-a) << '\n';
        }
    }
}
