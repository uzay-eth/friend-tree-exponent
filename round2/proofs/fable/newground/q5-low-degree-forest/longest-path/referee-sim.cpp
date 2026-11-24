#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <vector>

using Tree = std::vector<std::vector<int>>;

// Rejection sampler: exactly uniform bounded draw from the 64-bit engine.
int bounded(std::mt19937_64& rng, std::uint64_t b) {
    const auto threshold = (std::uint64_t{0} - b) % b;
    auto x = rng();
    while (x < threshold) x = rng();
    return static_cast<int>(x % b);
}

int largest(const Tree& t, int cutoff) {
    const auto keep = [&](int v) {
        return cutoff == 0 ? t[v].size() == 2 : t[v].size() <= static_cast<unsigned>(cutoff);
    };
    int best = 0;
    std::vector<unsigned char> seen(t.size());
    std::vector<int> stack;
    for (int v = 0; v < static_cast<int>(t.size()); ++v) {
        if (seen[v] || !keep(v)) continue;
        seen[v] = 1;
        stack.push_back(v);
        int size = 0;
        while (!stack.empty()) {
            const int x = stack.back();
            stack.pop_back();
            ++size;
            for (int y : t[x]) {
                if (!seen[y] && keep(y)) {
                    seen[y] = 1;
                    stack.push_back(y);
                }
            }
        }
        best = std::max(best, size);
    }
    return best;
}

struct Stats {
    std::map<int, int> hist;
    double sum{}, square{}, m2{}, m3{};
};

int main(int argc, char** argv) {
    const int reps = argc > 1 ? std::stoi(argv[1]) : 128;
    const std::uint64_t seed = argc > 2 ? std::stoull(argv[2]) : 4081123107ULL;
    constexpr std::array<int, 2> sizes{100000, 1000000};
    std::array<Stats, sizes.size()> stats;
    std::mt19937_64 rng(seed);
    for (int r = 0; r < reps; ++r) {
        Tree t{{1}, {0}};
        t.reserve(sizes.back());
        for (unsigned j = 0; j < sizes.size(); ++j) {
            for (int n = static_cast<int>(t.size()); n < sizes[j]; ++n) {
                const int donor = bounded(rng, n);
                const int recipient = t[donor][bounded(rng, t[donor].size())];
                t[recipient].push_back(n);
                t.push_back({recipient});
            }
            const int l = largest(t, 0), m2 = largest(t, 2), m3 = largest(t, 3);
            if (!(l <= m2 && m2 <= l + 2 && m2 <= m3)) return 2;
            auto& s = stats[j];
            ++s.hist[l];
            s.sum += l;
            s.square += l * l;
            s.m2 += m2;
            s.m3 += m3;
        }
    }
    std::cout << "reps=" << reps << " seed=" << seed << '\n';
    std::cout << std::fixed << std::setprecision(6);
    constexpr double mu = 0.56572366;
    for (unsigned j = 0; j < sizes.size(); ++j) {
        const int n = sizes[j];
        const double rhs = mu * std::log(n);
        double lo = 2, hi = 16;
        for (int i = 0; i < 60; ++i) {
            const double mid = (lo + hi) / 2;
            if (mid * (std::log(2 * mid) - 1) < rhs) lo = mid;
            else hi = mid;
        }
        const double u = (lo + hi) / 2;
        const auto& s = stats[j];
        const double mean = s.sum / reps;
        std::cout << n << " hist";
        for (const auto& [l, c] : s.hist) std::cout << ' ' << l << ':' << c;
        std::cout << " mean=" << mean
                  << " se=" << std::sqrt((s.square / reps - mean * mean) / (reps - 1))
                  << " u=" << u << " mean-u=" << mean - u
                  << " leading=" << rhs / std::log(std::log(n))
                  << " M2=" << s.m2 / reps << " M3=" << s.m3 / reps << '\n';
    }
}
