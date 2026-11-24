#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <random>
#include <vector>

// Each extant vertex rings at rate one and picks a uniform neighbour.
int main(int argc, char **argv) {
    if (argc != 4) return 2;
    const int reps = std::atoi(argv[1]);
    const double horizon = std::atof(argv[2]);
    std::mt19937_64 rng(std::strtoull(argv[3], nullptr, 10));
    std::exponential_distribution<double> exp(1.0);
    auto pick = [&](int n) { return std::uniform_int_distribution<int>(0, n-1)(rng); };
    for (int rep = 0; rep < reps; ++rep) {
        std::vector<int> parent{ -1, 0, 0 };
        std::vector<std::vector<int>> children(3);
        children[0] = {1, 2};
        double t = 0;
        while ((t += exp(rng) / parent.size()) <= horizon) {
            const int u = pick(parent.size());
            const int offset = (parent[u] >= 0);
            const int j = pick(children[u].size() + offset) - offset;
            const int v = (j < 0) ? parent[u] : children[u][j];
            children[v].push_back(parent.size());
            parent.push_back(v);
            children.emplace_back();
        }
        std::vector<int> index(parent.size(), -1);
        std::vector<int32_t> data;
        for (int u = 0; u < int(parent.size()); ++u) {
            if (u && children[u].empty()) continue;
            index[u] = data.size()/3;
            data.push_back(0);
            data.push_back(children[u].size() + (u != 0));
            data.push_back(u ? index[parent[u]] : -1);
        }
        for (int u = 1; u < int(data.size()/3); ++u) {
            ++data[3*data[3*u+2]];
            ++data[3*u];
        }
        int32_t n = data.size()/3;
        std::cout.write(reinterpret_cast<const char*>(&n), sizeof(n));
        std::cout.write(reinterpret_cast<const char*>(data.data()), data.size()*sizeof(int32_t));
    }
}
