// Direct jump-chain simulation of the isotropic-redirection tree, starting from T_3.
// Record the core size at fixed N, with one row per independently seeded realization.
#include <array>
#include <cstdint>
#include <iostream>
#include <random>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "usage: fixed_size FIRST_SEED RUNS\n";
        return 1;
    }
    const uint64_t first_seed = std::stoull(argv[1]);
    const int runs = std::stoi(argv[2]);
    if (runs < 1) return 1;
    constexpr std::array<int, 9> sizes{100, 316, 1000, 3162, 10000, 31623,
                                       100000, 316228, 1000000};
    std::cout << "seed";
    for (int size : sizes) std::cout << '\t' << size;
    std::cout << '\n';

    // A leaf stores -(its only neighbor + 1). A nonleaf stores the index of its adjacency list.
    // Every vertex can be selected directly and every edge is represented at both ends.
    std::vector<int> vertex(sizes.back());
    for (int r = 0; r < runs; ++r) {
        const uint64_t seed = first_seed + r;
        std::mt19937_64 rng(seed);
        std::vector<std::vector<int>> neighbors{{1, 2}};
        vertex[0] = 0;
        vertex[1] = vertex[2] = -1;
        int n = 3;
        std::cout << seed;
        for (int target : sizes) {
            while (n < target) {
                const int v = std::uniform_int_distribution<int>(0, n - 1)(rng);
                const int entry = vertex[v];
                int w;
                if (entry < 0) {
                    w = -entry - 1;
                } else {
                    const auto& adjacent = neighbors[entry];
                    const int j = std::uniform_int_distribution<int>(0, adjacent.size() - 1)(rng);
                    w = adjacent[j];
                }
                if (vertex[w] < 0) {
                    const int parent = -vertex[w] - 1;
                    vertex[w] = neighbors.size();
                    neighbors.push_back({parent, n});
                } else {
                    neighbors[vertex[w]].push_back(n);
                }
                vertex[n] = -w - 1;
                ++n;
            }
            std::cout << '\t' << neighbors.size();
        }
        std::cout << '\n';
    }
}
