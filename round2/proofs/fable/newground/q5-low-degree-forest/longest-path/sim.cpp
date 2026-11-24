#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>

struct Fenwick {
    std::vector<int> t{0};
    int prefix(int i) const {
        int s = 0;
        for (; i; i -= i & -i) s += t[i];
        return s;
    }
    void append(int value) {
        int i = t.size();
        t.push_back(value + prefix(i - 1) - prefix(i - (i & -i)));
    }
    void add(int core, int delta) {
        for (int i = core + 1; i < int(t.size()); i += i & -i) t[i] += delta;
    }
    int owner(int rank) const {
        int i = 0, step = 1;
        while (step < int(t.size())) step <<= 1;
        for (; step; step >>= 1) {
            int j = i + step;
            if (j < int(t.size()) && t[j] <= rank) {
                rank -= t[j];
                i = j;
            }
        }
        return i;
    }
};

struct Core {
    int parent, leaves;
    std::vector<int> children;
    int degree() const { return leaves + int(children.size()) + (parent >= 0); }
    int core_neighbor(int j) const {
        if (parent >= 0) {
            if (j == 0) return parent;
            --j;
        }
        return j < int(children.size()) ? children[j] : -1;
    }
};

struct Tree {
    int n = 3;
    std::vector<Core> cores{{-1, 2, {}}};
    Fenwick leaf_weights;
    Tree() { leaf_weights.append(2); }
    void step(int donor_rank, int neighbor_rank) {
        int c = cores.size(), to;
        if (donor_rank >= c) {
            to = leaf_weights.owner(donor_rank - c);
        } else {
            to = cores[donor_rank].core_neighbor(neighbor_rank);
            if (to < 0) {
                --cores[donor_rank].leaves;
                leaf_weights.add(donor_rank, -1);
                cores[donor_rank].children.push_back(c);
                cores.push_back(Core{donor_rank, 1, {}});
                leaf_weights.append(1);
                ++n;
                return;
            }
        }
        assert(to >= 0 && to < c);
        ++cores[to].leaves;
        leaf_weights.add(to, 1);
        ++n;
    }
    // Columns: degree-2 path (number of vertices), and cutoff volumes 2, 3, 4.
    std::array<int, 4> stats() const {
        std::vector<std::array<int, 4>> dp(cores.size());
        std::array<int, 4> out{0, 1, 1, 1};
        for (int i = int(cores.size()) - 1; i >= 0; --i) {
            const Core &v = cores[i];
            int d = v.degree();
            if (d == 2) {
                int m1 = 0, m2 = 0;
                for (int j : v.children) {
                    int x = dp[j][0];
                    if (x >= m1) { m2 = m1; m1 = x; }
                    else m2 = std::max(m2, x);
                }
                dp[i][0] = 1 + m1;
                out[0] = std::max(out[0], 1 + m1 + m2);
            }
            for (int q = 1; q < 4; ++q) if (d <= q + 1) {
                int s = 1 + v.leaves;
                for (int j : v.children) s += dp[j][q];
                dp[i][q] = s;
                out[q] = std::max(out[q], s);
            }
        }
        return out;
    }
};

// Explicit labelled-tree oracle: enumerate each possible donor and compare
// every recipient rank, then check degree cuts by ordinary-vertex DFS.
static void validate() {
    std::mt19937_64 gen(873625353);
    int checked = 0;
    for (int sample = 0; sample < 100; ++sample) {
        Tree t;
        std::vector<std::vector<int>> adj{{1,2},{0},{0}};
        std::vector<int> core_id{0};
        for (int n = 3; n <= 130; ++n) {
            const int c = core_id.size();
            std::vector<int> core_rank(n, -1);
            for (int i = 0; i < c; ++i) core_rank[core_id[i]] = i;
            std::vector<int> leaves;
            std::vector<std::vector<int>> neighbor(c);
            for (int i = 0; i < c; ++i) {
                const int v = core_id[i];
                if (t.cores[i].parent >= 0) neighbor[i].push_back(core_id[t.cores[i].parent]);
                for (int x : adj[v]) if (adj[x].size() > 1 &&
                    (neighbor[i].empty() || t.cores[i].parent < 0 ||
                    x != core_id[t.cores[i].parent])) neighbor[i].push_back(x);
                std::sort(neighbor[i].begin() + (t.cores[i].parent >= 0), neighbor[i].end(),
                          [&](int a, int b) { return core_rank[a] < core_rank[b]; });
                int l = 0;
                for (int x : adj[v]) if (adj[x].size() == 1) {
                    neighbor[i].push_back(x);
                    leaves.push_back(x);
                    assert(t.leaf_weights.owner(int(leaves.size()) - 1) == i);
                    ++l;
                }
                assert(l == t.cores[i].leaves);
                assert(int(adj[v].size()) == t.cores[i].degree());
                for (int r = 0; r < int(adj[v].size()); ++r) {
                    int got = t.cores[i].core_neighbor(r);
                    assert((got < 0) == (adj[neighbor[i][r]].size() == 1));
                    if (got >= 0) assert(core_id[got] == neighbor[i][r]);
                }
            }
            assert(t.leaf_weights.prefix(c) == n - c && int(leaves.size()) == n - c);
            std::array<int,4> expected{0,1,1,1};
            for (int q = 0; q < 4; ++q) {
                auto allowed = [&](int v) { return q ? adj[v].size() <= unsigned(q + 1)
                                                      : adj[v].size() == 2; };
                std::vector<bool> seen(n, false);
                for (int v = 0; v < n; ++v) if (!seen[v] && allowed(v)) {
                    std::vector<int> todo{v};
                    seen[v] = true;
                    for (size_t j = 0; j < todo.size(); ++j) {
                        for (int x : adj[todo[j]]) if (!seen[x] && allowed(x)) {
                            seen[x] = true;
                            todo.push_back(x);
                        }
                    }
                    expected[q] = std::max(expected[q], int(todo.size()));
                }
            }
            assert(t.stats() == expected);
            ++checked;
            int r = std::uniform_int_distribution<int>(0, n - 1)(gen);
            int k = r < c ? std::uniform_int_distribution<int>(0, int(neighbor[r].size()) - 1)(gen) : 0;
            int target = r < c ? neighbor[r][k] : adj[leaves[r - c]][0];
            if (adj[target].size() == 1) core_id.push_back(target);
            adj[target].push_back(n);
            adj.push_back({target});
            t.step(r, k);
        }
    }
    std::cout << "PASS: " << checked << " explicit labelled-tree states; donor/neighbor ranks, "
              << "leaf multiplicities and all four component scans agree\n";
}

int main(int argc, char **argv) {
    if (argc == 2 && std::string(argv[1]) == "check") { validate(); return 0; }
    if (argc != 4) return 2;
    int max_exp = std::stoi(argv[1]), rep = std::stoi(argv[2]);
    uint64_t seed = std::stoull(argv[3]);
    std::mt19937_64 gen(seed);
    std::cout << "run,n,cores,L,M2,M3,M4\n";
    for (int i = 0; i < rep; ++i) {
        Tree t;
        for (int obs = 100000, exponent = 5; exponent <= max_exp; ++exponent, obs *= 10) {
            while (t.n < obs) {
                int r = std::uniform_int_distribution<int>(0, t.n - 1)(gen);
                int k = r < int(t.cores.size()) ? std::uniform_int_distribution<int>(0, t.cores[r].degree() - 1)(gen) : 0;
                t.step(r, k);
            }
            auto s = t.stats();
            std::cout << i << ',' << t.n << ',' << t.cores.size();
            for (int x : s) std::cout << ',' << x;
            std::cout << '\n';
        }
    }
}
