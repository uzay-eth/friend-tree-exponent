// Retrieve the production weights requested on stdin. The checker still runs all rows.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using Key = std::array<int, 3>;
static std::map<Key, std::int64_t> requested;
static int largest_degree = 0;

static void observe_row(int degree, int core_degree, int state, std::int64_t weight) {
    if (degree > largest_degree) return;
    auto pos = requested.find({degree, core_degree, state});
    if (pos != requested.end()) pos->second = weight;
}

#define FRIEND_TREE_OBSERVE_ROW observe_row
#define FRIEND_TREE_NO_MAIN
#include "../banded-upper/check.cpp"


int main(int argc, char** argv) {
    if (argc != 3) return 1;
    try {
        Checker checker(argv[1], argv[2], 5657, 9983);
        std::vector<Key> order;
        int degree, core_degree;
        while (std::cin >> degree >> core_degree) {
            require(degree >= 2 && degree < checker.H, "finite degree required");
            auto& layer = checker.band(degree).layers.at(core_degree);
            Bag bag;
            for (int i = 0; i < core_degree; ++i) {
                int d, k, e, j;
                require(bool(std::cin >> d >> k >> e >> j), "incomplete label");
                bag.push_back(layer.label(d, k, e, j));
            }
            Key key{degree, core_degree, layer.index(bag)};
            requested[key] = 0;
            order.push_back(key);
            largest_degree = std::max(largest_degree, degree);
        }
        checker.run();
        for (auto key : order) {
            auto weight = requested.at(key);
            require(weight > 0, "missing checked weight");
            std::cout << "weight " << weight << '\n';
        }
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
