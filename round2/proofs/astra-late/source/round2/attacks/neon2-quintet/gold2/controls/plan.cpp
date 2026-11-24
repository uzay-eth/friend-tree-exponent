// Locate two different exact G3 envelope slots in the first env23 block.
#include "quintet_format.hpp"
#include <fstream>
#include <iostream>
#include <stdexcept>

static void require(bool ok) { if (!ok) throw std::runtime_error("control precondition"); }

int main(int argc, char** argv) {
    require(argc == 2);
    std::ifstream f(argv[1], std::ios::binary);
    QuintetHeader h{};
    require(bool(f.read(reinterpret_cast<char*>(&h), sizeof h)));
    uint64_t bytes = 0;
    require(q5_compute_expected_bytes(h, bytes) == 0);
    Q5Layout L;
    L.init(h);
    require(L.total_bytes == bytes);
    auto offset = [&](size_t i) { return sizeof h + 8*i; };
    auto word = [&](size_t i) {
        int64_t v;
        f.seekg(offset(i));
        require(bool(f.read(reinterpret_cast<char*>(&v), sizeof v)));
        return v;
    };
    const int x=L.G2.idx[L.M3.D[0]][L.M3.C[0]], up=L.G2.up(x);
    for (int i=0; i<L.G3.Nin_ex; ++i) {
        for (int j=i+1; j<L.G3.Nin_ex; ++j) {
            const int64_t m0=word(L.off_e23_m+i), m1=word(L.off_e23_m+j);
            const int64_t r0=word(L.off_e23_r+i), r1=word(L.off_e23_r+j);
            if (m0==m1) continue;
            const int w=L.G2.idx[L.G3.D[i]][L.G3.C[i]];
            const int z=L.G2.idx[L.G3.D[j]][L.G3.C[j]];
            const auto g = [&](int a,int b) { return word(L.off_g+L.g_off(a,b)); };
            const int64_t d0=g(up,w)-g(x,w), d1=g(up,z)-g(x,z);
            require(std::abs(d0-m0)<=r0 && std::abs(d1-m1)<=r1);
            require(std::abs(m0)<=Q5_MAX_E && std::abs(m1)<=Q5_MAX_E);
            require(0<=r0 && r0<=Q5_MAX_E && 0<=r1 && r1<=Q5_MAX_E);
            if (std::abs(d0-m1)<=r1 || std::abs(d1-m0)<=r0) continue;
            std::cerr << "env23 block 0 receiver " << L.M3.D[0] << ',' << L.M3.C[0]
                      << " slots " << i << ',' << j << " labels " << L.G3.D[i] << ','
                      << L.G3.C[i] << " and " << L.G3.D[j] << ',' << L.G3.C[j]
                      << " physical source deltas " << d0 << ',' << d1
                      << " swapped abs errors " << std::abs(d0-m1) << ',' << std::abs(d1-m0)
                      << " vs radii " << r1 << ',' << r0 << '\n';
            for (auto v : {std::tuple{L.off_e23_m+i,m0,m1}, std::tuple{L.off_e23_m+j,m1,m0},
                           std::tuple{L.off_e23_r+i,r0,r1}, std::tuple{L.off_e23_r+j,r1,r0}}) {
                const auto [k,old,replacement]=v;
                if (old==replacement) continue;
                std::cout << offset(k) << ' ' << old << ' ' << replacement << '\n';
            }
            return 0;
        }
    }
    throw std::runtime_error("no two-way incompatible swap found");
}
