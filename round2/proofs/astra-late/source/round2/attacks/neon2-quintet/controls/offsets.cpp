#include <fstream>
#include <iostream>
#include "quintet_format.hpp"

int main(int argc, char **argv) {
    if (argc != 2) return 1;
    QuintetHeader h{};
    std::ifstream in(argv[1], std::ios::binary);
    if (!in.read(reinterpret_cast<char *>(&h), sizeof(h))) return 2;
    uint64_t bytes{};
    if (q5_compute_expected_bytes(h, bytes)) return 2;
    Q5Layout l;
    l.init(h);
    // First M3 row = (2,2); fine sender (2,1) projects exactly to coarse (2,1).
    const int u = l.G2.idx[2][2], v = l.G2.idx[2][1], w = l.G3.idx[2][1];
    int ce, n, d;
    l.G3.proj_int(l.G2, v, ce, n, d);
    if (l.M3.idx[2][2] != 0 || ce != w || n != 1 || d != 1) return 3;
    std::cout << bytes;
    for (size_t i : {l.off_g + l.g_off(u,v), l.off_g + l.g_off(l.G2.up(u),v),
                     l.off_e23_m + w, l.off_e23_r + w})
        std::cout << ' ' << sizeof(h) + 8*i;
    std::cout << '\n';
}
