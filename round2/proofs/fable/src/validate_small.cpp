// validate_small.cpp -- mean C(t), S(t) at small t over many runs.  Usage: validate_small SEED RUNS T
#include "ftsim.hpp"
int main(int argc, char **argv) {
    uint64_t seed = strtoull(argv[1], 0, 10); long runs = atol(argv[2]); double T = atof(argv[3]);
    Tree tr; Rng g; double sC = 0, sC2 = 0, sS = 0, sN = 0;
    for (long r = 0; r < runs; r++) {
        g.seed2(seed, (uint64_t)r);
        tr.init(); tr.run_to(T, g);
        double C = tr.size(), S = 0, N = C;
        for (size_t u = 0; u < tr.size(); u++) { S += (double)tr.L[u] / (double)(tr.L[u] + tr.c[u]); N += tr.L[u]; }
        sC += C; sC2 += C * C; sS += S; sN += N;
    }
    double m = sC / runs, se = sqrt((sC2 / runs - m * m) / runs);
    printf("T=%g runs=%ld meanN=%.6f meanC=%.6f +- %.6f meanS=%.6f\n", T, runs, sN / runs, m, se, sS / runs);
}
