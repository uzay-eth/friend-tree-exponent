// ftsim.hpp -- continuous-time random friend tree, core-only representation.
//
// Model (continuous-time embedding of the random friend tree / isotropic redirection, r = 1):
// every vertex rings at rate 1 and gives a new leaf to a uniformly chosen neighbour.
// Only core vertices (degree >= 2) are stored.  For core vertex u: L[u] leaf neighbours, c[u] core
// neighbours, D = L + c.  The leaves of u each ring at rate 1 and always hit u, so between touches of u
// its leaf count is a pure-birth (Yule) process; it is advanced lazily and exactly in distribution
// (exact samplers throughout: geometric sums, Marsaglia-Tsang gamma, PTRS Poisson, large means split into parts):
//     L(t) = L(s) + NegBin(L(s), exp(-(t-s))).
// Core vertices ring at total rate C.  The ringing vertex u picks a uniform neighbour:
//     a leaf (prob L/D): that leaf becomes a core vertex with (L,c) = (1,1); u -> (L-1, c+1);
//     a core neighbour w (prob 1/D each): w gains a leaf.
// Start: the path on three vertices (one core vertex with L = 2, c = 0) at t = 0.
#pragma once
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>

struct Rng {
    uint64_t s[4];
    static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
    static uint64_t splitmix64(uint64_t &x) {
        uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
    }
    void seed(uint64_t x) { for (int i = 0; i < 4; i++) s[i] = splitmix64(x); }
    // independent stream for (stream id a, run index b): distinct (a,b) give distinct 64-bit pre-seeds unless a
    // 2^-64 hash collision occurs (the earlier a*1000003+b scheme made streams of neighbouring ids overlap).
    void seed2(uint64_t a, uint64_t b) {
        uint64_t x = a; uint64_t h = splitmix64(x);
        uint64_t y = h ^ (b * 0x9E3779B97F4A7C15ULL + 0x632BE59BD9B4E019ULL);
        for (int i = 0; i < 4; i++) s[i] = splitmix64(y);
    }
    inline uint64_t next() {  // xoshiro256**
        const uint64_t result = rotl(s[1] * 5, 7) * 9;
        const uint64_t t = s[1] << 17;
        s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t; s[3] = rotl(s[3], 45);
        return result;
    }
    inline double unif() { return ((next() >> 11) + 0.5) * 0x1.0p-53; }  // (0,1)
    inline uint64_t below(uint64_t n) {  // unbiased in [0,n), Lemire
        __uint128_t m = (__uint128_t)next() * n;
        uint64_t l = (uint64_t)m;
        if (l < n) {
            uint64_t t = (-n) % n;
            while (l < t) { m = (__uint128_t)next() * n; l = (uint64_t)m; }
        }
        return (uint64_t)(m >> 64);
    }
    inline double normal() {
        double u1 = unif(), u2 = unif();
        return sqrt(-2.0 * log(u1)) * cos(6.283185307179586476925 * u2);
    }
    inline double expo() { return -log(unif()); }
    double gamma_ge1(double a) {  // Marsaglia-Tsang, exact rejection sampler, shape a >= 1
        double d = a - 1.0 / 3.0, c = 1.0 / sqrt(9.0 * d);
        for (;;) {
            double x, v;
            do { x = normal(); v = 1.0 + c * x; } while (v <= 0.0);
            v = v * v * v;
            double u = unif();
            if (u < 1.0 - 0.0331 * (x * x) * (x * x)) return d * v;
            if (log(u) < 0.5 * x * x + d * (1.0 - v + log(v))) return d * v;
        }
    }
    uint64_t poisson(double lam) {
        if (lam < 10.0) {
            double el = exp(-lam), p = 1.0; uint64_t k = 0;
            for (;;) { p *= unif(); if (p <= el) return k; k++; }
        }
        if (lam > 1e8) {  // exact: a Poisson(lam) variable is a sum of m independent Poisson(lam/m) variables
            uint64_t m = (uint64_t)ceil(lam / 1e8), tot = 0; double part = lam / (double)m;
            for (uint64_t i = 0; i < m; i++) tot += poisson(part);
            return tot;
        }
        // PTRS (Hormann 1993), exact
        double slam = sqrt(lam), loglam = log(lam);
        double b = 0.931 + 2.53 * slam, a = -0.059 + 0.02483 * b;
        double invalpha = 1.1239 + 1.1328 / (b - 3.4), vr = 0.9277 - 3.6224 / (b - 2.0);
        for (;;) {
            double U = unif() - 0.5, V = unif();
            double us = 0.5 - fabs(U);
            double kf = floor((2.0 * a / us + b) * U + lam + 0.43);
            if (us >= 0.07 && V <= vr) return (uint64_t)kf;
            if (kf < 0 || (us < 0.013 && V > us)) continue;
            if (log(V) + log(invalpha) - log(a / (us * us) + b) <= -lam + kf * loglam - lgamma(kf + 1.0))
                return (uint64_t)kf;
        }
    }
    // Yule process: population after time tau starting from l >= 1.
    inline uint64_t yule(uint64_t l, double tau) {
        if (tau <= 0.0) return l;
        double q = -expm1(-tau);
        if (q <= 0.0) return l;
        if (l <= 12) {
            double lq = log(q);
            uint64_t add = 0;
            for (uint64_t i = 0; i < l; i++) add += (uint64_t)floor(log(unif()) / lq);
            return l + add;
        }
        double p = exp(-tau);
        double lam = gamma_ge1((double)l) * (q / p);
        return l + poisson(lam);
    }
};

static const uint32_t NONE = 0xffffffffu;

struct Tree {
    std::vector<uint64_t> L;
    std::vector<double> tl;
    std::vector<uint32_t> c, parent, first_child, next_sib;
    double t = 0.0;
    void init() {
        L.assign(1, 2); tl.assign(1, 0.0); c.assign(1, 0); parent.assign(1, NONE);
        first_child.assign(1, NONE); next_sib.assign(1, NONE); t = 0.0;
    }
    inline size_t size() const { return L.size(); }
    // advance to time tend (stops exactly at tend; all L brought up to date at tend)
    // model: 0 = isotropic redirection (the friend tree); 1 = PAN, 2 = DAN of Hartle-Krapivsky-Redner-Zhang
    // (a ringing core vertex with L > 0 always converts one of its leaves; with L = 0 it does nothing (PAN) or
    //  receives a leaf itself (DAN)).  PAN/DAN are used only to validate the estimator against exact exponents.
    // 3 = "IR-0": isotropic redirection in which a ring that picks a core neighbour is discarded (no core-to-core
    //  sends); used for the heuristic-closure comparison only.
    int model = 0;
    void run_to(double tend, Rng &g) {
        for (;;) {
            uint64_t C = L.size();
            double tev = t + g.expo() / (double)C;
            if (tev >= tend) break;
            t = tev;
            uint32_t u = (uint32_t)g.below(C);
            if (L[u]) L[u] = g.yule(L[u], t - tl[u]);
            tl[u] = t;
            uint64_t D = L[u] + c[u];
            uint64_t j = (model == 0 || model == 3) ? g.below(D) : (L[u] > 0 ? 0 : D);
            if ((model == 1 || model == 2) && L[u] == 0) { if (model == 2) L[u] = 1; continue; }
            if (model == 3 && j >= L[u]) continue;
            if (j < L[u]) {
                L[u]--; c[u]++;
                uint32_t w = (uint32_t)L.size();
                L.push_back(1); tl.push_back(t); c.push_back(1); parent.push_back(u);
                first_child.push_back(NONE); next_sib.push_back(first_child[u]);
                first_child[u] = w;
            } else {
                uint64_t k = j - L[u];
                uint32_t w;
                if (parent[u] != NONE) {
                    if (k == 0) w = parent[u];
                    else { w = first_child[u]; for (uint64_t i = 1; i < k; i++) w = next_sib[w]; }
                } else {
                    w = first_child[u]; for (uint64_t i = 0; i < k; i++) w = next_sib[w];
                }
                if (L[w]) L[w] = g.yule(L[w], t - tl[w]);
                tl[w] = t;
                L[w]++;
            }
        }
        // memorylessness: discarding the overshooting exponential is exact
        t = tend;
        uint64_t C = L.size();
        for (uint64_t u = 0; u < C; u++) {
            if (L[u]) L[u] = g.yule(L[u], tend - tl[u]);
            tl[u] = tend;
        }
    }
};
