// Exact continuous-time dynamics modulo RNG and floating-point arithmetic.
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <queue>
#include <random>
#include <vector>

constexpr int K = 25, M = 9;
constexpr std::array<double, 7> thresholds{.01, .05, .1, .2, .3, .5, .7};
struct Node {
    std::vector<int> core;
    std::vector<int> marked;
    int bulk = 0, label = -1;
    double last = 0;
    int degree() const { return bulk + int(marked.size() + core.size()); }
};
struct Stat {
    double sum = 0, ss = 0;
    void add(double x) { sum += x; ss += x*x; }
    double se(int n) const { return std::sqrt(std::max(0., (ss-sum*sum/n)/(n*(n-1.)))); }
};

struct Sim {
    std::mt19937_64 rng;
    std::vector<Node> nodes;
    using Event = std::pair<double,int>;
    std::priority_queue<Event,std::vector<Event>,std::greater<Event>> events;
    std::array<int,K> where;
    explicit Sim(uint64_t seed) : rng(seed) { nodes.reserve(4096); }
    double uniform() { return (double(rng() >> 12) + .5) * 0x1p-52; }
    int pick(int n) { return std::uniform_int_distribution<int>(0,n-1)(rng); }
    double exp() { return -std::log(uniform()); }

    void grow(int i, double t) {
        auto &v = nodes[i];
        const double dt = t-v.last;
        const int k = v.bulk + int(v.marked.size());
        if (dt > 0 && k > 0) {
            // Sum of k Yule families minus k; gamma-Poisson is negative binomial.
            // For k=1, use the geometric directly.
            if (k == 1) v.bulk += int(exp()/(-std::log(-std::expm1(-dt))));
            else {
                double g = std::gamma_distribution<double>(k,std::expm1(dt))(rng);
                v.bulk += std::poisson_distribution<int>(g)(rng);
            }
        }
        v.last = t;
    }

    void run(double stop, std::array<int,K>& deg, int& total) {
        std::array<std::vector<int>,K> adj;
        adj[0].push_back(1); adj[1].push_back(0);
        double t = 0;
        nodes.clear();
        for (int n=2; n<K; ++n) {
            t += exp()/n;
            if (t >= stop) {
                total = n;
                for (int x=0; x<K; ++x) deg[x]=int(adj[x].size());
                return;
            }
            int x = pick(n), y = adj[x][pick(int(adj[x].size()))];
            adj[y].push_back(n); adj[n].push_back(y);
        }
        where.fill(-1);
        for (int x=0; x<K; ++x) if (adj[x].size()>1) {
            where[x] = int(nodes.size());
            nodes.push_back(Node{{},{},0,x,t});
        }
        for (int x=0; x<K; ++x) if (where[x]>=0) {
            int i = where[x];
            for (int y:adj[x]) {
                if (where[y]>=0) nodes[i].core.push_back(where[y]);
                else nodes[i].marked.push_back(y);
            }
            events.emplace(t+exp(),i);
        }
        while (events.top().first < stop) {
            auto [time,i] = events.top(); events.pop();
            grow(i,time);
            int j = pick(nodes[i].degree());
            if (j < int(nodes[i].core.size())) {
                int y = nodes[i].core[j];
                grow(y,time);
                ++nodes[y].bulk;
            } else {
                j -= int(nodes[i].core.size());
                int label = -1;
                auto &m = nodes[i].marked;
                if (j < int(m.size())) {
                    label = m[j]; m[j] = m.back(); m.pop_back();
                } else --nodes[i].bulk;
                int y = int(nodes.size());
                if (label >= 0) where[label] = y;
                nodes[i].core.push_back(y);
                nodes.push_back(Node{{i},{},1,label,time});
                events.emplace(time+exp(),y);
            }
            events.emplace(time+exp(),i);
        }
        while (!events.empty()) events.pop();
        deg.fill(1);
        total = int(nodes.size());
        for (int i=0; i<int(nodes.size()); ++i) {
            grow(i,stop);
            const auto &v = nodes[i];
            total += v.bulk + int(v.marked.size());
            if (v.label >= 0) deg[v.label] = v.degree();
        }
    }
};

int main(int argc, char** argv) {
    int n = argc>1 ? std::stoi(argv[1]) : 200000;
    double t = argc>2 ? std::stod(argv[2]) : 13.;
    uint64_t seed = argc>3 ? std::stoull(argv[3]) : 330283;
    Sim sim(seed);
    std::array<std::array<Stat,M>,K> stats{};
    std::array<std::array<std::array<Stat,M>,K>,K> diffs{};
    std::array<Stat,K> deg50{}, deg100{};
    Stat sizes, cores;
    for (int r=0; r<n; ++r) {
        std::array<int,K> deg;
        int size;
        sim.run(t,deg,size);
        sizes.add(size); cores.add(double(sim.nodes.size()));
        std::array<std::array<double,M>,K> x{};
        for (int i=0; i<K; ++i) {
            x[i][0] = deg[i]>=20;
            deg50[i].add(deg[i]>=50); deg100[i].add(deg[i]>=100);
            x[i][1] = double(deg[i])/size;
            for (int j=0; j<7; ++j) x[i][j+2] = x[i][1]>thresholds[j];
            for (int j=0; j<M; ++j) stats[i][j].add(x[i][j]);
            for (int u=0; u<i; ++u) for (int j=0; j<M; ++j)
                diffs[u][i][j].add(x[u][j]-x[i][j]);
        }
    }
    std::cout<<std::setprecision(9)<<"n="<<n<<" t="<<t<<" seed="<<seed
             <<" mean_size="<<sizes.sum/n<<" se="<<sizes.se(n)<<" mean_cores="<<cores.sum/n<<'\n';
    std::cout<<"label,D>=20,SE,mean,SE,tail.01,SE,tail.05,SE,tail.10,SE,tail.20,SE,tail.30,SE,tail.50,SE,tail.70,SE,D>=50,SE,D>=100,SE\n";
    auto print = [&](Stat s) { std::cout<<','<<s.sum/n<<','<<s.se(n); };
    for (int u=0; u<K; ++u) {
        std::cout<<u+1;
        for (auto s:stats[u]) print(s);
        print(deg50[u]); print(deg100[u]); std::cout<<'\n';
    }
    std::cout<<"negative paired differences (u,v,metric[0..8],older-younger,SE,z)\n";
    for (int u=0;u<K;++u) for (int v=u+1;v<K;++v) for (int j=0;j<M;++j) {
        auto s=diffs[u][v][j];
        if (s.sum<0) {
            std::cout<<u+1<<','<<v+1<<','<<j; print(s);
            std::cout<<','<<s.sum/n/s.se(n)<<'\n';
        }
    }
}
