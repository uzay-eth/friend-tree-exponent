// Reuse the lead's exact tree encoding and signed-law propagation only.
#define main previous_main
#include "../../attacks/q3/signed_orbit_dp.cpp"
#undef main
struct Num {
    static constexpr int W=6;
    u64 a[W];
    Num(u64 x=0):a{x} {}
    bool operator==(const Num& b) const {return equal(a,a+W,b.a);}
    bool operator!=(const Num& b) const {return !(*this==b);}
    bool operator<(const Num& b) const {
        for(int j=W-1;j>=0;--j) if(a[j]!=b.a[j]) return a[j]<b.a[j];
        return false;
    }
    bool operator>=(const Num& b) const {return !(*this<b);}
    Num& operator+=(const Num& b) {
        u128 c=0;
        for(int j=0;j<W;++j) {c+=(u128)a[j]+b.a[j]; a[j]=c; c>>=64;}
        if(c || a[W-1]>>63) abort();
        return *this;
    }
    Num& operator-=(const Num& b) {
        u128 c=0;
        for(int j=0;j<W;++j) {c=(u128)a[j]-b.a[j]-(c>>127); a[j]=c;}
        if(c>>127) abort();
        return *this;
    }
    Num operator+(const Num& b) const {Num z=*this; return z+=b;}
    Num operator-(const Num& b) const {Num z=*this; return z-=b;}
    friend ostream& operator<<(ostream& o,const Num& x) {
        mpz_class z; mpz_import(z.get_mpz_t(),W,-1,8,0,0,x.a); return o<<z;
    }
};

struct Dinic {
    struct Arc { int to,next,idx; }; // -1 infinite, 0 unallocated zero, >0 capacity index
    vector<Arc> arcs;
    vector<Num> caps{0};
    vector<int> head, ptr, dist, fifo;
    Dinic(int n, size_t m):head(n,-1),ptr(n),dist(n),fifo(n) {arcs.reserve(m);}
    void add(int a,int b,const Num* c=nullptr) {
        int i=-1;
        if(c) {i=caps.size(); caps.push_back(*c);}
        arcs.push_back({b,head[a],i}); head[a]=arcs.size()-1;
        arcs.push_back({a,head[b],0}); head[b]=arcs.size()-1;
    }
    bool live(const Arc& e) const {return e.idx<0 || caps[e.idx]!=0;}
    void push(int j,const Num& x) {
        auto &a=arcs[j], &b=arcs[j^1];
        if(a.idx>=0) caps[a.idx]-=x;
        if(b.idx==0) {b.idx=caps.size(); caps.push_back(0);}
        if(b.idx>0) caps[b.idx]+=x;
    }
    bool levels(int s,int t) {
        fill(dist.begin(),dist.end(),-1); dist[s]=0;
        int h=0,k=1; fifo[0]=s;
        while(h<k) {
            int a=fifo[h++];
            if(dist[t]>=0 && dist[a]>=dist[t]) break;
            for(int j=head[a];j>=0;j=arcs[j].next) {
                const auto& e=arcs[j];
                if (dist[e.to]<0 && live(e)) {dist[e.to]=dist[a]+1; fifo[k++]=e.to;}
            }
        }
        return dist[t]>=0;
    }
    Num dfs(int a,int t,const Num& lim) {
        if(a==t) return lim;
        for(int &j=ptr[a];j>=0;j=arcs[j].next) {
            auto& e=arcs[j];
            if(dist[e.to]!=dist[a]+1 || !live(e)) continue;
            Num bound=e.idx<0 ? lim : min(lim,caps[e.idx]);
            Num x=dfs(e.to,t,bound);
            if(x!=0) {push(j,x); return x;}
        }
        return 0;
    }
    Num flow(int s,int t,const Num& want) {
        Num x=0;
        while(x<want && levels(s,t)) {
            ptr=head;
            while(x<want) {Num y=dfs(s,t,want-x); if(y==0) break; x+=y;}
        }
        return x;
    }
};

static Num capacity(const Big& a) {
    const Big b=bneg(a)?bnegate(a):a;
    if(b.l[LIMBS-1]>>63 || Num::W!=LIMBS) abort();
    Num x=0;
    copy_n(b.l,Num::W,x.a);
    return x;
}

static bool is_universal(const Tree& t) {
    if(t.deg[0]<2) return false;
    for(int j=0;j<t.deg[0];++j) if(t.deg[t.adj[0][j]]==1) return true;
    return false;
}

// At later label's birth, degree-two positive states are universal. First
// exhaust ordinary sources, then add a universal source for the exact deficit.
static bool solve(const Table& diff,int n,bool plant,const mpz_class& den) {
    vector<vector<u64>> shapes(n), leaf(n+1);
    unordered_map<u64,int> ids;
    Tree t;
    auto add=[&](u64 key,int s) {
        if(ids.emplace(key,ids.size()).second) shapes[s].push_back(key);
    };
    add(2,1);
    for(int s=1;s<n-1;++s) for(u64 key:shapes[s]) {
        decode(key,s,t);
        for(int v=0;v<s;++v) {
            t.adj[v][t.deg[v]++]=s; t.deg[s]=1; t.adj[s][0]=v;
            add(canon(t,0,-1).bits,s+1);
            --t.deg[v];
        }
    }
    ids.clear();
    auto leafcode=[](u64 key,int s) {return (1ULL<<(2*s+1))|(key<<1);};
    for(int s=1;s<n;++s) for(u64 key:shapes[s]) {
        u64 k=leafcode(key,s);
        ids.emplace(k,ids.size()); leaf[s+1].push_back(k);
    }
    int src=ids.size(), sink=src+1;
    Dinic f(sink+1,30ull*ids.size());
    vector<pair<int,int>> branches;
    for(int s=2;s<=n;++s) for(u64 key:leaf[s]) {
        decode(key,s,t);
        // Decode labels in preorder: the first neighbor of nonroot is parent.
        int depth[MAXN]={}; const int id=ids.at(key);
        branches.clear();
        for(int v=1;v<s;++v) {
            int p=t.adj[v][0], slot=0;
            depth[v]=depth[p]+1;
            if(p==0) continue;
            while(t.adj[p][slot]!=v) ++slot;
            swap(t.adj[p][slot],t.adj[p][--t.deg[p]]);
            int smaller=ids.at(canon(t,0,-1).bits);
            swap(t.adj[p][slot],t.adj[p][t.deg[p]++]);
            if((depth[p]%2==0)!=plant) branches.emplace_back(id,smaller);
            else branches.emplace_back(smaller,id);
        }
        sort(branches.begin(),branches.end());
        branches.erase(unique(branches.begin(),branches.end()),branches.end());
        for(auto [a,b]:branches) f.add(a,b);
    }
    Num want=0, sinkmass=0, universal=0;
    for(size_t i=0;i<diff.cap;++i) {
        if(!diff.keys[i]||biszero(diff.vals[i])) continue;
        u64 key=diff.keys[i]; Num mass=capacity(diff.vals[i]);
        decode(key,n,t);
        if(bneg(diff.vals[i])) {
            if(t.deg[0]!=1) abort();
            f.add(ids.at(key),sink,&mass); sinkmass+=mass;
        } else if(!plant && is_universal(t)) {
            universal+=mass;
        } else if(plant && is_universal(t)) {
            // Planted error: drop its permitted universal edges.
            want+=mass;
        } else {
            f.add(src,ids.at(key),&mass); want+=mass;
        }
    }
    if(want+universal!=sinkmass) abort();
    ids.clear(); ids.rehash(0);
    Num flow=f.flow(src,sink,want);
    cout<<"pair="<<n-1<<","<<n<<" leaf_shapes="<<leaf[n].size()
        <<" vertices="<<src<<" den_bits="<<mpz_sizeinbase(den.get_mpz_t(),2)
        <<" edges="<<f.arcs.size()/2<<" flow="<<flow<<"/"<<want
        <<" universal="<<universal<<" "<<(flow==want?"PASS":"FAIL")<<endl;
    return flow==want;
}

int main(int argc,char** argv) {
    if(argc<2) return 2;
    int stop=atoi(argv[1]), start=argc>2?atoi(argv[2]):3;
    bool plant=getenv("PLANT")!=nullptr;
    if(start<3 || stop>21 || start>stop) return 2;
    Table u; u.init(4); Big one=bzero(); one.l[0]=1; u.add(0b110100ULL,one);
    Table old;
    mpz_class den=1;
    for(int n=3;n<=stop;++n) {
        den*= (n-1)*lcm_upto(n-2);
        // Every path numerator by itself is <= Den. Even before
        // cancellation, signed additions have absolute sum <= 2 Den.
        if(mpz_sizeinbase(den.get_mpz_t(),2)+1 >= 64*LIMBS) abort();
        static const size_t rooted[]={0,1,1,2,4,9,20,48,115,286,719,1842,4766,12486,32973,87811,235381,634847,1721159,4688676,12826228,35221832};
        static const size_t free[]={0,1,1,1,2,3,6,11,23,47,106,235,551,1301,3159,7741,19320,48629,123867,317955,823065,2144505};
        Table q, un, diff;
        q.init(rooted[n-1]); propagate(u,n-1,q,1,false,true);
        if(n<stop) {un.init(free[n]); propagate(u,n-1,un,2,false,true);}
        if(q.total()!=den||q.cnt!=rooted[n-1]||(n<stop&&(un.total()!=den||un.cnt!=free[n]))) abort();
        u.clear();
        if(n>=start) {
            diff.init(q.cnt+rooted[n-2]);
            for(size_t i=0;i<q.cap;++i) if(q.keys[i]) diff.add(q.keys[i],bnegate(q.vals[i]));
            if(n==3) {diff.add(0b110100ULL,one); diff.add(0b111000ULL,one);}
            else propagate(old,n-1,diff,0);
            old.clear();
            if(n==stop) q.clear();
            if(diff.total()!=0) abort();
            if(!solve(diff,n,plant,den)) return 1;
        }
        old=std::move(q); u=std::move(un);
    }
}
