// Original mixed D18 row arithmetic, specialized to the oriented P4 chain.
// Included by solver.cc after the common integer types and rounding helpers.
constexpr int states=18*18;
struct Tables {
    const I *p4, *hv, *rp, *rc;
    const int *parent, *child, *middle, *triple, *state;
    I nc, n4, n3, scale, rscale, lambda, lscale, sign;
};
struct Charges {
    std::array<Tables,2> table;
    bool exclude_old;
    static constexpr I unit=153153000000000; // lcm(1e9, 2..17)
    static constexpr I grid=1000000000000000000LL;

    explicit Charges(const Tables* t,bool e): table{t[0],t[1]},exclude_old(e) {
        for (const auto& x:table)
            check(x.scale==100000000000 && x.rscale==1000000000 && x.lscale==1000000 &&
                  (x.sign==1 || x.sign==-1),"unsupported D18 table scales");
    }
    static int id(Census c) { return c.d*18+c.k; }
    static I tick(I x) {
        check(x>-(I(1)<<50) && x<(I(1)<<50),"original charge tick too large");
        return x;
    }
    static I p4(const Tables& t,const std::array<Census,4>& x) {
        auto b=t.middle[id(x[1])], c=t.middle[id(x[2])];
        if (b<0 || c<0) return 0;
        auto a=t.parent[id(x[0])], d=t.child[id(x[3])];
        check(a>=0 && d>=0,"invalid P4 axes");
        return tick(t.p4[((I(a)*t.n4+b)*t.n4+c)*t.nc+d]);
    }
    static I hv(const Tables& t,Census a,Census b,Census c) {
        int j=t.triple[id(b)];
        if (j<0) return 0;
        auto i=t.parent[id(a)], k=t.child[id(c)];
        check(i>=0 && k>=0,"invalid hV axes");
        return tick(t.hv[(I(i)*t.n3+j)*t.nc+k]);
    }
    static void pair(Value& v,int j,W num,I den) {
        v[j]=ratio(num,den,false); v[j+1]=ratio(num,den,true);
    }
    int operator()(const std::array<Census,4>& x,Reward& out) const {
        int outside=-3, size=0;
        for (auto [d,k]:x) { outside+=k; size+=d-k; }
        size+=2*outside;
        if (!(size>20 || (size>12 && outside>0))) return -1;
        if (exclude_old && x[0].d<=12 && x[1].d<=8 && x[2].d<=6 && x[3].d<=16 &&
            x[0].k==1 && x[1].k<=3 && x[2].k==1 && x[3].k==0) return -1;
        auto member=[&](const Tables& t) { return t.middle[id(x[1])]>=0 && t.middle[id(x[2])]>=0; };
        check(member(table[0])==member(table[1]),"mismatched D18 P4 membership");
        if (!member(table[0])) return -1;

        for (int j=0;j<4;j+=2) {
            const auto& t=table[j/2];
            I own=p4(t,x);
            W base=-W(t.lambda)*(unit/t.lscale)*own;
            std::array<I,5> slope{};
            std::array<int,5> counts{1};
            for (int i=0;i<4;++i) {
                auto [d,k]=x[i];
                auto z=x;
                ++z[i].d;
                I gain=p4(t,z)-own;
                I leaves=d-k-1;
                I gain_rate=leaves*unit+(i?unit/x[i-1].d:0)+(i<3?unit/x[i+1].d:0);
                base+=W(gain_rate)*gain;
                if (leaves) {
                    z=x; ++z[i].k;
                    // Extract assigns new embeddings here only at size five;
                    // the mixed witness has no size-five family.
                    base+=W(leaves*(unit/d))*(p4(t,z)-own);
                }
                if (!i) slope[0]=t.sign*gain;
                if (k>int(i<3)) { slope[i+1]=t.sign*gain; counts[i+1]=k-int(i<3); }
            }
            // Two connected deletions: incoming hV on either oriented triple.
            auto incoming=[&](Census a,Census b,Census c,Census from,bool parent) {
                if (t.triple[id(b)]<0) return W(0);
                I old=hv(t,a,b,c);
                int s=t.state[id(parent?a:c)];
                check(s>=0,"invalid transfer rate index");
                I prior=(parent?t.rp:t.rc)[s];
                check(prior>=0 && prior<=t.rscale,"invalid transfer reciprocal");
                if (parent) ++a.d; else ++c.d;
                return W(unit/from.d-prior*(unit/t.rscale))*(hv(t,a,b,c)-old);
            };
            base+=incoming(x[1],x[2],x[3],x[0],true);
            base+=incoming(x[0],x[1],x[2],x[3],false);
            base*=t.sign;
            I stretch=grid/t.scale;
            for (int i=0;i<5;++i) {
                I neg=std::min(I(0),slope[i]);
                base+=W(counts[i])*neg*(unit/2);
                pair(out.source[i],j,-W(neg)*stretch,2);
                pair(out.slope[i],j,W(slope[i])*stretch,1);
            }
            check(base>=0,"negative original pointwise minimum");
            pair(out.minimum,j,base*stretch,unit);
        }
        return 1;
    }
};
