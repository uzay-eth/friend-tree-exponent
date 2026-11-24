#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <string>
#include <algorithm>
#include <map>
using namespace std;
using I=__int128_t;
using U=uint64_t;
constexpr U S=10000000000ULL;
int m, d[12], a[12][12], qs[12], u;
U low[12][12], leaf0[12], leaf1[12][12], core0[12];
I mass, den=1, nodes, agg[8][4], orbmass, orbagg[2][4];
struct Root { I p[8]{}; U lo,hi; int d,l; };
map<unsigned,Root> roots;
unsigned shape(int v,int par) {
    unsigned b[11]; int k=0;
    for(int j=0;j<d[v];++j) {
        int x=a[v][j]; if(x!=par) b[k++]=shape(x,v);
    }
    for(int j=1;j<k;++j) for(int h=j;h && b[h-1]>b[h];--h) swap(b[h-1],b[h]);
    unsigned t=1;
    for(int j=0;j<k;++j) t=(t<<(32-__builtin_clz(b[j])))|b[j];
    return t<<1;
}
string str(I a) {
    string z;
    do { z+=char('0'+a%10); a/=10; } while(a);
    reverse(z.begin(),z.end()); return z;
}
void dfs(int n, I p) {
    if(n==m) {
        int l[12];
        for(int i=0;i<n;++i) {
            l[i]=0;
            for(int j=0;j<d[i];++j) l[i]+=(d[a[i][j]]==1);
        }
        ++nodes; mass+=p;
        U ub[12];
        for(int i=0;i<n;++i) {
            if(d[i]==1) {
                int w=a[i][0];
                ub[i]=(d[w]==l[w]?leaf0:leaf1[d[w]-l[w]])[l[w]];
            } else if(!l[i]) {
                ub[i]=0;
                for(int j=0;j<d[i];++j) ub[i]+=core0[l[a[i][j]]];
                ub[i]=min(ub[i],S);
            } else ub[i]=S;
        }
        for(int i=0;i<min(n,8);++i) {
            U h=low[l[i]][d[i]-l[i]];
            for(int j=0;j<d[i];++j) h=max(h,S-ub[a[i][j]]);
            agg[i][0]+=p*h;
            agg[i][1]+=p*ub[i];
            agg[i][2]+=p*max(l[i]-1,0);
            agg[i][3]+=p*d[i];
            auto &r=roots[shape(i,-1)];
            if(any_of(r.p,r.p+8,[](I x){ return x!=0; })) assert(r.lo==h&&r.hi==ub[i]);
            r.p[i]+=p; r.lo=h; r.hi=ub[i]; r.d=d[i]; r.l=l[i];
        }
        if(d[2]>=2 && a[2][0]==0 && a[3][0]==2) {
            // Canonical start: start with 1-2; attach 3 to 1, then 4 to 3.
            orbmass+=p;
            for(int j=0;j<2;++j) {
                int i=j?3:2;
                U h=low[l[i]][d[i]-l[i]];
                for(int jj=0;jj<d[i];++jj) h=max(h,S-ub[a[i][jj]]);
                orbagg[j][0]+=p*h;
                orbagg[j][1]+=p*ub[i];
                orbagg[j][2]+=p*max(l[i]-1,0);
                orbagg[j][3]+=p*d[i];
            }
        }
        return;
    }
    for(int i=0;i<n;++i) {
        int sum=0;
        for(int j=0;j<d[i];++j) sum+=qs[n]/d[a[i][j]];
        a[n][0]=i; d[n]=1; a[i][d[i]++]=n;
        dfs(n+1,p*sum);
        --d[i];
    }
}
int main(int argc,char**argv) {
    assert(argc==2); m=atoi(argv[1]); assert(3<=m&&m<=11);
    for(int i=0;i<=m;++i) for(int j=0;j<=m;++j) cin>>low[i][j];
    for(int i=0;i<=m;++i) cin>>leaf0[i];
    for(int i=0;i<=m;++i) for(int j=0;j<=m;++j) cin>>leaf1[i][j];
    for(int i=0;i<=m;++i) cin>>core0[i];
    assert(cin);
    for(int n=2;n<m;++n) {
        qs[n]=1; for(int k=1;k<n;++k) qs[n]=lcm(qs[n],k);
        den*=n*qs[n];
    }
    d[0]=d[1]=1; a[0][0]=1; a[1][0]=0; dfs(2,1);
    assert(den==mass); assert(nodes==[&](){I x=1;for(int n=2;n<m;++n)x*=n;return x;}());
    cout<<"m "<<m<<" nodes "<<str(nodes)<<" mass "<<str(mass)<<" hub_prob_den "<<str(den*S)<<" means_den "<<str(den*m)<<'\n';
    for(int i=0;i<min(m,8);++i) {
        cout<<i+1;
        for(auto x:agg[i]) cout<<' '<<str(x);
        cout<<'\n';
    }
    if(m>=4) {
        cerr<<"path m "<<m<<" den "<<str(orbmass*S)<<'\n';
        for(int i=0;i<2;++i) {
            cerr<<i<<' '<<str(orbagg[i][0])<<' '<<str(orbagg[i][1])<<'\n';
        }
    }
    cout<<"rooted shapes "<<roots.size()<<'\n';
    for(int i=0;i<min(m,8);++i) for(int j=i+1;j<min(m,8);++j) {
        I lo=0,hi=0,ml=0,mh=0;
        for(auto &[key,r]:roots) {
            auto z=r.p[i]-r.p[j];
            lo+=z*(z<0?r.hi:r.lo); hi+=z*(z<0?r.lo:r.hi);
            ml+=z*(z<0?r.d:max(0,r.l-1)); mh+=z*(z<0?max(0,r.l-1):r.d);
        }
        auto show=[](I x){return (x<0?"-":"")+str(x<0?-x:x);};
        cout<<i+1<<'-'<<j+1<<' '<<show(lo)<<' '<<show(hi)<<' '<<show(ml)<<' '<<show(mh)<<'\n';
    }
}
