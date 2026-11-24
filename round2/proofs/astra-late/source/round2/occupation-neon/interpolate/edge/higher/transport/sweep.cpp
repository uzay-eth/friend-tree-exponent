// Native skip-free Schur sweep. All arithmetic in the row proof is integer.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>
using I=std::int64_t;
using W=__int128_t;
constexpr I U=1000000000000LL,V=1000000000000000000LL;
struct S{int d,k;};
struct B{I lo=0,hi=0;};
struct Rate{I lo,hi; B gain;};
struct Toll {B a,b;};
I down(W a,I b){auto q=a/b; if(a<0&&a%b)--q;assert(q>=INT64_MIN&&q<=INT64_MAX);return I(q);}
I up(W a,I b){return -down(-a,b);}
B plus(B a,B b){return {a.lo+b.lo,a.hi+b.hi};}
B ratio(int n,int d){assert(d>0);return {down(W(n)*U,d),up(W(n)*U,d)};}
class Solver {
  std::array<std::vector<S>,5> dom;
  std::array<std::vector<B>,6> curr,next;
  std::array<size_t,6> len{};
  std::vector<Toll> tolls;
  using Vec=std::array<int,5>;
  I tested=0, sources=0, overflow=0; int apex,arity,sd,sk,paidmax,family,dmin,kmin;
  static constexpr I slo=565673000000LL,shi=565775000000LL;
  int idx(int i,S s)const{
    auto &p=dom[i];for(size_t j=0;j<p.size();++j)if(p[j].d==s.d&&p[j].k==s.k)return int(j);return -1;
  }
  size_t ix(int n,const Vec &x)const{
    size_t r=x[0];for(int i=1;i<n;++i)r=r*dom[i].size()+x[i];return r;
  }
  B get(int n,Vec x,bool gain=false)const{
    for(int i=0;i<n;++i)if(x[i]<0)return {};
    return (gain?next:curr)[n][ix(n,x)];
  }
  static void push(std::vector<Rate>&r,B b,B y){assert(0<=b.lo&&b.lo<=b.hi);if(b.hi)r.push_back({b.lo,b.hi,y});}
  B solve(B src,const std::vector<Rate>&r){
    B answer;
    for(bool hi:{false,true}){
      I d=slo;W a=W(hi?src.hi:src.lo)*U;
      std::vector<std::pair<I,I>> shifts{{0,shi-slo}};
      for(auto &e:r){I y=hi?e.gain.hi:e.gain.lo;a+=W(e.lo)*y;d+=e.lo;if(e.hi>e.lo)shifts.push_back({y,e.hi-e.lo});}
      std::sort(shifts.begin(),shifts.end());if(hi)std::reverse(shifts.begin(),shifts.end());
      for(auto [y,b]:shifts){if(hi?W(y)*d<=a:W(y)*d>=a)break;a+=W(y)*b;d+=b;}
      I v=hi?up(a,d):down(a,d);
      W check=W(hi?src.hi:src.lo)*U-W(v)*(slo+(((v<0)==hi)?shi-slo:0));
      for(auto &e:r){I z=(hi?e.gain.hi:e.gain.lo)-v;check+=W(z)*(((z>0)==hi)?e.hi:e.lo);}
      assert(hi?check<=0:check>=0); if(hi)answer.hi=v;else answer.lo=v;
    }
    assert(answer.lo<=answer.hi);return answer;
  }
  void work(int n,Vec x,int d){
    std::array<S,5> s{};s[0]={d,x[0]};for(int i=1;i<n;++i)s[i]=dom[i][x[i]];
    // P4: chain. Yb: first make the long arm, then the apex's short
    // arm, then its donor. Only that irreversible birth order is paid.
    const std::array<int,5> parent=family==0?std::array<int,5>{-1,0,1,2,3}:std::array<int,5>{-1,0,1,0,3};
    int promoter=n==5?-1:parent[n];
    // The paid shape can occur only after the required tagged promotions;
    // free physical promotions are always included below.
    for(int i=0;i<n;++i){int occ=0;for(int j=1;j<n;++j)occ+=parent[j]==i;if(s[i].d&&s[i].k<occ)return;}
    Toll toll{};
    if(n==5&&d>=dmin&&d<=paidmax&&x[0]>=kmin){Vec z=x;z[4]=0;toll=tolls[ix(5,z)/dom[4].size()];}
    std::vector<Rate> r;r.reserve(18);
    for(int i=0;i<n;++i){
      auto [di,k]=s[i];if(!di){
        assert(i==4&&k<=sk);++overflow;
        // For D>=sd+1 the live donor promotion has this interval.
        // Arity remains explicit, including the last retained promotion.
        auto z=x;z[i]=idx(i,{0,k+1});
        push(r,{ratio(std::max(0,sd-k),sd+1).lo,U},get(n,z));
        continue;
      }
      int children=0,L=di-k-1,external=0;assert(L>=0);
      Vec z=x;z[i]=i==0?x[0]:idx(i,{di+1,k});
      if(i==4&&di==sd)z[i]=idx(i,{0,k});
      B grown=(i==0&&d==apex)?B{}:get(n,z,i==0);
      B g{L*U,L*U};
      for(int j=0;j<n;++j){
        if(i==j)continue;
        if(parent[j]==i||parent[i]==j){
          if(parent[j]==i)++children;
          B rr=s[j].d?ratio(1,s[j].d):B{0,ratio(1,sd+1).hi};
          if(i==3&&j==4)push(r,rr,plus(grown,toll.b));else g=plus(g,rr);
        }
      }
      external=k-children+(i==0);assert(external>=0);
      push(r,plus(g,{0,external*U/2}),grown);
      if(L){
        z=x;z[i]=i==0?(k<arity?k+1:-1):idx(i,{di,k+1});
        B nextval=get(n,z);
        if(i==promoter&&z[i]>=0){
          z[n]=idx(n,{2,0});assert(z[n]>=0);B born=get(n+1,z);nextval=plus(nextval,born);
          sources+=(born.lo||born.hi);
        }
        push(r,ratio(L,di),nextval);
      }
    }
    curr[n][ix(n,x)]=solve(toll.a,r);++tested;
  }
  void visit(int n,int i,Vec x,int d){
    if(i==n){work(n,x,d);return;}
    int top=i?int(dom[i].size())-1:std::min(arity,d-1);
    for(int j=0;j<=top;++j){x[i]=i?j:top-j;visit(n,i+1,x,d);}
  }
public:
  void run(){
    int dm,km,dn,kn,dt,kt;
    std::cin>>apex>>arity>>dm>>km>>dn>>kn>>dt>>kt>>sd>>sk>>paidmax>>family;
    dmin=family==0?18:2;kmin=2;
    assert(std::cin&&(family==0||family==1)&&apex>=paidmax&&paidmax>=dmin&&sd>=2);
    std::array<int,5> dd{0,dm,dn,dt,sd},kk{0,km,kn,kt,sk};
    for(int i=1;i<5;++i){if(i==4)for(int k=sk;k>=0;--k)dom[i].push_back({0,k});
      for(int d=dd[i];d>=2;--d)for(int k=std::min(kk[i],d-1);k>=0;--k)dom[i].push_back({d,k});}
    len[1]=arity+1;
    for(int i=2;i<=5;++i)len[i]=len[i-1]*dom[i-1].size();
    for(int i=1;i<=5;++i){curr[i].resize(len[i]);next[i].resize(len[i]);}
    const size_t packing=len[5]/dom[4].size();
    // Reward tensors are streamed one apex degree at a time, descending.
    tolls.resize(packing);
    for(int d=apex;d>=2;--d){
      if(d>=dmin&&d<=paidmax){for(auto &q:tolls)std::cin>>q.a.lo>>q.a.hi>>q.b.lo>>q.b.hi;assert(std::cin);}
      for(int n=5;n>=1;--n)visit(n,0,{},d);
      if(d%16==0||d==2)std::cerr<<"completed apex degree "<<d<<", rows "<<tested<<"\n";
      std::swap(curr,next);
      for(auto &p:curr)std::fill(p.begin(),p.end(),B{});
    }
    std::cout<<next[1][0].lo<<" "<<next[1][0].hi<<" rows "<<tested<<" birth sources "<<sources<<" donor-tail rows "<<overflow<<" peak slab bytes "<<2*(len[1]+len[2]+len[3]+len[4]+len[5])*sizeof(B)<<"\n";
  }
};
int main(){std::ios::sync_with_stdio(false);std::cin.tie(nullptr);Solver{}.run();}
