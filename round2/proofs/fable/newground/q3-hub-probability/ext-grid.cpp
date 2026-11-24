// Exact directed rounding on a binary grid, with a certified missing-mass correction.
#define main exact_main
#include "dp.cpp"
#undef main
#include <memory>
using UW=__uint128_t;
using I=int64_t;
constexpr U grid=U(1)<<62;

struct GridLevel {
 struct Entry {size_t off;int nr;};
 unordered_map<U,Entry> rows;
 vector<U> order;
 I *x;size_t words=0,capacity,total=0;
 bool rooted;
 GridLevel(int n,bool root):rooted(root){rows.reserve(freecount[n]);order.reserve(freecount[n]);capacity=root?::rooted[n]:freecount[n];
  x=(I*)mmap(nullptr,capacity*sizeof(*x),PROT_READ|PROT_WRITE,MAP_PRIVATE|MAP_ANONYMOUS,-1,0);if(x==MAP_FAILED)abort();}
 ~GridLevel(){munmap(x,capacity*sizeof(*x));}
 I *add(U k,int n){auto [j,ins]=rows.try_emplace(k,Entry{words,n});assert(j->second.nr==n);if(ins){order.push_back(k);total+=n;words+=rooted?n:1;assert(words<=capacity);}return x+j->second.off;}
};

array<U,2> grow(GridLevel &old,GridLevel &next,int n,bool terminal){
 array<U,2> subtotal{};
 U lc=1;for(int d=1;d<n;++d)lc=lcm(lc,U(d));
 const U denom=n*lc;
 const size_t page=sysconf(_SC_PAGESIZE);size_t dropped=0;
 for(auto key:old.order){auto row=old.rows.at(key);Tree t(key,n);assert(t.canon()==key && t.nr==row.nr);
  size_t used=(row.off*sizeof(I)/page)*page;if(used>=dropped+(8UL<<20)){assert(madvise((char*)old.x+dropped,used-dropped,MADV_DONTNEED)==0);dropped=used;}
  int reps[M];copy(t.reps,t.reps+t.nr,reps);
  U rates[M]{};for(int i=0;i<n;++i){int z=old.rooted?i:reps[t.orb[i]];for(int j=0;j<t.d[i];++j)rates[z]+=lc/t.d[t.a[i][j]];}
  assert(accumulate(rates,rates+n,U(0))==denom);
  I *src=old.x+row.off;U sum=0;for(int o=0;o<(old.rooted?row.nr:1);++o){assert(src[o]>=0);sum+=src[o];}assert(sum<=grid);
  for(int w=0;w<n;++w){U rate=rates[w];if(!rate)continue;
   t.join(w,n);++t.n;U child=t.canon();I *a=next.add(child,t.nr);
   U birth=(UW(sum)*rate)/denom;
   if(!next.rooted)a[0]+=birth;else a[t.orb[n]]+=terminal?-I(birth):I(birth);
   subtotal[1]+=birth;
   if(terminal)for(int o=0;o<row.nr;++o){U z=(UW(U(src[o]))*rate)/denom;a[t.orb[reps[o]]]+=I(z);subtotal[0]+=z;}
   --t.n;t.split(w,n);
  }
 }
 assert(next.rows.size()==freecount[n+1] && next.total==::rooted[n+1]);
 assert(subtotal[0]<=grid && subtotal[1]<=grid);
 cerr<<"grid m="<<n+1<<" free="<<next.rows.size()<<" rooted="<<next.total<<" missing="<<(grid-subtotal[0])<<","<<(grid-subtotal[1])<<"\n";
 return subtotal;
}

void assess_grid(GridLevel &lev,int n,bool improved,array<U,2> total){
 CachedClassifier cache(improved);vector<array<U,4>> buckets;
 I mass=0;
 for(auto key:lev.order){auto row=lev.rows.at(key);Tree t(key,n);assert(t.canon()==key && t.nr==row.nr);Classifier cls(t,improved);
  int upp[M];for(int o=0;o<row.nr;++o){int v=t.reps[o];upp[v]=cache.upper(cls,v);}for(int v=0;v<n;++v)upp[v]=upp[t.reps[t.orb[v]]];
  for(int o=0;o<row.nr;++o){I v=lev.x[row.off+o];mass+=v;if(!v)continue;
   auto [l,h]=cache.limits(cls,t.reps[o],upp);assert(!cache.less(l,cache.zero) && !cache.less(h,l) && !cache.less(cache.one,h));
   if((int)buckets.size()<=max(l,h))buckets.resize(max(l,h)+1);
   if(v>0){buckets[l][0]+=v;buckets[h][1]+=v;}else{buckets[h][2]-=v;buckets[l][3]-=v;}
  }
 }
 assert(mass==I(total[0])-I(total[1]));
 Q lo=0,hi=0;for(size_t j=0;j<buckets.size();++j){auto b=buckets[j];lo+=cache.end[j].x*(I(b[0])-I(b[2]));hi+=cache.end[j].x*(I(b[1])-I(b[3]));}
 lo=(lo-(grid-total[1]))/grid;hi=(hi+(grid-total[0]))/grid;
 cout<<"m "<<n<<" improved "<<improved<<" free "<<lev.rows.size()<<" rooted "<<lev.total<<" denominator "<<grid<<"\n";
 cout.precision(12);cout<<"d "<<n-1<<" "<<lo.get_d()<<" "<<hi.get_d()<<" "<<lo<<" "<<hi<<endl;
}
int main(int argc,char **argv){int maxn=argc>1?atoi(argv[1]):10;assert(maxn>=4 && maxn<M);
 auto cur=make_unique<GridLevel>(2,false);cur->add(52,1)[0]=grid;array<U,2> totals{};
 for(int n=2;n<maxn;++n){auto next=make_unique<GridLevel>(n+1,n+1>=maxn-1);totals=grow(*cur,*next,n,n+1==maxn);cur=move(next);}
 assess_grid(*cur,maxn,false,totals);assess_grid(*cur,maxn,true,totals);
}
