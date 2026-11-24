#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <numeric>
#include <unordered_map>
#include <vector>
#include <gmpxx.h>
#include <sys/mman.h>
#include <unistd.h>
using namespace std;
using Z=mpz_class; using Q=mpq_class; using U=uint64_t; using W=mp_limb_t;
constexpr int M=25, G=8, K=40, cutoff=24;
const size_t rooted[]={0,1,1,2,4,9,20,48,115,286,719,1842,4766,12486,32973,87811,235381,634847,1721159,4688676,12826228,35221832,97055181,268282855,743724984};
const size_t freecount[]={0,1,1,1,2,3,6,11,23,47,106,235,551,1301,3159,7741,19320,48629,123867,317955,823065,2144505,5623756,14828074,39299897};
int width(U x){return bit_width(x);}
void append(U &x,U s){x=(x<<width(s))|s;}
struct Tree {
 int n=0, d[M]{}, a[M][M], nch[M+1], ch[M+1][M], leader[M+1], pos[M], ord[M], orb[M], reps[M], nr;
 U code[M+1];
 void join(int x,int y){a[x][d[x]++]=y;a[y][d[y]++]=x;}
 void split(int x,int y){assert(a[x][--d[x]]==y && a[y][--d[y]]==x);}
 int parse(U s,int &b,int p){assert((s>>b)&1); --b; int x=n++; if(p>=0)join(p,x);
   while((s>>b)&1)parse(s,b,x);
   --b; return x;
 }
 explicit Tree(U s,int size){int b=width(s)-1;if(b==2*size-1)parse(s,b,-1);
   else {assert(b==2*size+1);--b; int x=parse(s,b,-1); int y=parse(s,b,-1); join(x,y);--b;}
   assert(n==size && b==-1);
 }
 void enc(int x,int p){int k=0;for(int i=0;i<d[x];++i){int y=a[x][i];if(y==p)continue;enc(y,x);ch[x][k++]=y;} nch[x]=k; finish(x);}
 void finish(int x){int k=nch[x];auto c=ch[x]; for(int i=1;i<k;++i){int y=c[i],j=i;for(;j && code[y]<code[c[j-1]];--j)c[j]=c[j-1];c[j]=y;}
   U s=1;for(int i=0;i<k;++i)append(s,code[c[i]]);code[x]=s<<1;
 }
 void layout(int x,int &at,int rep){int start=at;if(x!=M){pos[x]=at;ord[at]=x;orb[x]=leader[at]=(rep==at?nr++:leader[rep]); if(rep==at)reps[nr-1]=x; ++at;}
   int off=x!=M; int first=off; U prev=0;
   for(int i=0;i<nch[x];++i){int y=ch[x][i]; if(code[y]!=prev){first=off;prev=code[y];}layout(y,at,rep+first);off+=width(code[y])/2;}
   assert(at-start==width(code[x])/2-(x==M));
 }
 U canon(){int deg[M];copy(d,d+n,deg);int queue[M]{},h=0,t=0;
   for(int i=0;i<n;++i)if(deg[i]==1)queue[t++]=i;
   while(n-h>2){int end=t;while(h<end){int x=queue[h++];for(int j=0;j<d[x];++j)if(--deg[a[x][j]]==1)queue[t++]=a[x][j];}}
   int x=queue[h], root=x;
   if(n-h==1)enc(x,-1);
   else{int y=queue[h+1];enc(x,y);enc(y,x);root=M;nch[M]=2;ch[M][0]=x;ch[M][1]=y;finish(M);}
   int at=0;nr=0;layout(root,at,0);assert(at==n);return code[root];
 }
};

struct Tables {
 Q B[M][M], F[M][K+2], FP[M][M], LP[M][M], H[M][K+2];
 Tables(){
  for(int c=0;c<M;++c){Q a=1,s=0;for(int k=0;k<40;++k){if(k+1<M)B[c][k+1]=s+a;s+=a;a/=c+k+1;}s+=a*Q(c+41,c+40);for(int k=1;k<M;++k)B[c][k]/=s;
   Q al[K+1],be[K+1];for(int k=1;k<=K;++k){int up=k*(k+c);Q den=up+k-(k-1)*be[k-1];al[k]=(1+(k-1)*al[k-1])/den;be[k]=up/den;}
   F[c][K+1]=Q(1,K);for(int k=K;k;--k){F[c][k]=al[k]+be[k]*F[c][k+1];assert(F[c][k]>=F[c][k+1] && F[c][k]<=1);}
  }
  // Improved leaf upper: preserve the parent degree across unmarked losses.
  for(int c=cutoff-1;c>=0;--c){
   for(int k=cutoff-c-1;k>=1;--k){int up=k*(k+c);
    FP[c][k]=((k==cutoff-c-1?F[c][k+1]:FP[c][k+1])*up+1+(k-1)*(k==1?Q(0):FP[c+1][k-1]))/(up+k);
    // Taking the minimum with the old table is valid, and helps preserve monotonicity.
    FP[c][k]=min(FP[c][k],F[c][k]);
    assert(FP[c][k]>=(k==cutoff-c-1?F[c][k+1]:FP[c][k+1]));
    assert((up+k)*FP[c][k]>=up*(k==cutoff-c-1?F[c][k+1]:FP[c][k+1])+1+(k-1)*(k==1?Q(0):FP[c+1][k-1]));
   }
  }
  // Improved lower: increase the frozen core degree after a loss, and continue with baseline at large degree.
  // At the all-leaf degree-cutoff star use zero as a conservative boundary value.
  const Q zero=0;
  for(int d=cutoff-1;d>=2;--d)for(int k=1;k<=d;++k){int c=d-k;const Q &birth=(k+1==cutoff?zero:d==cutoff-1?B[c][k+1]:LP[c][k+1]);const Q &loss=LP[c+1][k-1];
   Q next=(d*birth+loss)/(d+1);LP[c][k]=max(min(next,birth),B[c][k]);}
  // For a marked core target, count k unmarked leaves at its selecting neighbour.
  for(int c=1;c<M;++c){Q al[K+1],be[K+1];al[0]=1;
   for(int k=1;k<=K;++k){int up=k*(k+c);Q den=up+k+1-k*be[k-1];al[k]=(1+k*al[k-1])/den;be[k]=up/den;}
   H[c][K+1]=Q(1,K);for(int k=K;k;--k){H[c][k]=al[k]+be[k]*H[c][k+1];assert(H[c][k]>=H[c][k+1] && H[c][k]<=1);} H[c][0]=1;
  }
 }
} tabs;

struct Classifier {
 Tree &t; int l[M]{},c[M]; bool better;
 Classifier(Tree &tree,bool improve):t(tree),better(improve){for(int i=0;i<t.n;++i){for(int j=0;j<t.d[i];++j)l[i]+=t.d[t.a[i][j]]==1;c[i]=t.d[i]-l[i];}}
 Q upper(int i){if(t.d[i]==1){int w=t.a[i][0];return better?tabs.FP[c[w]][l[w]]:tabs.F[c[w]][l[w]];}
   if(l[i])return 1;
   Q q=0;for(int j=0;j<t.d[i];++j){int w=t.a[i][j];Q x(1,max(1,l[w]-1));if(better)x=min(x,tabs.H[c[w]][l[w]]);q+=x;if(q>=1)return 1;}return q;
 }
 pair<Q,Q> limits(int i){Q lo=0;if(t.d[i]>=2)lo=better?tabs.LP[c[i]][l[i]]:tabs.B[c[i]][l[i]];
   else if(better)lo=tabs.B[1][1]/(2*t.d[t.a[i][0]]);
   for(int j=0;j<t.d[i];++j){Q candidate=1-upper(t.a[i][j]);lo=max(lo,candidate);}return {lo,upper(i)};
 }
};

struct Row { size_t off; int nr; };
struct Level {
 unordered_map<U,Row> rows; vector<U> order; W *x=nullptr; int limbs,cols,masscols,shards,shard; bool unmarked; size_t total=0,words=0,capacity=0; Z den; vector<int> labels;
 Level(int n,Z denominator,int firstlabel=2,int lastlabel=M,int ns=1,int s=0):limbs(mpz_size(denominator.get_mpz_t())),shards(ns),shard(s),den(move(denominator)){
   // Before the first requested label, retain just the unrooted law.
   unmarked=n<firstlabel;if(unmarked)labels.push_back(n);else for(int j=max(2,firstlabel);j<=min(n,lastlabel);++j)labels.push_back(j);
   cols=labels.size();masscols=cols;assert(limbs<=G && n<M);auto budget=[&](size_t count){return ns==1?count:(count/ns)*101/100+1000;};rows.reserve(budget(freecount[n]));order.reserve(budget(freecount[n]));capacity=budget(unmarked?freecount[n]:rooted[n])*limbs*masscols;
   x=(W*)mmap(nullptr,capacity*sizeof(W),PROT_READ|PROT_WRITE,MAP_PRIVATE|MAP_ANONYMOUS,-1,0);if(x==MAP_FAILED)abort();
 }
 ~Level(){if(x)munmap(x,capacity*sizeof(W));}
 Level(const Level&)=delete;
 void swap(Level &other){rows.swap(other.rows);order.swap(other.order);std::swap(x,other.x);std::swap(limbs,other.limbs);std::swap(cols,other.cols);std::swap(masscols,other.masscols);std::swap(shards,other.shards);std::swap(shard,other.shard);std::swap(unmarked,other.unmarked);std::swap(total,other.total);std::swap(words,other.words);std::swap(capacity,other.capacity);den.swap(other.den);labels.swap(other.labels);}
 bool keeps(U key)const{if(shards==1)return true;key=(key^(key>>30))*0xbf58476d1ce4e5b9ULL;key=(key^(key>>27))*0x94d049bb133111ebULL;return (key^(key>>31))%shards==U(shard);}
 W* add(U key,int r){auto [j,inserted]=rows.try_emplace(key,Row{words,r});assert(j->second.nr==r);
   if(inserted){order.push_back(key);words+=size_t(unmarked?1:r)*masscols*limbs;assert(words<=capacity);total+=r;}
   return x+j->second.off;
 }
};
void propagate(const Level &old,Level &next,int n,U lc){
 const size_t page=sysconf(_SC_PAGESIZE);size_t dropped=0;
 for(U key:old.order){auto row=old.rows.at(key);
  size_t used=(row.off*sizeof(W)/page)*page;
  if(used>=dropped+(8UL<<20)){assert(madvise((char*)old.x+dropped,used-dropped,MADV_DONTNEED)==0);dropped=used;}
  Tree t(key,n);assert(t.canon()==key && t.nr==row.nr); int reps[M];copy(t.reps,t.reps+t.nr,reps);W unroot[G]{};
  W rates[M]{};for(int i=0;i<n;++i)for(int j=0;j<t.d[i];++j)rates[i]+=lc/t.d[t.a[i][j]];
  assert(accumulate(rates,rates+n,W(0))==n*lc);
  W const *src=old.x+row.off;
  for(int o=0;o<(old.unmarked?1:row.nr);++o)assert(!mpn_add_n(unroot,unroot,src+o*old.masscols*old.limbs,old.limbs));
  // Only the unrooted transfer may fold all attachment sites in the same orbit.
  bool unroot_only=old.unmarked;
  W group[M]{};if(unroot_only)for(int w=0;w<n;++w)group[reps[t.orb[w]]]+=rates[w];
  for(int w=0;w<n;++w){W rate=unroot_only?group[w]:rates[w];if(!rate)continue;t.join(w,n);++t.n;U child=t.canon();if(!next.keeps(child)){--t.n;t.split(w,n);continue;}W *dst=next.add(child,t.nr);
   auto add=[&](W *to,const W *from){W cy=mpn_addmul_1(to,from,old.limbs,rate);if(next.limbs>old.limbs){assert(next.limbs==old.limbs+1);assert(!__builtin_add_overflow(to[old.limbs],cy,&to[old.limbs]));}else assert(!cy);};
   if(next.labels.back()==n+1)add(dst+((next.unmarked?0:t.orb[n])*next.masscols+next.cols-1)*next.limbs,unroot);
   if(!unroot_only)for(int o=0;o<row.nr;++o){W *to=dst+t.orb[reps[o]]*next.masscols*next.limbs;
    for(int i=0;i<old.cols;++i)add(to+i*next.limbs,src+(o*old.masscols+i)*old.limbs);
   }
   --t.n;t.split(w,n);
  }
 }
}
struct Bucket {
 vector<W> x;
 Bucket(int size):x(size){}
};
struct CachedClassifier {
 struct Endpoint {Q x;U rank;int anti=-1;};
 map<Q,int> pool;
 vector<Endpoint> end;
 map<vector<int>,int> zeros;
 int one,zero,leaf_hi[M][M]{},solid_lo[M][M]{},leaf_lo[M]{};
 bool improved;
 int intern(const Q &x){auto [it,added]=pool.try_emplace(x,end.size());if(added){Z z=x.get_num()<<60;z/=x.get_den();end.push_back({x,z.get_ui()});}return it->second;}
 int complement(int i){if(end[i].anti<0){Q q=1-end[i].x;int j=intern(q);end[i].anti=j;end[j].anti=i;}return end[i].anti;}
 bool less(int x,int y){if(x==y)return false;if(end[x].rank!=end[y].rank)return end[x].rank<end[y].rank;return end[x].x<end[y].x;}
 CachedClassifier(bool better):improved(better){zero=intern(Q(0));one=intern(Q(1));
  for(int d=1;d<M;++d){leaf_lo[d]=improved?intern(tabs.B[1][1]/(2*d)):zero;
   for(int k=1;k<=d;++k){int c=d-k;leaf_hi[c][k]=intern(improved?tabs.FP[c][k]:tabs.F[c][k]);solid_lo[c][k]=intern(improved?tabs.LP[c][k]:tabs.B[c][k]);}}
 }
 int upper(Classifier &cls,int v){Tree &t=cls.t;if(t.d[v]==1){int w=t.a[v][0];return leaf_hi[cls.c[w]][cls.l[w]];}
  if(cls.l[v])return one;
  vector<int> key;key.reserve(t.d[v]);for(int j=0;j<t.d[v];++j){int w=t.a[v][j];if(cls.l[w]<=2 && !improved)return one;key.push_back((improved?cls.c[w]*M:0)+cls.l[w]);}sort(key.begin(),key.end());
  auto it=zeros.find(key);if(it!=zeros.end())return it->second;int id=intern(cls.upper(v));zeros.emplace(move(key),id);return id;
 }
 pair<int,int> limits(Classifier &cls,int v,const int *hi){Tree &t=cls.t;int id=t.d[v]==1?leaf_lo[t.d[t.a[v][0]]]:solid_lo[cls.c[v]][cls.l[v]];
  int best=hi[t.a[v][0]];for(int j=1;j<t.d[v];++j){int w=hi[t.a[v][j]];if(less(w,best))best=w;}
  int other=complement(best);if(less(id,other))id=other;
  return {id,hi[v]};
 }
};
Z readint(const W* x,int limbs){Z z;mpz_import(z.get_mpz_t(),limbs,-1,sizeof(W),0,0,x);return z;}
void assess(Level &lev,int n,bool improved,bool fractions){
 int C=lev.cols,L=lev.limbs,S=2*C-1;
 CachedClassifier cache(improved);vector<Bucket> buckets;vector<W> total(C*L), abs(L);
 for(auto key:lev.order){auto row=lev.rows.at(key);Tree t(key,n);assert(t.canon()==key && t.nr==row.nr);Classifier cls(t,improved);
  int upp[M];for(int o=0;o<row.nr;++o){int v=t.reps[o];upp[v]=cache.upper(cls,v);}for(int v=0;v<n;++v)upp[v]=upp[t.reps[t.orb[v]]];
  for(int o=0;o<row.nr;++o){auto [low,high]=cache.limits(cls,t.reps[o],upp);assert(!cache.less(low,cache.zero) && !cache.less(high,low) && !cache.less(cache.one,high));
   while((int)buckets.size()<=max(low,high))buckets.emplace_back((2*S+2*(C-1))*L);
   W *lo=buckets[low].x.data(), *hi=buckets[high].x.data();
   const W *src=lev.x+row.off+o*lev.masscols*L;
   auto plus=[&](W *a,const W *b){assert(!mpn_add_n(a,a,b,L));};
   for(int i=0;i<C;++i){const W *x=src+i*L;plus(total.data()+i*L,x);plus(lo+i*L,x);plus(hi+(S+i)*L,x);
    if(i+1==C || lev.labels[i]+1!=lev.labels[i+1])continue;
    int sign=mpn_cmp(x,x+L,L);int j=C+i;
    if(sign>0){mpn_sub_n(abs.data(),x,x+L,L);plus(lo+j*L,abs.data());plus(hi+(j+S)*L,abs.data());}
    if(sign<0){mpn_sub_n(abs.data(),x+L,x,L);plus(hi+(2*S+i)*L,abs.data());plus(lo+(2*S+C-1+i)*L,abs.data());}
   }
  }
 }
 // Signed sums for the differences require distinct positive and negative buckets.
 for(int i=0;i<C;++i){Z z=readint(total.data()+i*L,L);if(lev.shards==1)assert(z==lev.den);else cerr<<"mass "<<improved<<" "<<lev.labels[i]<<" "<<z<<" "<<lev.den<<"\n";}
 vector<Q> ans(2*S);
 for(size_t j=0;j<buckets.size();++j){const Q &q=cache.end[j].x;auto &b=buckets[j];for(int i=0;i<2*S;++i)ans[i]+=q*readint(b.x.data()+i*L,L);for(int i=0;i<C-1;++i){ans[C+i]-=q*readint(b.x.data()+(2*S+i)*L,L);ans[S+C+i]-=q*readint(b.x.data()+(2*S+C-1+i)*L,L);}}
 for(Q &a:ans)a/=lev.den;
 cout<<"m "<<n<<" improved "<<improved<<" free "<<lev.rows.size()<<" rooted "<<lev.total<<" denominator "<<lev.den<<"\n";
 for(int i=0;i<S;++i){if(i>=C && lev.labels[i-C]+1!=lev.labels[i-C+1])continue;
  cout<<(i<C?'p':'d')<<" "<<lev.labels[i<C?i:i-C]<<" "<<ans[i].get_d()<<" "<<ans[S+i].get_d();if(fractions)cout<<" "<<ans[i]<<" "<<ans[S+i];cout<<"\n";}cout<<flush;
}
int main(int argc,char **argv){int maxn=argc>1?atoi(argv[1]):10, first=argc>2?atoi(argv[2]):3, mode=argc>3?atoi(argv[3]):0, firstlabel=argc>4?atoi(argv[4]):2,lastlabel=argc>5?atoi(argv[5]):M,shards=argc>6?atoi(argv[6]):1,shard=argc>7?atoi(argv[7]):0;
 assert(maxn<M && maxn>=3 && firstlabel<=lastlabel && lastlabel>=2);cout.precision(12);Level cur(2,Z(1),firstlabel,lastlabel);assert(Tree(52,2).canon()==52);cur.add(52,1)[0]=1; // 110100 is the virtual root with two singleton children.
 assert(shards>0 && shard>=0 && shard<shards);
 for(int n=2;n<maxn;++n){U lc=1;for(int d=1;d<n;++d)lc=lcm(lc,U(d));Level next(n+1,cur.den*n*lc,firstlabel,lastlabel,n+1==maxn?shards:1,shard);propagate(cur,next,n,lc);
  if(next.shards==1)assert(next.rows.size()==freecount[n+1] && next.total==rooted[n+1]);
  cur.swap(next);
  if(n+1>=first && !cur.unmarked){if(mode!=2)assess(cur,n+1,false,true);if(mode)assess(cur,n+1,true,true);}
  cerr<<"completed m="<<n+1<<" rows="<<cur.rows.size()<<" rooted="<<cur.total<<" limbs="<<cur.limbs<<" words="<<cur.words<<"\n";
 }
 return 0;
}
