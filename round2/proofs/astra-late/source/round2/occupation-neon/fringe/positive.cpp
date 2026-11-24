// Positive directed grid fit; no upper output is required for killed charges.
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <vector>
#include <limits>
using I=int64_t;
using Z=__int128_t;
struct Edge{int j; I t;};
struct Row{I d; std::vector<Edge> e; std::vector<Z> h;};
int main(){
    std::ios::sync_with_stdio(false);std::cin.tie(nullptr);
    int n,m,p; if(!(std::cin>>n>>m>>p)||n<19||m<1||p<1)return 2;
    int start[19];for(auto& j:start)if(!(std::cin>>j)||j<0||j>=n)return 2;
    std::vector<Row> rows(n);
    for(int i=0;i<n;++i){
        auto& r=rows[i];int e;
        if(!(std::cin>>r.d>>e)||r.d<=0||e<0)return 2;
        for(int k=0;k<e;++k){int j;I t;if(!(std::cin>>j>>t)||j<0||j>=n||j==i||t<=0||t>r.d)return 2;r.e.push_back({j,t});}
        for(int k=0;k<m;++k){std::string s;if(!(std::cin>>s)||s.empty()||s.size()>28)return 2;
            Z x=0;for(char c:s){if(c<'0'||c>'9')return 2;x=10*x+c-'0';}r.h.push_back(x);}
    }
    std::vector<I> x(size_t(n)*m);
    for(int sweep=0;sweep<=p;++sweep)for(int i=0;i<n;++i){
        const auto& r=rows[i];for(int k=0;k<m;++k){
            Z h=r.h[k];for(auto [j,t]:r.e)h+=Z(t)*x[size_t(j)*m+k];
            auto& q=x[size_t(i)*m+k];
            if(h/r.d<q||h/r.d>std::numeric_limits<I>::max())return 3;
            if(sweep<p)q=h/r.d;
        }
    }
    for(int j:start){for(int k=0;k<m;++k){if(k)std::cout<<' ';std::cout<<x[size_t(j)*m+k];}std::cout<<'\n';}
    std::cerr<<"exact nonnegative finite residuals on "<<n<<" rows and "<<m<<" channels\n";
}
