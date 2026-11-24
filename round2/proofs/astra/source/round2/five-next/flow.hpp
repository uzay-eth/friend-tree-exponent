// Discounted scores on connected core sets in promotion order.
#pragma once
#include "../motif-high/snapshot.hpp"
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>

namespace next5 {
using I=mh::I;
using W=mh::Wide;
constexpr I S=I(1)<<40, sp=1133, sq=2000;
inline void need(bool good,const char* m){if(!good)throw std::runtime_error(m);}
inline I floorq(W n,I d){need(d>0,"nonpositive fixed-point divisor");W x=n/d-(n%d<0);
    need(x>=INT64_MIN&&x<=INT64_MAX,"fixed-point quotient overflow");return I(x);}
struct Cell {
    int d=0,c=0;bool high=false,tail=false;
    Cell gain()const{return {d+1,c,false,false};}
};
using Five=std::array<Cell,6>;
using Parent=std::array<int,6>;
struct Cut {int d,c,high_c;bool stop=false,high_only=false;};
using Cuts=std::vector<Cut>;
struct Axis {
    Cut cut;int children;
    std::vector<Cell> values;
    std::vector<int> gain,bud;
    std::vector<std::vector<int>> inner;
    std::vector<int> pooled;
    Axis(Cut c,int k):cut(c),children(k),inner(c.d+1,std::vector<int>(c.c+1,-1)),pooled(std::max(c.high_c,k+1)+1,-1){
        cut.high_c=std::max(c.high_c,k+1);
        for(int d=2;!c.high_only&&d<=c.d;++d)for(int j=k+1;j<=std::min(d,c.c);++j){
            inner[d][j]=int(values.size());values.push_back({d,j});
        }
        for(int j=k+1;j<=cut.high_c;++j){
            if(c.stop)continue;
            pooled[j]=int(values.size());values.push_back({c.d+1,j,true,j==cut.high_c});
        }
        for(auto a:values){
            gain.push_back(a.high?-1:index(a.gain()));
            if(a.tail){bud.push_back(-1);continue;}
            ++a.c;bud.push_back(index(a));
        }
    }
    int index(Cell x)const{
        if(x.c<=children)return -1;
        if(x.high&&(x.d<=cut.d||(x.tail&&x.c<cut.high_c)))return -1;
        if(x.high||x.d>cut.d)return pooled[std::min(x.c,cut.high_c)];
        return x.c<=cut.c&&x.c<=x.d?inner[x.d][x.c]:-1;
    }
};
struct Mesh {
    int n,id,donor_piece=-1;Parent par,k{},stride{},width{};
    std::vector<Axis> axes;
    std::vector<I> weights,births;
    int count=1;
    Mesh(Parent tree,int sh,const Cuts& cuts):n(int(cuts.size())),id(sh),par(tree){
        for(int j=1;j<n;++j)++k[par[j]];
        for(int j=0;j<n;++j)axes.emplace_back(cuts[j],k[j]);
        for(int j=n;j-->0;){
            width[j]=int(axes[j].values.size());stride[j]=count;
            need(int64_t(count)*width[j]<150000000,"mesh exceeds per-shape memory allowance");count*=width[j];
        }
        if(count)births.resize(count/width[n-1]);
    }
    I value(const Five& x)const{
        if(weights.empty())return 0;
        int ix=0;for(int j=0;j<n;++j){int v=axes[j].index(x[j]);if(v<0)return 0;ix+=v*stride[j];}
        return weights[ix];
    }
    I entry(const Five& x)const{
        if(births.empty())return 0;
        int ix=0;for(int j=0;j<n-1;++j){int v=axes[j].index(x[j]);if(v<0)return 0;ix+=v*stride[j];}
        return births[ix/width[n-1]];
    }
    I bounded_entry(const Five& x)const{
        if(births.empty())return 0;
        int ix=0;Parent uncertain{};int num=0,visits=1;
        for(int j=0;j<n-1;++j){
            int v=axes[j].index(x[j]);if(v>=0){ix+=v*stride[j];continue;}
            auto a=x[j];const auto& ax=axes[j];
            if(!a.high||!a.tail||a.d<=ax.cut.d||a.c>=ax.cut.high_c||ax.cut.stop)return 0;
            visits*=ax.cut.high_c-a.c+1;if(visits>96)return 0;uncertain[num++]=j;
        }
        if(!num)return births[ix/width[n-1]];
        I best=S;
        auto visit=[&](auto&& self,int k,int at)->void{
            if(k==num){best=std::min(best,births[at/width[n-1]]);return;}
            int j=uncertain[k];const auto& ax=axes[j];
            for(int c=x[j].c;c<=ax.cut.high_c&&best;++c){
                int v=ax.pooled[c];if(v<0){best=0;break;}self(self,k+1,at+v*stride[j]);
            }
        };
        visit(visit,0,ix);return best;
    }
    using Source=std::function<I(const Parent&,const Five&)>;
    using Entry=std::function<I(int,const Five&)>;
    void solve(const Source& force,const Entry& child,bool keep,
               const std::function<bool(const Parent&)>& active={}){
        weights.assign(count,0);
        for(int row=count;row-->0;){
            Parent pos{};Five x{};I q=2;
            for(int j=0;j<n;++j){pos[j]=(row/stride[j])%width[j];x[j]=axes[j].values[pos[j]];
                need(q<=INT64_MAX/x[j].d,"common denominator overflow");q*=x[j].d;}
            if(active&&!active(pos))continue;
            W num=W(sq)*q*(force?force(pos,x):0);I den=sp*q;
            std::array<I,12> lo{},hi{},to{};int rates=0;
            auto add=[&](I lower,I upper,I z){
                if(!upper)return;
                if(lower==upper){num+=W(sq)*lower*z;den+=sq*lower;}
                else{lo[rates]=lower;hi[rates]=upper;to[rates++]=z;}
            };
            for(int j=0;j<n;++j){
                auto a=x[j];
                auto dest=[&](int label){return label<0?I(0):weights[row+(label-pos[j])*stride[j]];};
                if(!a.high){
                    I leaf=a.d-a.c,gl=leaf*q,gh=gl+(a.c-k[j]-(j!=0))*q/2;
                    auto near=[&](int h){gh+=q/x[h].d;if(!x[h].high)gl+=q/x[h].d;};
                    if(j)near(par[j]);
                    for(int h=j+1;h<n;++h)if(par[h]==j)near(h);
                    add(gl,gh,dest(axes[j].gain[pos[j]]));
                    if(!leaf)continue;
                }
                if(a.tail)continue; // its own label stays fixed and omitted entries are nonnegative
                I dep=dest(axes[j].bud[pos[j]]);
                if(child){auto after=x;++after[j].c;after[n]={2,1};dep+=child(j,after);}
                I b=I(std::max(0,a.d-a.c))*(q/a.d);
                add(b,a.high?q:b,dep);
            }
            double guess=double(num)/den;
            for(int k=0;k<=rates;++k){double top=double(num),bottom=den;
                for(int j=0;j<rates;++j){double r=double(sq)*(to[j]<guess?hi[j]:lo[j]);top+=r*to[j];bottom+=r;}
                double next=top/bottom;if(std::abs(next-guess)<0.125){guess=next;break;}guess=next;
            }
            need(std::isfinite(guess)&&guess>=0&&guess<1e18,"unsafe Bellman proposal");
            I z=std::min(I(guess),S);
            auto residual=[&](I v){W r=num-W(den)*v;
                for(int j=0;j<rates;++j)r+=W(sq)*(to[j]-v)*(to[j]<v?hi[j]:lo[j]);
                return r;};
            W r=residual(z);if(r<0)z=std::max<I>(0,z+floorq(r,den));
            need(residual(z)>=0,"negative integer Bellman residual");weights[row]=z;
            if(!x[n-1].high&&x[n-1].d==2&&x[n-1].c==1)births[row/width[n-1]]=z;
        }
        if(!keep)std::vector<I>().swap(weights);
    }
};
inline std::array<std::vector<Parent>,6> shapes(){
    std::array<std::vector<Parent>,6> z;z[1].push_back({-1});
    for(int n=1;n<5;++n)for(auto p:z[n])for(int i=0;i<n;++i){p[n]=i;z[n+1].push_back(p);}
    return z;
}
inline I donor(I a,Cell z){
    return a<0?(z.high?0:(-a)/z.d):floorq(W(a)*(z.d-2),2*z.d);
}
} // namespace next5
