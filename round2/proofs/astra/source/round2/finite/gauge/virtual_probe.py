"""Floating diagnostics for virtual marks in the constant 1/2 residual gauge."""
import argparse

def apply(d,s,labels,K,u,extra='uniform',split=None):
    q=s-d
    actual=lambda D: 1/split if split and D>=split else .5
    cost=[actual(D)/(actual(D)+d)*u[max((i for i,cap in enumerate(labels) if cap<=D),default=0)] for D in range(K+1)]
    out=[]
    for i,tag in enumerate(labels):
        row=[cost[-1]/q]*(K+1)
        for D in range(K-1,1,-1):
            new=[0.]*(D+1)
            for c in range(D,0,-1):
                L=D-c; bud=L/D; f=cost[D]+(new[c+1] if c<D else 0)
                coeff={'uniform':1+(.5+d)/actual(tag),'intrinsic':1,'silent':0}[extra]
                e=(c-1)/2+coeff/tag
                new[c]=max((L*row[c]+bud*f)/(q+L+bud),((L+e)*row[c]+bud*f)/(q+L+e+bud))
            row=new
        out.append(row[1])
    return out

def estimate(d,s,labels,K,steps,extra,split=None):
    u=[1.]*len(labels)
    for _ in range(steps):
        v=apply(d,s,labels,K,u,extra,split)
        low=min(x/y for x,y in zip(v,u));hi=max(x/y for x,y in zip(v,u))
        u=[x/v[0] for x in v]
    return low,hi,u

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--labels',default='2,3,4,6,9,14,24')
    p.add_argument('--cuts',default='-0.08,-0.04,0,0.04,0.08')
    p.add_argument('--s',type=float,default=.565)
    p.add_argument('--K',type=int,default=36)
    p.add_argument('--n',type=int,default=12)
    p.add_argument('--extra',default='uniform')
    p.add_argument('--split',type=int)
    a=p.parse_args()
    for d in map(float,a.cuts.split(',')):
        lo,hi,u=estimate(d,a.s,list(map(int,a.labels.split(','))),a.K,a.n,a.extra,a.split)
        print('d',d,'operator ratio',lo,hi,'relative newborn profile',','.join(f'{x:.3f}' for x in u))
