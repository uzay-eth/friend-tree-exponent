"""Directed rational all-degree bounds for two positive physical F classes."""
from argparse import ArgumentParser
from fractions import Fraction as Q

G=10**20
LO=Q('.565676636772')
HI=Q('.565771871520')


def rnd(x,side):return Q(side*((side*x*G)//1),G)


def choose(v,side):return (min if side==1 else max)(v)


class Tables:
    def __init__(self,s,H,side,memory=True,mark=2,min_children=1):
        self.s,self.H,self.side,self.memory=s,H,side,memory
        self.lam=s+Q(1,mark)
        self.min_children=min_children

    def child(self,j):
        """∫(a−p) dt before the tagged child changes, jointly with root."""
        H,sgn=self.H,self.side
        a=Q(1,j);lam=self.lam+Q(3,2)
        tail=(a-Q(1,H))/(lam+Q(1,H)) if sgn==1 else a/lam
        vs=[None]*(H+1)
        vs[H]=[rnd(tail,sgn) if sgn<0 or k>=self.min_children else Q(0) for k in range(H)]
        for d in range(H-1,max(j,2)-1,-1):
            p=Q(1,d);vs[d]=[None]*d
            for k in range(d-1,0,-1):
                l=d-k-1;b=Q(l,d)
                x=vs[d+1][k]
                value=(a-p if k>=self.min_children else 0)+(b*vs[d][k+1] if l else 0)
                g0=l+(Q(1,2) if self.memory else 0)
                g1=l+Q(k,2)
                vs[d][k]=rnd(choose(((value+g*x)/(lam+p+b+g) for g in (g0,g1)),sgn),sgn)
        return vs

    def parent(self,fees,tail_rate,min_children=0):
        """Birth compensator, with a nonzero all-time high-degree boundary."""
        lam,H,sgn=self.lam,self.H,self.side
        v=[None]*(H+1);v[H]=[]
        for k in range(H):
            b=Q(H-k-1,H)
            eta=Q(1,6) if k==0 else Q(0)
            w=b/(lam+eta) if sgn==1 else 1/lam-(1-b)/(lam+1)
            v[H].append(rnd(tail_rate*w,sgn) if k>=min_children else Q(0))
        for d in range(H-1,1,-1):
            v[d]=[None]*d
            for k in range(d-1,-1,-1):
                l=d-k-1;b=Q(l,d);x=v[d+1][k]
                a=b*(fees(d,k)+v[d][k+1]) if l else 0
                v[d][k]=rnd(choose(((a+g*x)/(lam+b+g) for g in (l,l+Q(k,2))),sgn),sgn)
        return v[2][0]

    def classes(self):
        s,sgn=self.s,self.side
        step={j:self.child(j) for j in range(2,21)}
        immediate=lambda j:Q(1,2)/(s+Q(1,j)+Q(3,2))
        def fees(d,k):
            j=min(d,20)
            return immediate(j)*step[j][d][k+1]
        j=20
        tail=(Q(1,j)-Q(1,self.H))/(self.lam+Q(3,2)+Q(1,self.H)) if sgn==1 else Q(1,j)/(self.lam+Q(3,2))
        copy=self.parent(fees,tail*immediate(20),self.min_children-1 if sgn==1 else 0)
        # No old copy at all: exactly one nontrigger newborn, then trigger
        # before the new original particle's first intrinsic jump.
        direct=self.parent(lambda d,k: immediate(min(d,20)) if d>2 else 0,immediate(20))
        return copy,direct


def tests():
    # Same constant future reward in every state must survive all gain and
    # promotion actions.  The high-degree lower/upper boundary sandwiches it.
    # More concretely, compare a one-state virgin child row with its exact
    # exponential holding formula and check the copy source decomposition.
    s=Q(3,5);p=Q(1,7);a=Q(1,5);h=Q(1,2);r=Q(11,4)
    z=(a-p)/(s+Q(1,2)+Q(3,2)+p)
    assert (a-p+r*z)/(s+Q(1,2)+Q(3,2)+p+r)==z
    assert (a-p)*h/(s+a+Q(3,2)) + p*0 == (s+Q(1,2)+Q(3,2)+p)*z*h/(s+a+Q(3,2))
    for H in (24,32):
        v=[Tables(s,H,side).classes() for side in (1,-1)]
        assert all(0<l<u for l,u in zip(*v))
        if H==24:old=v
        else:
            assert all(old[0][i]<=v[0][i]<=v[1][i]<=old[1][i] for i in range(2))
    print('exact control and nested-boundary checks passed',flush=True)


def physical_tests():
    from pathlib import Path
    import sys
    sys.path.insert(0,str(Path(__file__).resolve().parents[4]/'finite'))
    import residual_entries as p
    p.MARK_CAP=20
    virgin=(3,1,())
    for L in range(1,5):
        old=(2,5,((3,3,()),))
        for ch in ((virgin,),(old,virgin)):
            t=(2,L,ch)
            _,d,_=p.statistics(t)
            rates=dict(p.intrinsic_moves(t))
            root_gain=rates[(2,L+1,ch)]
            assert root_gain==L+sum(Q(1,p.statistics(y)[1]) for y in ch)
            tagged_gain=tuple(sorted(p.plus(y) if y==virgin else y for y in ch))
            assert rates[(2,L,tagged_gain)]==1+Q(1,d)
            assert sum(r for u,r in rates.items() if u[1]==L and virgin not in u[2])==Q(3,2)+Q(1,d)
    # Fixed-degree paths with competing gain and kill: threshold occupation
    # is the first passage weight multiplied by the post-passage occupation.
    ds=range(3,9);s=Q(7,10);a=Q(1,3)
    g={d:Q(d-2) for d in ds};q={d:s+Q(3,2)+Q(1,d) for d in ds}
    reach={3:Q(1)}
    for d in ds:reach[d+1]=reach[d]*g[d]/(g[d]+q[d])
    occ={d:reach[d]/(g[d]+q[d]) for d in ds}
    original=sum(occ[d]*(a-Q(1,d)) for d in ds)
    assert original==sum(sum(occ[d] for d in ds if d>=m)/((m-1)*m) for m in range(4,9))
    for m in ds:
        assert sum(occ[d] for d in ds if d>=m)==reach[m]*sum(occ[d]/reach[m] for d in ds if d>=m)
    print('8 independent whole-fringe rate checks and 7 exact threshold identities passed',flush=True)


def direct_hierarchy(s,H,n,include_copy=False):
    """Arbitrary leading direct productions, then a disjoint terminal class."""
    labels=tuple(range(2,21))
    terminal={j:Q(1,2)/(s+Q(1,j)+Q(3,2)) for j in labels}
    def advance(w,j,sgn):
        t=Tables(s,H,sgn,mark=j)
        return t.parent(lambda d,k: w[min(d,20)] if d>2 else 0,w[20])
    ones=dict.fromkeys(labels,Q(1))
    alpha=advance(ones,20,-1)  # exact intrinsic process is monotone in root killing
    assert 0<alpha<1
    sources=[terminal,terminal]
    if include_copy:
        sources=[{j:terminal[j]+Tables(s,H,sgn,mark=j).classes()[0]
                  for j in labels} for sgn in (1,-1)]
    v=sources
    for _ in range(n):
        v=[{j:b[j]+advance(w,j,sgn) for j in labels} for sgn,w,b in zip((1,-1),v,sources)]
    # A posteriori order-preserving residual: bound the true vector above by
    # w + epsilon*1, where epsilon pays the entire next-generation residual.
    w=v[1]
    err=max(max(sources[1][j]+advance(w,j,-1)-w[j],0) for j in labels)/(1-alpha)
    assert v[0][2]<=w[2]+err<1
    print('full direct-reset hierarchy, includes virgin copy',include_copy,n,'alpha',float(alpha),
          'partial lo/hi',float(v[0][2]),float(w[2]),
          'all-generation upper',float(w[2]+err),'remainder budget',float(err),flush=True)
    return v[0][2],w[2]+err


if __name__=='__main__':
    p=ArgumentParser(description=__doc__)
    p.add_argument('--degree',type=int,default=48)
    p.add_argument('--hierarchy',type=int,default=0)
    p.add_argument('--extend-copy',action='store_true')
    p.add_argument('--broad',nargs='*',type=int,default=[])
    args=p.parse_args();H=args.degree;assert H>=21
    tests();physical_tests()
    for s in (LO,HI):
        plain=[];mem=[]
        for side in (1,-1):
            plain.append(Tables(s,H,side,False).classes())
            mem.append(Tables(s,H,side).classes())
        assert plain[0][0]<=mem[0][0]<=mem[1][0]<=plain[1][0]
        assert plain[0][1]==mem[0][1] and plain[1][1]==mem[1][1]
        immediate=Q(1,2)/(s+2)
        print('trial',str(s),'H',H,'immediate',float(immediate),flush=True)
        print('first copy raw',*[float(q[0]) for q in plain],
              'memory',*[float(q[0]) for q in mem],flush=True)
        print('one newborn, no copy',*[float(q[1]) for q in mem],flush=True)
        subtotal=[immediate+sum(q) for q in mem]
        print('full three-class subtotal',*[float(v) for v in subtotal],flush=True)
        assert 0<subtotal[0]<subtotal[1]<1
        if s==LO:
            print('critical complement minimum',float(1-subtotal[1]),flush=True)
            for k in args.broad:
                assert k>=2
                b=[Tables(s,H,side,min_children=k).classes()[0] for side in (1,-1)]
                assert 0<b[0]<=mem[0][0] and b[0]<b[1]<=mem[1][0]
                print('one virgin copy with crown k >=',k,*[float(z) for z in b],flush=True)
            if args.hierarchy:
                l,u=direct_hierarchy(s,24,args.hierarchy)
                print('all no-copy plus one-copy virgin class',float(l+mem[0][0]),float(u+mem[1][0]),
                      'critical complement >=',float(1-u-mem[1][0]),flush=True)
                if args.extend_copy:
                    l,u=direct_hierarchy(s,24,args.hierarchy,True)
                    print('critical complement after arbitrary resets then virgin copy >=',float(1-u),flush=True)
