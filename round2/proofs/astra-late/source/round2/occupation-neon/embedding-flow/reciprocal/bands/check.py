"""Exact point-atom replay of slab balances against physical jumps."""
from collections import defaultdict
from fractions import Fraction as Q
from math import prod

from bands import BandModel, Cell, basis


def make(local=True,threshold=3):
    return BandModel(threshold,1,1,inner=threshold-1,wide_roles=(0,),chain_arity=(1,1,0),local=local)


def atom(m,cells,ds,pressure, *, seam=True):
    coords=[]
    for c,d in zip(cells,ds):
        if c.degree:
            assert c.degree==d
            coords.append((Q(1),)); continue
        lo,hi=m.boundaries[c.slab:c.slab+2]
        y=Q(m.threshold,d)
        coords.append(tuple(w/(1+Q(1,d))**m.m for w in basis(m.n,(y-lo)/(hi-lo))))
    v={}
    for kind in ('mass',*m.ports(cells)):
        for ix,cols in m.blocks[cells,kind].items():
            v[cols[0]]=prod(a[j] for a,j in zip(coords,ix))*(1 if kind=='mass' else pressure[kind])
    rates=[]
    for i,(c,d) in enumerate(zip(cells,ds)):
        g=d-c.arity-1+sum((Q(1,ds[j]) for j in (i-1,i+1) if 0<=j<len(ds)),Q(0))
        g+=sum(z for (r,_),z in pressure.items() if r==i)
        rates.append(g)
        if (cells,i) in m.seams and d==17 and seam:
            coords0=coords[:i]+coords[i+1:]
            for ix,(a,b) in m.seams[cells,i].items():
                v[a]=prod(axis[j] for axis,j in zip(coords0,ix)); v[b]=g*v[a]
    return v,rates


def point(m,ds,ks):
    return tuple(Cell(d,k) if d<m.threshold else Cell(0,k,int(d<18)) for d,k in zip(ds,ks))


def verify(m,ds,ks):
    cells=point(m,ds,ks)
    p={k:Q(m.slots(cells,k),5) for k in m.ports(cells)}
    v,g=atom(m,cells,ds,p)
    expected=defaultdict(Q)
    def add(s,d,amount):
        for e in m.tests(tuple(c.label for c in s)):
            if (s,e) in m.local_rows:
                f=prod(Q(m.threshold,dj)**a if not c.degree else 1 for c,dj,a in zip(s,d,e))
                expected[s,e]+=amount*f
    add(cells,ds,m.low)
    for i,d in enumerate(ds):
        gain=ds[:i]+(d+1,)+ds[i+1:]
        add(cells,ds,g[i]); add(point(m,gain,ks),gain,-g[i])
        b=Q(d-ks[i]-1,d)
        promoted=ks[:i]+(ks[i]+1,)+ks[i+1:]
        add(cells,ds,b)
        add(point(m,ds,promoted),ds,-b)
        if len(ds)<3 and i==len(ds)-1:
            add(point(m,ds+(2,),promoted+(0,)),ds+(2,),-b)
    actual={key:sum(c*v.get(j,0) for j,c in row.items()) for key,row in m.local_rows.items()}
    assert {k:x for k,x in actual.items() if x}=={k:x for k,x in expected.items() if x},(ds,ks)
    for row,rhs in m.upper_rows[m.original_counts[2]:]:
        assert sum(c*v.get(j,0) for j,c in row.items())<=rhs,(ds,ks)
    return sum(bool(x) for x in actual.values())


def main():
    m=make()
    cases=[((2,),(0,)),((3,),(0,)),((16,),(1,)),((17,),(0,)),((18,),(1,)),
           ((17,2),(1,0)),((17,17),(1,0)),((18,17),(1,0)),((17,2),(1,1)),
           ((2,2,2),(1,1,0)),((16,2,2),(1,1,0)),((17,2,2),(1,1,0)),((18,2,2),(1,1,0))]
    count=sum(verify(m,ds,ks) for ds,ks in cases)
    # Remove true D=17 atom but keep its source mass: the discrete cut
    # must fail even though the simple alpha domination still passes.
    cells=point(m,(17,),(0,))
    p={(0,'mother'):Q(1,5)}
    v,_=atom(m,cells,(17,),p,seam=False)
    assert any(sum(c*v.get(j,0) for j,c in row.items())>rhs
               for row,rhs in m.upper_rows[m.original_counts[2]:])
    # Check all global balances are exactly sums of locals on this cut.
    sums=defaultdict(lambda:defaultdict(Q))
    for (cells,e),row in m.local_rows.items():
        r=sums[tuple(c.label for c in cells),e]
        for j,v in row.items(): r[j]+=v
    clean=lambda row:{j:z for j,z in row.items() if z}
    assert sums.keys()==m.balance.keys()
    assert all(clean(row)==clean(sums[k]) for k,row in m.balance.items())
    m4=make(threshold=4)
    extra=[((3,2),(1,0)),((17,3,2),(1,1,0)),((18,3,2),(1,1,0))]
    extra_rw=sum(verify(m4,d,k) for d,k in extra)
    # Reach both arity boundaries with a root-only side check; do not
    # build the enormous high-arity tensor family as a test fixture.
    single=BandModel(3,17,0,inner=2,wide_roles=(0,),chain_arity=(1,1,0),
                     project=False,port_transport=False,local=False)
    c16=(Cell(0,16,1),); c17=(Cell(0,17,1),)
    assert single.impossible(c17[0]) and not single.impossible(c16[0])
    start=len(single.equalities)
    single.add_seam(c16,0)
    p16={(0,'mother'):Q(1,3),(0,'children'):Q(4)}
    v16,_=atom(single,c16,(17,),p16)
    assert len(single.equalities)-start==single.n+1
    assert all(sum(c*v16.get(j,0) for j,c in row.items())==rhs for row,rhs in single.equalities[start:])
    print('literal physical atoms',len(cases),'nonzero local coefficients',count,
          'global rows recovered',len(sums),'planted missing seam rejected',
          'H4 atoms',len(extra),'H4 coefficients',extra_rw,'k16 atom/k17 impossible checked')


if __name__=='__main__':
    main()
