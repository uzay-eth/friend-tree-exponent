"""Full tree potential jumps and an independent local birth allocation."""
from collections import defaultdict
from charges import embedding


def verify(E,normal,c,cs,s,ev,vs,result):
    from check import tree, reach, drift
    Q,Z,J=E.Q,E.ZERO,E.J
    j=tree(1);b=tree(0,j);old=tree(4,b);fork=tree(0,b,b);big=tree(4,old,old,b)
    fixtures=list(reach(8))+[tree(2,tree(l,old,big,fork),tree(l,b,b)) for l in (4,12,15,16,17,19,24,32,40,64)]+[
        tree(2,tree(8,tree(0,tree(2,j)),tree(6,tree(2,j,j)))),tree(2,tree(11,tree(0,b),b,fork,old)),
        tree(1,tree(62,b),tree(65,old)),tree(2,tree(1,*([j]*19)))]
    margin,_,cv,tail=result;delta=abs(s-c.c.lam)
    def st(y):return sum(y),y[1]
    def census(y):
        d,k=st(y)
        return cv.get((d,k),tail[min(k,c.c.KF+1)])-cv[J]
    def occ(p,x):
        for f,ids in embedding(p):
            if ids[0]==0:continue
            if f=='g':tab=ev; ff='g'
            elif f in vs:tab=vs[f];ff=f
            else:continue
            yield ff,ids,tab.get(normal(ff,tuple(st(x[i]) for i in ids)),Z)
    def value(p,x):return -delta*cs.M/s+sum(census(y) for y in x[1:])+sum(v for _,_,v in occ(p,x))
    birth_count=defaultdict(int);min_slack=None
    for tr in fixtures:
        p,x=E.shape(tr,0);ch=[[j for j,k in enumerate(p) if k==i] for i in range(len(p))]
        a=value(p,x);whole=local=-s*a
        old={(f,ids):w for f,ids,w in occ(p,x)}
        for i,y in enumerate(x):
            l,k,h=y;b=Q(l,sum(y));g=l+sum(Q(1,sum(x[j])) for j in ch[i]+([p[i]] if i else []))
            up=E.label(x,i,E.up(y));whole+=g*(value(p,up)-a)
            local+=g*(sum(w-old[f,ids] for f,ids,w in occ(p,up))+(census(up[i])-census(y) if i else Z))
            if not b:continue
            post=E.label(x,i,E.bud(y))+(E.NB,);n=len(x);whole+=b*(value(p+(i,),post)-a)
            if i:local+=b*(census(post[i])-census(y))
            for f,ids,w in occ(p+(i,),post):
                if n not in ids:local+=b*(w-old[f,ids]);continue
                if not w:continue
                ancestor=tuple(j for j in ids if j!=n)
                if f=='g':assert ancestor==(i,)
                else:assert (('g' if f=='hV' else 'hV'),ancestor) in old
                birth_count[f]+=1;local+=b*w
        assert whole==local
        slack=c.sign*drift(c.c,tr)-delta*c.c.potential(ch,x)/c.c.S+whole-margin*(len(x)-1)
        assert slack>=0,(p,x,slack)
        min_slack=slack if min_slack is None else min(min_slack,slack)
    assert birth_count['g'] and birth_count['hV']
    print('exact full-tree fixtures',len(fixtures),'nonzero source matches',dict(birth_count),'minimum surplus',float(min_slack),flush=True)
