"""Independent exact contraction of charged inner–inner edge rows."""
from fractions import Fraction as Q
from itertools import product
import gzip
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent


def rows(data, omit_source=False):
    rates=data['Sr']; norm=data['S']; sign=data['sign']; lam=Q(data['lam'])
    out={}
    for a,b,h,x,y,i,v,pg,wg,pb,wb,hC,hV,dfp,dfw,rpp,rpc,rwp,rwc in data['rows']:
        p=(a,b,h); w=(x,y,i); d=sum(p); e=sum(w)
        assert 2<=d<18 and 2<=e<18 and b>=1 and h in (0,1) and i==1
        # own gain, included opposite vertex donor, centered external donors;
        # plus incoming transfers from the first-body rows of both endpoints.
        gain_p=(a+Q(1,e)+Q(h*rpp+(b-1)*rpc,rates))*(pg-v)
        gain_w=(x+Q(1,d)+Q(y*rwc,rates))*(wg-v)
        born=Q(a,d)*(pb+hC-v)+Q(x,e)*(wb+hV-v)
        source=(Q(1,e)-Q(rpc,rates))*dfp+(Q(1,d)-Q(rwp,rates))*dfw
        r=sign*(gain_p+gain_w+born+(0 if omit_source else source)-lam*v)/norm
        assert (omit_source or r>=0) and (p,w) not in out
        out[p,w]=r
    ps=[(d-k-h,k,h) for h in (0,1) for d in range(2,18) for k in range(1,d+1-h)]
    ws=[(d-k-1,k,1) for d in range(2,18) for k in range(d)]
    assert set(out)==set(product(ps,ws))
    return out


def main():
    with gzip.open(HERE/'pinned-fg-slabs.json.gz','rt') as f: pack=json.load(f)
    for mode in ('lower','upper'):
        out=rows(pack[mode]['g'])
        wrong=rows(pack[mode]['g'],True)
        assert min(wrong.values())<0
        print(mode,'planted omitted f→g donor source',sum(v<0 for v in wrong.values()),'negative rows')
        print(mode,'edge count',len(out),'full range',min(out.values()),max(out.values()))
        for title,accept in [
            ('all inner↦inner',lambda p,w: True),
            ('nonroot inner↦newborn',lambda p,w:p[2]==1 and w==(1,0,1)),
            ('nonroot degree<=4 ↦ newborn',lambda p,w:p[2]==1 and sum(p)<=4 and w==(1,0,1)),
            ('nonroot degree<=7 ↦ degree<=7',lambda p,w:p[2]==1 and sum(p)<=7 and sum(w)<=7),
            ('nonroot inner↦childless inner',lambda p,w:p[2]==1 and w[1]==0)]:
            rr=[(v,p,w) for (p,w),v in out.items() if accept(p,w)]
            mn,mp,mw=min(rr); mx,xp,xw=max(rr)
            print(f'{title}: n={len(rr)} min={mn} ({float(mn):.13g}) at {mp}->{mw}; max={mx} ({float(mx):.13g}) at {xp}->{xw}')
        for edge in [((0,1,1),(1,0,1)),((1,1,1),(1,0,1)),((2,1,1),(1,0,1)),
                     ((1,2,1),(1,0,1)),((1,1,1),(2,0,1))]:
            print('attainable',edge,out[edge],f'{float(out[edge]):.13g}')


if __name__=='__main__': main()
