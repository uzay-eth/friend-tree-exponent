#!/usr/bin/env python3
"""Gold2-only two-trial frontier, with actual marked births and signed LP bounds."""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from json import dumps, loads
from math import ceil
from mmap import mmap, ACCESS_READ
from pathlib import Path

HERE=Path(__file__).resolve().parent
spec=spec_from_file_location('gold2_pressure',HERE.parent/'replay.py')
pressure=module_from_spec(spec); spec.loader.exec_module(pressure)
life=pressure.life
R2=life.R2
GOLD=R2/'occupation-quintet-gold2'
TARGET=Q(241,426)
STEP=Q(2,10**6)
UNIT=10**6


def pinned(path,digest):
    data=path.read_bytes(); assert sha256(data).hexdigest()==digest,path
    return loads(data)


def inputs():
    receipt=pinned(life.tail.RECEIPT,'5e11f9c0b6eb3bd43e880be608c37d6f0b1fdbcc632b91bb5e8037db52e13284')
    joint=pinned(GOLD/'host/conditioned24.json','3f98385e0662b4a5aea71dcfef03cdbd66012de01f0083456e243592007de790')
    strong=pinned(GOLD/'host/h3-joint/receipt.json','ae6f5cde2ef710b4196c3c445a3786580e5ccc9dfd9206b2962f60628866c6ed')
    rat=life.load('gold2_release_floors',GOLD/'rationals/check.py')
    parts=(rat.normalized_floor(GOLD/'post-fy/threshold-lower.historical.log',
                                '39d0d02575cb05df74776466d787cf8f8ef90399003ee0709e2e7c0e9d5e9191',
                                r'^H3 normalized orientations .* total (\d+/\d+) '),
           rat.normalized_floor(GOLD/'edges/absorb12.log',
                                'eded191d098bb2c9296c5d7d2861a051bc303bce122953c518ebf1981d22d521',
                                r'^TOTAL (\d+/\d+) '),
           rat.normalized_floor(GOLD/'payer/README.md',
                                '0020b89ff3d30dcf654225d7f37626d62c8f21b60312b7b2b342036faac6904f',
                                r'Exact normalized value: `(\d+/\d+)`'))
    booked=sum(parts)
    a=receipt['lower']; lower=Q(a['nominal'])+(booked-a['err_per_core'])/a['m_max']
    upper=Q(receipt['upper']['certified_exact'])
    assert lower==Q(126905728962161948969119589499374909399,224328301684105796116575000000000000000)
    assert Q(99,175)<lower<lower+STEP<TARGET<upper-STEP<upper<Q(142,251)
    # Reuse the already cold-audited credits without repricing any one.
    audit=life.load('gold2_booked_credits',R2/'checks/gold2-host-h3-sol.py')
    credits={side:audit.credits(side) for side in ('lower','upper')}
    assert sum(credits['lower'][1].values())+credits['lower'][2]>=booked
    bounds=tuple(Q(strong['results'][side]['final_endpoint']) for side in ('lower','upper'))
    assert lower<bounds[0]<TARGET<bounds[1]<upper
    return receipt,joint['results'],bounds,credits


def price(lp,lo,hi,fn,up):
    sign=-1 if up else 1
    z=sign*UNIT*lp.lower(lambda k,x:sign*fn(x)/UNIT if k=='occ' else 0)
    return (max if up else min)(lo*z,hi*z)


def tagged_births(lp,limit,lo,hi):
    """Insert upper and lower physical next-generation degree flux constraints."""
    grid=10**12
    # Incoming mark is the degree of the real mother at the focal birth.
    # Her live degree can grow, so her send is in [0,1/source] thereafter.
    ports={}
    for source in range(2,limit+2):
        for target in range(2,limit+1):
            for up,disc in ((False,hi),(True,lo)):
                W={}
                for d in range(target,1,-1):
                    for c in range(d,0,-1):
                        L=d-c; b=Q(L,d)
                        A=b*((d==target)+W.get((d,c+1),0)); g=W.get((d+1,c),0)
                        h=(max if up else min)((A+(L+r)*g)/(disc+b+L+r)
                                               for r in (0,Q(1,source)+Q(c-1,2)))
                        W[d,c]=Q((ceil if up else lambda x:x//1)(h*grid),grid)
                ports[source,target,up]=W[2,1]
    tags=[(d,up) for d in range(2,limit+1) for up in (False,True)]
    assert all(0<=ports[s,d,False]<=ports[s,d,True]<=1/lo
               for s in range(2,limit+2) for d in range(2,limit+1))
    for tag in tags:
        lp.rows['tag',tag]=len(lp.rows)
    for i,(kind,col,x) in enumerate(lp.cols):
        if kind!='occ' or not x[2]: continue
        source=min(x[0],limit+1); b=x[2]
        lp.cols[i]=(kind,col+tuple((lp.rows['tag',(d,up)],z) for d,up in tags
                                   if (z:=b*((source==d)-ports[source,d,up]))),x)
    for tag in tags:
        # The exact dual repair has a 1/LO bound per kind, so each
        # inequality gets its own slack kind. Incoming birth mass is one.
        lp.col(('birth-tag-slack',tag),{('tag',tag):1 if tag[1] else -1})
    return ports


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--side',choices=('lower','upper'),required=True)
    p.add_argument('--bound',choices=('floor','ceiling'),required=True)
    p.add_argument('--mother',type=int,default=8)
    p.add_argument('--live-cut',type=int,default=700)
    p.add_argument('--tip',type=int,default=3)
    p.add_argument('--donors',type=int,default=2)
    p.add_argument('--live-arity',type=int,default=12)
    p.add_argument('--tag-limit',type=int,default=12)
    p.add_argument('--receipt',action='store_true',help='write exact output in this frontier directory')
    a=p.parse_args(); M,H,T,J,K=a.mother,a.live_cut,a.tip,a.donors,a.live_arity
    assert 2<=M<=24 and 2<=H<3999 and 2<=T<=4 and 1<=J<=2 and 1<=K<=30
    up=a.bound=='ceiling'; side=a.side; inward=1 if side=='lower' else -1
    dec=pressure.dec; receipt,rows,bounds,credits=inputs(); rec=receipt[side]; row=rows[side]
    base=bounds[side=='upper']; oldf,edges,h3=credits[side]; B=h3+sum(edges.values())
    assert Q(row['f_floor'])==oldf and not row['reduced_f']
    lo,hi=map(Q,(receipt['lower']['certified_exact'],receipt['upper']['certified_exact']))
    assert (row['payload_sha256'],row['nominal'],row['cap_error'],tuple(map(Q,row['prior'])))==(
        receipt[f'{side}_sha256'], rec['nominal'],rec['err_per_core'],(lo,hi))
    assert (life.tail.BASE/'hashes.log').read_text().splitlines()==['quintet_lower.bin: OK','quintet_upper.bin: OK']
    with (life.tail.BASE/f'quintet_{side}.bin').open('rb') as fp,mmap(fp.fileno(),0,access=ACCESS_READ) as raw:
        assert len(raw)==life.tail.SIZE
        t=life.tail.Snapshot(raw,rec); assert t.s==-inward and t.lam==Q(rec['nominal'])
        W={}; C=70 if up else 30
        _,n=life.price(t,3999 if up else 7999,C,lo if up else hi,rec['cap_f'],up,(H+1,W))
        roots,count,corners=pressure.live(t,rec['cap_f'],W,H,K,M,T,J,lo if up else hi,up)
    overflow=Q(ceil(10**9*130000*2**40/(1+lo)**C),10**9) if up else Q(0)
    model=life.load('gold2_frontier_lp',R2/'occupation-neon/pressure-memory/birth-law/replay.py')
    model.LO,model.HI=lo,hi; model.SAT=row['saturation']; lp=model.build(row['limit'],sided=True)
    cells={tuple(map(int,k.split(','))):tuple(map(int,v)) for k,v in row['cells'].items()}
    assert M<=row['limit']
    def mark(x):
        d,k,b,*_=x
        return b*roots.get((d,k+2),roots[0,0])
    plain=price(lp,lo,hi,mark,up)+overflow
    if a.tag_limit:
        assert 2<=a.tag_limit<=row['limit']
        lp=model.build(row['limit'],sided=True)
        tagged_births(lp,a.tag_limit,lo,hi)
    val=price(lp,lo,hi,mark,up)+overflow
    # Isolate the effect of reading more birth classes without changing the
    # live state solve: on a pooled class use the worst available root.
    coarse=(max if up else min)(v for (d,h),v in roots.items() if d==0 or d>5)
    coarseval=price(lp,lo,hi,lambda x:mark(x) if x[0]<=5 else x[2]*coarse,up)+overflow
    assert val<=coarseval if up else val>=coarseval
    print(side,a.bound,'birth degree cap',M,'live cut',H,'arity',K,'tip',T,'donors',J,flush=True)
    print('physical Bellman states/corners',count,corners,'broad',n,'sided LP rows/cols',len(lp.rows),len(lp.cols))
    print('original singleton ticks/core exact',val,'directed',dec(val,up),
          'signed gain from tagged generation flux constraints',dec((plain-val) if up else (val-plain),False),
          'pooled identical live replay at birth degree >5',dec(coarseval,up),
          'gain from denser birth marks',dec(abs(val-coarseval),False))
    print('non-singleton booked ticks/core',dec(B,False),'old singleton ticks/core',dec(oldf,False),
          'audited endpoint',dec(base,side=='upper',15),
          'distance to 241/426',dec(abs(TARGET-base),True,15),flush=True)
    E=rec['err_per_core']
    den_hi=Q(row['denominator'][1]); save=Q(row['saved_plus_reduced'])
    candidate=t.lam+inward*(B+max(val,oldf)+save-E)/den_hi
    record=dict(side=side,bound=a.bound,baseline=str(base),prior=[str(lo),str(hi)],
                nominal=str(t.lam),payload_sha256=receipt[f'{side}_sha256'],old_f_ticks=str(oldf),
                singleton_ticks=str(val),non_f_ticks=str(B),h3_ticks=str(h3),
                edge_ticks={k:str(v) for k,v in edges.items()},error_ticks=E,
                higher_cap_ticks=rec['higher_per_core'],higher_saved_floor_ticks=str(save),
                denominator_ticks=row['denominator'],options={k:v for k,v in vars(a).items() if k!='receipt'},
                physical_states=count,physical_corners=corners,rows=len(lp.rows),columns=len(lp.cols),trials={})
    # One finite candidate trial from the positive separate floor; a
    # uniform coupled margin can also be divided by the same M upper.
    for name,s in [('2ppm',base+inward*STEP),('241/426',TARGET)]+([] if up else [('candidate',candidate)]):
        delta=inward*(s-t.lam); assert delta>0
        def objective(x):
            d,k,*_=x; L,U,h=cells[d,k]
            m=max(L,rec['m_min']) if up else min(U,rec['m_max'])
            H=rec['higher_per_core']; assert 0<=h<=H
            return mark(x)-delta*m+(H if up else H-h)
        # Constant booked lower floor stays fixed in the obstruction test:
        # extra reward even inside an old mask requires a replacement proof.
        margin=B-E+price(lp,lo,hi,objective,up)+overflow
        den=Q(row['denominator'][0 if up else 1]); cap=rec['higher_per_core'] if up else Q(row['saved_plus_reduced'])
        g=delta*den+E-B-val-cap
        scalar=E+delta*rec['m_max']-B-val
        trial=record['trials'][name]=dict(point=str(s),shift=str(delta),margin_ticks=str(margin),
                                         residual_ticks=str(-margin),separate_sided_gap_ticks=str(g),
                                         original_m_max_gap_ticks=str(scalar))
        print(name,'trial',str(s),'coupled margin ticks/core',dec(margin,up),
              'residual gap bound',dec(-margin,not up), 'separate sided gap',dec(g,not up),
              'fixed original m_max gap',dec(scalar,not up),flush=True)
        if up and name=='2ppm':
            assert margin<0
            print('ceiling on supported inward move with fixed booked floors, full f and all cap savings',
                  dec(STEP+margin/den_hi,True,15),flush=True)
        if name=='candidate':
            endpoint=s+inward*max(0,margin)/den_hi
            trial['endpoint']=str(endpoint)
            print('certified coupled endpoint',dec(endpoint,side=='upper',15),'inward move',
                  dec(inward*(endpoint-base),False,15),flush=True)
    # Independent sufficient endpoint using just booked release ticks and f.
    if not up:
        imp=inward*(candidate-base)
        print('separate endpoint',dec(candidate,side=='upper',15),
              'inward move',dec(imp,False,15),flush=True)
    if a.receipt:
        assert (M,H,T,J,K,a.tag_limit)==(5,700,3,2,12,12)
        data=(dumps(record,sort_keys=True,indent=2)+'\n').encode()
        path=HERE/f'{side}-{a.bound}.json'; path.write_bytes(data)
        print('receipt',path.name,sha256(data).hexdigest(),flush=True)


if __name__=='__main__':main()
