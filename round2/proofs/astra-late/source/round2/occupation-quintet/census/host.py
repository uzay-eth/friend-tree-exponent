#!/usr/bin/env python3
"""Sparse host eligibility and g-only comparison using authenticated caps."""
from fractions import Fraction as Q
from functools import lru_cache
import json
from math import comb, ceil, floor
from pathlib import Path

GOLD = {'lower': (0,401,18785,99882299906,1808571640916),
        'upper': (-1053,946,19890,155236771435,1808494729937)}


class HostCaps:
    def __init__(self,s,receipt):
        from census import PINS
        if not s.authenticated: raise ValueError('authenticate the snapshot first')
        doc=json.loads(Path(receipt).read_bytes())[s.side]; v=doc['verify']
        if (v['sha256'],v['mode'],Q(v['P'],v['Q']),v['S'])!=(PINS[s.side][2],s.header['mode'],s.lam,s.scale):
            raise ValueError('receipt does not match witness')
        capf,capg,high,mn,mx=GOLD[s.side]
        assert (doc['cap_f'],v['caps']['g'],v['higher_per_core'],doc['m_min'],doc['m_max'])==(capf,capg,high,mn,mx)
        # Pins from the verified golden report. Using independently reported
        # maxima across a family may be looser than its ordered-block maxima.
        expected={'lower':(225,42,103,26,31,47),'upper':(408,41,139,27,32,45)}
        family=tuple(v['caps'][f] for f in ('h3','p4','s3','p5','b5','s4'))
        assert family==expected[s.side]
        self.s=s; self.family=family; self.high=high; self.comp=v['comp_fin']; self.alpha=v['alpha_ticks']
        self.min_weight=mn; self.max_weight=mx

    def multiplicities(self,d,c):
        h=self.s.header
        h3=comb(c,2) if d<h['d3m'] and 2<=c<=h['K3m'] else 0
        p4=0
        if d<h['d4i'] and 2<=c<=h['K4i']:
            near=min(h['d4i']-1,h['SP4']-d,h['K4i'])
            if near>=2:p4=comb(c,2)*(near-1)
        s3=comb(c,3) if d<h['dS3'] and 3<=c<=h['KS3'] else 0
        p5=0
        if d<h['d5i'] and 2<=c<=h['K5i']:
            # For fixed sum the product is maximal at balanced arities,
            # until both hit the per-slot degree/arity cap.
            n=min(h['d5i']-1,h['K5i']); s=min(2*n,h['SP5']-d)
            if s>=4:p5=comb(c,2)*(s//2-1)*((s+1)//2-1)
        b5=0
        if d<h['d5o'] and 3<=c<=h['K5o']:
            near=min(h['d5c']-1,h['K5c'],h['SB5']-d)
            if near>=2:b5=c*comb(c-1,2)*(near-1)
        s4=comb(c,4) if d<h['dS4'] and 4<=c<=h['KS4'] else 0
        return h3,p4,s3,p5,b5,s4

    def higher_ticks(self,d,c):
        """Valid hostwise outward cap; exact zero when no family is possible."""
        return min(self.high,sum(n*max(0,v) for n,v in zip(self.multiplicities(d,c),self.family)))

    @lru_cache(maxsize=512)
    def edge_column(self,i):
        values=[self.s.g(i,j) for j in range(self.s.inf+1)]
        return min(values),max(values)

    def comparison_ticks(self,d,c):
        """Local recharged lower / unrecharged upper ticks (cheap outside higher support)."""
        s=self.s; f=s.f(d,c)
        if not c:return floor(f),ceil(f)
        if any(self.multiplicities(d,c)):
            # A per-core GLOBAL min/max is not generally a local weighted
            # bound for new allocations. Use the bounded table weight instead.
            n=sum(self.multiplicities(d,c)); lo=-n*(16<<40); hi=-lo
        else:lo=hi=0
        if d<s.d2:
            i=s.inner[d,c] if c<=s.K2 else len(s.inner)+d-s.K2-1
            g0,g1=self.edge_column(i)
        else:
            j=s.tail_base+min(c,s.K2+1)-1
            t=Q(s.d2,d); a,b=self.edge_column(j); x,y=self.edge_column(s.inf)
            g0=(1-t)*x+t*a; g1=(1-t)*y+t*b
        lower=floor(f+lo+min(0,g0)+self.alpha*(c-2)); upper=ceil(f+hi+max(0,g1))
        # The two endpoints bound DIFFERENT allocations: the lower alone
        # includes the recharge. They need not be ordered at one host.
        return lower,upper
