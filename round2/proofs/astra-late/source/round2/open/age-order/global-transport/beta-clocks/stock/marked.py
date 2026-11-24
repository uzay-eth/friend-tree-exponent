"""Certify closed birth-scale boxes with the true siblings of each pending tip."""

from fractions import Fraction as Q
from functools import cache
from pathlib import Path
from runpy import run_path

HERE=Path(__file__).resolve().parent
A=run_path(str(HERE/'adverse.py'))
P=A['POS']
L=run_path(str(HERE/'sites.py'))
BOXES=[
    ('.50','.55',(1140,992,610,1945)),
    ('.55','.60',(1322,902,602,2100)),
    ('.60','.65',(1487,821,821,2014)),
    ('.65','.70',(1713,899,899,1831)),
    ('.70','.75',(1961,1052,801,1710)),
    ('.75','.80',(2150,1114,696,1732)),
    ('.80','.85',(2258,1192,589,1655)),
    ('.85','.90',(2526,1299,502,1729)),
    ('.90','1',(2047,1080,693,1661)),
    ('1','1.05',(2604,1352,381,1632)),
    ('1.05','1.15',(2370,1233,453,1590)),
    ('1.15','1.3',(2251,1169,500,1582)),
]


def check(lo,hi,tariff):
    middle,fav,host=P['prices'](lo,hi)
    prices=A['prices'](lo,hi)
    fallback=P['BASE']['price_caps'](lo)
    high=A['high_constant'](lo)
    @cache
    def upper(k,stock):
        u=fallback(k,lo)
        if k>=32:
            u=min(u,high/(k-1))
        return min(u,prices.get((k,stock),u),prices.get((k,0),u))
    cons=L['conditions'](middle,fav,host,upper,lo,high)
    tariffs=tuple(Q(i,100000) for i in tariff)
    assert L['certify'](cons,tariffs), (lo,hi)
    slack=min(sum(i*w for i,w in zip(row,tariffs))-rhs for row,rhs in cons.items())
    print(f'PASS [{lo},{hi}] for every m >= 21: M >= {middle}, U(2,0) <= {upper(2,0)}, minimum tariff slack {slack}',flush=True)


if __name__=='__main__':
    for i,(a,b,tariff) in enumerate(BOXES):
        assert i==0 or Q(a)==Q(BOXES[i-1][1])
        check(Q(a),Q(b),tariff)
    assert Q(BOXES[0][0])==Q(1,2) and Q(BOXES[-1][1])==Q(13,10)
    print('PASS closed [1/2,13/10] for all m >= 21; strict order for y >= 1/2 after the recorded splice.')
