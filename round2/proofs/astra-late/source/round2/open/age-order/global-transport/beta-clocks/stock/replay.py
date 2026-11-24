"""Verify the signed sharp-tail bounds at the original finite completion size."""

from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

HERE=Path(__file__).resolve().parent
K=run_path(str(HERE/'check.py'))
L=run_path(str(HERE/'ledger.py'))
BASE=K['BASE']
BOXES=[
    ('13/10','3/2',(1782,990,990,2231)),
    ('3/2','8/5',(2667,1367,786,2086)),
    ('8/5','19/10',(2449,1256,799,1992)),
    ('19/10','23/10',(2767,1410,597,1955)),
]


def main():
    for idx,(lo,hi,tariff) in enumerate(BOXES):
        left,right=Q(lo),Q(hi)
        assert idx==0 or Q(BOXES[idx-1][1])==left
        middle,fav,host=K['prices'](left,right)
        cap=BASE['price_caps'](left)
        constraints=L['prepare'](middle,fav,host,cap,left)
        assert L['OLD']['certify'](constraints,tuple(Q(w,100000) for w in tariff))
        print(f'PASS [{left},{right}], m >= 21. Serial {middle}, old-leaf host {host[2]}, favourable parent-2 {fav[2]}',flush=True)
    assert Q(BOXES[0][0])==Q(13,10) and Q(BOXES[-1][1])==Q(23,10)
    print('PASS strict physical adjacent tail order for b >= 20, (b+1)x >= 13/10 with beta-clock splice')


if __name__=='__main__':
    main()
