"""Upper bound on the homogeneous lower comparison, used to audit its lower signs."""
from fractions import Fraction as F
from math import ceil

from sector import BETA, Sector, UNIT, hy


if __name__=='__main__':
    s=F('.565615')
    g=F(2,3),F(1)
    # The target comparison omits integer rounding in its lower value tables.
    # Upward rounding here bounds that comparison from above. Tail moments
    # and the analytic direct barriers are kept at their original values.
    hy.rows.floor=ceil
    hy.tc.floor=ceil
    model=Sector(s,4,0,7)
    for a in (F(1,2),BETA):
        tl,_,_,mu=hy.ex.moments(s+a,hy.tc.STOP)
        assert (s+a+1)*tl>=2*mu
    points,C=model.payments(g)
    for tag,threshold in zip((2,16),g):
        _,upper,tail=model.occupation(tag,27)
        value=sum(ceil(m*point[0][1]) for m,point in zip(upper,points))/F(UNIT)
        value+=C*tail
        assert value<threshold
        print('unrounded homogeneous lower comparison, test (2/3,1), root',F(1,tag),
              'upper',ceil(value*10**12),'/ 1e12')
