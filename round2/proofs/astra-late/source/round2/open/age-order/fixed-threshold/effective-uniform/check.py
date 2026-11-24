"""Replay rational and symbolic dyadic cutoff budgets; no large integers/trees."""
from fractions import Fraction as Q
from pathlib import Path
import re


A=2**119
LOG_SPECTRAL_FLOOR=348+(515-1).bit_length()


def strip_pow2(expr):
    assert expr.startswith('2^{') and expr.endswith('}'), expr
    return expr[3:-1]


def check_height(expr):
    exponent=strip_pow2(strip_pow2(expr))
    if exponent.isdecimal():
        log_log_n=int(exponent)
    else:
        power=strip_pow2(exponent)
        assert power.isdecimal() and int(power)<=A.bit_length()+4, exponent
        log_log_n=2**int(power)
    assert log_log_n>=LOG_SPECTRAL_FLOOR, 'cutoff fails spectral floor'
    assert log_log_n==16*A, 'cutoff differs from the budgeted height'
    return log_log_n


def stated_height():
    readme=Path(__file__).with_name('README.md').read_text()
    (statement,)=re.findall(r'^n\\ge (\S+?)\\quad\\Longrightarrow', readme, re.M)
    log_log_n=check_height(statement)
    try:
        check_height('2^{2^{123}}')
    except AssertionError as err:
        assert str(err)=='cutoff fails spectral floor'
    else:
        raise AssertionError('smaller tower was accepted')
    return log_log_n


def large_window():
    H=512
    # Dyadic variance sum: 2^-1/2 < 3/4. Bernoulli second moments
    # use the 7/8 envelope for 2^-1/4, including the initial block.
    r=Q(7,8)
    env=6+r*(1+r)/(1-r)**3+10*r/(1-r)**2+25/(1-r)
    assert env < 2**11
    assert Q(1,8)+Q(9+2**11,H**2)<Q(1,7)
    assert Q(2*192,191)*Q(22,7)**2/6<4
    # First-block and geometric-tail sums used in the log estimate.
    r=Q(15,16)
    assert r**8>Q(1,2)  # 2^-1/8 < r
    assert 3*2*256+3*2*4/(1-r)<2**13
    assert 16+4/(1-Q(3,4))<=2**7
    assert 5*6*128+5*6*2/(1-Q(3,4))<2**14
    assert Q(2**13*H**2,2**39)+Q((2**7+2**14)*H**2,2**78)<1
    assert Q(6,7*3**6)>Q(1,2**10)
    assert 24*2**12 < 2**17  # K=2^(6*17)
    assert Q(113,200)-Q(1,2**21)-2*Q(1,128)>Q(1,2)


def compact_window(log_S0):
    B=A+100
    # In units of delta: retain enough inward slack even after rounding,
    # and remain inside the 16-delta stationary interval on both sides.
    center, rounding, size_shift, martingale=5,2,2,1
    assert center-rounding >= size_shift+martingale
    assert 2*(center+rounding+martingale)==16
    assert 4-rounding < size_shift+martingale
    # Bounds at exponent level: T+h+B < 2^(3B);
    # t < 2^(3B+6), S < 2^(3B+13).
    assert (256*B).bit_length()<B-4
    t_exp=3*B+6
    S_exp=3*B+13
    assert S_exp+4<5*B
    assert t_exp+2-8*B<-B  # core cap probability
    # Since t>=64(T+h+B), e^-t/y0 <= 2^-60(T+h+B).
    assert 32*B+3-60*B<-B
    assert 8*B+4-(16*B-1)<-B  # gamma density window
    # Stationary exclusion: log2 outside multiplier < 2^(6B),
    # and both negative tail exponents dominate this bound (use a
    # sharper multiplier upper bound 2^(5B+ceil(log2(6B))+1)).
    mult_exp=5*B+(6*B).bit_length()+1
    for negative_exp in (7*B-3,6*B-3):
        assert mult_exp+1<negative_exp and (B+16).bit_length()<negative_exp-1
    # log2 queried degree sum is O(B); the factorial error pays >= R-5.
    assert (48*B+50).bit_length()<5*B-1
    assert -B+12 < -A-4
    # S0=log2 N=2^(16A); use the checked log2 S0 without expanding S0.
    assert 348+10<log_S0
    assert 348+5*B+10<log_S0  # first spectral + radius requirement
    assert 12*B+12<log_S0    # last exponent using product of k and E
    # Since delta=2^-4t and K=2^102, log2 of the required n is
    # at most max(106+4t, 206+8t+B) < 2^(t_exp+4).
    assert (206+B).bit_length()<t_exp
    assert t_exp+4<log_S0


if __name__=='__main__':
    large_window()
    print('PASS rational likelihood, stationary flux, and beta separation; K=2^102')
    log_S0=stated_height()
    print('PASS symbolic height: N=2^(2^(2^123)); smaller tower rejected by the spectral floor')
    compact_window(log_S0)
    print('PASS symbolic cylinder and de-Poissonization budgets; A=2^119, log2(log2(N))=2^123')
