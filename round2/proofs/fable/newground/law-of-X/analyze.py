import sys
import numpy as np
from scipy import stats
from scipy.special import gamma

MU = .56572366
n1, c1, n2, c2, y = np.loadtxt(sys.argv[1], delimiter=',', skiprows=1).T
x1, x = c1 / n1**MU, c2 / n2**MU
z = x*y**MU
N = len(x)
rng = np.random.default_rng(789)
q = [0, .001, .01, .025, .05, .1, .25, .5, .75, .9, .95, .975, .99, .999, 1]

def desc(a):
    return np.array([a.mean(), a.var(ddof=1), stats.skew(a, bias=False), *np.quantile(a,q)])

def fit(a, grid, side='lower'):
    aa = np.sort(a)
    n = np.searchsorted(aa, grid, side='right')
    p = (n if side == 'lower' else (len(a)-n))/len(a)
    if np.any(p == 0):
        return np.full(3 if side == 'upper' else 1, np.nan)
    if side == 'lower':
        return np.array([stats.linregress(np.log(grid), np.log(p)).slope])
    return np.array([stats.linregress(np.log(grid), np.log(-np.log(p))).slope,
                     -stats.linregress(grid, np.log(p)).slope,
                     -stats.linregress(np.log(grid), np.log(p)).slope])

lower = np.array([.02,.03,.04,.06,.08,.1,.15,.2])
upper = np.array([1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8])

def compact(a):
    return np.r_[desc(a),fit(a,lower),fit(a,lower[3:]),fit(a,upper,'upper'),fit(a,upper[3:],'upper')]

print('N', N, 'sizes n1/n2: 1%, 50%, 99%', np.quantile(n1,[.01,.5,.99]), np.quantile(n2,[.01,.5,.99]))
print('n1 stats',desc(x1))
print('n2 stats',desc(x))
print('paired diff: mean, SE(mean), sd, RMSE, |delta| quantiles:',(x-x1).mean(), (x-x1).std(ddof=1)/N**.5,
      (x-x1).std(ddof=1), np.mean((x-x1)**2)**.5, np.quantile(abs(x-x1),[.5,.9,.99]))
print('Y mean,var,KS:',y.mean(),y.var(ddof=1),stats.kstest(y,'gamma',args=(2,)))
print('Z mean,var:',z.mean(),z.var(ddof=1),'independence predictions:', x.mean()*gamma(2+MU),
      np.mean(x*x)*gamma(2+2*MU)-(x.mean()*gamma(2+MU))**2)
for a,b in ((x1,y),(x,y),(z,y),(x1,x)):
    pr=stats.pearsonr(a,b)
    print('Pearson,CI; Spearman,p',pr.statistic,pr.confidence_interval(),stats.spearmanr(a,b))
print('quintiles of X by Y quintiles (N, mean, SE)')
for ix in np.array_split(np.argsort(y),5): print(len(ix), x[ix].mean(),x[ix].std(ddof=1)/len(ix)**.5)
bins=np.r_[np.arange(0,2.01,.1),np.inf]
counts,_=np.histogram(x,bins)
print('hist bin, counts, ascii bar (one * per 100)')
for i,v in enumerate(counts): print(f'{bins[i]:.1f}-{bins[i+1]:.1f}: {v:5} '+ '*'*round(v/100))
print('tail threshold, count early, count late, exact 95% binomial CI late')
for grid,op in ((lower,np.less_equal),(upper,np.greater)):
    for b in grid:
        k=int(op(x,b).sum()); print(b,op(x1,b).sum(),k,stats.binomtest(k,N).proportion_ci())
print('fits: lower all, lower .06+; upper Weibull/exponential/power slope all, then 1.4+')
print('early',compact(x1)[18:]); print('late',compact(x)[18:])
boots=[]
for j in range(500):
    ix=rng.integers(0,N,N)
    boots.append(compact(x[ix]))
print('late bootstrap 95% CI (rows: low, high)',np.nanquantile(boots,[.025,.975],axis=0))
