import numpy as np

mu = 0.56572366
small = np.loadtxt('small.csv', delimiter=',', skiprows=1, dtype=int)
print('small (c4,c5):', np.unique(small, axis=0, return_counts=True))
with open('samples.csv') as data:
    next(data)
    counts = np.loadtxt((row for row in data if row.endswith('\n')), delimiter=',', dtype=int)
x = counts / np.array([250_000., 1_000_000.])**mu
n = len(x)
assert np.all(counts[:, 1] >= counts[:, 0])
print(f'N={n}, mean={x.mean(axis=0)}, SD={x.std(axis=0, ddof=1)}, paired RMS={np.mean(np.diff(x,axis=1)**2)**0.5:.7f}')

def kernel(q, h, data):
    v = np.abs((data-q)/h)
    return 0.75/h * np.maximum(1-v*v, 0)

def stats(v):
    return v.mean(), v.std(ddof=1)/n**0.5

print('x    kappa=1     h=.01               h=.02               h=.04             paired hi-lo(h=.02)')
for q in [.01,.02,.04,.06,.08,.10,.15,.20,.30,.40,.50,.60,.70,.80,.90,1.,1.2,1.4,1.6,1.8,2.]:
    row = [f'{q:.2f}',f'{q**(1/mu-1)/mu:.5f}']
    for h in [.01,.02,.04]:
        est, se = stats(kernel(q,h,x[:,1])); row.append(f'{est:.5f} +/- {se:.5f}')
    est,se = stats(kernel(q,.02,x[:,1])-kernel(q,.02,x[:,0])); row.append(f'{est:+.5f} +/- {se:.5f}')
    print('   '.join(row))
print('left-bin comparison (series proxy fitted at kappa=1, bin mean exact)')
for a,b in [(0,.01),(.01,.02),(.02,.04),(.04,.06),(.06,.08),(.08,.10)]:
    obs = ((x[:,1]>a)&(x[:,1]<=b))/(b-a)
    est,se=stats(obs)
    print(f'({a:.2f},{b:.2f}] count={int(round(obs.sum()*(b-a))):4d} density={est:.5f} +/- {se:.5f}; ref={(b**(1/mu)-a**(1/mu))/(b-a):.5f}')
