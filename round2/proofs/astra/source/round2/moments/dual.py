"""Propose a sparse dual with HiGHS; verify its bound using rational arithmetic."""

from fractions import Fraction as F
from math import ceil, floor

from occupation import add, combine


DUAL_SCALE = 10**14


def verify(model, objective, witness):
    """Lower bound c.x over Ax=b, Gx>=h, and 0<=x<=cap (any dual is safe)."""
    scale = witness['dual_scale']
    if not isinstance(scale, int) or scale <= 0:
        raise ValueError('bad dual denominator')
    residual = dict(objective)
    bound = F()
    for name, rows in [('eq', model.equal), ('ge', model.atleast)]:
        seen = set()
        for i, n in witness[name]:
            if (type(i) is not int or type(n) is not int or i in seen or not 0 <= i < len(rows)
                    or (name == 'ge' and n < 0)):
                raise ValueError('bad dual coordinate')
            seen.add(i)
            y = F(n, scale)
            _, row, rhs = rows[i]
            bound += rhs*y
            for j, a in row.items():
                add(residual, j, -a*y)
    # This exact box correction also certifies rounded duals at a degenerate optimum.
    return bound + sum((min(0, z)*model.cap[j] for j, z in residual.items()), F())


class Solver:
    def __init__(self, model, warm=False):
        import numpy as np
        from scipy.sparse import coo_matrix
        self.model, self.np = model, np
        def matrix(rows, sign):
            i, j, data = [], [], []
            for idx, (_, coeff, _) in enumerate(rows):
                for col, val in coeff.items():
                    i.append(idx)
                    j.append(col)
                    data.append(sign*float(val))
            A = coo_matrix((data, (i, j)), shape=(len(rows), len(model.names))).tocsr()
            b = np.array([sign*float(rhs) for _, _, rhs in rows])
            return A, b
        self.eq, self.beq = matrix(model.equal, 1)
        self.ub, self.bub = matrix(model.atleast, -1)
        self.bounds = np.array([(0, u) for u in model.cap])
        self.warm = None
        if warm:
            try:
                self._warm_model()
            except (ImportError, AttributeError, TypeError, ValueError):
                self.warm = None

    def _warm_model(self):
        # The copy of HiGHS bundled with SciPy can keep a simplex basis when
        # Dinkelbach changes only the costs. The public linprog route is the fallback.
        from scipy.optimize._highspy import _core
        np=self.np
        h=_core._Highs()
        def check(status):
            if status==_core.HighsStatus.kError: raise ValueError('HiGHS model construction failed')
        for key,value in [('output_flag',False),('threads',1),('parallel','off'),('solver','simplex'),
                          ('primal_feasibility_tolerance',1e-9),('dual_feasibility_tolerance',1e-9)]:
            check(h.setOptionValue(key,value))
        n=len(self.model.names)
        b=np.asarray(self.bounds,dtype=float)
        check(h.addCols(n,np.zeros(n),b[:,0],b[:,1],0,np.zeros(n+1,dtype=np.int32),
                        np.empty(0,dtype=np.int32),np.empty(0)))
        for a,lo,hi in [(self.eq,self.beq,self.beq),(self.ub,np.full(len(self.bub),-np.inf),self.bub)]:
            check(h.addRows(a.shape[0],lo,hi,a.nnz,a.indptr,a.indices,a.data))
        self.warm=h,_core.HighsModelStatus.kOptimal,np.arange(n,dtype=np.int32)

    def optimize(self, obj):
        import warnings
        from scipy.optimize import OptimizeWarning, linprog
        weight = 10000
        c = self.np.zeros(len(self.model.names))
        for j, z in obj.items():
            c[j] = weight*float(z)
        if self.warm:
            h,optimal,indices=self.warm
            h.changeColsCost(len(c),indices,c)
            h.run()
            sol=h.getSolution()
            if h.getModelStatus()==optimal and sol.value_valid and sol.dual_valid:
                x=sol.col_value
                eq,ge=sol.row_dual[:len(self.beq)],sol.row_dual[len(self.beq):]
            else:
                self.warm=None
        if self.warm is None:
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', category=OptimizeWarning, message='Unrecognized options.*')
                options = {'threads': 1, 'parallel': False}
                for method, tolerance in [('highs-ds', 1e-9), ('highs-ipm', 1e-8), ('highs-ds', 1e-7)]:
                    options['dual_feasibility_tolerance'] = tolerance
                    options['primal_feasibility_tolerance'] = tolerance
                    res = linprog(c, self.ub, self.bub, self.eq, self.beq, bounds=self.bounds,
                                  method=method, options=options)
                    if res.success: break
            if not res.success: raise RuntimeError(res.message)
            x,eq,ge=res.x,res.eqlin.marginals,res.ineqlin.marginals
        def integers(mult, marginals):
            result = []
            for i, y in enumerate(marginals):
                n = round(mult*float(y)*DUAL_SCALE/weight)
                if mult < 0:
                    n = max(0, n)
                if n:
                    result.append([i, n])
            return result
        witness = {'dual_scale': DUAL_SCALE, 'eq': integers(1, eq), 'ge': integers(-1, ge)}
        return x, witness

    def ratio(self, numerator, denominator, direction, grid=10**9, initial=0):
        """direction +1 proves num>=q*den; direction -1 proves num<=q*den."""
        def value(row, x):
            return sum(float(v)*x[j] for j, v in row.items())
        q = F(initial)
        for _ in range(12):
            x, dual = self.optimize(combine((direction, numerator), (-direction*q, denominator)))
            den = value(denominator, x)
            if den <= 0:
                raise ValueError('ratio proposal has nonpositive denominator')
            rate = value(numerator, x)/den
            if abs(rate-float(q)) < 5e-12:
                break
            q = F(round(rate*10**14), 10**14)
        n = floor(rate*grid)-1 if direction == 1 else ceil(rate*grid)+1
        if direction == 1 and n <= 0 and all(v >= 0 for v in numerator.values()):
            return dict(dual_scale=DUAL_SCALE, eq=[], ge=[], bound='0', margin='0', direction=1)
        for _ in range(12):
            endpoint = F(n, grid)
            obj = combine((direction, numerator), (-direction*endpoint, denominator))
            margin = verify(self.model, obj, dual)
            if margin >= 0:
                dual.update(bound=str(endpoint), margin=str(margin), direction=direction)
                return dual
            # A dual for a nearby objective may price a different denominator.
            _, dual = self.optimize(obj)
            margin = verify(self.model, obj, dual)
            if margin >= 0:
                dual.update(bound=str(endpoint), margin=str(margin), direction=direction)
                return dual
            n -= direction*max(1, ceil(-float(margin)*grid/den))
        raise ArithmeticError('could not certify a nearby rational')
