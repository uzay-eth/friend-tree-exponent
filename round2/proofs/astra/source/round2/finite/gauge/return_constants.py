"""Integer-power and rational checks for the binary excess row."""
from fractions import Fraction as F

l,u = F('.565'),F('.5665')
anchor, row_radius = F('.54001'),F('5.435')
kA,kB,dA = F('.47131'),F('.049408'),F('.19316')
for s,R in ((F('.56499'),F('1.10131')),(l,F('1.10061')),(F('.565615'),F('1.0583')),
            (F('.565809'),F('1.0453')),(u,F(1))):
    power = (u-s)/(u-anchor)
    assert row_radius**power.numerator <= R**power.denominator
    excess = R-(R+dA)*kB/(R-kA)
    out = (10**7*excess).__ceil__()
    assert excess<(F('.999791') if s<l else F('.9991'))
    print('s',s,'radius upper',R,'excess high-mark column upper',out,'/ 10000000')
critical = 1-(1+dA)*kB/(1-kA)
print('at mu: excess high-mark column upper', (10**7*critical).__ceil__(), '/ 10000000')

# The excess itself is a positive time Laplace transform.
anchor_excess, critical_excess = F(5),F('.888495')
for s,ceiling in ((F('.56469'),F('.99983')),(F('.56470'),F('.99917')),
                  (F('.56475'),F('.99592')),(F('.56480'),F('.99268')),
                  (l,F('.97982')),(F('.565615'),F('.94130')),
                  (F('.565809'),F('.92947'))):
    power=(u-s)/(u-anchor)
    assert anchor_excess**power.numerator * critical_excess**(power.denominator-power.numerator) < ceiling**power.denominator
    print('s',s,'direct excess interpolation <',ceiling)
sstar = F('.56470')
hmin = F('.501218471')
w = F('.974')/(F(1,16)+F('.14'))
for ceiling,bound in ((F('.99917'),11600),(F('.99268'),1320)):
    assert 1+w/(1-ceiling)/hmin<bound
delta=sstar-F('.54')
assert delta*F('.4')*F('1.072734084')/(F(1,16)+F('.14'))<=hmin
assert F('.4')-F(1,16)+delta*(1+F('.4')/2)<sstar
print('active/dormant comparison factor < 11600 at .56470, < 1320 at .56480, size power',delta)
