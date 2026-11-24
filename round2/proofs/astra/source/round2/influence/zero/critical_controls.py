#!/usr/bin/env python3
"""Exact scalar checks for the signed critical theorem (no production replay)."""
from fractions import Fraction as F


l, u = F('.565658'), F('.565798')
assert F('.56564') < l < F('.5656581227') <= F('.5657979941') < u
assert F(50, 7) ** 4 < F('1.006') ** 1395
assert F('.978128') < F('.979')
assert F('.117') < min(F('.11793629'), F('.11910171'))
assert 1 / (1 - F('.8600511260')) < F('7.146')

P, M, m, delta = F(200, 21), F(50, 7), F('.117'), F('.021')
v = 1 + F('.006') / delta
A = P * (1 + M / (m * delta))
B = P * v * M / (m * delta)
inverse = 2 * u * (u + 2) ** 2 * (A + B)
critical_inverse = u * (u + 2) * P * M / m
assert inverse < 471522
assert critical_inverse < 845

print('PASS: bracket inclusion, positive-trigger comparisons, signed inverse norm < 7.146')
print('PASS: full-window inverse factor < 471522; critical inverse slope < 845')
