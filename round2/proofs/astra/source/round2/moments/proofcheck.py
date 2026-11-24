"""Second dual calculation with a common integer unit for the entire inequality."""

from fractions import Fraction
from math import lcm


def common_unit(numbers):
    unit = 1
    for z in numbers:
        unit = lcm(unit, z.denominator)
    return unit


class IntegerProof:
    def __init__(self, model):
        self.model = model
        unit = 1
        for _, row, rhs in [*model.equal, *model.atleast]:
            unit = lcm(unit, rhs.denominator, common_unit(row.values()))
        self.coeff_unit = unit
        self.cap_unit = common_unit(model.cap)

    def evaluate(self, objective, witness):
        model = self.model
        unit = lcm(self.coeff_unit, common_unit(objective.values()))
        yd = witness['dual_scale']
        if type(yd) is not int or yd <= 0:
            raise ValueError('invalid integer unit')
        def integer(z):
            return z.numerator*(unit//z.denominator)
        left = [0]*len(model.names)
        for j, z in objective.items():
            left[j] = yd*integer(z)
        constant = 0
        for tag, rows in [('eq', model.equal), ('ge', model.atleast)]:
            seen = set()
            for i, yi in witness[tag]:
                if (type(i) is not int or i in seen or not 0 <= i < len(rows)
                        or type(yi) is not int or (tag == 'ge' and yi < 0)):
                    raise ValueError('invalid dual index or sign')
                seen.add(i)
                _, row, rhs = rows[i]
                constant += yi*integer(rhs)
                for j, a in row.items():
                    left[j] -= yi*integer(a)
        cap_unit = self.cap_unit
        total = constant*cap_unit
        for residual, cap in zip(left, model.cap):
            if residual < 0:
                total += residual*cap.numerator*(cap_unit//cap.denominator)
        return Fraction(total, unit*yd*cap_unit)
