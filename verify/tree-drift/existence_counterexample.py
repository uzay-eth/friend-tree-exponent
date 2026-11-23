"""Enumerate the generator of the inside-core count twice at a fixed initial tree."""
from fractions import Fraction


def attach_leaf(tree, receiver):
    next_tree = [neighbors.copy() for neighbors in tree]
    next_tree[receiver].append(len(tree))
    next_tree.append([receiver])
    return next_tree


def free_tree():
    tree = [[1, 2], [0], [0]]
    for center in (1, 2):
        for _ in range(9):
            tree = attach_leaf(tree, center)
    return tree


def core_count(tree, inside):
    return sum(len(tree[v]) >= 2 for v in inside)


def derivatives(tree, inside):
    # Lf(T) sums f(T') - f(T) over ordered source/receiver pairs at rate 1/D_source.
    def generator(tree, inside, function):
        value = function(tree, inside)
        total = Fraction(0)
        for neighbors in tree:
            rate = Fraction(1, len(neighbors))
            for receiver in neighbors:
                next_inside = inside | {len(tree)} if receiver in inside else inside
                total += rate * (function(attach_leaf(tree, receiver), next_inside) - value)
        return total

    def slope(tree, inside):
        return generator(tree, inside, core_count)

    return core_count(tree, inside), slope(tree, inside), generator(tree, inside, slope)


tree = free_tree()
inside = set(range(len(tree)))
free = derivatives(tree, inside)
parent = len(tree)
tree = attach_leaf(tree, 0)
for _ in range(3):
    tree = attach_leaf(tree, parent)
actual = derivatives(tree, inside)
print("free count, slope, curvature:", *free)
print("actual count, slope, curvature:", *actual)
print("actual minus free curvature:", actual[2] - free[2])
assert free == (3, Fraction(9, 5), Fraction(1583, 1650))
assert actual == (3, Fraction(9, 5), Fraction(13229, 13200))
assert actual[2] - free[2] == Fraction(113, 2640)
