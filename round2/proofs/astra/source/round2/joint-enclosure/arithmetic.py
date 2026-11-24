"""Exact composition of the two signed ledgers; contains no solver calls."""

from fractions import Fraction as Q
import re

S = 2**40
PRIOR = (Q(113, 200), Q(1133, 2000))
SIDES = ("lower", "upper")
TRIAL = {"lower": Q(5656453, 10**7), "upper": Q(56581165, 10**8)}
ERROR = {"lower": 115075, "upper": 155608}
CEILING = {"lower": Q("2.024576014096"), "upper": Q("2.032667251118")}
OLD_Q4 = {"lower": 7402100, "upper": 7709975}
NEW_Q4 = {"lower": 36290118, "upper": 32703071}
PARENT = {"lower": 11253928, "upper": 8148713}  # upstream receiver names
CHILD = {"lower": 301604, "upper": 3226994}
MASK63 = {"lower": 15289028, "upper": 17458365}
LEGACY_Y = {"lower": 2654889165321, "upper": 2655887184018}
LEGACY = {"lower": Q(565662790716854568223, 1000008252270910000000),
          "upper": Q(1698046074520368376961, 3001152517940340000000)}
EXPECTED = (Q(1356995264017948174596364819, 2398901654969101844480000000),
            Q(2725370105716046104606536597, 4816977775848776007680000000))
TICKS = (565673570322, 565784239110)


def require(ok, message):
    if not ok:
        raise ValueError(message)


def one(pattern, text):
    matches = re.findall(pattern, text, flags=re.M)
    require(len(matches) == 1, f"expected exactly one row: {pattern}")
    return matches[0]


def read_ledgers(physical, split, ceilings):
    """Bind only the checked rows, including the source IDs and scale."""
    ledgers = {}
    for mode, side in ((2, "lower"), (1, "upper")):
        error = int(one(rf"^{side} S={S} cap=(\d+) min=\d+ max=\d+ base=\S+ "
                        r"universal arithmetic PASS$", physical))
        reserve = rf" rows [0-9 ]+ newborn reserve (\d+)/{S} at discount 1133/2000$"
        old = int(one(rf"^mode {mode} degree caps 192 48 16 9 6 core caps 30 14 9 6 5" +
                      reserve, physical))
        slots = [int(one(rf"^mode {mode} {kind} vertex slack degree caps 1399 1024 16 5 4 "
                         r"core caps 48 24 8 5 4" + reserve, physical)) for kind in ("parent", "child")]
        mask = int(one(rf"^{side} mask 63 at 1133/2000 rows 55259 510900 21172384 10109560 "
                       rf"positive sources(?: \d+){{4}} newborn (\d+)/{S}$", physical))
        archived, p0, p1, p2, p3, p4, total, increment = map(int, one(
            rf"^{side}\t" + "\t".join([r"(\d+)"] * 8) + "$", split))
        require(old == archived == OLD_Q4[side] and p0+p1+p2+p3+p4 == total and
                total-old == increment, f"{side}: quartet replacement accounting")
        require(slots == [PARENT[side], CHILD[side]] and mask == MASK63[side],
                f"{side}: changed reviewed credit")
        sign = 1 if side == "lower" else -1
        legacy = TRIAL[side] + sign*(old+sum(slots)+mask-Q(error)/PRIOR[0])/LEGACY_Y[side]
        require(legacy == LEGACY[side] == Q(one(rf"^{side} fraction ([0-9/]+) outward 0\.\d+$", physical)),
                f"{side}: legacy composition mismatch")
        ledgers[side] = {"sources": [("quartet-replacement", total), ("parent", slots[0]),
                                     ("child", slots[1]), ("mask63", mask)],
                         "trial": TRIAL[side], "loss": error, "ceiling": ceilings[side]}
    return ledgers


def enclosure(ledgers, prior=PRIOR):
    require(prior == PRIOR, "only the original [0.565, 0.5665] prior is allowed")
    lo, hi = prior
    ends = []
    for side in SIDES:
        row = ledgers[side]
        expected = {"quartet-replacement": NEW_Q4[side], "parent": PARENT[side],
                    "child": CHILD[side], "mask63": MASK63[side]}
        sources = row["sources"]
        require(len(sources) == len(expected) and dict(sources) == expected,
                f"{side}: overlapping, ineligible or changed source")
        lam, e, y = (row[key] for key in ("trial", "loss", "ceiling"))
        require(lam == TRIAL[side] and e == ERROR[side] and y == CEILING[side],
                f"{side}: wrong trial, residual or unsigned K7/144 ceiling")
        v = sum(value for _, value in sources)
        net = Q(v)-Q(e)/lo
        require(net > 0 and y > 0, f"{side}: no positive inward shift")
        ends.append(lam + (1 if side == "lower" else -1)*net/(S*y))
    require(lo < ends[0] <= ends[1] < hi and LEGACY["lower"] < ends[0] and
            ends[1] < LEGACY["upper"], "joint result not nested within prior and legacy")
    return tuple(ends)


def ticks(ends):
    scale = 10**12
    a, b = (q*scale for q in ends)
    return a.numerator//a.denominator, -(-b.numerator//b.denominator)


def crosscheck(ledgers, ends, reported):
    """Separate integer identities: cross-multiply, then check both adjacent ticks."""
    l, m = 113, 200
    for side, q, tick, sign in zip(SIDES, ends, reported, (1, -1), strict=True):
        w = ledgers[side]
        residual, total = w["loss"], sum(dict(w["sources"]).values())
        p, d = w["trial"].as_integer_ratio()
        a, b = w["ceiling"].as_integer_ratio()
        den = d*l*S*a
        num = p*l*S*a + sign*d*b*(l*total-m*residual)
        require(q.numerator*den == q.denominator*num, f"{side}: signed cross-product")
        point = num*10**12
        if sign == 1:
            require(tick*den <= point < (tick+1)*den, "lower tick not tight and outwards")
        else:
            require((tick-1)*den < point <= tick*den, "upper tick not tight and outwards")
    require(ends == EXPECTED and reported == TICKS, "frozen joint output differs")


def plants(ledgers, ends, reported):
    from copy import deepcopy

    def rejects(name, call, *args):
        try:
            call(*args)
        except ValueError:
            print(f"rejected {name}")
        else:
            raise ValueError(f"undetected error: {name}")

    for name, source in (("archived Q4 overlap", ("quartet-archived", OLD_Q4["lower"])),
                         ("duplicate mask63", ("mask63", MASK63["lower"])),
                         ("unreviewed mask64", ("mask64", 1))):
        other = deepcopy(ledgers)
        other["lower"]["sources"].append(source)
        rejects(name, enclosure, other)
    other = deepcopy(ledgers)
    other["lower"]["sources"][0] = ("quartet-replacement", NEW_Q4["lower"]+OLD_Q4["lower"])
    rejects("archived Q4 hidden in replacement", enclosure, other)
    rejects("narrowed growth prior", enclosure, ledgers, tuple(LEGACY[s] for s in SIDES))
    other = deepcopy(ledgers)
    other["upper"]["loss"] *= -1
    rejects("flipped residual", enclosure, other)
    wrong = (ends[0], TRIAL["upper"]+(TRIAL["upper"]-ends[1]))
    rejects("wrong upper shift sign", crosscheck, ledgers, wrong, reported)
    rejects("inward lower tick", crosscheck, ledgers, ends, (reported[0]+1, reported[1]))
    rejects("inward upper tick", crosscheck, ledgers, ends, (reported[0], reported[1]-1))
