#!/usr/bin/env python3
"""Audit explicit rational/decimal literals in the manuscript against pinned receipts."""
import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PAPER = 'paper/friend-tree-exponent.tex'
PIN = 'a0972085ad72ad5685ae7b424988648d1965e73c'
EARLY = 'round2/proofs/astra/source/round2/'
LATE = 'round2/proofs/astra-late/source/round2/'
AGE = LATE + 'open/age-order/'
ROUND2_PROOFS = 'round2/proofs/'
GOLD = ROUND2_PROOFS + 'neon2/quintet-cert/certificate_gold2.json'
# Historical receipt digest in round2/verify/gold/RESULTS.md.
GOLD_SHA256 = '5e11f9c0b6eb3bd43e880be608c37d6f0b1fdbcc632b91bb5e8037db52e13284'
MARGINALS = ROUND2_PROOFS + 'neon2/stationary-marginals/certificates_all.json'
EXCLUSION = LATE + 'occupation-quintet-gold2/rationals/repaired/README.md'
HUB_EXACT = ROUND2_PROOFS + 'fable/newground/q3-hub-probability/lead-hub_order_exact.py'
SIM = ROUND2_PROOFS + 'fable/CLAIM.md'

SOURCES = (
    'round2/verify/gold/RESULTS.md', GOLD, MARGINALS, EXCLUSION,
    'proofs/bracket-motifs.md', 'proofs/upper-bound-banded.md', 'simulations/README.md',
    'simulations/outputs/sim_summary.txt', 'simulations/route-exact/replay.txt',
    ROUND2_PROOFS + 'neon2/quintet-cert/PROOF.md', ROUND2_PROOFS + 'neon2/quintet-cert/README.md',
    'round2/LEDGER.md', 'round2/who-did-what.md',
    EARLY + 'fringe/renewal.md', EARLY + 'fringe/asymptotics.md',
    EARLY + 'standalone/claim-6.md', EARLY + 'fringe/local-certificates.md',
    EARLY + 'finite/one-port-upper.md', EARLY + 'fringe/finite-kernel.md',
    LATE + 'open/forest-age/inert/README.md',
    AGE + 'early/README.md', AGE + 'hub-ranking/README.md',
    AGE + 'hub-ranking/all-powers/README.md', AGE + 'hub-ranking/all-powers/positive/README.md',
    AGE + 'hub-ranking/tip-moments/check.py',
    ROUND2_PROOFS + 'fable/newground/q3-hub-probability/REPORT.md',
    ROUND2_PROOFS + 'fable/newground/stochastic-domination/COUPLING.md',
    SIM, ROUND2_PROOFS + 'neon2/q5/README.md',
    ROUND2_PROOFS + 'neon2/height-diameter/README.md', 'proofs/existence.md',
    AGE + 'global-transport/birth-window/README.md', AGE + 'global-transport/birth-window/check.py',
    AGE + 'global-transport/hostile/mixed-stock/README.md',
    AGE + 'global-transport/hostile/discharge/check.py',
    AGE + 'global-transport/hostile/interfaces/low-end.md',
    AGE + 'hub-ranking/mean-ledger/fractional/README.md',
    AGE + 'hub-ranking/mean-ledger/fractional/check.py',
    AGE + 'hub-ranking/fractional-middle/README.md',
    AGE + 'hub-ranking/fractional-middle/check.py',
    AGE + 'hub-ranking/mean-ledger/real/README.md',
    AGE + 'hub-ranking/mean-ledger/real/check.py',
    AGE + 'hub-ranking/mean-ledger/fractional/star_bound.py',
    PAPER,  # last fallback: the separately frozen round-one manuscript
)


@lru_cache(maxsize=None)
def receipt(path):
    out = subprocess.run(['git', 'show', f'{PIN}:{path}'], cwd=ROOT,
                         capture_output=True, text=True, check=True)
    return out.stdout


# Keep offsets stable while removing markup not rendered as claim text. A cite's
# locator, URL/label/path, source code, bibliography, and layout dimensions are
# identifiers or typesetting, not numerical assertions. The second argument of
# proofsrc/href remains visible and is scanned.
IGNORE = re.compile(
    r'\\(?P<cmd>proofsrc|roundsrc|href|url|path|label|(?:eq|auto|page)?ref|(?:text|paren)?cite|bibitem|'
    r'verb|includegraphics|setlength|addtolength|hspace|vspace)\*?(?![A-Za-z])')


def mask_tex(text):
    out = list(text)

    def erase(a, b):
        out[a:b] = ['\n' if c == '\n' else ' ' for c in text[a:b]]

    def argument(i, opening, closing):
        while i < len(out) and out[i].isspace():
            i += 1
        if i >= len(out) or out[i] != opening:
            return None
        depth = 1
        j = i + 1
        while j < len(out):
            if out[j] == '\\' and j + 1 < len(out) and out[j + 1] in (opening, closing, '\\'):
                j += 2
                continue
            if out[j] == closing:
                depth -= 1
                if depth == 0:
                    return j + 1
            elif out[j] == opening:
                depth += 1
            j += 1
        raise ValueError(f'unclosed TeX argument at line {text.count(chr(10), 0, i) + 1}')

    for match in re.finditer(r'(?<!\\)(?:\\\\)*(%[^\n]*)', text):
        erase(match.start(1), match.end())
    source = ''.join(out)
    for match in re.finditer(r'\\begin\{(verbatim\*?|lstlisting|minted)\}.*?\\end\{\1\}', source, re.S):
        erase(*match.span())
    source = ''.join(out)
    # Printed section/theorem numbers are locators, not decimals.
    for m in re.finditer(r'(?i)(?:\b(?:theorem|lemma|corollary|proposition|section|sec|fig|equation|eq)\.?|§)'
                         r'[\s~]+(?P<locator>\d+(?:\.\d+)+)(?:\([a-z]\))?', source):
        erase(*m.span('locator'))
    source = ''.join(out)
    intro = re.search(r'\\begin\{document\}', source)
    if intro:
        # Retain reader-facing front matter even when declared in the preamble.
        front = [(m.start(), argument(m.end(), '{', '}'))
                 for m in re.finditer(r'\\(?:title|author|date)\b', source[:intro.start()])]
        begin = 0
        for a, b in front:
            if b is not None:
                erase(begin, a)
                begin = b
        erase(begin, intro.end())
    for m in re.finditer(r'\\begin\{thebibliography\}.*?\\end\{thebibliography\}', ''.join(out), re.S):
        erase(*m.span())
    for m in IGNORE.finditer(''.join(out)):
        cmd = m['cmd']
        end = m.end()
        if cmd == 'verb':
            if end < len(out) and not out[end].isspace():
                close = text.find(out[end], end + 1)
                if close > end and '\n' not in text[end:close]:
                    end = close + 1
        else:
            while (opt := argument(end, '[', ']')) is not None:
                end = opt
            for _ in range(2 if cmd in {'setlength', 'addtolength'} else 1):
                end = argument(end, '{', '}') or end
        erase(m.start(), end)
    return ''.join(out)


# TeX accepts a single unbraced token per frac argument. Do not treat symbolic
# quotients as literals. Fractions in a compound expression still count.
GAP = r'(?:\s|\\[,;!: ]|~)*'
INT = r'[0-9]+(?:\{,\}[0-9]{3})*'
UNSIGNED_DEC = rf'(?:{INT}\.[0-9]+|\.[0-9]+)'
SIGNED_INT = rf'[+\-−]?{GAP}{INT}(?:{GAP}\^{GAP}(?:\{{[+\-−]?\d+\}}|\d))?'
FRAC_ARG = rf'(?:\{{{GAP}(?P<{{name}}>{SIGNED_INT}){GAP}\}}|(?P<{{name}}_single>[0-9]))'
NUM = FRAC_ARG.replace('{name}', 'num')
DEN = FRAC_ARG.replace('{name}', 'den')
TOKEN = re.compile(
    rf'(?P<texfrac>\\[dt]?frac{GAP}{NUM}{GAP}{DEN})'
    rf'|(?P<script>(?<![0-9.]){INT}{GAP}[/⁄]{GAP}(?P<sden>{SIGNED_INT})(?![\dA-Za-z]|\.\d))'
    rf'|(?P<decimal>(?<![0-9.]){UNSIGNED_DEC}(?![\d]|\.\d))'
    rf'(?:(?:{GAP}(?:\\times|\\cdot|×|\*){GAP}10{GAP}(?:\^{GAP}(?:\{{(?P<exp>[+\-−]?\d+)\}}|(?P<exp_single>\d))|(?P<sup>[⁺⁻⁰¹²³⁴⁵⁶⁷⁸⁹]+))'
    rf'|[eE](?P<eexp>[+\-−]?\d+)))?')
SIGN = re.compile(rf'(?P<sign>[+\-−]){GAP}$')
COMPARISON = re.compile(r'(?:<=|>=|≤|≥|<|>|=|\\(?:leq?|geq?|approx|sim|equiv))\s*(?:\{\s*)?$')


@dataclass(frozen=True)
class Literal:
    line: int
    spelling: str
    value: Fraction
    kind: str
    comparison: str
    start: int
    end: int


def integer(token):
    return int(re.sub(GAP, '', token.replace('{,}', '').replace('−', '-')))


def integer_power(token):
    base, sep, exp = token.partition('^')
    value = integer(base)
    return Fraction(value) ** integer(exp.strip().strip('{}')) if sep else Fraction(value)


def lex(text):
    for m in TOKEN.finditer(text):
        begin = m.start()
        prefix = text[max(0, begin - 100):begin]
        # Ignore the suffix of an identifier (b_0/50, x42.3), but admit
        # juxtaposition after TeX controls (e.g. a comparison with no space).
        if not m['texfrac'] and re.search(r'[A-Za-z_]$', prefix) and not re.search(r'\\[A-Za-z]+$', prefix):
            continue
        sign = SIGN.search(prefix)
        # A minus immediately following an operand is subtraction, not unary.
        # In particular `v-1/2` uses the literal +1/2; `x\ge-1/2` uses -1/2.
        unary = False
        if sign:
            before = prefix[:sign.start()].rstrip()
            unary = (not before or not re.search(r'[\w}\])]$', before)
                     or bool(COMPARISON.search(before)) or before.endswith((r'\pm', r'\mp')))
            if unary:
                begin -= len(prefix) - sign.start()
                prefix = prefix[:sign.start()]
        if m['texfrac']:
            value = integer_power(m['num'] or m['num_single']) / integer_power(m['den'] or m['den_single'])
            kind = 'fraction'
        elif m['script']:
            a = re.match(INT, m['script'])[0]
            value = integer(a) / integer_power(m['sden'])
            kind = 'fraction'
        else:
            value = Fraction(m['decimal'].replace('{,}', ''))
            exponent = m['exp'] or m['exp_single'] or m['eexp']
            if m['sup']:
                exponent = m['sup'].translate(str.maketrans('⁺⁻⁰¹²³⁴⁵⁶⁷⁸⁹', '+-0123456789'))
            if exponent:
                value *= Fraction(10) ** integer(exponent)
            kind = 'decimal'
        if unary and sign['sign'] in '-−':
            value = -value
        comparison = COMPARISON.search(prefix)
        yield Literal(text.count('\n', 0, begin) + 1, ' '.join(text[begin:m.end()].split()), value,
                      kind, comparison[0].strip() if comparison else '', begin, m.end())


def outward(exact, places, side):
    scale = 10**places
    ticks = exact * scale
    return Fraction(ticks.numerator // ticks.denominator if side == 'lower'
                    else -((-ticks.numerator) // ticks.denominator), scale)


def audited_intervals():
    gold = json.loads(receipt(GOLD))
    yield GOLD, 'certified_interval_dec', [gold[s]['certified_exact'] for s in ('lower', 'upper')], gold['certified_interval_dec']
    marginals = json.loads(receipt(MARGINALS))
    for key, row in marginals['observables'].items():
        yield (MARGINALS, f'observables.{key}.certified_interval_dec',
               row['certified_interval_exact'], row['certified_interval_dec'])


def derived_receipts():
    source = EARLY + 'fringe/renewal.md'
    text = receipt(source)
    weights = [re.search(r'h\(' + i + r'\)=([0-9/]+)', text) for i in ('0', '1', '2', '3', 'c')]
    params = [re.search(param + r'=([0-9]+/[0-9]+)', text) for param in ('d', 'q')]
    if not all(weights + params):
        raise ValueError('pinned renewal weight receipt no longer matches derivation')
    d, q = (Fraction(m[1]) for m in params)
    h = [Fraction(m[1]) for m in weights]
    h.append(h[-1])
    for c in range(1, 5):
        r, v = q*h[c], h[c+1]-h[c]
        a, b, e = (r-v)*d, r-v-1+c*d*v, c*(1+v)
        if c == 1:
            for v in (a, a*9 + b*3 + e, a*16 + b*4 + e):
                yield v, f'{source} (exact renewal quadratic, c=1)'
        else:
            yield b*b-4*a*e, f'{source} (exact renewal discriminant, c={c})'
    gold = json.loads(receipt(GOLD))
    lo, hi = (Fraction(gold[s]['certified_exact']) for s in ('lower', 'upper'))
    for places in (6, 7, 8, 9, 10, 11, 12):
        yield outward(lo, places, 'lower'), f'{GOLD}#lower.certified_exact (outward to {places} decimals)'
        yield outward(hi, places, 'upper'), f'{GOLD}#upper.certified_exact (outward to {places} decimals)'
    width = Fraction(gold['certified_width_ppm'])
    if (hi - lo) * 10**6 > width:
        raise ValueError(f'{GOLD}: recorded width rounds inward')
    yield width / 10**6, f'{GOLD}#certified_width_ppm (upper bound on exact width)'
    real = AGE + 'hub-ranking/mean-ledger/real/check.py'
    cutoff = int(re.search(r'^D = (\d+)$', receipt(real), re.M)[1])
    yield Fraction(1, cutoff + 1), f'{real} (degree-only cap U_p(d+1) <= 1/(d+1) at d = D)'
    sim = receipt(SIM)
    matches = list(re.finditer(r'^\(S[13]\).*? = (\d+\.\d+) \+-', sim, re.M))
    if len(matches) != 2:
        raise ValueError(f'{SIM}: the expected two simulation receipts are missing')
    for match in matches:
        yield Fraction(round(Fraction(match[1]) * 10**8), 10**8), (
            f'{SIM}:{sim.count(chr(10), 0, match.start()) + 1} (nearest 8 decimal places)')
    import os
    import sys
    checker = receipt(HUB_EXACT)
    # Replace only the display formatting so the pinned fixture emits the
    # underlying exact fractions instead of binary floating-point values.
    display = '{float(d):+.6f}, {float(e):+.6f}'
    if checker.count(display) != 1:
        raise ValueError('unexpected early-hub display format')
    checker = checker.replace(display, '{d}, {e}')
    result = subprocess.run([sys.executable, '-IB', '-c', checker, '9'],
                            env={**os.environ, 'NOEDGE': '1'}, check=True, capture_output=True, text=True)
    matches = re.findall(r'p_([2-7]) - p_([3-8]) in \[([0-9]+/[0-9]+),', result.stdout)
    if ([int(m[0]) for m in matches] != list(range(2, 8))
            or '40320 histories, 286 rooted shapes' not in result.stdout
            or 'certified: YES' not in result.stdout):
        raise ValueError('unexpected independent early-hub result')
    for label, _, endpoint in matches:
        yield outward(Fraction(endpoint), 4, 'lower'), (
            f'{HUB_EXACT} (NOEDGE=1, size 9, downward endpoint for p_{label} gap)')


def index_receipts():
    index = defaultdict(list)
    # Keep sign and endpoint direction exact. Never allow an arbitrary number
    # merely because it is on one side of an endpoint. Authenticate the
    # printed decimal receipt by recomputing its directed rounding.
    for path, key, exact, decimals in audited_intervals():
        for endpoint, shown, side in zip(exact, decimals, ('lower', 'upper'), strict=True):
            places = len(shown.partition('.')[2])
            if outward(Fraction(endpoint), places, side) != Fraction(shown):
                raise ValueError(f'{path}:{key}: bad outward {side} rounding: {shown}')
            index[Fraction(shown)].append(f'{path}#{key} ({side}, rounded outward from {endpoint})')
    for path in SOURCES:
        text = receipt(path)
        if path == PAPER:
            text = mask_tex(text)
        for lit in lex(text):
            index[lit.value].append(f'{path}:{lit.line}')
        if path.endswith('.py'):
            for m in re.finditer(r'\bQ\(\s*([+-]?\d+)\s*,\s*(\d+(?:\*\*\d+)?)\s*\)', text):
                value = Fraction(int(m[1])) / integer_power(m[2].replace('**', '^'))
                index[value].append(f'{path}:{text.count(chr(10), 0, m.start()) + 1} (exact Fraction constructor)')
    for value, source in derived_receipts():
        index[value].insert(0, source)
    return index


def reversed_bounds(tex, literals):
    """Reject a decreasing literal interval/comparison chain, even if both values are known."""
    sep = r'(?:\s|\$|\\[,;!: ]|~)*'
    edge = r'(?:\\(?:left|right)\s*)?'
    for a, b in zip(literals, literals[1:]):
        if a.value <= b.value:
            continue
        between = tex[a.end:b.start]
        bracketed = (re.fullmatch(sep + ',' + sep, between)
                     and re.search(r'\[' + sep + '$', tex[max(0, a.start-40):a.start])
                     and re.match(sep + edge + r'\]', tex[b.end:]))
        chain = re.fullmatch(sep + r'(?:≤|<|<=|\\leq?)' + sep + r'\\mu' + sep
                             + r'(?:≤|<|<=|\\leq?)' + sep, between)
        if bracketed or chain:
            yield a, b


def check_enclosure(tex_path):
    data = (ROOT / GOLD).read_bytes()
    digest = sha256(data).hexdigest()
    if digest != GOLD_SHA256:
        raise ValueError(f'Gold2 receipt SHA-256 mismatch: {digest}')
    gold = json.loads(data)
    exact = []
    for side, sign in (('lower', -1), ('upper', 1)):
        row = gold[side]
        endpoint = Fraction(row['nominal']) + sign * Fraction(row['err_per_core'], row['m_min'])
        if endpoint != Fraction(row['certified_exact']):
            raise ValueError(f'Gold2 receipt: recomputed {side} endpoint differs from certified_exact')
        exact.append(endpoint)
    decimals = gold['certified_interval_dec']
    for side, endpoint, shown in zip(('lower', 'upper'), exact, decimals, strict=True):
        if (not re.fullmatch(r'\d+\.\d{12}', shown) or Fraction(shown) != outward(endpoint, 12, side)
                or shown != gold[side]['certified_dec']):
            raise ValueError(f'Gold2 receipt: incorrect outward {side} decimal endpoint')

    tex = re.sub(r'(?<!\\)(?:\\\\)*(%[^\n]*)', lambda m: m[0][:m.start(1) - m.start()],
                 tex_path.read_text(encoding='utf-8'))
    theorem = re.findall(r'\\begin\{theorem\}(?:\[[^\]]*\])?\s*\\label\{th:gold2\}(.*?)\\end\{theorem\}', tex, re.S)
    abstract = re.findall(r'\\begin\{abstract\}(.*?)\\end\{abstract\}', tex, re.S)
    if len(theorem) != 1 or len(abstract) != 1:
        raise ValueError('expected one enclosure theorem and one abstract')
    theorem, abstract = mask_tex(theorem[0]), mask_tex(abstract[0])
    spacing = r'(?:\s|\\[,;!: \n]|~)*'
    chain = spacing + r'\\leq?' + spacing + r'\\mu' + spacing + r'\\leq?' + spacing

    def pair(body, kind, middle):
        tokens = list(lex(body))
        return [(a, b) for a, b in zip(tokens, tokens[1:])
                if a.kind == b.kind == kind and re.fullmatch(middle, body[a.end:b.start])]

    fractions = pair(theorem, 'fraction', chain)
    if len(fractions) != 1 or [lit.value for lit in fractions[0]] != exact:
        raise ValueError('Theorem 2 exact endpoints do not match the recomputed Gold2 enclosure')
    interval = [(a, b) for a, b in pair(theorem, 'decimal', GAP + ',' + GAP)
                if re.search(r'\[' + GAP + '$', theorem[:a.start]) and re.match(GAP + r'\]', theorem[b.end:])]
    for name, pairs in (('Theorem 2', interval), ('abstract', pair(abstract, 'decimal', chain))):
        if len(pairs) != 1 or [lit.spelling for lit in pairs[0]] != decimals:
            raise ValueError(f'{name} outward decimal endpoints do not match the Gold2 enclosure')

    print(f'Exact enclosure: [{exact[0]}, {exact[1]}]')
    print(f'Outward decimal enclosure: [{decimals[0]}, {decimals[1]}]')
    print('PASS: Gold2 receipt matches Theorem 2 and the abstract.')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--tex', type=Path, default=ROOT / PAPER, help='audit another TeX input (e.g. writer fragment)')
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument('--all', action='store_true', help='list every matched occurrence and a pinned source location')
    mode.add_argument('--enclosure', action='store_true', help='check the Gold2 enclosure offline without Git history')
    args = ap.parse_args()
    if args.enclosure:
        try:
            check_enclosure(args.tex)
        except (OSError, ValueError) as error:
            print(f'FAIL: {error}')
            return 1
        return 0
    text = mask_tex(args.tex.read_text(encoding='utf-8'))
    literals = list(lex(text))
    index = index_receipts()
    counts = Counter()
    missing = 0
    for lit in literals:
        if provenance := index.get(lit.value):
            counts[lit.kind] += 1
            if args.all:
                display = (lit.comparison + ' ' if lit.comparison else '') + lit.spelling
                print(f'{args.tex}:{lit.line}: {display} [{lit.value}] <- {provenance[0]}')
        else:
            missing += 1
            print(f'{args.tex}:{lit.line}: UNMATCHED {lit.kind} {lit.comparison} {lit.spelling} [{lit.value}]')
    reversals = list(reversed_bounds(text, literals))
    # The canonical manuscript's headline theorem must print both exact
    # endpoints. This catches deletion or replacement by a different known
    # bound, which a value-membership test alone cannot catch.
    gold_theorem = re.search(r'\\begin\{theorem\}[^\n]*(?:\n\s*)?\\label\{th:gold2\}(.*?)\\end\{theorem\}',
                             args.tex.read_text(encoding='utf-8'), re.S)
    gold_error = False
    if gold_theorem:
        exact = json.loads(receipt(GOLD))
        pair = [Fraction(exact[s]['certified_exact']) for s in ('lower', 'upper')]
        theorem_values = [lit.value for lit in lex(mask_tex(gold_theorem[1]))]
        gold_error = not all(v in theorem_values for v in pair)
        if gold_error:
            print(f'{args.tex}: th:gold2 must include both exact archived endpoints')
    for a, b in reversals:
        print(f'{args.tex}:{a.line}: REVERSED bound: {a.spelling} > {b.spelling}')
    checked = sum(counts.values())
    print(f'{checked}/{checked + missing} occurrences matched ({counts["fraction"]} fractions, '
          f'{counts["decimal"]} decimals); {missing} unmatched; {len(reversals)} reversed bounds; '
          f'pinned receipts at {PIN}')
    return int(bool(missing or reversals or gold_error or not checked))


if __name__ == '__main__':
    raise SystemExit(main())
