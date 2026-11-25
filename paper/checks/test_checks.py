#!/usr/bin/env python3
"""Focused positive and negative controls for the two editorial checks."""
from fractions import Fraction as Q
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest

import check_links as links
import check_literals as numbers


class LexingTests(unittest.TestCase):
    def test_literal_spellings(self):
        tex = r'$s\ge-0.551315$, $x\le0.56564$, $\frac{-97}{200}$, and $-3591/10^6$.' '\n' \
              r'$v-1/2$, $b_0/50$, $\tfrac9{5000}$, $a/3$, and $0.5657236606\pm6.2\times10^{-9}$.'
        self.assertEqual([lit.value for lit in numbers.lex(tex)], list(map(Q, (
            '-0.551315', '0.56564', '-97/200', '-3591/1000000', '1/2', '9/5000',
            '0.5657236606', '0.0000000062',
        ))))
        self.assertEqual([lit.comparison for lit in numbers.lex(tex)][:2], [r'\ge', r'\le'])

    def test_visible_text_not_paths(self):
        tex = r'''\usepackage[margin=1.05in]{geometry}
\title{Title with 7/8}
\begin{document}
% comment 1/100
\proofsrc{round2/0.00/99.md}{An interval $[.30,.50]$}.
\href{https://example.test/1.50}{Visible 1/2}. \label{eq:7.1}
\cite[Theorem 6.1]{C}
\begin{thebibliography}{1}
DOI: 42.4567
\end{thebibliography}
'''
        literals = list(numbers.lex(numbers.mask_tex(tex)))
        self.assertEqual([(lit.line, lit.value) for lit in literals],
                         [(2, Q('7/8')), (5, Q('.30')), (5, Q('.50')), (6, Q('1/2'))])

    def test_directional_rounding_with_negative_endpoints(self):
        self.assertEqual(numbers.outward(Q('-2/3'), 2, 'lower'), Q('-.67'))
        self.assertEqual(numbers.outward(Q('-2/3'), 2, 'upper'), Q('-.66'))

    def test_reversed_comparison_values(self):
        source = r'$0.565736986679\le\mu\le0.565713486333$ and $[0.565736986679,0.565713486333]$.'
        self.assertEqual(len(list(numbers.reversed_bounds(source, list(numbers.lex(source))))), 2)

    def test_unrecognized_inward_edit_fails(self):
        with TemporaryDirectory() as work:
            tex = Path(work) / 'mutated.tex'
            tex.write_text(r'Refined lower bound: $\mu\ge0.565714341934$.', encoding='utf-8')
            result = subprocess.run([sys.executable, '-B', numbers.__file__, '--tex', str(tex)],
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 1, result.stderr)
            self.assertIn('UNMATCHED decimal', result.stdout)
            self.assertIn('0.565714341934', result.stdout)


class EnclosureTests(unittest.TestCase):
    def test_offline_enclosure_rejects_one_changed_endpoint(self):
        # Recreate just the working-tree inputs, with no .git or Git on PATH
        # and with Python's site packages disabled.
        with TemporaryDirectory() as work:
            root = Path(work)
            script = root / 'paper/checks/check_literals.py'
            for rel in ('paper/checks/check_literals.py', numbers.PAPER, numbers.SUPPLEMENT, numbers.GOLD, numbers.REFINED_CHECK,
                        numbers.REFINED_RECEIPT, numbers.REFINED_CLASSES):
                dest = root / rel
                dest.parent.mkdir(parents=True, exist_ok=True)
                dest.write_bytes((numbers.ROOT / rel).read_bytes())
            command = [sys.executable, '-ISB', str(script), '--enclosure']
            options = {'env': {'PATH': ''}, 'cwd': root, 'capture_output': True, 'text': True}
            result = subprocess.run(command, **options)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertTrue(result.stdout.splitlines()[-1].startswith('PASS:'))

            paper = root / numbers.PAPER
            original = paper.read_text(encoding='utf-8')
            gold_path = root / numbers.GOLD
            import json
            del gold_path
            exact_lower = '3843726496355222669/6794465353700000000'
            numerator, denominator = exact_lower.split('/')
            printed = rf'\frac{{{numerator}}}{{{denominator}}}'
            changed = rf'\frac{{{int(numerator) + 1}}}{{{denominator}}}'
            # The first occurrence is the statement of Theorem ENC (th:gold2).
            self.assertLess(original.index(printed), original.index(r'\end{theorem}', original.index(r'\label{th:gold2}')))
            paper.write_text(original.replace(printed, changed, 1), encoding='utf-8')
            result = subprocess.run(command, **options)
            self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
            self.assertIn('Theorem ENC (th:gold2) exact endpoints do not match', result.stdout)
            self.assertNotIn('PASS:', result.stdout)

            # A supplement whose restated theorem differs from the paper's must fail.
            paper.write_text(original, encoding='utf-8')
            supplement = root / numbers.SUPPLEMENT
            supp_original = supplement.read_text(encoding='utf-8')
            supplement.write_text(supp_original.replace(printed, changed, 1), encoding='utf-8')
            result = subprocess.run(command, **options)
            self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
            self.assertIn('does not restate the enclosure theorem', result.stdout)
            supplement.write_text(supp_original, encoding='utf-8')

            # A changed per-class file must also fail, even with the manuscript restored.
            classes = root / numbers.REFINED_CLASSES
            rows = classes.read_bytes()
            classes.write_bytes(rows.replace(b'\n2 135889307074 1941\n', b'\n2 135889307074 1940\n', 1))
            result = subprocess.run(command, **options)
            self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
            self.assertIn('refined lower endpoint', result.stdout)


class TranscriptionTests(unittest.TestCase):
    def test_changed_formula_or_word_fails(self):
        script = numbers.ROOT / 'paper/checks/check_transcriptions.py'
        original = (numbers.ROOT / numbers.SUPPLEMENT).read_text(encoding='utf-8')
        ok = subprocess.run([sys.executable, '-B', str(script)], capture_output=True, text=True)
        self.assertEqual(ok.returncode, 0, ok.stdout)
        start = original.index('% BEGIN TRANSCRIPTION round3/q3-birth-order/PROOF.md')
        formula = original.index('\\frac{1}{d(u)', start)
        word = original.index('symmetric leaves', start)
        for mutated, message in ((original[:formula] + '\\frac{2}{d(u)' + original[formula + len('\\frac{1}{d(u)'):], 'formulas differ'),
                                 (original[:word] + 'symmetrical leaves' + original[word + len('symmetric leaves'):], 'prose differs')):
            with TemporaryDirectory() as work:
                tex = Path(work) / 'supplement.tex'
                tex.write_text(mutated, encoding='utf-8')
                result = subprocess.run([sys.executable, '-B', str(script), '--tex', str(tex)], capture_output=True, text=True)
                self.assertEqual(result.returncode, 1, result.stdout)
                self.assertIn(message, result.stdout)


class LinkTests(unittest.TestCase):
    def test_relative_paths_anchors_and_references(self):
        docs = {
            'x/a.md': '''[valid](../guide%20file.md#repeat-1) [missing][M]

[M]: <../guide file.md#absent> "Quoted title"

`[not a link](never)`

```md
[nor this](also-never)
```

<a href='../guide%20file.md#explicit'>html</a>
![bad image](absent.png)
''',
            'guide file.md': '# Repeat\n## Repeat\n<a id="explicit"></a>',
        }
        issues, count, local = links.scan(set(docs), docs.items())
        self.assertEqual((count, local), (2, 4))
        self.assertEqual([(row[2], row[3]) for row in issues],
                         [('../guide%20file.md#absent', 'reference definition'), ('absent.png', 'inline')])

    def test_pinned_proof_macro_paths(self):
        macro = (r'\newcommand{\proofsrc}[2]{\href{' + links.REPO_URL + 'blob/' + links.PIN + '/#1}{#2}}')
        text = macro + '\n' + r'\proofsrc{paper/README.md}{A good source}' '\n' \
            + r'\proofsrc{missing.md}{A bad source}'
        issues, count = links.paper_links(text, {'paper/README.md'}, {})
        self.assertEqual(count, 2)
        self.assertEqual([row[2] for row in issues], ['missing.md'])
        issues, _ = links.paper_links(text.replace(links.PIN, 'main'), {'paper/README.md'}, {})
        self.assertTrue(any('pinned GitHub blob URL' in row[-1] for row in issues))
        issues, _ = links.paper_links(text + '\n' + r'\\proofsrc{paper/README.md}{bad escape}', {'paper/README.md'}, {})
        self.assertTrue(any('single backslash' in row[-1] for row in issues))

    def test_order_macro_paths(self):
        base = (r'\newcommand{\proofsrc}[2]{\href{' + links.REPO_URL + 'blob/' + links.PIN + '/#1}{#2}}' '\n'
                r'\proofsrc{paper/README.md}{A pinned source}' '\n')
        macro = r'\newcommand{\ordersrc}[2]{\href{' + links.REPO_URL + 'blob/' + links.ORDER_PIN + '/#1}{#2}}'
        calls = '\n' + r'\ordersrc{round3/README.md}{good}' '\n' + r'\ordersrc{round3/absent.md}{bad}'
        extra = {'ordersrc': (links.ORDER_PIN, {'round3/README.md'}, {})}
        issues, count = links.paper_links(base + macro + calls, {'paper/README.md'}, {}, extra)
        self.assertEqual(count, 3)
        self.assertEqual([row[2] for row in issues], ['round3/absent.md'])
        # A path that exists only at the later pin is not accepted against the older pin.
        issues, _ = links.paper_links(base + r'\proofsrc{round3/README.md}{wrong pin}', {'paper/README.md'}, {}, extra)
        self.assertEqual([row[2] for row in issues], ['round3/README.md'])
        # The macro must expand to its own pinned commit.
        issues, _ = links.paper_links(base + macro.replace(links.ORDER_PIN, 'main') + calls,
                                      {'paper/README.md'}, {}, extra)
        self.assertTrue(any('ordersrc must expand' in row[-1] for row in issues))
        # Without the later pin configured, its expansion is an unpinned repository link.
        issues, _ = links.paper_links(base + macro + calls, {'paper/README.md'}, {})
        self.assertTrue(any('must use proofsrc or ordersrc' in row[-1] for row in issues))

    def test_other_pinned_markdown_links(self):
        good = links.REPO_URL + 'blob/' + links.ORDER_PIN + '/round3/README.md'
        bad = links.REPO_URL + 'blob/' + links.ORDER_PIN + '/round3/absent.md'
        unknown = links.REPO_URL + 'blob/' + '0' * 40 + '/anything.md'
        docs = {'a.md': f'[good]({good}) [bad]({bad}) [unknown]({unknown})'}
        others = {links.ORDER_PIN: ({'round3/README.md'}, {})}
        issues, _, _ = links.scan({'a.md'}, docs.items(), set(), {}, others)
        self.assertEqual([row[2] for row in issues], [bad])

    def test_order_macro_is_masked_for_literals(self):
        literals = list(numbers.lex(numbers.mask_tex(r'\ordersrc{round3/q2-a/1.50.md}{Visible 1/2}')))
        self.assertEqual([lit.value for lit in literals], [Q('1/2')])


if __name__ == '__main__':
    unittest.main()
