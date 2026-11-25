#!/usr/bin/env python3
"""Check that the single-file manuscript of commit 6fee4c4 is carried into the paper and supplement.

Every non-blank line of that manuscript must appear verbatim as a line of the new paper or of the
supplement, after the reference renames listed in RENAMES. The exceptions are the lines that were
rewritten on purpose; they are listed in REWRITTEN with the reason. The script prints a summary, and
with --list it prints every exception. Run from the repository root (needs Git history):
  python3 -B paper/checks/check_proof_parity.py
"""
import argparse, re, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASE_COMMIT = '6fee4c4d2464878ad06fff03c749c6fbede9baad'
NEW = [ROOT / 'paper' / 'friend-tree-exponent.tex', ROOT / 'paper' / 'friend-tree-exponent-supplement.tex']

SUPPSEC = ['sec:exist', 'sec:core-details', 'app:hub-certificates', 'app:stoch-order', 'app:refined', 'app:replay']


def renames(s):
    """The reference renames applied when lines moved into the supplement (see the assembly notes)."""
    s = re.sub(r'Appendices~\\ref\{(%s)\}(,? and~)\\ref\{(%s)\}' % ('|'.join(SUPPSEC), '|'.join(SUPPSEC)), r'Sections~\\ref{\1}\2\\ref{\3}', s)
    s = re.sub(r'Appendices~\\ref\{app:nonholonomic\} and~\\ref\{app:refined\}', r'Appendix~\\ref{M-app:nonholonomic} of the paper and Section~\\ref{app:refined}', s)
    s = re.sub(r'Appendix~\\ref\{(%s)\}' % '|'.join(SUPPSEC), r'Section~\\ref{\1}', s)
    s = s.replace('Appendix~\\ref{app:nonholonomic}', 'Appendix~\\ref{M-app:nonholonomic} of the paper')
    for lab in ['th:nonholonomic', 'prop:window-bound', 'lem:window-law', 'prop:nearstar', 'lem:nearstar']:
        s = s.replace('\\ref{%s}' % lab, '\\ref{M-%s}' % lab)
    s = s.replace('\\ref{sec:exponent}', '\\ref{sec:enc}')
    s = s.replace('\\subsection{Degree counts at growing degree}', '\\subsection{Degree counts at growing degree}\\label{sec:further-degrees}')
    return s


def main_renames(s):
    """Renames applied to lines that moved into the paper's appendix or statements."""
    s = s.replace('the light checker of Appendix~\\ref{app:replay} repeats this', 'the light checker described in \\Sref{sec:replay} repeats this')
    s = s.replace('The light checker described in Appendix~\\ref{app:replay}, run', 'The light checker described in \\Sref{sec:replay}, run')
    s = s.replace('reproduction operator $K_s$ of \\eqref{eq:core-kernel}', 'reproduction operator $K_s$ of \\eqref{S-eq:core-kernel} in the supplement')
    return s


# Baseline line ranges rewritten on purpose, with the reason.
REWRITTEN = [
    ((1, 34), 'preamble, front matter and abstract: rewritten for the two documents'),
    ((35, 65), 'introduction: rewritten (one result list, worked example); prior-work paragraphs 50-52 kept verbatim in the supplement'),
    ((66, 67), 'section heading and notation paragraph: notation rewritten in the paper (Section 2)'),
    ((134, 135), 'section heading'),
    ((294, 296), 'section heading'),
    ((299, 299), 'subsection heading'),
    ((333, 333), 'subsection heading'),
    ((409, 410), 'subsection heading (forest section now headed in the supplement)'),
    ((431, 436), 'section heading and definitions of the arithmetic section: rewritten in the paper'),
    ((440, 449), 'mechanism paragraph and section heading: rewritten in the paper; the scope remark is kept'),
    ((490, 501), 'open questions and provenance headings: open questions rewritten in the paper'),
    ((507, 511), 'appendix headings'),
    ((653, 657), 'appendix heading'),
    ((778, 781), 'appendix heading'),
    ((889, 890), 'subsection heading'),
    ((900, 904), 'appendix heading and the sentence describing the appendix as an outline (the complete proof is now reproduced)'),
    ((946, 948), 'appendix heading (kept in the paper)'),
    ((1107, 1108), 'appendix heading'),
    ((1163, 1164), 'appendix heading'),
    ((1189, 1191), 'bibliography opening'),
    ((1217, 1219), 'bibliography closing'),
]


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--list', action='store_true', help='print every line that is not carried verbatim')
    args = ap.parse_args()
    base = subprocess.run(['git', '-C', str(ROOT), 'show', f'{BASE_COMMIT}:paper/friend-tree-exponent.tex'],
                          check=True, capture_output=True, text=True).stdout.split('\n')
    new_lines = set()
    for p in NEW:
        new_lines.update(p.read_text(encoding='utf-8').split('\n'))
    carried = rewritten = missing = 0
    for n, line in enumerate(base, 1):
        if not line.strip():
            continue
        if line in new_lines or renames(line) in new_lines or main_renames(line) in new_lines:
            carried += 1
            continue
        reason = next((r for (a, b), r in REWRITTEN if a <= n <= b), None)
        if reason:
            rewritten += 1
            if args.list:
                print(f'rewritten  {n}: {reason}: {line[:90]}')
        else:
            missing += 1
            print(f'MISSING    {n}: {line[:140]}')
    print(f'{carried} baseline lines carried verbatim (after reference renames), {rewritten} in documented rewritten ranges, {missing} missing')
    print('RESULT: PASS' if missing == 0 else 'RESULT: FAIL')
    return 0 if missing == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
