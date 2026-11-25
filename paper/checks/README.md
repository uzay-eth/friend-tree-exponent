# Editorial checks for the paper and supplement

These checks compare the printed documents with pinned sources. They are separate from mathematical verification.

## Quick enclosure check

Run from the repository root:

```sh
python3 -B paper/checks/check_literals.py --enclosure
```

This uses only the Python standard library and works in a shallow clone or a source copy without Git. It does the following:

- It authenticates the working-tree Gold2 receipt against its historical SHA-256. It then recomputes the receipt's two exact endpoints as nominal rate minus or plus residual divided by the minimum weight.
- It runs the authenticated light checker in `verify/refined-lower/` on its pinned receipt and per-class file. That checker recomputes the refined lower endpoint class by class.
- It checks that Theorem ENC (`th:gold2`) prints the refined lower and the Gold2 upper endpoint exactly, and that the theorem and the abstract print the outward decimals.
- It checks that the supplement restates Theorem ENC exactly as the paper does.

It does not rerun the drift certificate or check the proof.

## Full editorial checks

Run these from the repository root:

```sh
python3 -B paper/checks/check_literals.py --all
python3 -B paper/checks/check_links.py
python3 -B paper/checks/check_transcriptions.py
python3 -B paper/checks/check_proof_parity.py
python3 -B paper/checks/test_checks.py -v
```

The link checker and the test suite use `markdown-it-py` 3.0.0. If it is absent, install it in a private virtual environment and run the commands with that environment's Python:

```sh
work=$(mktemp -d)
python3 -m venv "$work/venv"
"$work/venv/bin/python" -m pip install 'markdown-it-py==3.0.0'
```

### Literals (`check_literals.py`)

It scans both TeX files and lists each occurrence with its exact rational value and a location in the pinned archive.

- **What it matches:** slash and TeX fractions, signed decimals, and decimal scientific notation.
- **What it skips:** plain integers, symbolic quotients, comments, link destinations, citation arguments, the bibliography and typesetting metadata. It also skips the coordinates of the paper's figure, table column widths, and the transcribed proofs in the supplement, which `check_transcriptions.py` checks instead.
- **Where the values come from:**
  - Decimal endpoints are checked against the exact Gold2 and stationary-marginal receipts by directed `Fraction` rounding.
  - Several constants are derived from the pinned renewal, simulation and independent early-hub receipts.
  - The constants of Appendix PREC come from the authenticated `verify/mean-arithmetic/check_window.py`, and from an exact recomputation of the residue calculation that asserts its final values.
  - The refined-endpoint constants come from the authenticated files in `verify/refined-lower/`.
  - The promotion probabilities of the paper's Figure 1 are recomputed exactly.
- **Options:**
  - `--tex FILE` (repeatable) audits other TeX inputs instead.
  - It rejects reversed square-bracket pairs and ascending comparison chains around `\mu`.

Matching a value is otherwise a lexical check. It does not prove that the value's use, inequality direction or attribution in the documents is correct.

### Links (`check_links.py`)

It checks relative paths, Markdown headings and explicit anchors in tracked Markdown, including inline links, images, reference definitions and HTML links.

- **Pinned links:** archive-pinned GitHub links in Markdown are resolved against the pinned Git tree. In both TeX files, `\proofsrc` calls are resolved against the tree of `a0972085ad72ad5685ae7b424988648d1965e73c`, and `\ordersrc` calls (the stochastic-order sources) against the tree of `62081cad01112c42c2011ac3eed11aa3aa665e09`.
- **Other trees:** Markdown links pinned to that commit, to `c250a1f2fdd53ced91b7b485e94ac8a16106b9c7` or to `56f2d807ba63d19e9cb0e56d7678ae50843d3f6a` are checked against those trees. Links to other commits are treated as external.
- **Failures:** any other link into this repository from a TeX file fails. A relative link from one PDF to the other must name a file next to the TeX source.
- **Limits:** it does not fetch external URLs, and it does not diagnose a missing fragment in a non-Markdown file.
- **Options:**
  - `--include path/to/new.md` checks a new Markdown document before staging.
  - `--allow-existing` fails only on regressions from the pinned baseline. Omit it for the release gate.

### Transcriptions (`check_transcriptions.py`)

The supplement reproduces four reviewed Markdown proofs between `% BEGIN TRANSCRIPTION` and `% END TRANSCRIPTION` markers. Each marker records the source path, line range, SHA-256 and commit.

For each region the script does three things:

- it checks the source digest;
- it requires every formula of the source lines to appear in the TeX region in order and byte for byte;
- it requires the words of the prose to agree once markup is removed.

A negative control in the test suite changes one formula and one word.

### Proof parity (`check_proof_parity.py`)

This check needs Git history. It reads the earlier single-file manuscript at commit `6fee4c4d2464878ad06fff03c749c6fbede9baad`. It then requires every non-blank line to appear verbatim in the paper or the supplement, after the listed reference renames. The only exceptions are the documented ranges of rewritten front matter, introduction, open questions and headings. `--list` prints every exception with its reason.

The pins are `a0972085ad72ad5685ae7b424988648d1965e73c` for `\proofsrc` and `62081cad01112c42c2011ac3eed11aa3aa665e09` for `\ordersrc`. See the [supporting-results verifier guide](../../round2/verify/README.md) for replay requirements and scope.
