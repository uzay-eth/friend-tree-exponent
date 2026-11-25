# Editorial checks for the paper

## Quick enclosure check

Run from the repository root:

```sh
python3 -B paper/checks/check_literals.py --enclosure
```

This uses only the Python standard library and works in a shallow clone or a source copy without Git. It authenticates the working-tree Gold2 receipt against its historical SHA-256 and recomputes its two exact endpoints as nominal rate minus or plus residual divided by the minimum weight. It then runs the authenticated light checker in `verify/refined-lower/` on its pinned receipt and per-class file, which recomputes the refined lower endpoint class by class. It checks that Theorem 2 prints the refined lower and the Gold2 upper endpoint exactly, and that Theorem 2 and the abstract print their outward decimals. It does not rerun the drift certificate or check the proof.

## Full editorial checks

Run these from the repository root:

```sh
python3 -B paper/checks/check_literals.py --all
python3 -B paper/checks/check_links.py
python3 -B paper/checks/test_checks.py -v
```

The link checker and test suite use `markdown-it-py` 3.0.0. If it is absent, install it in a private virtual environment and run the commands with that environment’s Python:

```sh
work=$(mktemp -d)
python3 -m venv "$work/venv"
"$work/venv/bin/python" -m pip install 'markdown-it-py==3.0.0'
```

- **Literals:** `--all` lists each occurrence with its exact rational value and a location in the pinned archive. The scanner matches slash and TeX fractions, signed decimals, and decimal scientific notation. It skips plain integers, symbolic quotients, comments, link destinations, citation arguments, bibliography and typesetting metadata. It checks decimal endpoints against the exact Gold2 and stationary-marginal receipts by directed `Fraction` rounding. It also derives several constants from the pinned renewal, simulation and independent early-hub receipts, takes the constants of Appendix E from the authenticated `verify/mean-arithmetic/check_window.py` and an exact recomputation of the residue calculation that asserts its final values, and takes the refined-endpoint constants from the authenticated files in `verify/refined-lower/`. `--tex` accepts a draft fragment. It rejects reversed square-bracket pairs and ascending comparison chains around `\mu`. Matching a value is otherwise a lexical check, not proof that its use, inequality direction or attribution in the manuscript is correct.
- **Links:** checks relative paths, Markdown headings and explicit anchors in tracked Markdown; includes inline links, images, reference definitions and HTML links. Also resolves archive-pinned GitHub links in Markdown and `\proofsrc` calls in the manuscript against the pinned Git tree, and `\ordersrc` calls (the stochastic-order sources) against the tree of commit `62081cad01112c42c2011ac3eed11aa3aa665e09`. Markdown links pinned to that commit, to `c250a1f2fdd53ced91b7b485e94ac8a16106b9c7` or to `56f2d807ba63d19e9cb0e56d7678ae50843d3f6a` are checked against those trees; links to other commits are treated as external. Any other link into this repository from the manuscript fails. It does not fetch external URLs. A missing non-Markdown fragment is not diagnosed, because file formats use different fragment rules. Use `--include path/to/new.md` to check an explicitly provided new Markdown document before staging it. `--allow-existing` reports broken links also found in the pinned baseline, but fails only on regressions. Omit that flag for the release gate. `--tex` can point at assembled draft TeX; the checker requires at least one pinned source call.

The pin is `a0972085ad72ad5685ae7b424988648d1965e73c`; the pin for `\ordersrc` is `62081cad01112c42c2011ac3eed11aa3aa665e09`. The literal scanner masks the path argument of `\ordersrc` as it does for `\proofsrc`. These are editorial checks, separate from mathematical verification. See the [supporting-results verifier guide](../../round2/verify/README.md) for replay requirements and scope.
