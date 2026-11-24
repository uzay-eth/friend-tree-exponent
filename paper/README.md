# Building the paper

[Read the paper (PDF)](friend-tree-exponent.pdf) · [LaTeX source](friend-tree-exponent.tex) · [Submission sheet](submission-sheet.md) · [Round-one paper (archived PDF)](https://github.com/uzay-eth/friend-tree-exponent/blob/a0972085ad72ad5685ae7b424988648d1965e73c/paper/friend-tree-exponent.pdf)

The paper builds with pdfLaTeX from a single source file. It has an embedded bibliography and needs the standard LaTeX, AMS, Latin Modern, geometry, mathtools, booktabs, xcolor and hyperref packages. No external figures, shell escape or network access are needed for the build.

From the repository root, unpack the [arXiv bundle](friend-tree-exponent-arxiv.tar.gz) into an empty directory and run

```sh
build=$(mktemp -d)
tar -xzf paper/friend-tree-exponent-arxiv.tar.gz -C "$build"
cd "$build"
export SOURCE_DATE_EPOCH=1795478400 FORCE_SOURCE_DATE=1 TZ=UTC
pdflatex -halt-on-error -interaction=nonstopmode friend-tree-exponent.tex &&
pdflatex -halt-on-error -interaction=nonstopmode friend-tree-exponent.tex &&
pdflatex -halt-on-error -interaction=nonstopmode friend-tree-exponent.tex
```

The same commands work with a copy of the checked-in TeX file in place of the archive. The [theorem map](../docs/proof-map.md) links each result to its argument and checker. The [editorial checks](checks/README.md) compare the printed numbers with pinned records and validate repository links.
