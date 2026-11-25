# Building the paper

[Read the Round 3 paper (PDF)](friend-tree-exponent.pdf) · [LaTeX source](friend-tree-exponent.tex) · [Submission sheet](submission-sheet.md) · [Round 2 paper (archived PDF)](https://github.com/uzay-eth/friend-tree-exponent/blob/c250a1f2fdd53ced91b7b485e94ac8a16106b9c7/paper/friend-tree-exponent.pdf) · [Round-one paper (archived PDF)](https://github.com/uzay-eth/friend-tree-exponent/blob/a0972085ad72ad5685ae7b424988648d1965e73c/paper/friend-tree-exponent.pdf)

This is the Round 3 version of 25 November 2026. It adds Theorem 12 (strict stochastic order of the limiting degree shares by birth label, at every threshold), Appendix D (its proof route, finite checks and pinned sources) and a supplementary mean result for growing exact degrees. The Round 2 version stays available at commit `c250a1f`.

The paper builds with pdfLaTeX from a single source file. It has an embedded bibliography and needs the standard LaTeX, AMS, Latin Modern, geometry, mathtools, booktabs, xcolor and hyperref packages. No external figures, shell escape or network access are needed for the build.

From the repository root, unpack the [arXiv bundle](friend-tree-exponent-arxiv.tar.gz) into an empty directory and run

```sh
build=$(mktemp -d)
tar -xzf paper/friend-tree-exponent-arxiv.tar.gz -C "$build"
cd "$build"
export SOURCE_DATE_EPOCH=1795564800 FORCE_SOURCE_DATE=1 TZ=UTC
pdflatex -halt-on-error -interaction=nonstopmode friend-tree-exponent.tex &&
pdflatex -halt-on-error -interaction=nonstopmode friend-tree-exponent.tex &&
pdflatex -halt-on-error -interaction=nonstopmode friend-tree-exponent.tex
```

The same commands work with a copy of the checked-in TeX file in place of the archive. The bundle contains only `friend-tree-exponent.tex`, byte-identical to the checked-in file. It was made from the repository root with GNU tar 1.35 and gzip 1.12:

```sh
src=$(mktemp -d)
cp paper/friend-tree-exponent.tex "$src/"
(cd "$src" && tar --format=gnu --sort=name --mtime=@1795564800 --owner=0 --group=0 --numeric-owner \
  --mode='u=rw,go=r' -cf - friend-tree-exponent.tex | gzip -n -6) > paper/friend-tree-exponent-arxiv.tar.gz
```

The checked-in files were built this way with pdfTeX 1.40.25 (TeX Live 2023/Debian), three passes from the unpacked bundle in an empty directory. Two clean builds gave identical bytes. SHA-256: TeX `f87a300608376e609e2065096a5eb7e592bf2b93801d25b57b0d1044e26b3311`, bundle `1be04e3d7cab32549118803029a3de89e963d489401cf2ce1db821953583c317`, PDF `3b2d692d6008d5556d06fecf782dc7ca45d77f9584ca62b4281743364f3d20e1` (44 pages). Another TeX distribution can produce a different but equivalent PDF.

The [theorem map](../docs/proof-map.md) links each result to its argument and checker. The [editorial checks](checks/README.md) compare the printed numbers with pinned records and validate repository links.
