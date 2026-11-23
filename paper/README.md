# Building the paper

From the repository root, run

```bash
cd paper
SOURCE_DATE_EPOCH=1795392000 tectonic friend-tree-exponent.tex
```

This writes [friend-tree-exponent.pdf](friend-tree-exponent.pdf) beside the [LaTeX source](friend-tree-exponent.tex). The release uses Tectonic 0.17.0. The first build downloads the TeX packages it needs. The figure comes from [the simulation data and plotting script](../simulations/README.md) and is included from its checked-in location at `../simulations/local_slopes.png`.

The [arXiv source archive](friend-tree-exponent-arxiv.tar.gz) contains the source and figure at its top level. From the repository root, the following commands build it in an empty directory with TeX Live.

```bash
build_dir=$(mktemp -d)
tar -xzf paper/friend-tree-exponent-arxiv.tar.gz -C "$build_dir"
cd "$build_dir"
for pass in 1 2 3; do
  pdflatex -interaction=nonstopmode -halt-on-error friend-tree-exponent.tex || exit 1
done
```

The same Tectonic command also works with the two files side by side. The bibliography is embedded in the source. The [submission sheet](submission-sheet.md) gives the title, author, plain abstract and other metadata.
