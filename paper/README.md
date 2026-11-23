# Building the paper

From the repository root, run

```bash
cd paper
SOURCE_DATE_EPOCH=1795392000 tectonic friend-tree-exponent.tex
```

This writes [friend-tree-exponent.pdf](friend-tree-exponent.pdf) beside the [LaTeX source](friend-tree-exponent.tex). The release uses Tectonic 0.17.0. The first build downloads the TeX packages it needs. The figure comes from [the simulation data and plotting script](../simulations/README.md) and is included from its checked-in location at `../simulations/local_slopes.png`.

For a two-file source bundle, place the LaTeX file and `local_slopes.png` together in one directory and run the same Tectonic command there. The bibliography is part of the LaTeX source. Building the paper uses the saved figure and requires neither a simulation run nor a certificate search.
