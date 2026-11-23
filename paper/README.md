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

## Replaying the headline certificates offline

From the repository root, these commands regenerate the two motif arrays and verify them and the independent degree-band construction. They require GCC with OpenMP, G++ with C++17 and undefined-behaviour sanitizer support, GMP development headers and gzip. The full seven-suite run also installs pinned Python packages.

```bash
set -e
out=$(mktemp -d)
gcc -O3 -Wall -Wextra -fopenmp verify/motif-bracket/cert_build.c -lm -o "$out/build"
gcc -O3 -fopenmp -fsanitize=undefined -fno-sanitize-recover=all \
  verify/motif-bracket/cert_verify_gmp.c -lgmp -o "$out/verify"
"$out/build" 1 1133 2000 "$out/upper.bin"
"$out/build" 2 113 200 "$out/lower.bin"
sha256sum "$out/upper.bin" "$out/lower.bin"
"$out/verify" "$out/upper.bin"
"$out/verify" "$out/lower.bin"
gzip -dc verify/banded-upper/recipe.txt.gz > "$out/recipe.txt"
g++ -std=c++17 -O2 -fsanitize=undefined -fno-sanitize-recover=all \
  verify/banded-upper/check.cpp -o "$out/banded"
"$out/banded" verify/banded-upper/geometry.txt "$out/recipe.txt" 5657 9983
```

Appendices A and B give the reference digests, header requirements, integer bounds and expected counts. The seven-suite command also runs the separately written degree-band verifier.
