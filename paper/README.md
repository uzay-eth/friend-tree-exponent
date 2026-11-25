# Building the paper

[Paper (PDF)](friend-tree-exponent.pdf) · [LaTeX source](friend-tree-exponent.tex) · [Submission sheet](submission-sheet.md) · [Statement map](../docs/proof-map.md) · [Editorial checks](checks/README.md)

The paper is dated 25 November 2026. It builds with pdfLaTeX from a single source file with an embedded bibliography, and needs the standard LaTeX, AMS, Latin Modern, geometry, mathtools, booktabs, xcolor and hyperref packages. No external figures, shell escape or network access are needed for the build.

From the repository root, unpack the [source archive](friend-tree-exponent-arxiv.tar.gz) into an empty directory and run

```sh
build=$(mktemp -d)
tar -xzf paper/friend-tree-exponent-arxiv.tar.gz -C "$build"
cd "$build"
export SOURCE_DATE_EPOCH=1795564800 FORCE_SOURCE_DATE=1 TZ=UTC
pdflatex -halt-on-error -interaction=nonstopmode friend-tree-exponent.tex &&
pdflatex -halt-on-error -interaction=nonstopmode friend-tree-exponent.tex &&
pdflatex -halt-on-error -interaction=nonstopmode friend-tree-exponent.tex
```

The same commands work with a copy of the checked-in TeX file in place of the archive. The archive contains only `friend-tree-exponent.tex`, byte-identical to the checked-in file. It is made from the repository root with GNU tar and gzip:

```sh
src=$(mktemp -d)
cp paper/friend-tree-exponent.tex "$src/"
(cd "$src" && tar --format=gnu --sort=name --mtime=@1795564800 --owner=0 --group=0 --numeric-owner \
  --mode='u=rw,go=r' -cf - friend-tree-exponent.tex | gzip -n -6) > paper/friend-tree-exponent-arxiv.tar.gz
```

The checked-in files were built this way with pdfTeX 1.40.25 (TeX Live 2023/Debian), three passes in an empty directory, once from the TeX file and once from the unpacked archive; the two PDFs were byte-identical. SHA-256: TeX `401446f8a50817ef3be8ecbda8b1b9492cccb7a1a13429333920dd56b7fc1cf5`, archive `813b916e19c150b8462052dc4986cdea8e44bdf5e44f90c8d04e02551fc002a4`, PDF `76c2170ca020ef8540372529d52afc6266c5b1eaa8c810f1fdfb1d0f302e537c` (55 pages). The archive was made with GNU tar 1.35 and gzip 1.12. Another TeX distribution can produce a different but equivalent PDF.

The [statement map](../docs/proof-map.md) links each result to its argument and checker. The [editorial checks](checks/README.md) compare the printed numbers with pinned records and validate repository links.
