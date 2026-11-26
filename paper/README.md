# Building the paper and its supplement

[Paper (PDF)](friend-tree-exponent.pdf) · [Supplement (PDF)](friend-tree-exponent-supplement.pdf) · [Source bundle](friend-tree-exponent-source.tar.gz) · [Statement map](../docs/proof-map.md) · [Editorial checks](checks/README.md)

There are two documents, both dated 25 November 2026:

- **The paper,** [`friend-tree-exponent.pdf`](friend-tree-exponent.pdf). It states the results, gives the arguments behind them, and proves Theorem PREC in full in its Appendix PREC.
- **The supplement,** [`friend-tree-exponent-supplement.pdf`](friend-tree-exponent-supplement.pdf). It has the same title and author. It contains the complete proofs of every other result, in three parts: growth, the birth-order module, and further material.

The documents cite each other by result tags such as Theorem ENC, by section tags such as Supplement §ORD, and by page. Subsections are lettered (§ENC.d) and lemmas, propositions and equations are numbered within their tagged section (Lemma ENC.1), so the two kinds of number cannot be confused. Their cross-links are relative links to the other PDF, so keep both files in the same folder.

Both build with pdfLaTeX from two TeX files with embedded bibliographies. They need the standard LaTeX, AMS, Latin Modern, geometry, mathtools, booktabs, array, graphicx, xcolor and hyperref packages; `xr-hyper` is part of hyperref. No external figures, shell escape or network access are needed. The documents read each other's `.aux` files, so [`build.sh`](build.sh) runs both through three rounds in an empty directory, with a fixed `SOURCE_DATE_EPOCH`. It fails on any undefined reference.

From the repository root:

```sh
build=$(mktemp -d)
tar -xzf paper/friend-tree-exponent-source.tar.gz -C "$build"
bash "$build/build.sh" "$build/out"
```

The bundle contains exactly `build.sh`, `friend-tree-exponent-supplement.tex` and `friend-tree-exponent.tex`, byte-identical to the checked-in files. Running `bash paper/build.sh OUT_DIR` on the checked-in files gives the same result. The bundle is made from the repository root with GNU tar and gzip:

```sh
src=$(mktemp -d)
cp paper/build.sh paper/friend-tree-exponent.tex paper/friend-tree-exponent-supplement.tex "$src/"
(cd "$src" && tar --format=gnu --sort=name --mtime=@1795564800 --owner=0 --group=0 --numeric-owner \
  --mode='u=rw,go=r' -cf - build.sh friend-tree-exponent-supplement.tex friend-tree-exponent.tex | gzip -n -6) \
  > paper/friend-tree-exponent-source.tar.gz
```

The checked-in PDFs were built this way with pdfTeX 1.40.25 (TeX Live 2023/Debian), GNU tar 1.35 and gzip 1.12. They were built once from the checked-in files and once from the unpacked bundle, and the PDFs were byte-identical. Another TeX distribution can produce a different but equivalent PDF.

The bundle is a source bundle for these two documents. It has not been tested as an arXiv upload. arXiv processes one main document per submission, so distributing the supplement there would need a different layout, for example a separate submission or an ancillary file.

The supplement reproduces four reviewed Markdown proofs (the three- and five-vertex reductions and the complete stochastic-order proof). They were converted mechanically. Every formula keeps its content; the only changes inside formulas are documented layout commands that scale a wide display or break it across lines. `python3 -B paper/checks/check_transcriptions.py` compares each formula (after removing those layout commands), the words of the prose and the list numbers with the pinned source; three framed review notes in the stochastic-order proof are reported as additions. `python3 -B paper/checks/check_proof_parity.py` checks that every line of the earlier single-file manuscript (commit `6fee4c4`) reappears in the paper or the supplement, apart from a listed set of rewritten introductions and headings; it reports separately the lines carried verbatim, those changed only by reference renames, and those with a documented wording edit.

## Earlier numbering

Before this edition the manuscript was a single 55-page file with numbered theorems (commit `6fee4c4`). Some code comments and historical notes in `verify/` still use that numbering. The correspondence is:

| Old | Tag | Old | Tag |
|---|---|---|---|
| Theorem 1 | EXP | Theorem 9 | MASS |
| Theorem 2 | ENC | Theorem 10 | HUB |
| Theorem 3 | KER | Theorem 11 | MOM |
| Theorem 4 | AMP | Theorem 12 | ORD |
| Theorem 5 | DEG | Theorem 13 | EARLY |
| Theorem 6 | SPA | Theorem 14 | FOR |
| Theorem 7 | DRF | Theorem 15 | PREC |
| Theorem 8 | DEN | | |

| Old | Now |
|---|---|
| Appendix A | supplement §EXP |
| Appendix B (B.1–B.3) | supplement §CORE.d (Lemmas CORE.1–CORE.3) |
| Appendix B (B.4) | supplement §SPA (Lemma SPA.1) |
| Appendix C (C.1–C.2) | supplement §CERT |
| Appendix C (C.3) | supplement §EARLY |
| Appendix D | supplement §ORD |
| Appendix E (E.1–E.4) | paper Appendix PREC (Proposition PREC.1, Lemmas PREC.2–PREC.3, Proposition PREC.4) |
| Appendix F | supplement §ENC.d |
| Lemma F.1 | Lemma ENC.1 |
| Appendix G | supplement §VER |
| §6 | supplement §NOTES |
| §8 | supplement §PROV |
