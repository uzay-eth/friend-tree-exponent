# The mean core size is not P-recursive

This folder checks the finite steps of Theorem PREC of the [paper](../../paper/friend-tree-exponent.pdf): the sequence $a_n=\mathbb E C_n$ satisfies no linear recurrence with polynomial coefficients, so its generating function is not D-finite. The proof is in Appendix PREC of the paper. It says nothing about whether the exponent $\mu$ is rational.

```sh
python3 -I verify/mean-arithmetic/check_window.py 17 1200
```

Standard library only; about 16 seconds and 33 MB. It ends with `RESULT: PASS`.

- **Part A.** The exact law of $T_n$ over unlabelled trees for $n\le17$: total probability one, the prime support of the denominators of $a_n$, the five size-seven probabilities of Lemma PREC.3, the near-star recurrences against the exhaustive law, and the window law of Lemma PREC.2 directly from the exhaustive law at $p=5,11,13$ (at $p=7$, $\Lambda_7=0$ and the law shows the predicted extra cancellation).
- **Part B.** The near-star recurrences run forward for all 190 primes $17\le p\le1200$: condition (C1) and the residues of Proposition PREC.4, and, for the 98 primes $p\equiv2\pmod3$ among them, $x_a(p)\equiv-1/10$ and $\Lambda_p=-4$.
- **Part C.** The identities for the broom solution, the telescopers $R_a$ and $R_c$, the Gosper polynomial $G$, the constants, and the Wilson step.

These are checks of the finite statements, not a substitute for the proof.

**Provenance.** The proof was first written as research notes, not part of the repository: the window bound, window law and near-star chain (SHA-256 `619969f1b929ac1c8d75fc8764954f6dc6f42988ab59f767d458e1b0b4764df2`), the near-star residues for $p\equiv2\pmod3$ (`23a8b1fc537bab71f96abe3132c0cf7442d24406e986ae422aa7f3d33aac0138`) and an alternative residue derivation (`0fc654a9c2fcd1771b5698d7929150da6baa7af5c963c255501d2c8b872bcbf6`). Appendix PREC is a self-contained write-up with the same scope. A GPT-family referee checked the full chain, including the classification, the grouped transitions, the Wilson step and the recurrence obstruction, and found that it holds. No human has refereed it.

**Numbering.** The checker's comments and output use the numbering of the earlier single-file manuscript (commit `6fee4c4`): Theorem 15 is now Theorem PREC, and Appendix E with its Proposition E.1, Lemmas E.2–E.3 and Proposition E.4 is now Appendix PREC of the paper with Proposition PREC.1, Lemmas PREC.2–PREC.3 and Proposition PREC.4. The checker bytes are unchanged, so the digest pinned in `paper/checks/check_literals.py` still applies.
