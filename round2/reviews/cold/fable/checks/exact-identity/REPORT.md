RESULT: IDENTITY HOLDS

# Exact check of the rival fringe identity

All arithmetic is `fractions.Fraction`; the run is single-threaded, offline, and
uses only the Python standard library. The only code file is
`identity_check.py`.

## Command

```
cd checks/exact-identity
python3 identity_check.py
```

## Real output

```
reachable marked trees by vertex count: {3: 1, 4: 2, 5: 6, 6: 18, 7: 61, 8: 208, 9: 743}
total trees: 1039
checks (tree,f): 5195
mismatches: 0
RESULT: IDENTITY HOLDS
all-f coefficient check: mismatching trees 0 / 1039
error 1 [global]: failing checks 5160 / 5195
  first failing f=f1 tree=O(L1)[1/2|L1|]
error 1 [rb_only]: failing checks 4010 / 5195
  first failing f=f1 tree=O(L1)[1/2|L0|[1/2|L1|]]
error 2 [global]: failing checks 3646 / 5195
  first failing f=f2 tree=O(L1)[1/2|L0|[1/2|L1|]]
error 2 [rb_only]: failing checks 2653 / 5195
  first failing f=f2 tree=O(L1)[1/2|L0|[1/2|L0|[1/2|L1|]]]
error 3 [global]: failing checks 2078 / 5195
  first failing f=f4 tree=O(L2)
error 3 [rb_only]: failing checks 2064 / 5195
  first failing f=f4 tree=O(L1)[1/2|L1|]
```

## What the numbers mean

- Reachable marked trees from the 3-path with at most 9 vertices: 1039 in total
  (1, 2, 6, 18, 61, 208, 743 for sizes 3, 4, 5, 6, 7, 8, 9). States are
  deduplicated by the full marked rooted-tree structure, so marks are kept.
- The five required test functions give 5 x 1039 = 5195 exact checks. There are
  0 mismatches; the identity holds on every enumerated tree for all five.
- Extra check: for a fixed tree, LHS - RHS is linear in `f`. Comparing the
  coefficient of every fringe key is therefore equivalent to checking the
  identity for every function `f` on the fringes that occur, not only the five
  chosen ones. That check also passes on all 1039 trees (0 mismatches).

Tree notation used in the output: `O(Lk)` is the centre with `k` leaves;
`[m|Lk|...]` is a core with mark `m`, `k` leaves, and nested core children.

## Planted errors

Two readings of the error instructions are reported, because the centre
boundary terms are written with `p_o` / `J_{p_o}` rather than with the `R`/`B`
notation.

- `rb_only`: the replacement is made only inside `R` and `B`; the centre
  boundary terms stay correct. This is the literal reading of (i) (only
  `p_S` occurs) and (iii) (only `L0` and `B` are named).
- `global`: the same replacement is also made in the analogous centre boundary
  terms (`p_o`, `J_{p_o}`).

| planted error | fails, `rb_only` | fails, `global` |
|---|---|---|
| (i) replace `(a_Y - p_S)` by `a_Y` | 4010 / 5195 | 5160 / 5195 |
| (ii) replace `f(Y+)` by `f(Y)` | 2653 / 5195 | 3646 / 5195 |
| (iii) newborn mark `1/(D+1)` in `L0` and `B` | 2064 / 5195 | 2078 / 5195 |

The first failing tree for each case is printed in the output above.

## What was implemented

- `gen_events` enumerates every event of the whole tree: each leaf ringing
  (rate 1 each, giving its core a leaf), and each core ringing (rate 1) with
  every neighbour choice: a leaf (promotion with mark `1/D_u`), a core child or
  the parent (that vertex gets a leaf), and for the centre no parent choice.
- `build_fringe` builds `S_v` for every core `v != o`, with the parent slot
  counted in the degree.
- `internal_events` implements `L0` directly on a fringe, including the child
  root ringing and picking its parent inside `S`; the root picking the slot
  contributes nothing and is omitted.
- `rhs` implements `R`, `B`, and the centre boundary terms exactly as stated.
- `lhs` applies every whole-tree event and recomputes `Sigma_f`.
- The five functions are exactly those requested, including the md5-based
  pseudo-random `f5` on the canonical string with marks and leaf counts.
