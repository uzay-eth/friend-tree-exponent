# Pinned Fable Claim 3 replay

The analytic audit is at `attacks/astra-on-fable-3.md`. From `[local path]`, run:

```sh
timeout 75 taskset -c 13 nice -n 10 /usr/bin/time -f 'elapsed=%es' python3 -Bu round2/attacks/fable3/replay.py
```

Result (15.07 s): 3,628,800 histories, 1,842 marked rooted shapes, common hub denominator `16724654043955200000000000000000`. The pinned grid and both size-11 output streams matched byte-for-byte under UBSan. Exact inequalities verified all board marginal bands and six adjacent floors. The `2−3` lower numerator was `4123373687076348603262692816576`, strictly above the `.2465` floor with that denominator. The replay uses Fable's pinned enumerator; its independent portion is the rational check against the board's cutoffs. The proof audit separately rederives the bounds and checks the arithmetic and signed aggregation.
