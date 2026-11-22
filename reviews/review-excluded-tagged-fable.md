VERDICT: GAP
Scope, in my own words: this is a partial review that the lead stopped early. GAP is about the strict claim (an upper exponent below 17/30): no positive reserve for the 17/30 table has been produced or checked, and I did not replay that table. In the parts I did review (the construction, the centering identity, the reserve algebra, the arithmetic, and a tree-level test of the small capped pilot table) I found no error.

# Cold review of the tagged upper certificate (`bounds/TAG.md`, `RESERVE.md`, `tag.cpp`)

Auditor: Claude Fable 5.2 (claude-code harness), child thread of the Astra lead, 2026-11-22.
The lead called an early stop while this review was under way. This file records what was actually
checked, how, and what was not.

## What this review supports, and what it does not

| Statement | Status after this review |
|---|---|
| The construction (uniform mean over pooled neighbours, entry and stop rules, tag dynamics, centering, reserve shift) is mathematically sound | Checked by reading and hand derivation; no error found |
| The centering identity is justified under the correct filtration | **Yes**, see §1 |
| A tagged table that passes `tag.cpp` satisfies the drift inequality on real trees | Tested on the small **capped** pilot table only: 124,361 vertex checks, 0 violations (§5) |
| Banded core-degree labels, `--joint`, `--four-stop`, tags at c = 4 | Read only. Not exercised by any run of mine |
| `PASS tagged upper 17/30 … 170377363 rows` (`tag-threshold.log`) | **Not replayed by me.** The log was also not produced by the hashed source (§6, item 2) |
| A strict upper exponent below 17/30 | **Unverified** by the author (their message) and by me. Without a positive reserve the certificate gives limsup ≤ 17/30, which does not exclude 17/30 |

Files (sha256):

    f5666e2e8da4e652f5341b938702d523078e0d75729ab7fbf0b48ac8fd8cfc8f  tag.cpp      (final; read in full)
    5cfe9dd44f1a0c239874ec22bc3f5d9f0d90c66c0c6cee8e4278cb987e98cd9a  tag.cpp      (earlier snapshot; used for my pilot dump)
    a34422b8df1e3fc1f7888e9567486cec46300721684a603d9c650f6c9129530e  transfer.cpp (included by tag.cpp; read in full)
    324a610486f059d199d7cdc08da8b153bf85e8b1e3a5cf50ecea192a9a6b7467  RESERVE.md
    db2f84d1a7d128118282fdf5fcdf985d8a13b888c072b4fbcf7dd92e0b97a8ff  TRANSFER.md
    fb2304d94b9356d6f0a4e3e29e0d711e9e080c2eb9a40daccd59f3fbf5287343  TAG.md       (final; I read the earlier 94ed5478… and the diff)

The two `tag.cpp` versions differ only in the helper `weighted`: the final one adds a 64-bit fast path
`product/d + (product%d > 0)`. For d > 0 this is the ceiling of a signed quotient, the same value as
the old `divide(x*n, d, true)`. I checked this by reading, not by a run.

Not read: any other audit report, the board, `check-tag-events.cpp`, `SENSITIVITY.md`,
`tag-upper.seed`. I read the one-line `tag-events.log` but did not examine or rerun that test.

## 1. Is the centering identity justified under the correct filtration? Yes.

The correct filtration is the physical one, F_n = history of the tree. No internal randomness is
needed, and I recommend stating the potential without it:

* For a core vertex u whose tag is on (2 ≤ c ≤ C, pool nonempty, D below the stop), every pooled
  neighbour j has a label t_j that it would carry if it were the tag: its entry label, then updated by
  its own gains and conversions. t_j is a function of the physical history of j since entry. The
  reservoir rule never selects a member later than its entry, so this label does not depend on any coin.
* Define pot(u) = (1/p) Σ_{j in pool} F(x_u, t_j), and the untagged weight when the tag is off.
  Φ_n = Σ_u pot(u) is F_n-measurable. This is exactly "the conditional average over the internal
  sampling given the physical history": the coins are independent of the physical process, so given the
  whole physical history the tag is uniform on the pool.
* The reservoir mixture in `destination()` is the statement that the uniform mean over a pool of size
  p+m equals the mean over the old members j of [p/(p+m)·F(x′, t′_j) + (1/(p+m))·Σ_entrants F(x′, t_e)].
  So n·E[Φ_{n+1} − Φ_n | F_n] = Σ_u mean_j (row-j expression with the **true** rates), exactly.
* In that expression the pooled credit Σ_i q_i and the constant b(x) are F_n-measurable numbers that do
  not depend on j, and mean_j q_j = (1/p) Σ_i q_i. Hence
  mean_j[(Σ_i q_i) a_j] = mean_j[(Σ_i q_i)(a_j − b) + p q_j b]. This is algebra on fixed numbers under a
  uniform average. It needs nothing beyond uniform weights 1/p and b independent of j.

Consequences that the write-up should say explicitly:

1. A single row inequality is **not** a bound on the conditional drift given the tag. Only the mean over
   j is. The proof therefore bounds E[Φ_{n+1} | F_n] for the averaged potential, never a drift conditional
   on the enlarged history with the coins. TAG.md says "averaged potential", which is right; the phrase
   "probabilities of sampling a tag" invites the wrong reading.
2. Bounding the two centered pieces separately (Σ_i q_i in [0, Σ caps] for the residual, q_j from the tag
   for the constant) ignores that q_j is one of the q_i. That is a relaxation, so it is valid.
3. Edge transfers do not depend on tags and pooled neighbours carry none, so the cancellation of
   TRANSFER.md is untouched by the averaging.

## 2. Uniform sampling, entry, drop and refresh transitions (read against `tag.cpp`)

No error found. Points I checked:

* **Members never leave the pool.** Cutoffs h_c are nonincreasing in c and degrees never decrease.
* **Entrants.** At a central conversion: every exact neighbour with d ≥ h_{c+1}, and the newborn when
  h_{c+1} = 2. At a neighbour's gain: that neighbour when d+1 = h_c. `check_mix` enforces
  keep + entrants = destination pool size for every transition.
* **Switching on.** A source row with pooled members is inactive only for c = 1 (or c > C, where the
  destination is inactive too). So the label `top` given to "previously pooled" members is only used for
  a neighbour pooled in layer 1, whose degree is at least h_1. The overflow guarantee min(T, h_1) holds
  for all three ways of reaching `top` (degree reached T, entered with d ≥ T, pooled in layer 1).
* **Stopping is permanent and consistent.** `mix` is built per row, not per L. This is safe because the
  stop is nonincreasing in c (U for c ≤ 3, min(U, four-stop) for c ≥ 4): a stopped source always has a
  stopped destination. The code relies on this silently; a parameter set with a larger stop at higher c
  would break it. A one-line guard would make it explicit.
* **Tag range of a row.** Tags run from degree `lower` (h_c, or the single stored mark when p = 1) to
  `top`. A kept tag always lies in the range because an exact tag degree is the true degree, which is at
  least any valid mark. `at()` throws otherwise.
* **Tag dynamics.** Gain rate in [d−k₊+1/D, d−k₋+1/D+(k₋−1)/2], conversion rate in
  [(d−k₊)/d, (d−k₋)/d]: correct (maximise d−j+(j−1)/2 over j). Capped mode: a saturated label never
  moves on conversion. Banded mode: bounds K+1, then b → 2b−2, then T; a conversion moves to the next
  bound, valid because the true core degree rises by one and 2b−2 ≥ b+1 for b ≥ 3; the terminal bound
  T > d is vacuous. Entry from a saturated coarse label uses the old degree as upper bound, valid also
  when the trigger is a gain.
* **`--joint`.** Σ_{i≠j} q_i(a−b) + q_j(a+(p−1)b) is the same sum. Dropping the cap of the largest mark
  bounds Σ_{i≠j} q_i whichever member is tagged; dropping an ordinary cap is allowed when the exact tag
  degree is below the largest mark, since then the tagged member cannot be the one holding it.

## 3. Reserve shift (RESERVE.md)

Correct, by hand:

* Upper. Rows are verified with provisional rewards ≥ actual + e, and the tail coefficient against the
  provisional pooled reward, so the actual summed drift is ≤ λΦ − e·Σ_u q_u. With Ψ = Φ + e·C and
  A·C = Σ_u q_u: AΨ ≤ λΦ ≤ [λM/(M+e)]Ψ because Φ ≤ M·C. Prefactor (f₀+e)/(inf f+e).
* Lower means. Ψ = Φ − e·C gives AΨ ≥ δΦ ≥ [δM/(M−e)]Ψ, also when Ψ < 0. Needs f₀ − e > 0.
* I recomputed the quoted numbers in exact fractions (run): 0.563·17611596558/17605492535 =
  4957664431077/8802746267500 with prefactor 15424346584/17605492535; 0.568·17605633804/17607347966 =
  625000000042/1100459247875; and the pilot value 0.5706·17525411848/17527651385 =
  6250000000293/10954782115625. All three agree with RESERVE.md and TAG.md.
* `certify_tags` uses M = max(all cells, ⌈A(1+λ)⌉) and e = min(provisional − actual) ≥ 0, as required.

The shift is only as good as e > 0 for the table in question. For the 17/30 table that number has not
been produced.

## 4. Exact arithmetic (read only)

* `slack`: each term ⌈(target − z)·rate⌉ with the high endpoint when target > z, compared with ⌊λz⌋.
  A nonnegative result implies the real inequality.
* Reservoir means are rounded up and multiply nonnegative rates; the constant ⌈b·p/d⌉ is rounded up; for
  an overflow tag it is ⌈b·p/d_top⌉ if b > 0 and 0 otherwise, the maximum over q_j in [0, 1/d_top].
* The final pass rechecks every cell with the provisional rewards, checks positivity, the tail
  condition A·λ(1+λ) ≥ provisional pooled reward, the slot/transfer dictionary, and the cell count.
* I did not find an overflow path: weights ≤ 16·10^10, at most 512 summands, products in `__int128`,
  UBSan in the build line.

## 5. Tree-level test on the pilot table (the only runs I made)

Table: `./tag 5706 10000 56 20 3 3 56 24 8 6 4 3 3` (capped tag, K = 3, tags at c = 2, 3, no joint
split, no early stop), dumped after `PASS` by a patched copy. The patched build printed the same lines
as `tag-pilot.log`.

My script `tag_tree_check.py` grows real trees and keeps, with its own implementation of the written
rules, the coarse state of TRANSFER.md (including the c = 1 witness) and the would-be tag label of every
pooled neighbour. After every move it checks all of them against the true tree, including the rate
intervals of TAG.md. On sampled trees it checks, in exact rationals and with the true rates, for every
core vertex u:

    Σ_events rate·Δpot(u) + q_u·(provisional newborn reward) + r_u·Σ_x g(u→x) − Σ_w r_w·g(w→u)  ≤  λ·pot(u)

with pot(u) the uniform mean of §1 and g the dumped transfers; it also asserts that both ends of every
edge see the same transfer.

| run | trees checked | vertex checks | tagged | pool > 1 | D ≥ H | violations |
|---|---|---|---|---|---|---|
| mixed, seed 1 | 429 | 5,883 | 1,034 | 61 | 0 | 0 |
| mixed, seed 11 | 4,112 | 77,059 | 13,312 | 1,420 | 0 | 0 |
| natural, seed 12 | 5,597 | 41,419 | 11,302 | 2,304 | 1,793 | 0 |

The largest ratio (amended local drift)/pot is 0.570600000, the certified exponent, so the pilot table
is tight on real trees. Negative controls, all caught: the potential taken at one member instead of the
uniform mean (first violation at tree size 24 in the second trajectory); the tag label not updated when
the tagged neighbour gains (first violation at tree size 33 in the first trajectory); and the same
table tested against 0.5705 (3,772 violations in 30 trajectories).

**Limits.** 43 distinct tag labels were visited. No tree hit a stopped tag (U = H in the pilot). This
run says nothing about banded labels, `--joint`, `--four-stop` or c = 4. My script implements the banded
and four-stop rules but they were never run against a table. I did not mutate `tag.cpp`.

## 6. Concrete gaps and open items

1. **Strictness is open.** The only 17/30 result on file is non-strict. To exclude 17/30 one needs a
   positive reserve on the 17/30 table (or a pass at a smaller exponent) from the frozen source and seed.
2. **Provenance of `tag-threshold.log`.** The hashed `tag.cpp` prints a line "verified at provisional
   newborn rewards; minimum reserve …" before `PASS`. The log has no such line, so it came from a
   different version of the checker. The 17/30 pass needs a replay with f5666e2e… and `tag-upper.seed`
   (0bcd2a54…) before it is quoted as verified by that source.
3. **Write-up.** Define the potential as the uniform mean of §1 and say that each member's label is a
   function of the physical history. State the monotone-stop assumption of §2.
4. **Inherited base.** Everything rests on the transfer model in `transfer.cpp` (witness of
   TRANSFER.md, edge transfers). I read it and my pilot test exercises it on trees (witness labels,
   transfer symmetry at both ends of every edge), but I did not audit the transfer certificates
   themselves (0.563 / 0.568), nor rerun `reproduce-transfer.sh`.
5. **Shared-machine note.** `reproduce-tag.sh` pins the first permitted CPU at nice 15. On this server
   that CPU is contended; a pinned baseline run of mine got 5–16 % of a core. Use `BOUND_CPU`.

## 7. Reproducing my runs

Scratch, not part of the submission: `attempts/strogatz-astra/branch-work/tag-scratch/`
(`make_tag_dump.py`, `tag_tree_check.py`, `tag_tree_check.log`).

    python3 make_tag_dump.py ../../bounds/tag.cpp
    g++ -std=c++17 -O2 -fsanitize=undefined -fno-sanitize-recover=all -I../../bounds tag_dump.cpp -o tag_dump
    TAG_DUMP=pilot.dump ./tag_dump 5706 10000 56 20 3 3 56 24 8 6 4 3 3          # 2.5 s, 11 MB
    python3 tag_tree_check.py pilot.dump 24,8,6,4,3,3 512 120 160 11 mixed       # not timed; all five runs finished within a few minutes in parallel
    python3 tag_tree_check.py pilot.dump 24,8,6,4,3,3 512 120 200 12 natural
    python3 tag_tree_check.py pilot.dump 24,8,6,4,3,3 512 30 120 5 mixed first_member     # must report violations
    python3 tag_tree_check.py pilot.dump 24,8,6,4,3,3 512 30 120 5 mixed none 5705/10000  # must report violations

`tag_dump.cpp` includes `transfer.cpp`, hence the `-I`. My own build was made beside copies of the
sources; I have not rerun this exact command line with `-I`.
