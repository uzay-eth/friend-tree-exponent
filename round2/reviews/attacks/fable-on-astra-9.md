VERDICT: HOLDS

# Fable on Astra claim 9 (certified bracket 0.5656581227 <= mu <= 0.5657979941, rung 3)

**Scope of the verdict:** The Claude-family lead rebuilt the lifetime tables from Astra’s per-core ceilings, checked the endpoint arithmetic and finite oracle tests, and read the combination, occupation and founder arguments. The fringe tests exercise selected finite trees and faults with Astra’s oracle; they are not a full row scan. The witness potentials, per-core ceilings and founder estimates are inherited. Replay of the full certificate and separate implementations of the component caps and donor tables were checked only by GPT-family children.

Attacker: the Fable lead (Claude Fable 5.2, a different family from the GPT claimant). Lead scripts and working logs are in Fable folder `round2/attacks/astra9/`; sections 8 and 9 credit two same-family GPT-6.2 Sol Ultrafast children. The earlier proof bracket here is the rung-3 result, distinct from the later board Fact with a tighter bracket.

## 1. What the claim consists of

Trial exponents lambda_L = 0.5656453, lambda_U = 0.56581165 with frozen potentials Phi (two files of 611,658,376
bytes; their sha256 equal Astra folder `round2/motif-cert/certificate.json` (74,226,396 base rows per witness before envelopes)). Three inequalities in ticks of 2^-40 per core:
(1) the signed pointwise row inequality sigma(lambda Phi - L Phi) >= F - e C - E with residual cap e and nonnegative
physical credits F; (4) an occupation bound F >= (u - L) T + v B giving a credit v per newborn; (6) a lifetime
bound Phi <= (l - L) Y_tot + y B + E_1 giving a ceiling y. In Laplace-normalised form at s -> mu these give
sigma(lambda - mu) >= (v - e/l)/y with l = 113/200.

## 2. Combination argument, read by hand

With [Z]_s = int e^{-st} E Z / int e^{-st} E C: [L Z]_s = s [Z]_s - Z(0)/int e^{-st} E C, the denominator diverges as
s decreases to mu, so [B]_s -> mu and [L Phi]_s - mu [Phi]_s -> 0 along a subsequence with [Phi]_s -> m.
(1) and (4): sigma(lambda - mu) m >= (u - mu) lim [T]_s + mu v - e >= mu v - e, using u >= mu and T >= 0.
(6): m <= (l - mu) lim [Y_tot]_s + mu y <= mu y, using l <= mu and Y >= 0. Hence, when mu v > e,
sigma(lambda - mu) >= (mu v - e)/(mu y) >= (v - e/l)/y. Signs and the direction of each use of l <= mu <= u are
consistent with the argument; the finiteness of the transforms for s > mu needs E T, E Y_tot = O(e^{mu t}), which follows from bounded
scores per set and the core-degree moment bound they cite. No flaw found.

In the occupation step each connected set of at most five cores is followed in time; a promotion gives the old
set in its new state plus the enlarged set, and every new set arises exactly once by deleting its youngest member;
gain rates of members are enclosed by [known part, known part + 1/2 per unobserved core slot]. Summing their row
inequality (3) over all sets gives (4) provided the sources are distinct physical credits. Read with the question
"is anything counted twice": the credit (2) is the gap between the maximised donor input used in the rows and
the actual input; the per-core charge of the higher families uses positive parts with worst-case multiplicities
(so the negative parts are unused and may be credited); for an edge, the signed row value plus the donor credits and the separate row credit is at
most max(0, R_max), below the positive edge cap. No double counting found. This is reading, not recomputation.

## 3. Recomputed exactly with own code

`python3 -B lifetime_own.py` (exact integers and rationals):

    upper: own lifetime ceiling y = 2655887184018 (claimed 2655887184018) EQUAL; states 82169; all rows re-verified with integers
       first calculation: lambda - (v - e/l)/y = 1698065802472818376961/3001152517940340000000 = 0.565804567519  (EQUAL to the stated fraction)
       with the additional credit 17458365: 1698046074520368376961/3001152517940340000000 = 0.565797994061  (EQUAL to the stated fraction)
       outward ten-digit endpoint 0.5657979941 (claimed 0.5657979941)
    lower: own lifetime ceiling y = 2654889165321 (claimed 2654889165321) EQUAL; states 82169; all rows re-verified with integers
       first calculation: lambda + (v - e/l)/y = 1696971095548923704669/3000024756812730000000 = 0.565652363933  (EQUAL to the stated fraction)
       with the additional credit 15289028: 565662790716854568223/1000008252270910000000 = 0.565658122753  (EQUAL to the stated fraction)
       outward ten-digit endpoint 0.5656581227 (claimed 0.5656581227)
    planted error 'one-endpoint' (upper): y = 2242943518112, failed endpoint inequalities: 82169
    planted error 'round-down' (upper): y = 2655887183965, failed endpoint inequalities: 82169

The lifetime table was rebuilt from the text of their section 3 (recursion over (D, c), both end points of the
gain-rate interval [D - c, D - c + c/2], outer constant ceil(M/l)) from their per-core ceiling file, and the row
inequality was then re-verified with integers in a separate pass. The per-core ceiling m(D, c) itself is theirs.

## 4. The pointwise witness inequality on new inputs

Their own audit uses 206 small enumerated trees. I generated random planted fringes with an independent simulation
of the friend-tree dynamics (root with an inert parent of infinite degree; `random_fringes.py`) and passed each to
their exact rational tree oracle, which checks on the whole tree: the transfer identity, sigma(L Phi - lambda Phi)
plus both physical credits <= cores x e, the comparison bounds, and every realised local row against the corner
caps returned by their compiled exact probes. `./run_all.sh`:

    upper: 12000 trees of 5..40 cores: 497790 actual events, 813428 local rows, 259722 physical row types -> PASS
           400 trees of 60..200 cores:  97545 actual events, 157894 local rows,  64791 physical row types -> PASS
           3000 trees of 1..12 cores (young fringes) -> PASS
    lower: the same three samples -> PASS, PASS, PASS
    15400 distinct trees, each tested against both witnesses; 617055 event terms per witness
    largest leaf count 5034280, largest number of children 15

In the separate slack run reported by the lead (seed 21, 4,000 trees of 5–40 cores per witness), the largest `(drift + credits)/(cores × e)` was +0.625 for the lower witness (62 positive trees) and −0.007 for the upper witness. These are sample maxima; log in Fable folder `round2/attacks/astra9/expect_test_seed21.log` (command on its first line).

## 5. Planted errors

* Theirs, on my inputs: `--plant fy` fails at the 4th tree (transfer identity), `--plant double-gap` fails
  (a credited row exceeds its corners).
* Mine: trial exponent moved towards the reference estimate 0.56572366. Upper witness with lambda lowered by 2e-5, 5e-5, to
  0.56572366 and by 2e-4: the whole-tree bound fails at the 5th, 2nd, 2nd, 2nd tree. Lower witness raised by
  2e-5, 5e-5, to 0.56572366, by 2e-4: fails at the 5th, 2nd, 2nd, 1st tree. These results establish sensitivity to the specified changes on these inputs, including the 2e-5 changes to the trial exponents.
* Lifetime table built with one end point of the gain interval, or rounded down: rejected (section 3).

## 6. Finite-sample measurements of the certificate’s quantities

On 160 simulated fringes of 60–219 cores per witness, with potential, drift and baseline credit evaluated by Astra’s oracle; see Fable folder `round2/attacks/astra9/expect_test_seed7.log` (command on line 1):

    upper: m = 1.2401e12 ± 3e9 (μy = 1.5025e12); G = 4.54e7 ± 0.35e7 (baseline μv = 1.08e7)
           sum(L Phi)/sum(Phi) = 0.565722 ± 0.0000036
    lower: m = 1.2396e12 ± 3e9 (μy = 1.5019e12); G = 4.45e7 ± 0.32e7 (baseline μv = 1.07e7)
           sum(L Phi)/sum(Phi) = 0.565725 ± 0.0000036

Here μy is about 21% above the measured mean m, and measured baseline G is about 4.2 times baseline μv. Relative to μv with the extension (upper 2.07e7, lower 1.94e7), measured baseline G is about 2.2 times as large; G does not include the extension’s additional sources. The two ratios are consistent with the reference estimate 0.56572366. Substituting the sample means gives illustrative endpoints about 0.5656811 and 0.5657751, not certified bounds or estimates of the Laplace-limit means.

## 7. Against the independent number, and what rung 3 needs

The reference estimate 0.56572366 has a reported statistical error below 1e-8 in selected Fable estimator windows; this is not a full error bound.
It lies inside the bracket, 6.55e-5 above the lower and 7.43e-5 below the upper end point; also inside the base
bracket without credits, [0.5656449338, 0.5658122517].

For the window [0.565615, 0.565809] of rung 3: the lower end point needs no credit at all (the base value is
inside). The upper end point needs v >= 7,313,514 ticks; the claim has 36,544,047, the first calculation
19,085,682, and each of the two largest families alone (8,148,713 and 7,709,975) suffices.

## 8. Inside the GPT family: Sol Ultrafast child (replay and recomputation)

Child report (Fable folder `round2/children/astra9-gpt/REPORT.md` and tracked `round2/checks/astra9-replay-gpt-family.md`; first line
`RESULT: REPRODUCED AND RECOMPUTED`):

* Second-party replay of the claimant's full script (Astra folder `round2/occupation-tail/replay.py`, using its `--snapshots` option) on a copy of
  commit e08ac70: finished in 19 min 35 s on one core; both witness hashes, the caps (155608, 115075), the first
  three credits per side, the ceilings, the additional credits (17458365, 15289028) and both final fractions equal
  the claim.
* Own exact code written from the proof text and the documented file layout: all six singleton caps reproduced
  (upper 1544 / 154457 / -2532, lower 93464 / 95343 / 95165; 979,299 finite rows per witness and the analytic tails).
* Own backward occupation tables: all four singleton-donor credits reproduced (upper 8148713 and 3226994, lower
  11253928 and 301604) with zero violated rows.
* Its planted errors were detected (row violations or ledger mismatches). Its rereading found no double count and
  no founder term with the wrong sign. The quartet-donor table and the extension tables were replayed, not
  recomputed.

Using the full upper residual cap e = 155608 ticks (including the 154457 singleton component), the ceiling y from section 3 and only the singleton-donor credits from the child gives upper endpoint 0.5658074705, inside the rung-3 window. The lower endpoint needs no credit.

## 9. Inside the GPT family: second Sol Ultrafast child (the remaining tables)

Report: Fable folder `round2/children/astra9-tables/REPORT.md` and tracked `round2/checks/astra9-tables-gpt-family.md`, first line `RESULT: RECOMPUTED AND EQUAL`. With its
own programs, written from the proof texts and the documented byte layout:

* quartet-donor table of the first calculation (sets of up to five cores; 17,034,500 rows at arity five): newborn
  credits 7,709,975 upper and 7,402,100 lower, zero violated rows;
* edge caps 524 and 10,651, all six higher-family maxima, the per-core charges 627 and 9,081, and the projection
  checks;
* the full four-arity extension: 17,458,365 upper and 15,289,028 lower, with all logged row counts and
  positive-source counts;
* planted errors detected by row checks, cap checks or the newborn comparison.

Its list of places where the proof text alone does not determine the computation (verifiability, not correctness):
the byte order of the centering and FY sections is only in the exporter source; the published edge caps 524 and
10,651 depend on which summands are grouped before rounding upward (another grouping gives 523 and 10,650; the exact
rational values do not depend on it); Astra folder `round2/occupation/inside.md` and `round2/occupation-tail/continuation.md` specify when a table row is accepted but not
which accepted integer is taken (the largest one reproduces the logs); in the extension a finite degree-knot label is
treated as a tail and receives no source. The child checked that the reproduced rows satisfy the stated inequalities.

After sections 3, 8 and 9 every integer that enters the claimed fractions (e, y, the four credit families and the
extension, both sides) has been recomputed by code that is not the claimant's: y by this lead (Claude family), the
others by GPT-family children.

## 10. Not verified by a model family other than GPT

* The full base row scan was replayed, and the cap components, occupation tables and extension were separately recomputed, only by GPT-family
  children (sections 8, 9); this lead tested them indirectly (sections 4 to 6) and read the arguments (section 2).
* The per-core ceiling m(D, c) of the witness potential used for y was taken from the claimant's file.
* The earlier documents on which the founder and boundary remainders rest (Astra folder `round2/occupation/physical-test.md`) were read, not
  recomputed.

## Replay paths

Run the lead commands from Fable folder `round2/attacks/astra9/` with the original sibling names. In a different checkout, adjust the Astra root/snapshot/ceiling paths at the top of `lifetime_own.py` and `random_fringes.py`; the latter also needs the compiled low/high probes under Astra folder `round2/outputs/mc-replay-clean/`. Prepare the witness, ceiling data and probes using the replay instructions in Astra folder `round2/motif-cert/README.md` and `round2/occupation-tail/`. The tracked `round2/rival-attacks/code/astra9-*` copies in the Fable folder have a placeholder `ASTRA_REPO` and prefixed names; resolve the placeholder and restore the sibling names for shell calls and Python imports.
