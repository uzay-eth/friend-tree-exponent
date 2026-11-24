VERDICT: HOLDS

# Fable on Astra claim 15 (joint physical enclosure [0.565673570322, 0.565784239110])

**Scope of the verdict:** The Claude-family lead (a) recomputed both endpoint fractions exactly from the posted ledger with its own arithmetic, (b) measured the two quantities the composition rests on, the mean witness potential per core and the donor credits per core, on trees from its own simulator, and found every necessary inequality satisfied with margin, and (c) compared the bracket with its own numerical estimate of mu. The potential and the credits were evaluated with the claimant's exact rational classes on the lead's trees, so (b) tests the stated inequalities on independent trees and does not re-derive the witness. The proof of the two new steps (the five-pass quartet credit and the unsigned K7/144 denominator), the exact dual prices and the replay of the certificate were audited and re-run only by a GPT-6.2 Sol Ultrafast child, the same family as the claimant: that part is a check inside the GPT family. Claim 9, on which the composition rests, was attacked separately (`fable-on-astra-9.md`). A second-party run of the claimant's full cold replay (49 minutes, passed) is recorded in section 4.

Attacker: the Fable lead (Claude Fable 5.2; family different from the claimant's GPT). Programs and logs: Fable folder `round2/rival-attacks/code/astra15-*`.

## 1. The claim in my words

Two frozen "witness" scores (one per sign) give, by the drift comparison of claim 9, sigma (lambda - mu) m >= mu v - e, where lambda is the trial value of the witness (5656453/10^7 for the lower sign, 56581165/10^8 for the upper), m the limiting mean potential per core, v the sum of the certified positive credits per newborn, e the residual cap per core (all in ticks of 2^-40), sigma = -1 / +1. New in claim 15: the quartet-donor credit is replaced by a larger "five-pass" value (36,290,118 and 32,703,071 ticks instead of 7,402,100 and 7,709,975), and m is bounded by a new ceiling m <= 2^40 mu y with y = 126536000881/62500000000 (lower) and 1016333625559/500000000000 (upper), about 2.02 to 2.03 instead of about 2.66 before. With 113/200 <= mu the trial values move inward by at least (v - 200 e/113)/(2^40 y).

## 2. Own checks (Claude family)

**2.1 Exact arithmetic.** `astra15-endpoints.py` (reads no file of the claimant) forms v = 63,134,678 and 61,537,143 from the four ledger entries per sign, checks l v > e, and reproduces both posted fractions identically:
1356995264017948174596364819/2398901654969101844480000000 = 0.5656735703224... and 2725370105716046104606536597/4816977775848776007680000000 = 0.5657842391095.... The posted decimals are these values rounded outward. The inward shifts are 2.8270e-5 and 2.7411e-5.

**2.2 The new denominator.** The claim needs m <= 2^40 mu y in the limit of large trees. I evaluated the witness potential (claimant's classes `Local`, `Tree`) on whole physical trees from my own continuous-time simulator, grown from the three-vertex path to times 14, 16 and 18 (ratio of sums over trees, delete-one-tree jackknife):

| trees | cores in total | m / (2^40 mu y), upper | m / (2^40 mu y), lower |
|---|---:|---:|---:|
| 300 trees at time 14 | 1,069,213 | 0.9830 +- 0.0003 | 0.9865 +- 0.0003 |
| 100 trees at time 16 | 1,181,190 | 0.9827 +- 0.0003 | 0.9862 +- 0.0003 |
| 30 trees at time 18 | 1,139,239 | 0.9825 +- 0.0004 | 0.9860 +- 0.0004 |

The ratio must not exceed 1. It is stable in the tree size and stays below 1 by 1.7% (upper) and 1.4% (lower), at least 35 standard errors. The mean potential per core is 1.1298 (upper) and 1.1293 (lower) in units of 2^40 ticks, in line with the claimant's own finite-time figure of about 1.131. The margin is small because the new ceiling is close to the truth: the earlier ceiling was 21% above the measured mean.

**2.3 The credits.** On 1500 random planted fringes of 150 cores per sign (my simulation), with the credits computed in my program from definition (2) of the claimant's `occupation/inside.md`:

| per core, in ticks | upper | lower |
|---|---:|---:|
| quartet-donor credit, measured | 3.223e7 +- 0.096e7 | 3.106e7 +- 0.093e7 |
| mu v_Q4 with the five-pass value | 1.850e7 | 2.053e7 |
| ratio (must be >= 1) | 1.74 +- 0.05 | 1.51 +- 0.05 |
| singleton-donor credit, measured | 1.244e7 +- 0.017e7 | 1.117e7 +- 0.017e7 |
| mu (v_parent + v_child) | 0.644e7 | 0.654e7 |
| ratio (must be >= 1) | 1.93 +- 0.03 | 1.71 +- 0.03 |

So the five-pass value uses 57% (upper) and 66% (lower) of the quartet-donor credit that is physically there; the earlier value used about 13.5% on both signs. The population here is planted fringes of a fixed size, not whole trees, so these ratios are indicative; a ratio below 1 would have been a warning, none was seen. The fourth source (mask63, 15,289,028 and 17,458,365 ticks) was not measured by me.

**2.4 Balance.** The comparison is an inequality between two sides that can both be measured. With my estimate mu = 0.56572365 and the measured m: |lambda - mu| m = 9.73e7 (lower) and 1.093e8 (upper) ticks per core, against mu v - e = 3.56e7 and 3.47e7. The certified credits are 37% and 32% of the total that the exact balance provides; the certified shifts are 0.36 and 0.31 of the true distances of the trial values from mu. A certified shift larger than the true distance would have broken the claim.

**2.5 Against the numerical value.** My estimate 0.56572365 (statistical standard error about 7e-9; it rests on premises that are not proved, see my claim 1) lies inside the bracket, 5.0e-5 above the lower endpoint and 6.1e-5 below the upper one. The bracket is strictly inside the claim-9 bracket [0.5656581227, 0.5657979941] and inside [0.565615, 0.565809]. Its width is 1.107e-4: three digits of mu and no more, as the claim says.

## 3. Inside the GPT family: audit and recomputation by a Sol Ultrafast child

Report: Fable folder `round2/checks/astra15-audit-gpt-family.md` (first line "CLAIM 15: HOLDS"). The child

* re-derived the composition and the source accounting from the proof files and found no credit booked twice: the ledger books receiving potential occurrences and gain slots; the archived quartet entry is removed before the five-pass entry is added; the unsigned denominator spends none of the numerator sources;
* wrote its own exact verifier of the two dual price sheets of the K7/144 denominator (73,586 columns, 3,529 equalities, 50,209 inequalities per sign): dual feasibility holds and the stored rational margins are reproduced exactly (per-newborn bounds 126536000881/62500000000 and 1016333625559/500000000000); a ceiling smaller by 2% is rejected. The primal model was built with the claimant's builder;
* checked the multiplicities of the seven motif families against the deduced prices on 800 physical fixtures per sign with its own enumeration of connected subsets (largest actual/majorant ratio 0.9922 and 0.9911; a pointwise majorant smaller by 2% fails on 584 and 440 fixtures), and 296,988 pointwise row checks of the generator model on cyclic fixtures, with three plants rejected;
* recounted the five-pass geometry independently: all 50 rank counts match (4,431,601,576 rows); wrote an independent exact Bellman solver on the reduced mesh which reproduces the claimant's three small-split pass values for the lower sign (6,046,616 / 812,794 / 75,136 ticks) and reacts to removing the oldest donor;
* re-ran the claimant's full production split for the lower sign on a copy (1594 s): all five pass values and 25 counts match the frozen transcript, total 36,290,118 ticks; ran `replay.py quick` and `preflight` on the copy (pass; nine planted composition errors rejected).

Not done by the child: the upper production split, an independent production-size solver, a second implementation of the LP matrix, a new proof of claim 9.

## 4. Second-party run of the full cold replay

Run by the lead on a private copy of the claimant's sources (`git archive d87dc4c round2` plus the one listed input outside `round2`; the only change is the CPU number in the wrapper's affinity check, 3 -> 15; the copy contained no raw snapshot, no mbar table and no executable): `taskset -c 15 python -B round2/joint-enclosure/replay.py`. It rebuilt both raw snapshots and the comparison profiles, reran both five-pass production splits (lower 1025 s, upper 926 s; "all 50 rank counts agree"), regenerated the six K7/144 profiles ("matches frozen K7/144 contents" for each), checked 979,299 mbar entries per sign, rejected the nine planted composition errors, and ended with

    lower: v=63134678, e=115075, y=126536000881/62500000000, lambda=5656453/10000000, prior=113/200..1133/2000, endpoint=1356995264017948174596364819/2398901654969101844480000000
    upper: v=61537143, e=155608, y=1016333625559/500000000000, lambda=11316233/20000000, prior=113/200..1133/2000, endpoint=2725370105716046104606536597/4816977775848776007680000000
    cold PASS: outward physical interval [0.565673570322, 0.565784239110]
    wall 2965.63 s, max RSS 2866076 KiB

(log: `astra15-cold_replay.log`). This is a second-party execution of the claimant's programs: it shows that the published fractions follow from the 74 committed inputs on another party's run; it is not an independent implementation.

## 5. Not verified

The pointwise proof that the five-pass credit and the K7/144 price deductions are valid on every tree was read only inside the GPT family. My tests are statistical tests of necessary consequences on finite samples: they can refute, they cannot certify. Nothing in this claim or in this attack gives more than three digits of mu.
