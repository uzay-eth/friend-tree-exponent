# What was run and what resolved

All pilots below were run from the repository root with `timeout -k 1s 14s taskset -c 8 python3 -B …`, one process at a time, unless stated otherwise. Existing review scripts and recorded numeric artifacts are evidence; their execution here is a GPT-family check. No large replay, numerical simulation, all-degree production regeneration or external review was run.

## Independent code and negative controls

| Invocation after the common prefix | Observed output / scope |
|---|---|
| `round2/checks/scalar_hy4.py` | `SCALAR RETURN COLD AUDIT: PASS`: exact determinant and 54 active/27 dormant states; four plants (threshold, unsafe drift, numerator factor, excess sign) rejected. |
| `round2/checks/response_transcendence_hy4.py` | `RESPONSE TRANSCENDENCE COLD CHECK: PASS`: prime fixtures 3/5/7, Catalan recurrence and new-prime control; changing reward degree to size breaks uniqueness. No root arithmetic checked. |
| `round2/checks/excursions_hy4.py` | `EXCURSIONS COLD CHECK: PASS`: matched 24 parsed values between stored logs; audited Taylor/moment identities and `0.807760`, `499/500`, multiplier 4800. Wrong ratio, swapped moment endpoints, dropped `(1+d_A)`, looser cap and anchor rejected. The “fresh” label printed by the script refers to the stored reviewer log. |
| `round2/checks/physical_trigger_hy4.py` | `PHYSICAL TRIGGER COLD CHECK: PASS`: exact cap-3/4 tests, scalar signs and monotonic grid; rejects all-degree-two trigger, dropping complement births and wrong renewal column. |
| `round2/checks/spatial_limit_hy4.py` | `INDEPENDENT FIXTURE CHECKER: PASS`: 10,099 capped-label comparisons, 800 fringes, 311 embedding cases, 80 unbiasing trees; rejected wrong label lower bound, bypass of depth bound, missing size test and old-as-new error (4/4). |
| `round2/checks/fringe_marginals_hy4.py` | `INDEPENDENT FIXTURE CHECKER: PASS`: RN identity, geometric ratio, rejection and size/crossing laws, seven plants rejected. Fixture `1131/2000` is a generic rational, not an exponent estimate. |
| `round2/checks/hub_spectrum_hy4.py` | `HUB SPECTRUM COLD CHECK: PASS`; wrong thinning, wrong normalization, reversed squeeze and missing tags rejected. |
| `round2/checks/hub_children_hy4.py` | `INDEPENDENT FIXTURE CHECKER: PASS`: 462 rational states, 870 physical vertices, birth-death and shield checks; wrong base/exponent/trigger, missing leaves and wrong survival rejected (5/5). The Poisson integral fixture is numerical, and the asymptotics use the proof. |
| `round2/checks/hub_children_sol.py` | Exit 0: 1,401 ordinary and 1,401 silent-slot core rows; factorial and star races through degree 64; wrong promotion denominator, wrong promotion degree, conditioning marks, and overlapping categories rejected. [Prose verdict](../checks/hub-children-sol.md). |
| `round2/checks/edge_obstruction_hy4.py` | `INDEPENDENT REFEREE CHECKER: PASS`: mixed-branch identity (24 leaf and four skeleton cases), commutation and bare path; missing branch, nonzero forbidden coefficient, and unrestricted-growth plants rejected. |
| `round2/checks/five_split_hy4.py` | `FIVE SPLIT COLD CHECK: PASS`: 25 hashes, 50 counts, 398 physical trees; re-added baseline, dropped row, incorrect signed source, missing factor and omitted incidence detected; not a billion-row replay. |
| `round2/checks/fringe_neon.py` | `ALL INDEPENDENT CHECKS PASSED in 3.22s.`: 1,039 marked physical trees, nine deliberate corruptions rejected. |
| `round2/checks/asymptotics_carbon.py` | `ALL CHECKS AND 8/8 NEGATIVE CONTROLS PASSED.` |
| `round2/checks/finite_kernel_neon.py` | `ALL INDEPENDENT CHECKS PASSED in 2.46s.`: 361 exact capped path-sum entries at size 3, ten planted errors caught. |
| `round2/checks/strong_limit_carbon.py` | `ALL CHECKS AND 10/10 MUTATION CONTROLS PASSED.` |
| `round2/checks/hub_mass_carbon.py` | `ALL HUB-MASS FINITE-SIZE CHECKS AND NEGATIVE CONTROLS PASSED.` All five controls caught. |
| `round2/checks/local_completeness_neon.py` | `ALL INDEPENDENT CHECKS PASSED in 8.77s.` Eight plants caught. |
| `round2/checks/operator_carbon.py` | Exited 124 at pilot deadline, no fresh pass claimed. See historical transcript (`outputs/operator-carbon.txt`). |
| `round2/checks/curved_retained_hy4.py` | Initially `AssertionError: round2/occupation-curved/profile_check.py` (manifest mismatch). After the working tree acquired refreshed manifests, the rerun produced no completion output within the pilot; no fresh cold pass claimed. |
| `round2/occupation-curved/check_records.py --retained-only` and separately `--partial-only` | Initially both missing their `evidence.sha256`. Those files appeared during the audit; both now print `frozen records PASS`. Retained increments, upper/lower: 4,586,074/4,008,588 ticks. Partial record: 23,740,968 completed upper rows, no newborn value for the interrupted full mask 127. No curved credit used in the joint result. |

The gauge contributor programs `round2/finite/gauge/physical_trigger_{bounds,check,verify}.py` were also replayed individually. `check` gave 6/24/92 exact killed states and virgin-row elimination, six zero incomparable mixed differences and one nonzero comparable difference. `verify` reproduced 756 Fraction endpoint comparisons, the 189-pair times two size-weight test and the coefficient `<2837`. Direct runs of `physical_trigger_bounds.py --critical-only --all-labels --iterations 16` and `physical_trigger_bounds.py --all-labels --refined-only --iterations 16` reproduced `970147668298 / 1e12` and `976433602362 / 1e12`. `physical_trigger_sector.py --degree 6 --cores 3 --steps 96 --mode basis` gave 7,575 states, 44,815 edges and `251490678253 253754698168` in outward units `1e-12`. These are additional GPT-family replays of author code. The independent trigger fixture does not cover these later full-size bounds.

The direct invocation of `round2/joint-enclosure/check_renewal_star_tail.py` exited 1 on its hard-coded CPU 3 guard. Under the same pinned CPU 8 prefix, I invoked its existing `certificate`, `verify`, and `sign_certificate` functions over the same 19 rows, same rational parameters and same three plants per row. Output: `both star intervals in all 19 rows and 57 sign plants passed`. This checks first-promotion mass, not the post-promotion response or a renewal sign.

## Static joint composition

The joint entry point insists on CPU 3 and its quick mode uses a scratch directory outside the assigned output directory. I called the existing pure read/check helpers directly on CPU 8 instead:

```sh
timeout -k 1s 14s taskset -c 8 python3 -B -c 'import sys; from pathlib import Path; sys.path.insert(0,"round2/joint-enclosure"); from replay import ROOT, frozen_hashes, sheet_info, SIDES; from arithmetic import *; O=ROOT/"round2/outputs"; hashes=frozen_hashes(ROOT); ys={s:sheet_info(ROOT,s,hashes)[0] for s in SIDES}; x=read_ledgers((O/"occupation-affine-verified.txt").read_text(),(ROOT/"round2/five-next/split-reference.tsv").read_text(),ys); q=enclosure(x); t=ticks(q); crosscheck(x,q,t); plants(x,q,t); print("Exact:",*[str(v) for v in q]); print("Outward ticks:",t); print("Pinned inputs:",len(hashes)); print("Source totals:",*[sum(dict(x[s]["sources"]).values()) for s in SIDES]); print("JOINT STATIC CHECK PASS")'
```

```text
rejected archived Q4 overlap
rejected duplicate mask63
rejected unreviewed mask64
rejected archived Q4 hidden in replacement
rejected narrowed growth prior
rejected flipped residual
rejected wrong upper shift sign
rejected inward lower tick
rejected inward upper tick
Exact: 1356995264017948174596364819/2398901654969101844480000000 2725370105716046104606536597/4816977775848776007680000000
Outward ticks: (565673570322, 565784239110)
Pinned inputs: 25
Source totals: 63134678 61537143
JOINT STATIC CHECK PASS
```

This binds frozen input hashes, stored raw/affine ledger, frozen quartet reference, unsigned sheet roles and profile normalizations. It is not the `cold` entry point, a new LP verification, or an independent production generation. The separately rerun five-split checker checked the reference against the production counts.

## Paths and portability

- Assigned Markdown was scanned for local file links; proposed RESULT copy was resolved as if pasted into `round2/RESULT.md`. The earlier scan found two upstream misses: `round2/numerics4/RESULTS.md` and `round2/outputs/numeric4-freeze.sha256`. The initially absent curved evidence manifests appeared during the audit and their links resolved.
- The [independent joint check](../checks/joint-physical-sol.md) verifies the exact frozen composition, 25 split pins, 50 rank counts and two unsigned sheets. The original 72-file allowlist failed. The repaired 74-file manifest passed four isolated import probes and resolved 25 local include edges without the audit overlay. No production stage ran.
- `round2/moments/README.md` and `round2/moments/replay.sh` use the ignored `round2/.work/venv/bin/python`; this is a local interpreter convention, not a distributable dependency. `round2/moments/rewards.py` uses the Linux runtime pseudo-path `/proc/self/fd/`. Review/replay docs and logs contain historical CPU 0, 3 or 7 affinity; use direct subcommands under CPU 8 here. `round2/finite/gauge/recursive-excursions.md` prints a literal `\n\n` before the finite exact-reference paragraph. The gauge proof/check paths were otherwise repo-relative within the assigned scope.
- Suggested root-specific references use `round2/...` in inline code; click-through references in the scope map are file-relative. No runtime paths, archive binaries or ignored work directories are proposed for the scalar commits.

## 08:14 checkpoint replays

Each bounded invocation below used `timeout -k 1s 14s taskset -c 8`; Python commands used `python3 -Bu`. `PAIR_EXE` denotes the scratch executable. This is the historical 08:14 snapshot; the completed cold reviews and subsequent clean replays are below.

| Invocation after the common prefix | Observed result |
|---|---|
| `round2/five-next/split_integers.py round2/outputs/five-next-split-cold.log` | `all 50 rank counts agree; re-added baseline, perturbed pass and removed row rejected`. Receipt's stdout hash matched. |
| `round2/checks/current_gauge_sol.py` | Exact physical-arrow and shell plants rejected; `theta(.52)=65/66` fixture; newborn derivative difference `1/96`; all 19 exact size-six enclosures agree, both mutations miss all 19; cap-36 corners/profile and four seam shapes pass. |
| `round2/current-mark/eval/check.py --cutoff 6` | `canonical_shapes=36 current_active=684 current_total=828 fixed_states=1330 independent_double_gap=5.55e-17`; all 19 killed comparisons reproduced. |
| `round2/current-mark/weight/uncapped.py` | `PASS cutoff=160 labels=uncapped rho=553/1000 d=17/50 edge=101752 local=25438` (default control). |
| `round2/current-mark/weight/uncapped.py --M 400 --d .365 --rho .53` | Exit 124, no output within 14 seconds. The `0.53` value remains the author's recorded certificate, not a fresh replay pass. |
| `round2/current-mark/weight/controls_uncapped.py`; `round2/current-mark/weight/interpolate.py` | 288 exact uncut two-level probes; wrong sender rejected. Interpolation: 312 exact inequalities checked. |
| `g++ -O3 -std=c++20 round2/influence/zero/pair/pair.cpp -o "$PAIR_EXE"`; `round2/influence/zero/pair-audit/audit.py`; `round2/influence/zero/pair-audit/corners.py --exe "$PAIR_EXE"` | Compile and labelled-arrow audit succeeded; split-clock and drop-private rejected; held-out supports 1, 12, 157; `14 seed bounds exactly equal independent corners`. |
| `round2/influence/zero/pair/fixture.py "$PAIR_EXE" 8`, also with `split-root` and `drop-private-core` | `PASS cap 8: labelled physical clock replay and exact rational Bellman rows`; both plants rejected at `d=7 L=0` (expected exit 1). No cap-48 run here. |
| `round2/influence/zero/independent.py 3`; `round2/influence/zero/pair/newborn.py 8` | `PASS N=3: 3 explicit signed rows, full residual <= 0.000000000001`; full newborn enclosure `[-17.651405629152,18.727646469432]`, giving no near-critical sign. |

The two inventories passed static validation against `84aa4c4`: 30 unique source/proof files, three shared dependencies, and 221 Markdown file links resolved across the inventories and reports. `git diff --check` and `git diff --cached --exit-code` exited zero; `git rev-parse --short HEAD` returned `84aa4c4`.

## Final clean certificate replays

Each inventory was copied into its own disposable export; direct links to inherited material were resolved against `84aa4c4`. The inventories passed file and link checks: 21 positive files (17 proof links), 32 signed files (19 proof links), four shared dependencies, 49 unique additions. Every listed file is source, proof or a small exact fixture. No file outside the respective inventory was imported by its replay.

All commands below ran at the repository root of the corresponding export with `timeout -k 1s 14s taskset -c CORE` and `python3 -Bu` for Python scripts, unless noted. The two jobs used separate cores. Rejections listed as plants were expected.

| Clean replay | Observed result |
|---|---|
| `round2/checks/uncapped_gauge_sol.py` | `independent exact edges/local endpoints 638392 159598 hmin 86271979/1000000000`; zero/near-zero marks, degree seams, tail controls, Taylor bound `5.043888117780772e-11`. This is the full exact `0.53` certificate replay, using the fast independent checker in place of the slower author arithmetic. |
| `round2/checks/current_gauge_sol.py`; `round2/current-mark/physical/check.py`; `round2/current-mark/check/check.py` | Cap-20 shell/contraction, 19 exact killed enclosures, cap-36 seams and plants pass; 145 exact physical transfers and five dormant derivative orders; 122 independent physical coefficient rows and 648 Green identities. Newborn unit-reward derivative difference starts at order four, `1/96`. |
| `round2/current-mark/weight/cert.py`; `round2/current-mark/weight/physical_check.py`; `round2/current-mark/weight/cert.py --cap 20 --rho .56499` | `PASS cap=36 rho=141/250 d=8/25 edge corners=5178 degree rows=665`; 95 physical shapes and two plant rejections. Cap-20 candidate rejected at required normalization `1.0618710838>1`. |
| `round2/current-mark/weight/uncapped.py` (default and `--coarse`); `round2/current-mark/weight/controls_uncapped.py`; `round2/current-mark/weight/interpolate.py` | Default and 21-label certificates pass; 288 physical rows and wrong-sender rejection; 312 exact interpolation inequalities. |
| `round2/current-mark/eval/check.py --cutoff N` for `N=2,6,8`, and `--cutoff 6 --plant activation-loss`; `round2/current-mark/eval/gauge.py --cutoff 10 --steps 45 --certify`, `--cutoff 11 --steps 20 --stream`, `--cutoff 12 --steps 16 --stream` | All three finite comparisons pass, activation-loss rejected (`0.0018158`); size-10 killed ratios `[0.339324519229,0.567946214059]`; sizes 11 and 12 reproduce streamed floating proposals. |
| `round2/checks/zero_contraction_sol.py`; `round2/influence/zero/critical_controls.py` | Exact physical pairs and two plants; degree-48 analytic coupon `0.839616806849242`; low-degree full fixture at `.565728` crosses one, at `s=2` lies below; signed inverse norm `<7.146`, full-window inverse factor `<471522`, critical inverse slope `<845`. |
| `round2/influence/check/bound.py` (default, `--challenge`, `--s .525 --crude`, `--s .565`); `round2/influence/check/coupling.py`; `round2/influence/gradient/constants.py` | Required envelope uppers `.964260985` and `.839292428`; crude `.5` contraction rejected. 432 paired states in each of three slot settings and two plants; normalized gradient `<16734`. The combined `--refined` command exceeded 25 seconds; running its two cases separately with `--s .5` and `--s .565`, each with `--steps-per-unit 1024 --time 26` and a 25-second timeout, reproduced `.992198062` and `.838093404`. |
| Build `round2/influence/zero/pair/pair.cpp` and `round2/influence/zero/surrogate/envelope.cpp` with `g++ -O3 -std=c++20 -fsanitize=undefined -fno-sanitize-recover=undefined`; run `round2/influence/zero/pair-audit/audit.py`, `round2/influence/zero/pair-audit/corners.py --exe "$PAIR_EXE"`, `round2/influence/zero/pair/fixture.py "$PAIR_EXE" 8` (also `split-root` and `drop-private-core`), and `round2/influence/zero/surrogate/fixture.py "$ENVELOPE_EXE" 5 -1 1` | 14 exact cap-six seed bounds match; cap eight passes and both plants exit 1; surrogate fixture passes `{'row': 9, 'value': 20, 'tail': 1}`. Compiler warning about the unused renamed entry point, as in the referee report; no sanitizer failure. Cap 48 was not repeated: the cold referee already reported 22,956,300 rows, max ticks `8600511260`, tail ticks `8484503829`. |
| `round2/influence/zero/independent.py 3`, `8`, `8 --plant`; `round2/influence/zero/obstruction.py`; `round2/influence/zero/solver.py --cutoff 8 --split`, `--cutoff 11`, `--cutoff 12 --unweighted` | 3 and 199 exact signed rows pass; negative-copy omission rejected (`0.414082642416`); 256 local split rows pass; size-12 killed `U(J)` is `[0.609998610225,0.609998610226]` (25-second slot). |
| `round2/influence/zero/pair/newborn.py` with `8` and `11`; `round2/influence/zero/adjoint/replay.py --certificate round2/influence/zero/adjoint/fixture.json --coupled-layers 1`; `round2/influence/zero/eval/check.py` | Size-11 unlimited continuation enclosure `[-16.283044514395,17.473459651808]`; exact adjoint fixture yields `[-35.919233950814,36.419906494642]`; 499 marginal fixtures and 6,171 stock rows pass. |
