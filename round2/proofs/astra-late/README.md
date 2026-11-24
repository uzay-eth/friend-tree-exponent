# Astra (GPT family): Claims 17–28

**Proof author: Astra (GPT family). Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. The source subdirectory preserves the original relative layout. [SOURCES.tsv](SOURCES.tsv) records original commits and Git blobs, editorial changes, canonical review destinations, and exclusions. The [public review index](../../reviews/README.md) holds the consolidated cold reviews. Clean hash-bound originals are stored as text data under `authenticated/`; see the replay section for how they are restored. Editorial copies do not substitute for originals in an authenticated replay. Finite fixtures check the listed identities and certificates; the analytic limits depend on the corresponding proofs and [earlier Astra claims](../astra/README.md).

| Claim | Astra proof and scope | Checks and public reviews |
|---|---|---|
| 17 | [Amplitude law](source/round2/open/limit-law/README.md), [left tail](source/round2/open/limit-law/left/README.md), [prefactor](source/round2/open/limit-law/prefactor/README.md), [positive density a.e.](source/round2/open/limit-law/density/README.md), [exact right-tail rate](source/round2/open/limit-law/right/README.md). No numerical rate or prefactor is certified. [Translation cusp](source/round2/open/limit-law/smoothing/cusp/README.md) additionally assumes repaired Fable Claim 6. | Arrow/smoothing and urn fixtures; [Gemini checks](independent/audit_claim17.py) for equivalence and right rate under the cited growth inputs. [Fable](../../reviews/attacks/fable-on-astra-17.md), [Neon](../../reviews/attacks/neon-on-astra-17.md). |
| 18 | [Birth-label comparison](source/round2/open/age-order/README.md), [late hub probability](source/round2/open/age-order/late/README.md), [adjacent gap](source/round2/open/age-order/adjacent/README.md), [near-one equivalent](source/round2/open/age-order/near-one-constant/README.md), [compact-interior equivalent](source/round2/open/age-order/fixed-threshold/uniform/README.md). The conditional size-204 reversal concerns a completed placement, not the pre-placement comparison in 26. | Endpoint/arrow fixtures; [Fable](../../reviews/attacks/fable-on-astra-18.md) distinguishes Claude checks of symmetry and reversal from same-family asymptotic audits. |
| 19 | [Fixed degree cutoff](source/round2/open/low-forest/README.md): almost-sure finite-type TV limits for core sampling and once-per-component counting; ordinary-vertex isolation; tail estimates at fixed cutoff. | [Gemini fixtures](independent/audit_claims_19_20.py); [Fable](../../reviews/attacks/fable-on-astra-19.md), [Neon](../../reviews/attacks/neon-on-astra-19.md). |
| 20 | [Sequential count, TV, and log-volume laws](source/round2/open/forest-age/README.md), [inert seed](source/round2/open/forest-age/inert/README.md), [fixed shapes](source/round2/open/forest-age/fixed-shapes/README.md), [marked cutoff repair](source/round2/open/low-forest/marked-repair/README.md). The repair does not restore the original log-power survivor TV range or bound (3). | Boundary-flux, marked repair and [Gemini fixtures](independent/audit_claims_19_20.py); [Fable base/addenda](../../reviews/attacks/fable-on-astra-20.md), [Fable gap](../../reviews/attacks/fable-on-astra-20b.md), [Neon base](../../reviews/attacks/neon-on-astra-20.md). |
| 21 | [D18 high101](source/round2/occupation-neon/high-hub-limit/release/README.md) on [Yt](source/round2/occupation-neon/pressure-memory/fork-lifetime/yt-release/README.md) and [repaired census](source/round2/occupation-neon/interpolate/compose/README.md); [D20 baseline](source/round2/occupation-neon20/renewal/README.md), [signed corrector](source/round2/occupation-neon20/corrector/README.md) and [older-arm alternative](source/round2/occupation-neon20/renewal/older/README.md). D18 and D20 price distinct original witnesses; the D20 alternatives are not composed. | Original payloads required for production; [Fable](../../reviews/attacks/fable-on-astra-21.md) documents Claude archive replay and endpoint arithmetic, and same-family physical pricing audits. |
| 22 | [Within-tree marked spectrum](source/round2/open/age-order/spectrum/quenched/README.md), [sampling law](source/round2/open/age-order/spectrum/README.md), [adjacent moments](source/round2/open/age-order/moments/README.md). Fixed positive thresholds, no rate or hub-independence result. | Backward forecasts, signed moments and [Gemini fixtures](independent/audit_claims_22_23.py); [Fable](../../reviews/attacks/fable-on-astra-22.md), [Neon](../../reviews/attacks/neon-on-astra-22.md). |
| 23 | [Sublinear fixed shapes](source/round2/open/forest-age/fixed-shapes/diagonal/README.md) and [proportional stable Poisson limit](source/round2/open/forest-age/macroscopic/README.md). General sublinear convergence is in probability; the slow-cutoff counterexample precludes a general a.s. law. | Residence and size-clock controls; [Gemini fixtures](independent/audit_claims_22_23.py); [Fable](../../reviews/attacks/fable-on-astra-23.md), [Neon](../../reviews/attacks/neon-on-astra-23.md). |
| 24 | [Eventual uniform tail order](source/round2/open/age-order/eventual-domination/README.md), [effective tower cutoff](source/round2/open/age-order/fixed-threshold/effective-uniform/README.md). One late label cutoff for all thresholds, conditional on inherited spectral constants; does not close the intervening label range at general thresholds. | Exact placement/cylinder fixtures; [Fable](../../reviews/attacks/fable-on-astra-24.md) consistency check, with same-family full audit. |
| 25 | [Gold1 composition](source/round2/occupation-quintet/conditioned-joint/composition/README.md) gives `[0.565704609764900, 0.565745628170000]`; [previous stationary ledger](source/round2/occupation-quintet/stationary/README.md). Credits belong only to these Gold1 witnesses. | 613 equations and 5,554 columns per side in the archived calculation; [Fable](../../reviews/attacks/fable-on-astra-25.md) replay covers earlier freezes, not the final composition. |
| 26 | [All-tree hubhood](source/round2/open/age-order/global-transport/hostile/RELEASE.md) and [founder/telescoping synthesis](source/round2/open/age-order/hub-ranking/README.md). [One-tenth proof](source/round2/open/age-order/global-transport/quarter/README.md) and [separate label-20 check](source/round2/open/age-order/global-transport/quarter/b20.py) with [cold scope](../../reviews/cold/astra/open/age-order/hub-ranking/beta-tests/physical-tail/poisson/quarter.md) and the Fable early-label splice. The external hubhood attacks do not review that positive-threshold band. | [Gemini reimplementation](independent/audit_claim26.py); [Fable](../../reviews/attacks/fable-on-astra-26.md), [Neon](../../reviews/attacks/neon-on-astra-26.md), [early Fable transports](../../reviews/attacks/astra-on-fable-5.md). |
| 27 | [Every positive real power](source/round2/open/age-order/hub-ranking/all-powers/positive/README.md): conditional lower bound, founder seam, earlier-history and stopping-size comparisons. This does not imply order of sharp tails. | Joins and beta telescopes; [Gemini reimplementation](independent/audit_claim27.py); [Fable](../../reviews/attacks/fable-on-astra-27.md), [Neon](../../reviews/attacks/neon-on-astra-27.md). |
| 28 | [Repaired four-reward rational exclusion](source/round2/occupation-quintet-gold2/rationals/repaired/README.md); [separate corrected Gold2 frontier](source/round2/occupation-quintet-gold2/host/frontier/README.md) gives `0.565714467904 < μ < 0.565735974375`. The old three-reward and D24 comparisons are withdrawn; `241/426` remains in the interval. | Archived-input check and [Claude exact recomputation](independent/astra28-lifetime_recompute.py) both reuse authenticated local cells and reward floors. [Fable](../../reviews/attacks/fable-on-astra-28.md) independently covers the stationary and Farey reduction and finite ceiling, not the four reward floors or full frontier. |

## Claim-specific review scope

The [Astra ledger](../../LEDGER.md#astra) identifies reviewed revisions, distinguishes review families from inherited dependencies, and records the limits of each verdict.

### Claim 17

GPT proof audits and scoped Claude and Gemini reviews; Claim 6, inert limits, and attachment influence; the cusp also uses repaired Fable Claim 6.

Arrow/smoothing and urn fixtures; [Gemini checks](independent/audit_claim17.py) for equivalence and right rate under the cited growth inputs. [Fable](../../reviews/attacks/fable-on-astra-17.md), [Neon](../../reviews/attacks/neon-on-astra-17.md).

### Claim 18

GPT asymptotic audits and scoped Claude founder and reversal checks; earlier degree, hub, and spatial results.

Endpoint/arrow fixtures; [Fable](../../reviews/attacks/fable-on-astra-18.md) distinguishes Claude checks of symmetry and reversal from same-family asymptotic audits.

### Claim 19

GPT, Claude, and Gemini within their recorded scopes; Claim 10 and the certified growth lower bound.

[Gemini fixtures](independent/audit_claims_19_20.py); [Fable](../../reviews/attacks/fable-on-astra-19.md), [Neon](../../reviews/attacks/neon-on-astra-19.md).

### Claim 20

GPT, Claude, and Gemini within their recorded scopes; Claims 3 and 10; Gemini reviewed the base sequential law, with separate scopes for addenda.

Boundary-flux, marked repair and [Gemini fixtures](independent/audit_claims_19_20.py); [Fable base/addenda](../../reviews/attacks/fable-on-astra-20.md), [Fable gap](../../reviews/attacks/fable-on-astra-20b.md), [Neon base](../../reviews/attacks/neon-on-astra-20.md).

### Claim 21

GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation.

Original payloads required for production; [Fable](../../reviews/attacks/fable-on-astra-21.md) documents Claude archive replay and endpoint arithmetic, and same-family physical pricing audits.

### Claim 22

GPT, Claude, and Gemini within their recorded scopes; Claims 6, 10, 11, and 18.

Backward forecasts, signed moments and [Gemini fixtures](independent/audit_claims_22_23.py); [Fable](../../reviews/attacks/fable-on-astra-22.md), [Neon](../../reviews/attacks/neon-on-astra-22.md).

### Claim 23

GPT, Claude, and Gemini within their recorded scopes; earlier growth, spatial, spectrum, birth-label, and fixed-shape results.

Residence and size-clock controls; [Gemini fixtures](independent/audit_claims_22_23.py); [Fable](../../reviews/attacks/fable-on-astra-23.md), [Neon](../../reviews/attacks/neon-on-astra-23.md).

### Claim 24

GPT proof audits and scoped Claude consistency checks; earlier spectral constants and, for the named early-label splice, Fable transports.

Exact placement/cylinder fixtures; [Fable](../../reviews/attacks/fable-on-astra-24.md) consistency check, with same-family full audit.

### Claim 25

GPT proof audits and Claude replay of earlier freezes; original Neon-2 Gold1 witnesses and inherited stationary results.

613 equations and 5,554 columns per side in the archived calculation; [Fable](../../reviews/attacks/fable-on-astra-25.md) replay covers earlier freezes, not the final composition.

### Claim 26

GPT, Claude, and Gemini within their recorded scopes for hubhood, with GPT-only audits of the accepted positive-threshold proof; threshold addenda use Claim 24 or pinned Fable finite transports where stated.

[Gemini reimplementation](independent/audit_claim26.py); [Fable](../../reviews/attacks/fable-on-astra-26.md), [Neon](../../reviews/attacks/neon-on-astra-26.md), [early Fable transports](../../reviews/attacks/astra-on-fable-5.md).

### Claim 27

GPT proof audits, scoped Claude derivations and replay, and Gemini proof review and independent checks; prior two-birth identities and finite box certificates as cited.

Joins and beta telescopes; [Gemini reimplementation](independent/audit_claim27.py); [Fable](../../reviews/attacks/fable-on-astra-27.md), [Neon](../../reviews/attacks/neon-on-astra-27.md).

### Claim 28

GPT physical-floor and frontier audits and scoped Claude conditional reduction and arithmetic; original Neon-2 Gold2 witnesses and the four reward floors.

Archived-input check and [Claude exact recomputation](independent/astra28-lifetime_recompute.py) both reuse authenticated local cells and reward floors. [Fable](../../reviews/attacks/fable-on-astra-28.md) independently covers the stationary and Farey reduction and finite ceiling, not the four reward floors or full frontier.

## Replay from a release clone

`bash round2/verify/astra-late/run.sh` runs the modest checks, including the five unchanged Gemini programs and the Claude lifetime check. `quick` skips the Gemini programs; `independent` runs only those programs. These modes use Bash and Python 3.10 or newer with the standard library; choose the interpreter with `PYTHON`. The heavy original replay tools use the upstream versions of their inputs and may require NumPy and SciPy, GCC/G++ with OpenMP, or GMP.

For a checker with all of its other inputs already available, `python3 -IB round2/verify/astra-late/with-originals.py CHECKER.py [ARGS...]` restores the clean archived Markdown inputs at their original paths in a fresh temporary copy and launches the unchanged checker there. `CHECKER.py` is relative to the original Astra source root. Pass absolute paths for any arguments that name external inputs. The wrapper authenticates the archived text before restoring it; the checker still applies its own pins. It does not reconstruct omitted originals or raw certificates.

`bash round2/verify/astra-late/run.sh threshold` runs the Poisson overlap, the nine exact label-20 boxes and the producer’s closed-box sweep. For shards, from `round2/proofs/astra-late/source`, run `python3 -IB round2/open/age-order/global-transport/quarter/check.py --start 21 --stop 22` (stop exclusive), then the other ranges in the proof. The cold check reran only the `21→22` shard and the separate label-20 extension; the producer ran all other shards. The separate one-sixteenth producer record had no completed cold review and supplies no additional all-label result.

The modest Claim 28 replay checks SHA-256s of the original Gold2 receipt, local ceiling table and compact reward records, reuses unmodified Astra cold and writer arithmetic, and compares the two authenticated frontier receipts. It does not rebuild original physical reward floors, full Gold2 drift rows, or the frontier LP. The separate [arithmetic sensitivity checker](source/round2/checks/arithmetic-transform-sol.py) enumerates bounded heights on its stated certified intervals. It does not identify the root.

The modest Claim 21 replay authenticates the two saved high101 price receipts, reruns finite seam/registration fixtures and computes the two exact D18 high101 endpoints from frozen D18/Yt constants and the saved new prices. It does not read raw D18 rows, regenerate prices, or replay D20.

## Original numerical replay dependencies

These commands describe the historical releases, run from the root of the corresponding original Astra archive. Obtain the exact original payloads or regenerate them using the pinned [Neon](../neon/README.md) and [Neon-2](../neon2/README.md) instructions. Authenticate both files before a heavy replay; rebuilding an input only substitutes for the archived original when its SHA-256 matches. Original production and proof pins are different from release provenance pins.

- **D18 and D20 (Claim 21):** D18 pin `cec3855cc52ab4265ce3e36b8a7a38f1ae257f77`, 3,495,091,968 bytes per side; D20 pin `51f4cfadac675fdf1cef84d44935eebb5456aabb`, 11,190,881,056 bytes per side. D18 ledger: `python3 -EBu round2/occupation-neon/high-hub-limit/release/verify.py "$D18_LOWER" "$D18_UPPER"`. D20 strongest alternative, each sign: `python3 -Bu round2/occupation-neon20/corrector/replay.py --side lower --mode preserve --edge-degree 16 --arity 5 --stationary 40 --fixtures --exact-den` (repeat for `upper`, with original payload and source preparation in the corrector proof). The independent D18 archive manifest requires `round2/checks/high-hub-sol.md` at SHA-256 `04309b76fe60f8a6466037f28fc34beca9246accb594fd5a99c7dc1677f7b29f`; the [canonical cold review](../../reviews/cold/astra/checks/high-hub-sol.md) is an editorial copy. Some archived D18/Yt and D20 programs also pin original files omitted for machine-local paths. The raw GMP production scans are separate from stationary repricing.
- **Gold1 (Claim 25):** original Neon-2 pin `b04cc50c1e391e0420eb7b006ac781d484c68543`, 1,588,521,616 bytes per side. With the original `quintet_{lower,upper}.bin` and receipt `gold.json` (SHA-256 `d26d22a31d8f18790f517c451708a35b2721cbb4809613d1316c05df5071faeb`) in `GOLD1_WORK`, run `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -IBu round2/occupation-quintet/conditioned-joint/composition/check.py "$GOLD1_WORK"`; `--solve` regenerates the saved duals. This needs the upstream `round2/attacks/neon2-quintet/physical/check.py` via `occupation-quintet/stationary/weights.py`; the physical reader with machine-local paths is not imported. The historical cold main also has machine-local source assumptions. Diagnostic H3 dumps are not required for the exact composition.
- **Gold2 (Claim 28):** original Neon-2 pin `b435b4c824e2ba5835f58b85e67a27de7d6a0d49`, 6,344,902,168 bytes per side. Lower `41565d49c10809eb1f732794822359410af3ed804577476563bef17f27946413`, upper `f1dcd7ad3e6df82a985eb3dd6fd793b5a3e38c78368d94e4b135ab174ab92876`. Use two threads for the original verifier (recorded 322.14 s, peak 12,414,080 KiB). The frozen `rationals/repaired/check.py` main consumes `lifetime/README.md` at SHA-256 `17ad57b53161e61517a8ea776fa536abec39a76f967a1c57fe044b14fb418a5c`; the public documentation copy is sanitized and does not substitute. `--rebuild-host`, `checks/gold2-rational-repair-sol.py` and `host/frontier/check.py --live-sibling` hard-bind `host/reader.py` (blob `796d43b52d9f472513c77823e862b26fc580bcc4`, SHA-256 `6be2e8aa5835ad0f1087179af5d84664d1176ed3ab071ca64f58af8d628fc910`), omitted because it embeds a machine-local input path; it expects both authentication receipts and an exact detached Neon-2 source checkout. The related tail/physical reader files with the same class of path are likewise excluded. The compact sibling and upper checks additionally require the original source/inputs pinned in `checks/gold2-frontier-sol.py` and `receipt-live-sibling.json`.
- **Six-marginal alternative (Claim 28):** `checks/gold2-no-h3-sol.py` uses a source tar at SHA-256 `2f2c2c1ee81577255aa4a8834b7c97850e7445f92008351bc5d2862e35cf7177` of Neon-2 pin `6c8afec5cd1b7f98bbcd7c6e7848c0d631a50f7d` and Astra's `attacks/neon2-marginals/replay.json` at SHA-256 `a3565adcb81ca54db62112c74a7c59b6056387e501691fc4c1badd14008891a6`. The tar and the receipt with machine-local identifiers are excluded. The recorded alternative uses only six lower mass bounds and no H3 reward, clears `99/175` conditional on those bounds, and does not narrow either final endpoint.
- **Density cusp (Claim 17):** the original `checks/limit-cusp-sol.py` pins three documents including `attacks/fable6-density/README.md` at blob `b2dabc52e9e2bebc293c2708b29ae9d166867550`. The original audit text embeds private paths and is excluded; the modest replay invokes only the unchanged scale-budget and model translation fixtures, without the historical source-pin step. The [public Fable audit](../../reviews/attacks/astra-on-fable-6.md) records its scope.

The archived exploratory program `occupation-neon/arms/replay.py` has an indentation error at line 127 in its pinned source. It supplies no passing replay or endpoint in this release.
