VERDICT: CORRECT

# Independent Cold Audit of the Multi-Band Upper Certificate (`lambda = 5657/9983 < 17/30`)

Independent audit of commit `f515157` ([BOUND.md](../attempts/strogatz-astra/banded/BOUND.md), [check.cpp](../attempts/strogatz-astra/banded/check.cpp), [geometry.txt](../attempts/strogatz-astra/banded/geometry.txt), `recipe.txt.gz`, and [reproduce.sh](../attempts/strogatz-astra/banded/reproduce.sh)).

**Isolation protocol:** Our standalone C++17 verifier ([verify_banded.cpp](strogatz-astra-banded/verify_banded.cpp)) and exact-rational Python audit suite ([test_banded_audit.py](strogatz-astra-banded/test_banded_audit.py)) were designed, implemented, compiled with UBSan, tested on fixtures and mutation negative controls, and run to completion on the full 226,545,617-state production input at `5657/9983` **strictly from the mathematical description in `BOUND.md` (and baseline `TRANSFER.md`), without reading `check.cpp` or any proposal source**. Only after our full independent `PASS` was in hand did we inspect commit `f515157` and compare `check.cpp` and the `BOUND.md` diff line by line against our model.

**Summary:** The certificate at $\lambda = 5657/9983 = 17/30 - 1/299490 \approx 0.566663328 < 17/30$ is **mathematically sound and independently verified**. For every $n \ge 3$,
$$\mathbb{E}[C_n] \le \frac{1539738113358}{3908512428} \prod_{j=3}^{n-1}\!\left(1 + \frac{5657}{9983\,j}\right) = O\!\left(n^{5657/9983}\right),$$
establishing $\limsup_{n\to\infty} \frac{\log \mathbb{E}[C_n]}{\log n} \le \frac{5657}{9983} < \frac{17}{30}$ (and, by Markov + Borel–Cantelli on dyadic sizes and monotonicity of $C_n$, almost surely $\limsup_{n\to\infty} \frac{\log C_n}{\log n} \le \frac{5657}{9983} < \frac{17}{30}$).

---

## 1. Independent Full Verification Results (`verify_banded.cpp`)

Built with `g++-13 -std=c++17 -O2 -Wall -Wextra -fsanitize=undefined -fno-sanitize-recover=all` (zero warnings, zero UBSan errors), pinned to a single low-priority CPU (`nice -n 15`).

### 1.1 Production certificate at `lambda = 5657/9983`
Verified bit-for-bit equality of `attempts/strogatz-astra/banded/geometry.txt` and `gzip -dc attempts/strogatz-astra/banded/recipe.txt.gz` against the input files tested (`sha256sum` of `recipe.txt.gz` = `584bb93c2797138ae3de54cb85c7d8e487bb7d7cc5320c021b5a64ec0cf0f60a`, matching `BOUND.md` line 9).

- **Geometry:** `7` central-degree bands, `H = 2400`, constant `h_1 = 620`, `12,697` `c = 1` labels, `226,545,617` finite states:
  - Band 6 `[1000, 2400)`: `66,455,200` states (`45.61s`)
  - Band 5 `[500, 1000)`: `21,291,500` states (`12.56s`)
  - Band 4 `[200, 500)`: `22,446,000` states (`16.33s`)
  - Band 3 `[80, 200)`: `17,550,960` states (`12.91s`)
  - Band 2 `[32, 80)`: `22,401,552` states (`17.54s`)
  - Band 1 `[16, 32)`: `23,199,024` states (`19.12s`)
  - Band 0 `[2, 16)`: `53,201,381` states (`46.84s`)
- **Recipe:** `scale = 10^12`, `threshold = 8`, `63,309` mutually observable signed edge transfers (all mirror-checked across both endpoints and all witness variants; zero on pooled/ambiguous-saturated slots), `18` nonnegative arrow entries with `a_star = 398,437,500,000`.
- **Analytic tail (`D >= 2400`):** `A = 1,380,739,398,522`, `B = 383,976,882,193`, checked in exact 128-bit arithmetic against both `s_trial = 10^12` and `s_actual <= s_trial`:
  - $9983(B + a_*) \le 5657 A$ (exact slack `6,235 > 0`),
  - $9983 s_{\text{trial}} \le 5657(A + B)$ (exact slack `4,755 > 0`).
- **Newborn consistency:** `min(trial - actual) = 4,341,454` integer units ($4.34 \times 10^{-6}$ at worst `c = 1` index `286`: `actual = 1,015,179,606,101`, `trial = 1,015,183,947,555`).
- **Weights & prefactor (in units of $10^{-12}$):**
  - Finite weight range `[3,908,512,428, 1,764,716,218,642]`, tail `[1,380,739,398,522, 1,764,716,280,715]`, global `[3,908,512,428, 1,764,716,280,715]`.
  - Root ($T_3$ centre) weight `1,539,738,113,358`.
  - Exact $T_3$ scaled drift / root weight `872,508,277,036 / 1,539,738,113,358 = 0.566660180369` ($30\,\text{drift}_3 - 17\,\Phi_3 = -299,616,006 < 0$).
  - Minimum finite row slack `0` (non-strict inequality satisfied everywhere with upward-rounded drift and downward-rounded $\lambda z$).
  - **Result:** `PASS upper 5657/9983`, product prefactor `1539738113358/3908512428` (`= 393.944791`, Gamma prefactor `220.11`).

### 1.2 Side-by-side comparison with production `check.cpp` (`f515157`)
| Metric | Production `check.cpp` (`BOUND.md`) | Independent `verify_banded.cpp` | Notes |
|---|---:|---:|---|
| Finite states | `226,545,617` | `226,545,617` | Exact match |
| Signed edge transfers | `63,309` | `63,309` | Exact match |
| Tail `A + B` | `1,764,716,280,715` | `1,764,716,280,715` | Exact match |
| Min newborn reserve | `4,340,525` | `4,341,454` | `verify_banded` solves minimal `z*` without `check.cpp`'s `+up((terms+2)*q,p)` pad |
| Global minimum `m` | `3,908,512,529` | `3,908,512,428` | Agrees to 8 significant digits |
| Root weight $\Phi_3$ | `1,539,738,114,735` | `1,539,738,113,358` | Agrees to 9 significant digits |
| Verdict at `5657/9983` | `PASS` | `PASS` | Independent agreement |

### 1.3 Smoke calibration at `3/5` (`fixtures/bands-smoke.txt`, `fixtures/smoke-recipe.txt`)
- Verified `180,354` finite states + `D >= 48` analytic tail in `0.13s` (also tested via `reproduce_audit.sh` through `gzip -dc`), `min(trial - actual) = 55,885,959,070`, `PASS upper 3/5`, prefactor `1458767210441/75502835165`.

---

## 2. Exact-Tree Brute-Force Audit & Mutation Negative Controls

### 2.1 Independent exact-rational tree verification (`test_banded_audit.py`)
Implemented the full multi-band history-dependent state machine from scratch in Python `fractions.Fraction` on a 3-band geometry (`H = 22, T = 6`, cutoffs `[18,12,8,5,4]`, `[18,9,6,4,3]`, `[18,6,4,3,2]`) with random rational weights, signed edge transfers, nonnegative arrows, and analytic tail $(A, B)$:
- **All 436 exhaustive RFT histories ($n = 3..7$):** Exact rational equality between brute-force $n\mathbb{E}[\Phi(T_{n+1})-\Phi(T_n)\mid\mathcal{H}_n]$ and summed amended row drifts, exact zero net edge-transfer sum, exact zero net arrow-transfer sum, and exact upper-bound enclosure on **436/436 histories**.
- **50 random & adversarial multi-band trees:** Explicitly exercised band-boundary gains ($D = 9 \to 10$, $D = 15 \to 16$ preserving `c = 1` witness), entry into tail ($D = 21 \to 22 \to 24$), saturated $k = K_3 = 5$ sender arrows ($d = 7, k_{\text{true}} = 7 > 5$), and pooled receiver arrows with both $b < T$ and $b \ge T$. **All 50/50 passed identically.**
- **3 exact-tree negative controls (all caught):**
  1. Asymmetric signed edge transfer (`Edge cancellation failed: net=2997733/55440`).
  2. Omitting consistent-type max over $k' \ge K_c$ on saturated sender arrow (`Arrow max failed: true=500 > max=0`).
  3. Unsafely subtracting arrow on pooled neighbour with mark $b < T$ (`Finite row upper bound failed at u=0: 1297869/200 > 1232797/200`).

### 2.2 Black-box CLI & production mutation negative controls (all 16 rejected)
- **11 smoke-geometry mutations (`test_banded_audit.py`):** stricter $\lambda = 56/100$, understated reward, unobservable edge, negative arrow, arrow at $D \ge H$, missing reward, varying $h_1$, witness cutoff $h_3 < 6$, inflated arrow breaking bound, non-monotone bands, and corrupted edge breaking bound — **11/11 rejected**.
- **5 production `recipe2.txt` (`bands-2.txt`) mutations:** unobservable edge (`edge 50 2 90 2 1000`), negative arrow (`arrow 12 2 -1`), tail arrow (`arrow 2400 2 100`), missing reward (`reward 100`), and ambiguous saturated edge (`edge 20 3 10 5 1000`) — **5/5 rejected** immediately with non-zero exit code.

---

## 3. Line-by-Line Audit of Commit `f515157` (`check.cpp`, `BOUND.md`, `reproduce.sh`)

- **`check.cpp` (380 lines):** Audited every component against `BOUND.md` and our independent model:
  - `up(W n, I d)` (`n/d + (n%d > 0)` with signed `__int128_t`) is exact ceiling for all signed `n` and `d > 0`; `W(p)*value/q` is exact floor for positive `p, value, q`.
  - Geometry validation enforces `h_1` constant across bands, nonincreasing cutoffs in `c` and across bands, and `h_3 >= 6` on every band (ensuring exact witness initialization for $c \in \{1,2,3\}$).
  - Mutual edge observability enforces strict `d < u.cap && k < u.kcap && D < w.cap && c < w.kcap`, matching the exact `w.k < l.kcap` lookup in `solve()` (lines 279–282).
  - Arrow receiver handling (`adjusted = target - (d >= threshold ? arrow : 0)`, `lo = exact ? count : 0`, `hi = count`) and `receiver_max()` (consistent-type max over pooled $d \ge b$, saturated $k \ge K_c$, and exact $k$) match `BOUND.md` identically.
  - All 4 move-event rate bounds (`GAIN` witnessed/unwitnessed, `CONVERT`, `FAR_GAIN`, `FAR_CONVERT`), witness preservation on `c = 1` cross-band gain, witness swap to `(2, 1)` on neighbour conversion, and `newborn()` initialization match `BOUND.md` identically.
  - Sliding-window `weights[D%2]` with explicit runtime guard `row.degrees[D%2] == D` prevents reading any uncomputed successor; `W` (`__int128_t`) + `narrow()` + UBSan prevent overflow.
  - Tail inequalities `W(B + amax)*q <= W(p)*A` and `W(s)*q <= W(p)*(A+B)` are checked explicitly in 128-bit arithmetic, and final `gap >= 0` confirms `s_actual <= s_trial` and every actual newborn weight $\le$ its trial reward.
- **`BOUND.md` diff:** History-dependent filtration & tie-breaking (line 19), event exhaustion (line 23), explicit per-term upper bound (line 39), and the optional strict reserve improvement $\lambda M/(M+e)$ (line 51, verified algebraically via $\Psi_n = \Phi_n + e C_n$ and $\Phi_n \le M C_n$) are all mathematically sound.
- **Resolved packaging note:** In `f515157`, `BOUND.md` linked to `verified.log` before that transcript was tracked; commit `de98dc7` added [verified.log](../attempts/strogatz-astra/banded/verified.log), closing the broken relative link on a clean checkout.

---

## 4. Portable Independent Reproduction

All independent audit files live in [strogatz-astra-banded/](strogatz-astra-banded/) with zero hardcoded server paths:
- `verify_banded.cpp` — standalone C++17 UBSan exact verifier
- `test_banded_audit.py` — exact-rational tree & mutation negative-control suite
- `fixtures/bands-smoke.txt`, `fixtures/smoke-recipe.txt` — self-contained smoke calibration fixtures
- `reproduce_audit.sh` — one-command wrapper compiling with `g++-13` (UBSan), running smoke calibration + `test_banded_audit.py` + the full `5657/9983` production verification on `../../attempts/strogatz-astra/banded/geometry.txt` and `recipe.txt.gz` (failing if either production file is absent).

---

## 5. Final Packaging & Diff Sign-Off (`de98dc7`, `f7be3e0`)

**PACKAGING VERDICT: PASS**

Audited commits `de98dc7` and `f7be3e0` (and the two pending wording additions in the working tree on [RESULT.md](../attempts/strogatz-astra/RESULT.md) and [BOUNDS.md](../attempts/strogatz-astra/bounds/BOUNDS.md)):

1. **Report statements, finite-$n$ constants, and limits ([RESULT.md](../attempts/strogatz-astra/RESULT.md)):**
   - Exact rational arithmetic verified: $17/30 - 5657/9983 = 1/299490$ and $5657/9983 \approx 0.566663327657$.
   - Finite-$n$ upper bound for every $n \ge 3$,
     $$\mathbb{E}[C_n] \le \frac{1539738114735}{3908512529}\prod_{j=3}^{n-1}\!\left(1 + \frac{5657}{9983\,j}\right) = O\!\left(n^{5657/9983}\right),$$
     and independent prefactor $1539738113358/3908512428 \approx 393.944791457626$, together with the lower bound $(15430450607/17611596558)\prod_{j=3}^{n-1}(1+0.563/j)$ ([transfer.log](../attempts/strogatz-astra/bounds/transfer.log)), the three-vertex obstruction $557549077 > 0$ at $17/30$ for the $0.568$ weights, and the existence / almost-sure growth statements ($\mathbb{E}[C_n] = n^{\mu+o(1)}$ and $C_n = n^{\mu+o(1)}$ a.s. with $563/1000 \le \mu \le 5657/9983 < 17/30$), match every underlying proof and transcript verbatim.
   - Packaging metadata verified: `recipe.txt.gz` is `426,026` bytes (`416.041 KiB`, ceiling `417 KiB`) with SHA-256 `584bb93c2797138ae3de54cb85c7d8e487bb7d7cc5320c021b5a64ec0cf0f60a`; `de98dc7` added the 10-line explicit per-term bound paragraph and transcript block to [BOUND.md](../attempts/strogatz-astra/banded/BOUND.md), committed [verified.log](../attempts/strogatz-astra/banded/verified.log), and packaged our `checks/strogatz-astra-banded/` suite unmodified.
   - Pending edits checked: `RESULT.md` cleanly notes the stopped sampled-neighbour alternative (with its provisional numerical bounds excluded) and Carbon's two `HOLDS` attacks (Claude Opus 5.5); [BOUNDS.md](../attempts/strogatz-astra/bounds/BOUNDS.md) cleanly justifies the two $(2,2)$ witness origins (a degree-$2$ core had degree $2$ throughout its core lifetime and converted its leaf while of type $(2,1)$) and conditions on the full growth history $\mathcal{F}_n$; and the one-line `"${CXX:-g++-13}"` pin in [reproduce.sh](../attempts/strogatz-astra/bounds/reproduce.sh) and [reproduce-transfer.sh](../attempts/strogatz-astra/bounds/reproduce-transfer.sh) resolves to the identical `g++-13` executable as `g++` here while remaining overridable via `CXX`.
2. **Relative links & portability:**
   - All `16` relative markdown links in [RESULT.md](../attempts/strogatz-astra/RESULT.md) and both relative links in [BOUND.md](../attempts/strogatz-astra/banded/BOUND.md) resolve to git-tracked files at `f7be3e0`.
   - Zero hardcoded absolute filesystem paths or local URI schemes in `RESULT.md`, `attempts/strogatz-astra/banded/`, `checks/strogatz-astra-banded/`, or this review.
3. **Both reproduction entry points verified live from repository root (and cross-checked against clean-clone logs):**
   - `bash attempts/strogatz-astra/banded/reproduce.sh` (pinned low-priority single core, GCC 13.3.0 UBSan + gzip 1.12): exited `0`, output byte-for-byte identical to [verified.log](../attempts/strogatz-astra/banded/verified.log), `outputs/band-clean.log`, and `RESULT.md` lines 40–42 (`PASS upper 5657/9983`), with scratch directory cleaned up on exit.
   - `bash checks/strogatz-astra-banded/reproduce_audit.sh` (pinned low-priority single core, GCC 13.3.0 UBSan + Python stdlib + gzip 1.12): exited `0`, passed smoke calibration (`180,354` states at `3/5`), all `436` exhaustive histories + `50` random/adversarial trees + `3` exact-tree negative controls + `11` CLI negative controls, and the full `5657/9983` production check (`226,545,617` states + analytic $D \ge 2400$ tail in `167.87s`), emitting `PASS upper 5657/9983 prefactor 1539738113358/3908512428 (=393.944791457626)` verbatim matching `RESULT.md` line 52 and `outputs/independent-clean.log`.
