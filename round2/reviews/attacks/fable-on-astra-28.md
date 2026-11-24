VERDICT: HOLDS (the four reward floors and the route-A frontier are checked inside the GPT family only; see the scope)

# Fable on Astra claim 28 (corrected 16:53 UTC: if mu is rational, its reduced denominator is at least 426; lower bound above 99/175 from four stationary reward floors on Neon-2's gold2 witness)

**Scope of the verdict.**
- Attacker: the Fable lead, Claude Opus 5.5.
- Pins: c488e00 for the four-reward certificate (route B) and e0708ef for the corrected two-sided interval (route A).
- The upper end comes from Neon-2's original gold2 certificate, which I attacked separately: HOLDS (`fable-on-neon2-6.md` section 6).
- Done by me (sections 1-3):
  - the Farey step;
  - the stationary drift argument;
  - the finite lifetime ceiling, recomputed;
  - the convexity rule for simultaneous siblings;
  - the replay of route B's checker.
- Not done by me: the four reward floors themselves. These are the singleton lifetime price, the two Fenchel-Young edge prices and the actual-sibling price. They were checked by Astra's own cold reviews and by a GPT-6.2 Sol Ultrafast child of this lead (section 4). Both checks are inside the GPT family.

## 1. Farey step and endpoints (lead)

- 142·175 − 99·251 = 1. So any reduced a/b strictly between 99/175 and 142/251 satisfies b = 251(175a − 99b) + 175(142b − 251a) ≥ 251 + 175 = 426.
- The gold2 upper end 0.565736986679 lies below 142/251 = 0.565737052, by 6.5e-8.
- The smallest denominator inside the corrected interval (0.565714467904, 0.565735974375) is 241/426, found by search in exact rationals.
- So the denominator statement follows once mu > 99/175 is established.

## 2. Route B, recomputed by the lead with independent code

The program is `astra28-lifetime_recompute.py`. Its inputs are the authenticated receipt (sha256 5e11f9c0...) and the 299-entry local upper table u(D, c).

**Lifetime recursion.**
- It recomputes the stopped discounted lifetime V(D, c) of B = m_max − u(D, c) on 2 ≤ D ≤ cut. The gain rate ranges over [L, L + c/2], the promotion rate is L/D, the discount is the gold2 upper end, and exits pay 0.
- It then forms M = ceil(m_max − lo·V(2,1)), the endpoint lambda + r/M, and the margin r − (99/175 − lambda) M.
- The generator inequality B + g(V(D+1,c) − V) + (L/D)(V(D,c+1) − V) − s V ≥ 0 holds at 61 gain rates in each interval, not only at the two ends.

**Results.** These are identical to Astra's table:

| retained states | V(2,1) | M | lower endpoint | margin at 99/175 (ticks/core) |
|---:|---:|---:|---:|---:|
| 20 | 619269501408 | 1458274085072 | 0.565714288946 | 4,712.6 |
| 170 | 680528936614 | 1423618796411 | 0.565714303770 | 25,703.8 |
| 299 | 686326891945 | 1420338814887 | 0.565714305210 | 27,690.5 |

**Sensitivity.**
- All four reward floors (207,990; 395,261; 52,341; 251,635; total 907,227) are essential. Removing any one makes the margin negative; without the smallest it is −24,650.5 at 299 states.
- The total slack is 27,690.5 ticks/core, about 3.1% of the floors.
- Route B's lower end exceeds 99/175 by only 1.95e-8. Route A's corrected interval clears it by 1.82e-7.
- So the conclusion is exactly as sound as the four floors.

**Drift argument, checked.** Write Q(t) = e^{−mu t} E P(T_t) for the lower potential.
- The original certificate gives, pointwise, a signed drift of at least −19,217 ticks per core.
- The four rewards are stationary averages of disjoint nonnegative terms of that drift. The two FY terms have disjoint payer degrees; the singleton and the paths are different terms.
- In the Cesàro limit this gives (mu − lambda) lim Q ≥ r c_E, with c_E = lim e^{−mu t} E C > 0.
- Together with Q ≤ M c_E, this gives mu ≥ lambda + r/M.
- The ceiling uses π B ≥ mu V(2,1) + (s − mu) π V ≥ lo V(2,1). This needs V ≥ 0, s ≥ mu and B ≥ 0, all of which hold.

**Simultaneous siblings, checked.** F(z) = [z − P_0]_+ − [−P_0]_+ is convex with F(0) = 0, hence superadditive on z ≥ 0. With P ≤ B − Σ z_i this gives Σ F(z_i) ≤ F(Σ z_i) ≤ [−P]_+. So the sibling and upstream payments never exceed the negative part of the actual row, and P + Σ F(z_i) ≤ max(P, 0) stays under the original production cap.

## 3. Replay (lead)

On a clean archive of c488e00, `check.py` passed (0.2 s): "PASS four integer reward floors, stationary potential bounds and exact Farey exclusion through denominator 425". The `--rebuild-host` replay, which rebuilds the 299 local ceilings from the 6.3 GB gold2 lower payload, was run by the child (section 4).

## 4. GPT-family child and Astra's cold reviews

Report: `round2/checks/astra28-audit-gpt-family.md`, "A-lower (`e0708ef`): HOLDS. B (`c488e00`): HOLDS. Denominator statement: HOLDS". Its programs and logs are archived as `round2/rival-attacks/code/astra28-child-*`.

**Hashes and route B.**
- Both original gold2 payloads matched their full SHA-256 hashes, and were read in place.
- Route B: both README commands passed, including `PASS rebuilt 299 local ceilings from the authenticated Gold2 lower payload`.

**Reward floors.**
- The singleton, small-payer FY, wide-payer FY and sibling producers were replayed from the authenticated payloads. Each matched its published amount (singleton 207990.215716345).
- Line audit of the four floors:
  - The two FY streams use disjoint payer degrees: at most 5, and 6 to 34.
  - Each extra maternal child is registered in exactly one birth order.
  - The sibling reward is rounded down before clock scaling.
  - The convexity rule prevents overspending the H3 row.
  - The corrected mother-80-100 seam is disjoint from the compact 19-79 payments.

**Route A.**
- `host/frontier/check.py --live-sibling --replay-singleton` reproduced the receipt byte for byte (SHA-256 7a6ddda0...) and the lower end 0.565714467904.
- At that end its fixed margin is 27.02 ticks/core. At 99/175 the restricted margin is about 247,577 ticks/core.
- No use of the withdrawn H3 price remains in either route. A cold helper parses the historical H3 amount and then discards it.

**Plants.** All three were rejected:
- the small-payer FY floor +10%;
- the withdrawn extra rate-denominator division, which makes the physical residual −1/2;
- a doubled divisor on the 20-state saving, which makes the margin negative.

**Independent recomputation.** The child's own code reproduces the lifetime table of section 2 exactly, as mine does.

**Limits stated by the child.** It inherited the original gold2 production verification and the physical stable-age theorem, and did not audit route A's upper half (which the denominator statement does not use).

Astra's own cold reviews (`gold2-rational-repair-sol`, `gold2-frontier-sol`, `gold2-lower-sibling-sol`, `gold2-singleton-repair-sol`) are also GPT-family.

## 5. Conclusion

**Denominator.** If mu is rational, its reduced denominator is at least 426. I found no error in this conclusion or in its reduction. That reduction is:
- the Farey step;
- the stationary drift argument;
- the finite lifetime ceiling, recomputed independently;
- the rule for simultaneous siblings.

**What the cross-family check does not reach.** The four reward floors are each essential (the total slack is 3.1%). They and route A's frontier have been checked only inside the GPT family, by Astra's cold reviews and by this lead's GPT child, which also replayed them from the authenticated payloads.

**How this lead states the result.**
- Proved with a cross-family check: denominator at least 175 (from gold2).
- Denominator at least 426: holds on the evidence above; its reward floors carry a GPT-family check only.

**Not claimed here.** Nothing about rationality, irrationality or transcendence.
