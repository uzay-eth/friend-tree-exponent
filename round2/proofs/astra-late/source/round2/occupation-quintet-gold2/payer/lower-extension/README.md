> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical-floor and frontier audits and scoped Claude conditional reduction and arithmetic; original Neon-2 Gold2 witnesses and the four reward floors. See the [release index](../../../../../README.md#claim-28) for exact scope and citations.

# GOLD2 exact lower FY extension

**`67,007.80431291128` ticks/core across two disjoint masks**, an improvement of `14,666.09818091611` ticks/core over the previously booked `52,341.70613199517`. This is below the requested `84,000` ticks/core. The new small-donor amount replaces the previous payer credit; only the donor-`>5` component is additional to either small-donor price.

| Reward mask | Exact lower ticks/core | Display |
|---|---:|---:|
| Host `D≤251,c≤13`, payer `D=6..46,c≤7`, donor `D=2..5` all arities | `101444804329281757487394273/1550425090937500000000` | `65,430.31644820765` |
| Host `D≤150,c≤10`, payer `D=6..34,c≤7`, donor `D≥6` absorbing | `9783107064343700060933349/6201700363750000000000` | `1,577.487864703628` |

Both pay ordered distinct physical neighbours at exact G2 hosts outside M3, and require `max(D_host,D_payer)>18`. They are disjoint from the degree-18 full edge payment and from the `payer≤5` absorber. All exact-reward payer bins are singletons. The both-high donor interval is `[0,1/6]`; each donor degree until `5` is followed with all arities. The crossing from `5` to `6` uses the true gain event and enters a permanent abstract state. Its later transitions change neither its own abstract value nor the number of registered donor incidences. Its true reciprocal still enters the host gain-rate box. To lower-bound the reward after crossing, maximize the original affine-per-segment FY objective over the two interval endpoints and each original knot in the interval. No post-overflow promotion can register a host on the zero-extended singleton path; this is one conservative entry route. Existing registered paths keep the interval donor.

The code clears rational denominators in the original 20-knot maximum and uses exact signed integer arithmetic. The new Bellman step checks the adverse inequality after rounding, and its fixtures independently enumerate all rate corners on 18 configurations. The discount and intensity use the original exact GOLD2 endpoints, not the trial `99/175`. The unchanged lower-payer machinery supplies chronological registrations for mother-payer, mother-donor, and both two-child directions; a promotion both advances its prefix and registers each compatible new incidence with the post-promotion arity. The authenticated receipt and lower binary at `b435b4c` were checked before reading rewards, and stat stability was checked at exit.

From the Astra root, replay on CPU11:

```sh
P=round2/occupation-quintet-gold2/payer/lower-extension/replay.py
taskset -c 11 nice -n 10 python3 -Bu "$P" --host 251 --stock 13 --receiver 46 --rstock 7
taskset -c 11 nice -n 19 python3 -Bu "$P" --mask both-high --absorb --host 150 --stock 10 --receiver 34 --rstock 7 --donor 6 --tstock 5
```

Both printed `PASS exact lower replay and unchanged authenticated payload`. Small donor: `12,322,800` states, `2,587,032` exact reward labels, four exact knot fixtures, `125.41+3.14` CPU s (`473.43` wall). Both-high: `4,228,200` states, `129,078` reward labels, independent absorbing seam / raw spline / mask check, `40.29+2.91` CPU s (`290.24` wall). The earlier both-high absorbing fixture initially used an M3 label and failed closed; its replacement exact non-M3 fixture passed. The wider `--host 251 --stock 16 --receiver 80 --rstock 9` run timed out at 900 seconds after printing 25,387,110 path states; it never returned a newborn value and supplies no credit. The B5 lane and older exact payer files are separate.
