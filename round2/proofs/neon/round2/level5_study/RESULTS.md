# Level-5 prototype benchmark results

**Authors:** Gemini 4 Neon (lead) and GPT-6.2 Sol (prototype and checks). **Archive pin:** `c2498cc482bfef861644092d33f727cbbc93b280`. **Review and dependencies:** development checks and Level-4 baseline are detailed in the [study report](REPORT.md); these tables report numerical benchmarks.

| D0 | Case | λ lo | λ hi | Width | Reduction | Wall lo/hi (s) | Peak RSS lo/hi (KiB) | Max last-two-sweep shift |
|---:|:---|---:|---:|---:|---:|:---|:---|---:|
| 12 | base | 0.56565788754263469 | 0.56579149436263099 | 0.00013360681999630142 | 0.000% | 215.76/168.31 | 661528/661656 | 2.6e-10 |
| 12 | center | 0.56568394486865603 | 0.56576480694192277 | 8.0862073266740886e-05 | 39.478% | 173.29/147.32 | 661780/661652 | 9.4e-11 |
| 12 | d5-6 | 0.56567495571581106 | 0.56577823581536657 | 0.00010328009955551387 | 22.698% | 174.15/156.25 | 661652/661532 | 1.5e-10 |
| 12 | d5-7 | 0.56567645783386855 | 0.56577698989866121 | 0.00010053206479265775 | 24.755% | 182.01/185.27 | 661652/661652 | 1.4e-10 |
| 12 | d5-8 | 0.56567762241469222 | 0.56577601440211078 | 9.8391987418566806e-05 | 26.357% | 224.72/199.80 | 661652/661780 | 1.4e-10 |
| 12 | h-center | 0.56566754344890324 | 0.5657804535430605 | 0.00011291009415725917 | 15.491% | 198.41/203.95 | 661652/661780 | 2.4e-10 |
| 12 | l5-6 | 0.56566050397627932 | 0.56578850648188217 | 0.00012800250560285331 | 4.195% | 198.95/202.56 | 674024/674276 | 1.5e-09 |
| 12 | l5-6-no-inner-slack | 0.56566218554896741 | 0.56578692013676202 | 0.000124734587794606 | 6.641% | 175.48/195.55 | 674024/674024 | 1.5e-09 |
| 12 | l5-6-no-mixed-slack | 0.5656846387582406 | 0.56576972204110343 | 8.5083282862830067e-05 | 36.318% | 130.38/133.92 | 674280/674024 | 5.8e-10 |

`sample_build`: `TIME wall=284.64 user=278.81 sys=1.07 peak_KiB=667092 exit=0`.

`sample_verify`: `TIME wall=110.42 user=109.46 sys=0.45 peak_KiB=647424 exit=0` — `=== VERDICT: PASS (total_fails = 0, exact P/Q = 565655/1000000) ===`.

GMP benchmark: 82,492,914 printed row slots; 109.91 CPU seconds; 1332.4 ns per printed slot (includes startup and comparison).

Approximate CPU hours per endpoint by original table size:
- D0=20 production: 0.52 h
- D0=22 production caps: 0.80 h
- D0=24 base: 0.65 h
- D0=24 large: 1.24 h
- D0=25 large: 1.50 h

Real clipped L5 (D5=6): width change 4.195% narrower versus D0=12 baseline (below the 30% decision threshold).
- Suppressing all-inner residual: width falls by 3.2679178082473115e-06 (2.553% of the real L5 width; non-rigorous sensitivity).
- Suppressing mixed residual: width falls by 4.2919222740023244e-05 (33.530% of the real L5 width; non-rigorous sensitivity).
