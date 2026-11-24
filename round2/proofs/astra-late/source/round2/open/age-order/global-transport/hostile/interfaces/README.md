> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT, Claude, and Gemini within their recorded scopes for hubhood, with GPT-only audits of the accepted positive-threshold proof; threshold addenda use Claim 24 or pinned Fable finite transports where stated. See the [release index](../../../../../../../README.md#claim-26) for exact scope and citations.

# Discharging low components through high boundary ports

At hubhood, the mixed-stock result extends to trees with the following low components. In the old core tree call total degrees 2–9 low and degrees at least 10 high. Each connected low component may be:

* a path of cores with exactly two core slots each, joining two high ports (any length and any live leaf stocks);
* a path from one high port to a low terminal, with exactly two core slots at each nonterminal;
* one low core with any number of slots, all joining high cores.

High cores may have arbitrary stock and branching, and include high terminals. This is an arbitrary-size conditional comparison, but does not cover all branching low components. No negative physical score was found.

## Charges and port bounds

Use the [mixed-stock credit](../mixed-stock/README.md) `C(d,ds)` after paying old leaf edges and serial tips. Adjacent low profiles with coefficients `s,t`, outgoing first-send bounds `U,V` and favourable leaf lower bounds `L_d,L_e` pay the shared signed debit

```
max(0, (s−t)*(U−L_e)_+, (t−s)*(V−L_d)_+).
```

Here both stocks and both incident reciprocals are read in the profiles. At a high port the high side's outgoing exposure is jointly bounded by `s_high U_high < .012`. Indeed for `10≤d≤41` verify `(d+l)U(d+1,l)/[2d(d+1)]<.012` for every stock `0≤l<d`; above that range use `1/[2(d+1)]`. Drop the favourable high tip and the coefficient subtractions at the boundary only. The resulting cap for the mixed edge is `max(s_low U_low, .012*(1−L_low))`, because `U_high≤1`. This deliberately relaxes each high endpoint separately but does not multiply a maximum high coefficient by an independently maximized high send bound.

Each high port contributes its actual reciprocal to the incident low coefficient. Conditional on all low degrees and stocks, the low credits are concave functions of these reciprocals and all debit bounds are convex: their signed slopes are constant on the profile once stock is fixed. Thus the component lower charge is concave jointly on the box of high reciprocals `[0,1/10]`. It suffices to check its corners, independently at every port. Replacing every reciprocal by the same endpoint need not suffice on a multi-vertex component.

## Path transfer and junctions

Use the state `(previous degree, current degree, next degree)`, with `0,10` denoting the high reciprocal endpoints. The 800 states encode all allowed low degrees and hence all stocks when the current vertex has two slots. Transitions match the two shared degrees and pay the exact paired low–low debit above, rounded upward; each entering vertex earns its credit rounded downward. The initial labels pay the left high port; for a connector the terminal read pays the right high port. For a pendant path the terminal read instead adds the one-slot low terminal credit and its paired edge debit. A singleton connector is covered by a state with two high endpoints; a one-slot singleton low terminal uses the [existing terminal inequality](../mixed-stock/README.md).

The shortest-path relaxation stabilizes at step two, and the verified fixed labels are lower bounds for **every finite path**, by induction along the path. On scale `10^10`, the worst connector label is `−341943716` and the least pendant label is `930252329`. Thus every connector costs less than `b=1/22`, and every pendant pays surplus greater than `a=.00634`. The separate bare-path check computes the corresponding short chains directly as a control. A singleton internal junction with `c≥2` high ports costs less than `(c−1)b`: its worst excess above this charge is `0.0228065105…`. All its reciprocals enter just their sum, so two endpoint profiles suffice for each `(d,c)`.

## Reassembly

Delete pendant low components and low singletons with one port after crediting `a` for each. Treat a connector or singleton junction with `k≥2` high ports as a tree connection of `k−1` links, each charged at most `b`. Direct high–high links are charged at most `b` as before. Replacing each such component by a tree on its distinct high neighbours yields a tree on the high vertices: contraction preserves connectivity and any duplicate link or cycle would imply a cycle in the original core tree. Each high port counts once toward the original high core's census, even when a junction produces multiple links in the replacement. Equivalently, keep the junction as a zero-credit auxiliary vertex with census `k` and assign its charge `(k−1)b` to the `k−1` net links it adds. For the high internal credits the needed census identity is the original one: the total number of terminals is `2+Σ_internal(c−2)`. All nonterminal low path vertices have `c=2`; a low junction with `k>2` adds `k−2` to that sum and only helps. Hence there are at least `2+Σ_high_internal(c−2)` terminal payments of `a`, including high terminals, whose incident path/link is covered by the `b` account. For a high terminal its large-stock credit is greater than `b+a`, by the previous scalar bound. Counting at most `B−1` net links among high cores when there are `B` high cores, give each high vertex one `b` credit (subtract `a(c−2)` only at internal vertices), and pay the terminal credits separately. The resulting numerator is greater than `b+2a`. All-low and single-core cases use the earlier certificates.

## Replay and unpaid cases

```
python3 -IB round2/open/age-order/global-transport/hostile/interfaces/check.py
python3 -IB round2/open/age-order/global-transport/hostile/interfaces/check.py 9
```

The default runs the exact port scalar checks, the 800-state path transfer and 72 singleton junction endpoints. The second command is an exploratory check of a sufficient **terminal-child template** on general low profiles with high slots. It checks 350,056 oriented local rows and returns negative margins for that template (high `−691769`, low `−13917006`, scale `10^8`). These are not upper bounds on the physical score. In particular, the template fails even on a pendant profile already paid by the component transfer. A branching low component with several low–low edges and high ports still needs a different subtree potential or a physical upper certificate to resolve its actual sign.
