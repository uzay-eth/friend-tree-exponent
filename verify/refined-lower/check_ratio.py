#!/usr/bin/env python3
"""Light, fail-closed check of the refined lower endpoint of Theorem 2 (paper, Appendix F). Standard library only.

It does NOT replay the 13 GB witness (see replay_full.sh). It checks, in exact rational arithmetic:
  1. the receipt and the per-class file classes.txt are the pinned bytes (SHA-256);
  2. the receipt is a LOWER (mode 2) receipt for the pinned payload digest, ok = true, and its
     verify record is internally consistent (alpha, caps, higher_per_core, m_min, E_tot);
  3. the per-class file is well formed: exactly the rows c = 1..DF-1, each "c minlo hmax", no duplicates;
     the per-c data reproduce the receipt's min_W (at the receipt alpha) and higher_per_core;
  4. the classwise aggregation (Lemma F.1) recomputed from scratch gives exactly kappa = 1907/67944653537 and the endpoint
     X = 3843726496355222669/6794465353700000000 (and the pinned combiner gives the receipt's
     certified_exact, which is below 99/175);
  5. with the unchanged Gold2 upper endpoint U from the repository Gold2 receipt (pinned by SHA-256):
     99/175 < X, U < 142/251, no fraction with denominator <= 425 lies in [X, U], 241/426 does,
     outward 12-digit decimals, and the width reduction.
`--self-test` also feeds planted defects (changed c=2 cap, missing / duplicate / malformed row,
wrong mode, wrong digest) and requires each one to be rejected.
Usage: python3 -I check_ratio.py [--self-test] [--aux FRESH_CLASS_FILE]
Exit status 0 only if every check passes.
"""
import hashlib, json, math, os, sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
RECEIPT = os.path.join(HERE, "receipt.json")
AUX = os.path.join(HERE, "classes.txt")
GOLD2 = os.path.join(ROOT, "round2", "proofs", "neon2", "quintet-cert", "certificate_gold2.json")
PIN = {
    RECEIPT: "f46b5f8dd7b73bfec3d91fd80c83b9fa4e1d96f326a7f4581175ed41fb44dc17",
    AUX: "278c423e76650feff8766405c4ab766152fa1df2aed18cffb93d5bab6db6a612",
    GOLD2: "5e11f9c0b6eb3bd43e880be608c37d6f0b1fdbcc632b91bb5e8037db52e13284",
}
PAYLOAD_SHA = "80e3b80a9f962b458854828abab5b00e2806a6c1ce709835b5e884b0afd4ad28"
X_EXPECTED = F(3843726496355222669, 6794465353700000000)
KAPPA_EXPECTED = F(1907, 67944653537)
U_EXPECTED = F(127993535256830043, 226242120050000000)
L_GOLD2 = F(35083855339718511, 62017003637500000)


class Reject(Exception):
    pass


def need(ok, msg):
    if not ok:
        raise Reject(msg)


def parse_aux(text, DF):
    rows = {}
    for i, line in enumerate(text.split("\n")):
        if line == "" and i == len(text.split("\n")) - 1:
            continue
        parts = line.split(" ")
        need(len(parts) == 3 and all(p.lstrip("-").isdigit() and p != "-" for p in parts),
             f"per-class line {i + 1} malformed: {line!r}")
        c, mn, hm = map(int, parts)
        need(c not in rows, f"per-class duplicate row c={c}")
        need(hm >= 0, f"per-class negative hcap at c={c}")
        rows[c] = (mn, hm)
    need(set(rows) == set(range(1, DF)), f"per-class rows are not exactly c = 1..{DF - 1}")
    return rows


def lemma_r(rec, rows, alpha):
    """Return (kappa, classes) for the classwise aggregation (Lemma F.1) at comparability shift alpha (ticks)."""
    v = rec["verify"]; caps = v["caps"]; cf = v["comp_fin"]; t = rec["tail"]
    e_f = max(0, caps["f_in"], caps["f_tab"], t["cap_ticks"])
    e_g = max(0, caps["g"])
    cl = [("C=1", e_f, min(cf["min_W_c0"], t["w_c0_min_ticks"]))]
    for c in sorted(rows):
        mn, hm = rows[c]
        cl.append((f"finite c={c}", e_f + e_g + hm, mn + alpha * (c - 2)))
    cl.append(("tail c=1", e_f + e_g, t["w_c1_min_ticks"] + cf["g_far_min"] - alpha))
    cl.append(("tail c=2", e_f + e_g, t["w_c2_min_ticks"] + cf["g_far_min"]))
    cl.append(("tail c>=3", e_f + e_g, t["w_min_ticks"] + cf["g_far_min"] + alpha))
    need(all(m > 0 for _, _, m in cl), "a class has non-positive mass")
    need(all(e >= 0 for _, e, _ in cl), "a class has negative charge")
    return max(F(e, m) for _, e, m in cl), cl


def check(receipt_bytes, aux_bytes, gold2_bytes, out=print):
    need(hashlib.sha256(receipt_bytes).hexdigest() == PIN[RECEIPT], "receipt digest differs from the pin")
    need(hashlib.sha256(aux_bytes).hexdigest() == PIN[AUX], "per-class file digest differs from the pin")
    need(hashlib.sha256(gold2_bytes).hexdigest() == PIN[GOLD2], "Gold2 certificate digest differs from the pin")
    full = json.loads(receipt_bytes)
    rec = full["lower"]; v = rec["verify"]; t = rec["tail"]
    need(full["lower_sha256"] == PAYLOAD_SHA and v["sha256"] == PAYLOAD_SHA, "receipt is not for the pinned payload")
    need(rec["ok"] is True and rec["mode"] == 2 and v["mode"] == 2 and t["mode"] == 2, "not an accepted LOWER receipt")
    need(v["format"] == "FT2_U5V1", "format")
    lam = F(v["P"], v["Q"])
    need(rec["nominal"] == f"{v['P']}/{v['Q']}" and lam == F(56571437, 10**8), "nominal")
    need(v["S"] == 2**40 and v["S_r"] == 2**30 and v["alpha_ticks"] % (v["S"] // v["S_r"]) == 0, "scales")
    alpha = v["alpha_ticks"]
    need(0 <= alpha <= v["S"], "alpha outside the validated range [0, 1]")
    DF = t["DF"]
    need(DF == 4096, "DF")
    rows = parse_aux(aux_bytes.decode("ascii"), DF)
    # consistency of the per-c dump with the trusted (unmodified-verifier) receipt fields
    need(min(mn + alpha * (c - 2) for c, (mn, _) in rows.items()) == v["comp_fin"]["min_W"], "per-class minima vs receipt min_W")
    need(max(hm for _, hm in rows.values()) == v["higher_per_core"] == rec["higher_per_core"], "per-class caps vs receipt higher_per_core")
    e_f = max(0, v["caps"]["f_in"], v["caps"]["f_tab"], t["cap_ticks"])
    E = max(0, e_f + max(0, v["caps"]["g"]) + v["higher_per_core"])
    need(E == rec["err_per_core"], "E_tot")
    kappa, cl = lemma_r(rec, rows, alpha)
    need(min(m for _, _, m in cl) == rec["m_min"], "class masses vs receipt m_min")
    pinned = lam - F(E, rec["m_min"])
    need(pinned == F(rec["certified_exact"]), "pinned combiner endpoint vs receipt certified_exact")
    need(pinned < F(99, 175), "pinned combiner unexpectedly clears 99/175")
    X = lam - kappa
    binding = max(cl, key=lambda r: F(r[1], r[2]))
    need(kappa == KAPPA_EXPECTED and X == X_EXPECTED, f"classwise endpoint differs: {X}")
    out(f"receipt, per-class file, Gold2 digests pinned; LOWER receipt ok, mode 2, payload {PAYLOAD_SHA[:8]}...")
    out(f"per-class file: {len(rows)} rows c=1..{DF - 1}, consistent with receipt min_W={v['comp_fin']['min_W']}, higher_per_core={v['higher_per_core']}")
    out(f"pinned combiner: {pinned} = {float(pinned):.12f} < 99/175 (does not exclude 99/175)")
    out(f"classwise aggregation (Lemma F.1): kappa = {kappa}, binding class {binding[0]} (e, m) = ({binding[1]}, {binding[2]})")
    out(f"X = {X}")
    g2 = json.loads(gold2_bytes)
    U = F(g2["upper"]["certified_exact"]); L2 = F(g2["lower"]["certified_exact"])
    need(g2["upper"]["ok"] is True and g2["upper"]["mode"] == 1 and U == U_EXPECTED and L2 == L_GOLD2, "Gold2 endpoints")
    need(F(99, 175) < X < U < F(142, 251), "ordering 99/175 < X < U < 142/251")
    need(175 * 142 - 99 * 251 == 1, "Farey neighbours")
    small = [F(a, b) for b in range(1, 426) for a in range(math.ceil(X * b), math.floor(U * b) + 1) if X <= F(a, b) <= U]
    need(not small, f"a fraction with denominator <= 425 lies in [X, U]: {small[:3]}")
    need(X <= F(241, 426) <= U, "241/426 in [X, U]")
    lo = F(math.floor(X * 10**12), 10**12); hi = F(math.ceil(U * 10**12), 10**12)
    need((lo, hi) == (F(565714341933, 10**12), F(565736986679, 10**12)), "outward decimals")
    red = 1 - (U - X) / (U - L2)
    need(F(364, 10000) <= red < F(365, 10000), "width reduction")
    out(f"99/175 < X: margin {X - F(99, 175)}; U < 142/251: margin {F(142, 251) - U}")
    out(f"outward [{float(lo):.12f}, {float(hi):.12f}]; width {float(U - X):.4e} vs Gold2 {float(U - L2):.4e} (reduction {float(red) * 100:.2f}%)")
    out("no fraction with denominator <= 425 in [X, U]; 241/426 is in [X, U] (not asserted to be mu)")
    return X


def self_test(rb, ab, gb):
    """Planted defects: each must be rejected. Digest pins are bypassed here so that the
    later checks are exercised as well."""
    saved = dict(PIN)
    bad = 0
    def variant(name, r=rb, a=ab, keep_pins=False):
        nonlocal bad
        if not keep_pins:
            PIN[RECEIPT] = hashlib.sha256(r).hexdigest(); PIN[AUX] = hashlib.sha256(a).hexdigest()
        try:
            check(r, a, gb, out=lambda *_: None)
            print(f"  self-test {name}: NOT rejected"); bad += 1
        except (Reject, ValueError, KeyError, json.JSONDecodeError) as e:
            print(f"  self-test {name}: rejected ({e})")
        PIN.update(saved)
    lines = ab.decode().split("\n")
    def aux_with(f):
        return "\n".join(f(list(lines))).encode()
    variant("c=2 hcap set to 0", a=aux_with(lambda L: [("2 135889307074 0" if x.startswith("2 ") else x) for x in L]))
    variant("row c=7 removed", a=aux_with(lambda L: [x for x in L if not x.startswith("7 ")]))
    variant("row c=7 duplicated", a=aux_with(lambda L: L[:7] + [L[6]] + L[7:]))
    variant("malformed row", a=aux_with(lambda L: [("5 12x 3" if x.startswith("5 ") else x) for x in L]))
    variant("upper-side receipt", r=rb.replace(b'"mode": 2', b'"mode": 1', 1))
    variant("ok = false", r=rb.replace(b'"ok": true', b'"ok": false', 1))
    variant("other payload digest", r=rb.replace(PAYLOAD_SHA.encode(), b"0" * 64))
    variant("per-class bytes not the pinned ones", a=ab.replace(b"\n2 135889307074 1941", b"\n2 135889307074 1940"), keep_pins=True)
    return bad


if __name__ == "__main__":
    # --aux PATH checks a freshly dumped per-class file instead of the stored one; it must be the pinned bytes.
    aux_path = sys.argv[sys.argv.index("--aux") + 1] if "--aux" in sys.argv else AUX
    rb, gb = open(RECEIPT, "rb").read(), open(GOLD2, "rb").read()
    ab = open(aux_path, "rb").read()
    try:
        check(rb, ab, gb)
    except Reject as e:
        print("FAIL:", e); sys.exit(1)
    if "--self-test" in sys.argv:
        if self_test(rb, ab, gb):
            print("FAIL: a planted defect was accepted"); sys.exit(1)
    print("RESULT: PASS")
