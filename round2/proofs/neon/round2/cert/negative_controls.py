#!/usr/bin/env python3
"""
Negative-control test suite for cert_b3f_verify_gmp.
Verifies that the GMP verifier REJECTS (non-zero exit code):
  1. Impossible exponent override (--override-pq) in the tighter direction.
  2. Corrupted magic number, invalid mode, invalid KT < D0 - 1, truncated file.
  3. Non-zero g[INF, INF].
  4. Single-entry adversarial perturbations in EVERY table:
     f, g, hV, hC, qP4, qYt, qYb, qS3, c2_f, c3_f, Yp_ref, Yw_ref.
"""
import sys, os, struct, subprocess, tempfile, shutil

def run_verifier(verifier, path, extra=None):
    cmd = [verifier, path, "--threads", "4", "--tree-check", "0"] + (extra or [])
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    return r.returncode, r.stdout

def main():
    if len(sys.argv) < 3:
        print("Usage: negative_controls.py <verifier_bin> <cert.bin>"); sys.exit(2)
    verifier, cert_path = sys.argv[1], sys.argv[2]

    with open(cert_path, "rb") as fp:
        raw = bytearray(fp.read())
    (magic, mode, version, P, Q, S, S_r,
     LF, KF, D0, KT, DP4, KP4, DYT, KYT, DYB, KYB, DS3, KS3) = struct.unpack("<Qiiqqqqiiiiiiiiiiii", raw[:96])

    # Compute exact byte offsets of every table (must match cert_b3f_format.h)
    NLF, NKF = LF + 1, KF + 1
    NV = 2 * NKF * NLF + 2
    NT = 2 * (KT + 2)
    INNER_N = sum(D - h + 1 for h in (1, 0) for D in range(2, D0))
    INF_ID = INNER_N
    NG = INNER_N + 1 + 3 * NT
    npars = 1 + sum(1 for h in (1, 0) for D in range(2, D0) for k in range(1, D - h + 1)) + 3 * 2 * (KT + 1)
    nchs  = 1 + sum(D for D in range(2, D0)) + 3 * (KT + 2)
    NHU   = sum(1 for D in range(2, D0) for k in range(1, min(KT, D - 1) + 1)) + 3 * KT
    NCU   = sum(1 for h in (1, 0) for D in range(2, D0) for k in range(2, min(KT, D - h) + 1)) + 3 * 2 * (KT - 1)
    NP4   = sum(max(0, DP4 - k - 1) for k in range(1, KP4 + 1))
    NYT   = sum(max(0, DYT - k - 1) for k in range(2, KYT + 1))
    NYBU  = sum(max(0, DYB - k - h) for h in (1, 0) for k in range(2, KYB + 1))
    NYBW  = sum(max(0, DYB - k - 1) for k in range(1, KYB + 1))
    NS3   = sum(max(0, DS3 - k - h) for h in (1, 0) for k in range(3, KS3 + 1))

    lens = [
        ("rpar",   INNER_N), ("rch",    INNER_N), ("trpar",  NT), ("trch",   NT),
        ("c2_f",   NT),      ("c3_f",   NT),      ("Yp_ref", NT * nchs), ("Yw_ref", npars * (KT + 2)),
        ("f",      NV),      ("g",      NG * NG), ("hV",     npars * NHU * nchs),
        ("hC",     NCU * nchs * (nchs + 1) // 2),
        ("qP4",    npars * NP4 * NP4 * nchs),
        ("qYt",    npars * NYT * nchs * (nchs + 1) // 2),
        ("qYb",    NYBU * NYBW * nchs * nchs),
        ("qS3",    NS3 * nchs * (nchs + 1) * (nchs + 2) // 6),
    ]
    offs = {}
    pos = 96
    for name, n in lens:
        offs[name] = (pos, n)
        pos += 8 * n
    assert pos == len(raw), f"Offset calculation mismatch: {pos} != {len(raw)}"

    tmpdir = tempfile.mkdtemp(prefix="b3f_neg_")
    failures = 0
    total = 0

    def expect_reject(label, data=None, extra=None):
        nonlocal failures, total
        total += 1
        p = cert_path
        if data is not None:
            p = os.path.join(tmpdir, f"mut_{total}.bin")
            with open(p, "wb") as fp:
                fp.write(data)
        rc, out = run_verifier(verifier, p, extra)
        if rc != 1:
            print(f"  [UNEXPECTED RC={rc}!] {label}")
            if out.strip():
                print("   ", out.strip().splitlines()[-1])
            failures += 1
        else:
            last = out.strip().splitlines()[-1] if out.strip() else f"rc={rc}"
            print(f"  [REJECTED OK] {label} -> {last}")
        if data is not None:
            os.remove(p)

    # 1. Impossible exponent override (tighter than true mu ~ 0.5657237)
    bad_P = (P * 99) // 100 if mode == 1 else (P * 101) // 100
    expect_reject(f"override-pq impossible ({bad_P}/{Q})", extra=["--override-pq", str(bad_P), str(Q)])

    # 2. Header corruption & generic integer-range guard tests
    m = bytearray(raw); m[0] ^= 0xFF
    expect_reject("corrupted magic", data=m)
    m = bytearray(raw); struct.pack_into("<i", m, 8, 3)
    expect_reject("invalid mode=3", data=m)
    m = bytearray(raw); struct.pack_into("<i", m, 52, D0 - 2)
    expect_reject("invalid KF < D0", data=m)
    m = bytearray(raw); struct.pack_into("<i", m, 60, D0 - 2)
    expect_reject("invalid KT < D0-1", data=m)
    m = bytearray(raw); struct.pack_into("<q", m, 40, (1 << 63) - 1)
    expect_reject("out-of-range S_r = INT64_MAX", data=m)
    m = bytearray(raw); struct.pack_into("<i", m, 48, (1 << 31) - 1)
    expect_reject("out-of-range LF = INT32_MAX", data=m)
    m = bytearray(raw); struct.pack_into("<i", m, 60, (1 << 31) - 1)
    expect_reject("out-of-range KT = INT32_MAX", data=m)
    m = bytearray(raw); f_off, _ = offs["f"]; struct.pack_into("<q", m, f_off + 16, 2 * 10**15)
    expect_reject("out-of-range table entry f[2] = 2e15", data=m)
    expect_reject("truncated file (-8 bytes)", data=raw[:-8])

    # 3. Non-zero g[INF, INF]
    m = bytearray(raw)
    g_off, _ = offs["g"]
    struct.pack_into("<q", m, g_off + 8 * (INF_ID * NG + INF_ID), 1)
    expect_reject("nonzero g[INF, INF] = 1", data=m)

    # 4. Single-entry perturbation in EVERY table (delta large enough to break inequality on an active interior entry)
    delta = -S if mode == 1 else +S  # breaks upper (decreases v so -lam*v increases) / breaks lower (increases v)
    targets = [
        ("f",      2),              # inner f(2,0,1)
        ("g",      0),              # inner g(0,0) -> wait: p=0 has k=0 (not a parent!), use an active (p,w)!
        ("hV",     0),
        ("hC",     0),
        ("qP4",    0),
        ("qYt",    0),
        ("qYb",    0),
        ("qS3",    0),
        ("c2_f",   1),              # k=1, h=1
        ("c3_f",   1),
        ("Yp_ref", 1 * nchs + 0),   # kap=1, h=1, w=chs[0]
        ("Yw_ref", 0 * (KT + 2) + 1), # p=pars[0], kap=1
    ]
    # Find an active (p, w) index in g: first inner parent (k=1, L=0, h=1 -> D=2) and first inner child (L=1, k=0, h=1 -> D=2)
    # In lut order h=1, D=2: k=0 (L=1) is idx 0, k=1 (L=0) is idx 1 -> (p=1, w=0) is valid!
    targets[1] = ("g", 1 * NG + 0)

    for tname, idx in targets:
        base, n = offs[tname]
        assert 0 <= idx < n
        # Perturb by +delta and/or -delta (for c2_f/c3_f/Yp_ref/Yw_ref, sign of transfer flips across entries, so huge +10*S definitely breaks at least one row!)
        bump = (10 * S) if tname in ("c2_f", "c3_f", "Yp_ref", "Yw_ref") else delta
        m = bytearray(raw)
        old = struct.unpack_from("<q", m, base + 8 * idx)[0]
        struct.pack_into("<q", m, base + 8 * idx, old + bump)
        expect_reject(f"perturb {tname}[{idx}] by {bump:+d}", data=m)

    shutil.rmtree(tmpdir)
    print(f"Negative controls summary ({cert_path}): {total - failures}/{total} rejected cleanly (failures={failures})")
    sys.exit(0 if failures == 0 else 1)

if __name__ == "__main__":
    main()
