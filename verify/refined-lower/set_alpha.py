#!/usr/bin/env python3
"""Patch the LOWER-side comparability shift alpha_lo_num (header int64 at byte offset 48) in place.

alpha enters only the comparability accounting S*Phi = 2 alpha S + sum_u (W(u) + alpha S (c_u - 2))
(quintet-cert PROOF.md section 4); Phi, every row and every envelope are unchanged.  The unmodified
verifier recomputes everything from the patched bytes, so a patched file is just a new proposal.
Usage: set_alpha.py FILE NUM   (alpha = NUM / 2^30, 0 <= NUM <= 2^30; mode must be 2)
"""
import struct
import sys

path, num = sys.argv[1], int(sys.argv[2])
assert 0 <= num <= (1 << 30)
with open(path, "r+b") as fp:
    hdr = fp.read(512)
    mode = struct.unpack_from("<q", hdr, 8)[0]
    assert mode == 2, "alpha shift is only allowed on the LOWER side"
    old = struct.unpack_from("<q", hdr, 48)[0]
    fp.seek(48)
    fp.write(struct.pack("<q", num))
print(f"{path}: alpha_lo_num {old} -> {num}")
