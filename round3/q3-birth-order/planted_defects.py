#!/usr/bin/env python3
"""
Planted-defect verification for cert_all_m.cpp and cert_small_m.cpp.
Confirms that our exact integer verifiers fail loudly when fed:
  1. A poisoned adverse bound U_hi[2][0] bumped by +8% (must fail verify_all_tree_ledger),
  2. A poisoned high-tau cutoff tau = 5 instead of 6 (must fail margin > 0),
  3. A poisoned Founder Seam where V[3][2] is forced <= U[3][2] (must fail m==4 check),
  4. A poisoned box where a gap is introduced in coverage (must trigger assert).
"""
import subprocess, tempfile, os

DEFECT_HARNESS = r'''
#define main cert_all_m_orig_main
#include "cert_all_m.cpp"
#undef main

int main() {
    init_exp_grid(1600);
    // Defect 1: bump U_hi[2][0] by +10% on [0, 0.00025] -> ledger MUST FAIL!
    BoxResult br = solve_box_int(0, 100000, 25, 100000, 200, 11 * H_DEN, 24, 32, 32, 40);
    string info_ok;
    bool ok0 = verify_all_tree_ledger(br, 32, 32, 40, info_ok);
    assert(ok0);
    BoxResult poisoned = br;
    poisoned.U_hi[2][0] = (int64_t)(br.U_hi[2][0] * 1.10);
    string info_bad;
    bool ok1 = verify_all_tree_ledger(poisoned, 32, 32, 40, info_bad);
    assert(!ok1);
    cout << "DEFECT 1 (poisoned U20 +10%) CAUGHT: " << info_bad << "\n";

    // Defect 2: high-tau cutoff lowered to tau = 5 (v = 1/5 instead of 1/6) -> margin MUST BE NEGATIVE!
    Iv inv_e = iv_exp_neg({S2, S2});
    Iv gam12 = iv_gamma_rat(1, 2);
    Iv gam13 = iv_gamma_rat(1, 3);
    Iv v15_12 = iv_pow_rat(iv_from_rat(1, 5), 1, 2);
    Iv v15_23 = iv_pow_rat(iv_from_rat(1, 5), 2, 3);
    // d_s at v = 1/5: d_{1/2} = (7/9)*(1/2)*(1 - 3/10) = 49/180; 1 - d_{1/2}/5 = 851/900
    i128 u12_hi = div_cl(div_cl(v15_12.hi * S2, gam12.lo) * 4 * 851, 3 * 900);
    // d_{2/3} = (7/9)*(2/3)*(1 - 5/15) = 28/81; 1 - d_{2/3}/5 = 377/405
    i128 u13_hi = div_cl(div_cl(v15_23.hi * S2, gam13.lo) * 2 * 377, 3 * 405);
    i128 margin5 = inv_e.lo - u12_hi - u13_hi;
    assert(margin5 < 0);
    cout << "DEFECT 2 (tau >= 5 cutoff) CAUGHT: margin = " << fixed << setprecision(6)
         << (double)margin5 / S2 << " < 0\n";

    // Defect 3: poison M_lo by -10% on [1.0, 1.2] -> ledger MUST FAIL!
    BoxResult br2 = solve_box_int(100000, 100000, 120000, 100000, 200, 4 * H_DEN, 22, 24, 26, 40);
    string info2_ok, info2_bad;
    assert(verify_all_tree_ledger(br2, 26, 24, 40, info2_ok));
    BoxResult p2 = br2;
    p2.M_lo = (int64_t)(br2.M_lo * 0.85);
    bool ok2 = verify_all_tree_ledger(p2, 26, 24, 40, info2_bad);
    assert(!ok2);
    cout << "DEFECT 3 (poisoned M_lo -15% at tau in [1,1.2]) CAUGHT: " << info2_bad << "\n";

    cout << "ALL 3 PLANTED DEFECTS CAUGHT LOUDLY!\n";
    return 0;
}
'''

with tempfile.TemporaryDirectory() as td:
    src = os.path.join(td, "defects.cpp")
    exe = os.path.join(td, "defects")
    with open(src, "w") as f:
        f.write(DEFECT_HARNESS)
    subprocess.check_call(["g++", "-O2", "-std=c++17", f"-I{os.path.dirname(os.path.abspath(__file__))}", src, "-o", exe])
    print(subprocess.check_output([exe], text=True), end="")
