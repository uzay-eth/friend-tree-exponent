#!/usr/bin/env python3
"""Analyze simulation replicates for IR, DAN, and PAN models across grids, modes, and fit variants."""
import csv
import math
import os
import sys
import numpy as np
from scipy.optimize import curve_fit

GRID_17 = np.array([
    100, 178, 316, 562,
    1000, 1778, 3162, 5623,
    10000, 17783, 31623, 56234,
    100000, 177828, 316228, 562341,
    1000000
], dtype=np.float64)

GRIDS = {
    "5pt_decade": np.array([0, 4, 8, 12, 16]),
    "9pt_halfdec": np.array([0, 2, 4, 6, 8, 10, 12, 14, 16]),
    "17pt_qtrdec": np.arange(17),
}

REF_MU = {
    "IR": 0.565723,
    "DAN": 0.771192302407,
    "PAN": 0.549734972504,
}

PAPER_TARGET = {
    "IR": (">=", 0.569),
    "DAN": (">=", 0.773),
    "PAN": ("<=", 0.545),
}


def fit_all_variants(Ns, Cs):
    """Return dict of fitted exponent b for (i) linear NLS, (ii) log-log OLS, (iii) linear NLS + c."""
    logN = np.log(Ns)
    logC = np.log(Cs)
    # (ii) unweighted OLS on log<C> vs log N
    b_log, log_a = np.polyfit(logN, logC, 1)
    a_init = math.exp(log_a)

    # (i) unweighted nonlinear least squares in linear space: a * N^b
    def power_law(x, a, b):
        return a * np.power(x, b)

    popt_nls, _ = curve_fit(
        power_law, Ns, Cs,
        p0=[a_init, b_log],
        maxfev=50000,
        ftol=1e-12, xtol=1e-12, gtol=1e-12,
    )
    b_nls = float(popt_nls[1])

    # (iii) unweighted nonlinear least squares with additive constant: a * N^b + c
    def power_law_c(x, a, b, c):
        return a * np.power(x, b) + c

    try:
        popt_nlsc, _ = curve_fit(
            power_law_c, Ns, Cs,
            p0=[popt_nls[0], popt_nls[1], 0.0],
            maxfev=100000,
            ftol=1e-12, xtol=1e-12, gtol=1e-12,
        )
        b_nlsc = float(popt_nlsc[1])
    except Exception:
        b_nlsc = float("nan")

    return {
        "linear_NLS": b_nls,
        "loglog_OLS": float(b_log),
        "linear_NLS_plus_c": b_nlsc,
    }


def load_reps(path):
    with open(path, "r") as f:
        rows = list(csv.DictReader(f))
    n_reps = len(rows)
    trees = int(rows[0]["trees"])
    means = np.zeros((n_reps, len(GRID_17)), dtype=np.float64)
    for r, row in enumerate(rows):
        for g, N in enumerate(GRID_17):
            means[r, g] = float(row[f"mean_{int(N)}"])
    return means, trees


def load_pooled_cov(path):
    with open(path, "r") as f:
        rows = list(csv.DictReader(f))
    n_g = len(rows)
    total_trees = int(rows[0]["total_trees"])
    mean = np.array([float(r["mean"]) for r in rows])
    sd = np.array([float(r["sd"]) for r in rows])
    cv = np.array([float(r["cv"]) for r in rows])
    cov = np.zeros((n_g, n_g), dtype=np.float64)
    for i, r in enumerate(rows):
        for j, N in enumerate(GRID_17):
            cov[i, j] = float(r[f"cov_{int(N)}"])
    return total_trees, mean, sd, cv, cov


def compute_local_slopes(Ns_idx, total_trees, mean, cov):
    """Compute local slopes between adjacent indices in Ns_idx with delta-method SE."""
    out = []
    for k in range(len(Ns_idx) - 1):
        i, j = Ns_idx[k], Ns_idx[k + 1]
        Ni, Nj = GRID_17[i], GRID_17[j]
        dlogN = math.log(Nj) - math.log(Ni)
        mi, mj = mean[i], mean[j]
        slope = (math.log(mj) - math.log(mi)) / dlogN
        # Delta method: Var(log(mj) - log(mi)) = (1/M) * (Var_j/mj^2 + Var_i/mi^2 - 2*Cov_ij/(mi*mj))
        var_diff = (cov[j, j] / (mj * mj) + cov[i, i] / (mi * mi) - 2.0 * cov[i, j] / (mi * mj)) / total_trees
        se = math.sqrt(max(0.0, var_diff)) / dlogN
        out.append((int(Ni), int(Nj), slope, se))
    return out


def analyze_dir(data_dir):
    models = ["IR", "DAN", "PAN"]
    modes = ["coupled", "indep"]

    print("=" * 96)
    print("FIT REPLICATE SUMMARY (2400 realizations per replicate)")
    print("=" * 96)

    for model in models:
        ref = REF_MU[model]
        op, thresh = PAPER_TARGET[model]
        thresh_str = f"b {op} {thresh:.3f}"

        for mode in modes:
            rep_path = os.path.join(data_dir, f"{model.lower()}_{mode}_reps.csv")
            if not os.path.exists(rep_path):
                print(f"[MISSING] {rep_path}")
                continue
            means, trees = load_reps(rep_path)
            R = means.shape[0]
            print(f"\n### Model: {model} | Mode: {mode} | Replicates: R={R} x T={trees} | Ref mu={ref:.6f} | Target: {thresh_str}")
            header = f"{'Grid':12s} {'Fit Method':18s} {'Mean(b)':>10s} {'SD(b)':>9s} {'Bias':>10s} {'SE(mean)':>9s} {'Frac(' + thresh_str + ')':>18s} {'z_target':>9s} {'Min(b)':>9s} {'Max(b)':>9s}"
            print(header)
            print("-" * len(header))

            for gname, gidx in GRIDS.items():
                Ns = GRID_17[gidx]
                bs = {"linear_NLS": [], "loglog_OLS": [], "linear_NLS_plus_c": []}
                for r in range(R):
                    res = fit_all_variants(Ns, means[r, gidx])
                    for k, v in res.items():
                        bs[k].append(v)

                for fname in ["linear_NLS", "loglog_OLS", "linear_NLS_plus_c"]:
                    arr = np.array(bs[fname], dtype=np.float64)
                    valid = arr[np.isfinite(arr)]
                    n_ok = len(valid)
                    m_b = float(np.mean(valid))
                    s_b = float(np.std(valid, ddof=1))
                    se_m = s_b / math.sqrt(n_ok)
                    bias = m_b - ref
                    if op == ">=":
                        hits = int(np.sum(valid >= thresh))
                    else:
                        hits = int(np.sum(valid <= thresh))
                    frac = hits / n_ok
                    z_tgt = (thresh - m_b) / s_b if s_b > 0 else float("nan")
                    print(
                        f"{gname:12s} {fname:18s} {m_b:10.6f} {s_b:9.6f} {bias:+10.6f} {se_m:9.6f} "
                        f"{hits:3d}/{n_ok:<3d} ({100*frac:5.1f}%) {z_tgt:+9.2f} {np.min(valid):9.6f} {np.max(valid):9.6f}"
                    )

                # Also fit the POOLED mean curve across all R replicates to measure pure finite-size bias cleanly
                pooled_mean = np.mean(means[:, gidx], axis=0)
                pooled_res = fit_all_variants(Ns, pooled_mean)
                print(
                    f"  [Pooled {R*trees:,d} trees {gname}] "
                    f"linear_NLS={pooled_res['linear_NLS']:.6f} (bias {pooled_res['linear_NLS']-ref:+.6f}), "
                    f"loglog_OLS={pooled_res['loglog_OLS']:.6f} (bias {pooled_res['loglog_OLS']-ref:+.6f}), "
                    f"NLS+c={pooled_res['linear_NLS_plus_c']:.6f} (bias {pooled_res['linear_NLS_plus_c']-ref:+.6f})"
                )

    print("\n" + "=" * 96)
    print("POOLED STATISTICS: CV(C) AT N=10^6 AND LOCAL SLOPES (ALL COUPLED TREES POOLED)")
    print("=" * 96)

    for model in models:
        ref = REF_MU[model]
        cov_path = os.path.join(data_dir, f"{model.lower()}_coupled_pooled_cov.csv")
        if not os.path.exists(cov_path):
            continue
        total_trees, mean, sd, cv, cov = load_pooled_cov(cov_path)

        # Check also independent run at N=10^6 if available and pool all independent trees at N=10^6
        indep_path = os.path.join(data_dir, f"{model.lower()}_indep_pooled_cov.csv")
        print(f"\n### Model: {model} (coupled pooled M={total_trees:,d} trees, ref mu={ref:.6f})")
        print(f"  At N=10^6 (coupled): mean(C) = {mean[-1]:.4f} +- {sd[-1]/math.sqrt(total_trees):.4f}, "
              f"sd(C) = {sd[-1]:.4f}, CV = {cv[-1]:.6f}")
        if os.path.exists(indep_path):
            M_i, m_i, sd_i, cv_i, _ = load_pooled_cov(indep_path)
            print(f"  At N=10^6 (indep):   mean(C) = {m_i[-1]:.4f} +- {sd_i[-1]/math.sqrt(M_i):.4f}, "
                  f"sd(C) = {sd_i[-1]:.4f}, CV = {cv_i[-1]:.6f}")
            # Combined grand pool at N=10^6 (coupled + indep are disjoint seeds)
            M_all = total_trees + M_i
            m_all = (total_trees * mean[-1] + M_i * m_i[-1]) / M_all
            # Exact pooled variance
            ss_c = (total_trees - 1) * sd[-1]**2 + total_trees * (mean[-1] - m_all)**2
            ss_i = (M_i - 1) * sd_i[-1]**2 + M_i * (m_i[-1] - m_all)**2
            sd_all = math.sqrt((ss_c + ss_i) / (M_all - 1))
            print(f"  At N=10^6 (ALL M={M_all:,d}): mean(C) = {m_all:.4f} +- {sd_all/math.sqrt(M_all):.4f}, "
                  f"sd(C) = {sd_all:.4f}, CV = {sd_all/m_all:.6f}")

        print("\n  Pooled <C(N)> and CV(N) across 17-point grid (coupled):")
        for i in range(len(GRID_17)):
            se_i = sd[i] / math.sqrt(total_trees)
            print(f"    N={int(GRID_17[i]):7d}: <C>={mean[i]:12.4f} +- {se_i:8.4f}  (CV={cv[i]:.5f})")

        for gname, gidx in [("5pt_decade", GRIDS["5pt_decade"]),
                            ("9pt_halfdec", GRIDS["9pt_halfdec"]),
                            ("17pt_qtrdec", GRIDS["17pt_qtrdec"])]:
            print(f"\n  Local slopes ({gname}, coupled delta-method SE):")
            slopes = compute_local_slopes(gidx, total_trees, mean, cov)
            for Ni, Nj, sl, se in slopes:
                print(f"    [{Ni:7d} -> {Nj:7d}]: slope = {sl:.6f} +- {se:.6f}  (diff from ref: {sl - ref:+.6f})")


if __name__ == "__main__":
    data_dir = sys.argv[1] if len(sys.argv) > 1 else "fit569"
    analyze_dir(data_dir)
