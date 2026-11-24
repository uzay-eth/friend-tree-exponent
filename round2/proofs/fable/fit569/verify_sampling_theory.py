#!/usr/bin/env python3
"""Verify empirical replicate SD(b) and hit probabilities in independent-grid mode
against exact linearized (delta-method) error propagation and 100,000 parametric Gaussian draws
around the high-precision 192,000-tree coupled pooled mean <C(N)> and CV(N).
"""
import math
import numpy as np
from scipy.optimize import curve_fit
from analyze import GRID_17, GRIDS, REF_MU, PAPER_TARGET, load_pooled_cov, fit_all_variants

DATA_DIR = "fit569"
T = 2400
N_MC = 20000


def exact_linearized_sd(Ns, Cs, rel_se, fit_type):
    """Exact first-order delta-method SD(b) for independent observations with SD_i = rel_se_i * Cs_i."""
    sigma = rel_se * Cs
    if fit_type == "loglog_OLS":
        x = np.log(Ns)
        xc = x - np.mean(x)
        w = xc / np.sum(xc * xc)
        # d(b)/d(log C_i) = w_i; Var(log C_i) = rel_se_i^2
        return float(math.sqrt(np.sum((w * rel_se) ** 2)))
    elif fit_type == "linear_NLS":
        # Fit a*N^b to exact Cs, then J = [N^b, a*N^b*log(N)], Cov(theta) = (J^T J)^-1 (J^T diag(sigma^2) J) (J^T J)^-1
        logN, logC = np.log(Ns), np.log(Cs)
        b0, la0 = np.polyfit(logN, logC, 1)
        popt, _ = curve_fit(lambda x, a, b: a * np.power(x, b), Ns, Cs, p0=[math.exp(la0), b0], maxfev=50000)
        a, b = popt
        f = a * np.power(Ns, b)
        J = np.column_stack([np.power(Ns, b), f * np.log(Ns)])
        JtJ_inv = np.linalg.inv(J.T @ J)
        M = J.T @ (np.diag(sigma ** 2)) @ J
        cov = JtJ_inv @ M @ JtJ_inv
        return float(math.sqrt(cov[1, 1])), float(b)
    elif fit_type == "linear_NLS_plus_c":
        logN, logC = np.log(Ns), np.log(Cs)
        b0, la0 = np.polyfit(logN, logC, 1)
        popt0, _ = curve_fit(lambda x, a, b: a * np.power(x, b), Ns, Cs, p0=[math.exp(la0), b0], maxfev=50000)
        popt, _ = curve_fit(lambda x, a, b, c: a * np.power(x, b) + c, Ns, Cs, p0=[popt0[0], popt0[1], 0.0], maxfev=100000)
        a, b, c = popt
        f = a * np.power(Ns, b)
        J = np.column_stack([np.power(Ns, b), f * np.log(Ns), np.ones_like(Ns)])
        JtJ_inv = np.linalg.inv(J.T @ J)
        M = J.T @ (np.diag(sigma ** 2)) @ J
        cov = JtJ_inv @ M @ JtJ_inv
        return float(math.sqrt(cov[1, 1])), float(b)


def main():
    from scipy.stats import norm
    print("Exact linearized (delta-method) predictions for INDEPENDENT 2400-tree replicates")
    print("(using 192,000-tree coupled pooled <C(N)> and CV(N) as ground truth):")
    print()
    for model in ["IR", "DAN", "PAN"]:
        _, mean, sd, cv, _ = load_pooled_cov(f"{DATA_DIR}/{model.lower()}_coupled_pooled_cov.csv")
        rel_se_all = cv / math.sqrt(T)
        ref = REF_MU[model]
        op, thresh = PAPER_TARGET[model]
        print(f"=== {model} (ref={ref:.6f}, target b {op} {thresh:.3f}, rel_SE(10^6)={100*rel_se_all[-1]:.3f}%) ===")
        for gname, gidx in GRIDS.items():
            Ns = GRID_17[gidx]
            Cs = mean[gidx]
            rse = rel_se_all[gidx]
            sd_nls, b_nls = exact_linearized_sd(Ns, Cs, rse, "linear_NLS")
            sd_nlsc, b_nlsc = exact_linearized_sd(Ns, Cs, rse, "linear_NLS_plus_c")
            b_log, _ = np.polyfit(np.log(Ns), np.log(Cs), 1)
            sd_log = exact_linearized_sd(Ns, Cs, rse, "loglog_OLS")
            for fname, b0, s0 in [("linear_NLS", b_nls, sd_nls),
                                  ("loglog_OLS", b_log, sd_log),
                                  ("linear_NLS_plus_c", b_nlsc, sd_nlsc)]:
                z = (thresh - b0) / s0
                p_hit = (1.0 - norm.cdf(z)) if op == ">=" else norm.cdf(z)
                print(f"  {gname:12s} {fname:18s}: b_true_fit={b0:.6f} (bias {b0-ref:+.6f}), "
                      f"SD_theory={s0:.6f}, z={z:+5.2f}, P_hit={100*p_hit:5.1f}%")
        print()


if __name__ == "__main__":
    main()
