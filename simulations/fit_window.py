#!/usr/bin/env python3
"""Finite-window fits of mean core counts at fixed sizes and fixed times."""
import argparse
import csv
import gzip
import math
from pathlib import Path

from local_slopes import checkpoints

HERE = Path(__file__).resolve().parent


def jack_se(values):
    center = sum(values) / len(values)
    return math.sqrt((len(values) - 1) / len(values) * sum((v - center)**2 for v in values))


def power_fit(log_n, means):
    # For a fixed exponent b, the least-squares amplitude is (x.y)/(x.x), x=N**b.
    def residual(b):
        x = [math.exp(log_size * b) for log_size in log_n]
        a = sum(v*y for v, y in zip(x, means)) / sum(v*v for v in x)
        return sum((y - a*v)**2 for y, v in zip(means, x))

    left, right = 0., 1.
    ratio = (math.sqrt(5.) - 1.) / 2.
    lo, hi = right - ratio * (right-left), left + ratio * (right-left)
    low_loss, high_loss = residual(lo), residual(hi)
    while right-left > 1e-12:
        if low_loss <= high_loss:
            right, hi, high_loss = hi, lo, low_loss
            lo = right - ratio * (right-left)
            low_loss = residual(lo)
        else:
            left, lo, low_loss = lo, hi, high_loss
            hi = left + ratio * (right-left)
            high_loss = residual(hi)
    return (lo+hi) / 2


def log_fit(log_n, means):
    center = sum(log_n) / len(log_n)
    return (sum((x - center)*math.log(y) for x, y in zip(log_n, means))
            / sum((x-center)**2 for x in log_n))


def read(path):
    with gzip.open(path, "rt") as src:
        reader = csv.reader(src, delimiter="\t")
        name, *sizes = next(reader)
        sizes = list(map(int, sizes))
        seeds = set()
        rows = []
        for seed, *counts in reader:
            values = list(map(int, counts))
            if seed in seeds or len(values) != len(sizes):
                raise ValueError("duplicate seed or incomplete trajectory")
            seeds.add(seed)
            rows.append(values)
    if name != "seed" or len(rows) < 2:
        raise ValueError("expected at least two seeded trajectories")
    return sizes, rows


def fit_summary(sizes, rows):
    n = len(rows)
    totals = [sum(col) for col in zip(*rows)]
    means = [s/n for s in totals]
    log_n = [math.log(v) for v in sizes]
    leave_one = [[(s-v)/(n-1) for s, v in zip(totals, row)] for row in rows]
    for name, fit in [("Unweighted nonlinear least squares on mean C", power_fit),
                      ("Unweighted least squares on log(mean C)", log_fit)]:
        print(f"{name}: {fit(log_n, means):.7f} +- {jack_se([fit(log_n,m) for m in leave_one]):.7f}")
    return means, log_n, leave_one


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=HERE/"fixed-size"/"core_counts.tsv.gz")
    parser.add_argument("--table", type=Path, required=True)
    args = parser.parse_args()
    sizes, rows = read(args.data)
    n = len(rows)
    print(f"Fixed N: {n} trajectories, {len(sizes)} fixed sizes from {sizes[0]} to {sizes[-1]}")
    means, log_n, leave_one = fit_summary(sizes, rows)
    with args.table.open("w", newline="") as dest:
        writer = csv.writer(dest, delimiter="\t", lineterminator="\n")
        writer.writerow(["N", "runs", "mean_core", "mean_core_se", "adjacent_slope", "adjacent_slope_se"])
        for i, size in enumerate(sizes):
            mean_se = jack_se([m[i] for m in leave_one])
            fields = [size, n, f"{means[i]:.9f}", f"{mean_se:.9f}"]
            if i:
                width = log_n[i] - log_n[i-1]
                slope = math.log(means[i]/means[i-1])/width
                se = jack_se([math.log(m[i]/m[i-1])/width for m in leave_one])
                fields.extend([f"{slope:.9f}", f"{se:.9f}"])
            else:
                fields.extend(["NA", "NA"])
            writer.writerow(fields)

    # Match the published size range using the model's exact mean population, 3*exp(t).
    by_time = checkpoints()
    times = sorted(t for t in by_time if 100 <= 3*math.exp(t) <= 1e6)
    counts = [[row[1] for row in by_time[t]] for t in times]
    if len({len(col) for col in counts}) != 1:
        raise ValueError("CT fitting window must have the same trajectories throughout")
    ct_rows = list(zip(*counts))
    print(f"Fixed time: {len(ct_rows)} trajectories, t={times[0]}..{times[-1]}, x=3 exp(t)")
    fit_summary([3*math.exp(t) for t in times], ct_rows)


if __name__ == "__main__":
    main()
