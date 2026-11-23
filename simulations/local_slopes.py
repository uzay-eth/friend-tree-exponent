#!/usr/bin/env python3
"""Build finite-time slope estimates from the three archived CT batches."""
import argparse
import csv
import gzip
import math
from collections import defaultdict
from pathlib import Path

from analyze_ct import jack

ROOT = Path(__file__).resolve().parent
BATCHES = {
    "ct_t28_seedA.txt.gz": (527, 28),
    "ct_t28_seedB.txt.gz": (294, 28),
    "ct_t31_seedC.txt.gz": (91, 31),
}


def checkpoints(batches=None):
    if batches is None:
        batches = [(ROOT / "outputs" / name, runs, tmax)
                   for name, (runs, tmax) in BATCHES.items()]
    by_time = defaultdict(list)
    for filename, runs, tmax in batches:
        opener = gzip.open if filename.suffix == ".gz" else open
        with opener(filename, "rt") as source:
            if source.readline().split() != "# run t n C S F Cv A Z".split():
                raise ValueError(f"{filename}: unexpected column names")
            for run in range(runs):
                previous_core = 1  # deterministic initial T_3
                for t in range(1, tmax + 1):
                    fields = source.readline().split()
                    if len(fields) != 9 or tuple(map(int, fields[:2])) != (run, t):
                        raise ValueError(f"{filename}: expected run {run}, checkpoint {t}")
                    n, core = map(int, fields[2:4])
                    if core - previous_core != int(fields[6]):
                        raise ValueError(f"{filename}: core increment differs from conversions")
                    by_time[t].append((n, core, float(fields[4]), previous_core))
                    previous_core = core
            if source.read().strip():
                raise ValueError(f"{filename}: unexpected trailing records")
    return by_time


def paired_secant(core, previous, width):
    runs = len(core)
    c, p = sum(core), sum(previous)
    estimate = math.log(c / p) / width
    leave_one = [math.log((c - ci) / (p - pi)) / width for ci, pi in zip(core, previous)]
    center = sum(leave_one) / runs
    error = math.sqrt((runs - 1) / runs * sum((v - center)**2 for v in leave_one))
    return estimate, error


def estimates(by_time):
    rows = []
    for t, sample in sorted(by_time.items()):
        n, core, propensity, previous = zip(*sample)
        runs = len(core)
        slope, slope_se = jack(propensity, core)
        # Same trajectories on both sides, including the change in cohort at t=29.
        secant, secant_se = paired_secant(core, previous, 1)
        rows.append(dict(t=t, runs=runs, mean_n=sum(n) / runs, mean_core=sum(core) / runs,
                         slope=slope, slope_se=slope_se, secant=secant, secant_se=secant_se))
    return rows


def plot(rows, filename):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FormatStrFormatter, LogLocator

    fig = plt.figure(figsize=(11.2, 7.2), layout="constrained")
    grid = fig.add_gridspec(2, 2, height_ratios=[1.2, 1])
    full = fig.add_subplot(grid[0, :])
    small = fig.add_subplot(grid[1, 0])
    late = fig.add_subplot(grid[1, 1])
    fig.suptitle("Growth of the mean nonleaf core at fixed time", fontsize=17)

    for ax, selected in [(full, rows), (small, rows[3:13]), (late, rows[20:])]:
        ax.set_xscale("log")
        ax.grid(axis="y", color="#dee2e6", linewidth=0.7)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(axis="both", labelsize=9)
        ax.set_xlabel("Mean number of vertices at checkpoint")
        ax.set_ylabel("Local slope")
        for count, color in [(912, "#15597d"), (91, "#bf581f")]:
            group = [r for r in selected if r["runs"] == count]
            if group:
                ax.errorbar([r["mean_n"] for r in group], [r["slope"] for r in group],
                            yerr=[r["slope_se"] for r in group], fmt="o-", markersize=4,
                            linewidth=1.1, elinewidth=1, capsize=2, color=color,
                            label=f"{count} trajectories")

    full.set_xlim(5, 1.4e14)
    small.set_xlim(70, 1.8e6)
    late.set_xlim(2.5e9, 1.3e14)
    full.set_title("All 31 checkpoints, with one jackknife standard error", loc="left")
    full.xaxis.set_major_locator(LogLocator(base=10, numticks=16))
    full.legend(frameon=False)
    small.set_title("Sizes around the published fit", loc="left")
    small.hlines(0.569, 1e2, 1e6, ls="--", color="#6f4917", lw=1.4,
                 label="Hartle et al., fixed-size fit")
    small.legend(frameon=False, fontsize=9)
    small.yaxis.set_major_formatter(FormatStrFormatter("%.3f"))
    late.set_title("Late checkpoints, expanded scale", loc="left")
    late.yaxis.set_major_formatter(FormatStrFormatter("%.5f"))
    fig.savefig(filename, dpi=165, metadata={"Software": "local_slopes.py"})
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", type=Path, required=True)
    parser.add_argument("--figure", type=Path, help="PNG plot (requires matplotlib)")
    parser.add_argument("--sample", nargs=3, metavar=("FILE", "RUNS", "TMAX"),
                        help="use a fresh friend_ct sample with unit checkpoint spacing")
    args = parser.parse_args()
    if args.sample:
        filename, runs, tmax = args.sample
        runs, tmax = int(runs), int(tmax)
        if runs < 2 or tmax < 1:
            parser.error("a sample requires at least two runs and one checkpoint")
        if args.figure:
            parser.error("--figure renders the three-panel comparison of archived batches")
        by_time = checkpoints([(Path(filename), runs, tmax)])
    else:
        by_time = checkpoints()
    rows = estimates(by_time)
    with args.table.open("w", newline="") as dest:
        writer = csv.DictWriter(dest, fieldnames=rows[0].keys(), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: format(v, ".12g") if isinstance(v, float) else v for k, v in row.items()})
    if args.figure:
        plot(rows, args.figure)
    if args.sample:
        last = rows[-1]
        print(f"{len(rows)} checkpoints from {last['runs']} runs; "
              f"final instantaneous slope {last['slope']:.7f} +- {last['slope_se']:.7f}")
        return
    print(f"{len(rows)} checkpoints, t={rows[0]['t']}..{rows[-1]['t']}; "
          f"{rows[0]['runs']} runs through t=28, {rows[-1]['runs']} through t=31")

    # First and last checkpoints with sample mean n in [100, 1,000,000].
    inside = [r["t"] for r in rows if 100 <= r["mean_n"] <= 1e6]
    t1, t2 = inside[0], inside[-1]
    left = [row[1] for row in by_time[t1]]
    right = [row[1] for row in by_time[t2]]
    if len(left) != len(right):
        raise ValueError("different cohorts at comparison endpoints")
    change, error = paired_secant(right, left, t2 - t1)
    print(f"CT endpoint slope, t={t1}..{t2}: {change:.7f} +- {error:.7f} ({len(left)} runs)")


if __name__ == "__main__":
    main()
