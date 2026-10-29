#!/usr/bin/env python3
"""Factorium Benchmark Visualization Script (Phase 2).

Generates comparative analysis plots from benchmark CSV data using pandas and matplotlib:
  1. Runtime versus input bit length (log scale).
  2. Fermat iterations/runtime versus factor distance |q - p|.
  3. Pollard Rho iteration distributions across seeds.
  4. Pollard p-1 success and runtime versus smoothness bound B.
  5. Success and timeout rates by algorithm and dataset category.

Outputs plots as high-resolution PNG images in benchmarks/plots/.
"""

import os
import subprocess
import sys
import numpy as np

# Use Agg backend for headless environments
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

PLOTS_DIR = os.path.join(os.path.dirname(__file__), "..", "benchmarks", "plots")
RAW_CSV_PATH = "bench_raw.csv"
SUMMARY_CSV_PATH = "bench_summary.csv"


def ensure_data():
    """Ensure raw and summary benchmark CSV files exist; generate if missing."""
    if not os.path.exists(RAW_CSV_PATH) or not os.path.exists(SUMMARY_CSV_PATH):
        bin_path = os.path.join(".", "build", "factorium_bench")
        if not os.path.exists(bin_path):
            bin_path = "factorium_bench"
        print(f"Benchmark CSV not found. Running {bin_path}...")
        cmd = [bin_path, "--full", "--raw-out", RAW_CSV_PATH, "--summary-out", SUMMARY_CSV_PATH]
        subprocess.run(cmd, check=True)


def plot_runtime_vs_bit_length(raw_df):
    """Plot 1: Runtime vs input bit length for balanced scaling semiprimes."""
    fig, ax = plt.subplots(figsize=(8, 5))
    df_scale = raw_df[raw_df["category"] == "bit_length_scaling"].copy()
    df_scale["runtime_us"] = df_scale["runtime_ns"] / 1000.0

    colors = {"fermat": "#1f77b4", "rho": "#ff7f0e", "pminus1": "#2ca02c"}
    markers = {"fermat": "o", "rho": "s", "pminus1": "^"}

    for algo in ["fermat", "rho", "pminus1"]:
        sub = df_scale[df_scale["algorithm"] == algo]
        if sub.empty:
            continue
        # Group by bit length to get median runtime
        grouped = sub.groupby("bit_length")["runtime_us"].median().reset_index()
        ax.plot(
            grouped["bit_length"],
            grouped["runtime_us"],
            label=f"{algo.capitalize()}",
            color=colors.get(algo, "gray"),
            marker=markers.get(algo, "o"),
            linewidth=2,
            markersize=6,
        )

    ax.set_yscale("log")
    ax.set_xlabel("Input Bit Length (bits)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Median Runtime (μs, log scale)", fontsize=11, fontweight="bold")
    ax.set_title("Runtime Scaling vs. Input Bit Length", fontsize=12, fontweight="bold")
    ax.grid(True, which="both", ls="--", alpha=0.5)
    ax.legend(title="Algorithm", fontsize=10)
    plt.tight_layout()
    out_file = os.path.join(PLOTS_DIR, "runtime_vs_bit_length.png")
    plt.savefig(out_file, dpi=200)
    plt.close()
    print(f"Generated: {out_file}")


def plot_fermat_distance(raw_df):
    """Plot 2: Fermat iterations and runtime vs factor distance |q - p|."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    df_close = raw_df[(raw_df["category"] == "close_semiprimes") & (raw_df["algorithm"] == "fermat")].copy()
    df_close["delta"] = (df_close["ground_truth_q"] - df_close["ground_truth_p"]).abs()
    df_close["runtime_us"] = df_close["runtime_ns"] / 1000.0
    df_close = df_close.sort_values(by="delta")

    # Left: Iterations vs Factor Distance
    ax1.plot(df_close["delta"], df_close["iterations"], marker="o", color="#1f77b4", linewidth=2)
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel("Factor Distance |q - p| (log scale)", fontsize=11, fontweight="bold")
    ax1.set_ylabel("Fermat Iterations (log scale)", fontsize=11, fontweight="bold")
    ax1.set_title("Fermat Iterations vs. Factor Distance", fontsize=12, fontweight="bold")
    ax1.grid(True, which="both", ls="--", alpha=0.5)

    # Theoretical reference line: iterations ~ delta / 2
    x_theory = np.linspace(df_close["delta"].min(), df_close["delta"].max(), 50)
    ax1.plot(x_theory, np.maximum(1, x_theory / 2), ls=":", color="red", label="Theory (Δ/2)")
    ax1.legend()

    # Right: Runtime vs Factor Distance
    ax2.plot(df_close["delta"], df_close["runtime_us"], marker="s", color="#d62728", linewidth=2)
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel("Factor Distance |q - p| (log scale)", fontsize=11, fontweight="bold")
    ax2.set_ylabel("Runtime (μs, log scale)", fontsize=11, fontweight="bold")
    ax2.set_title("Fermat Runtime vs. Factor Distance", fontsize=12, fontweight="bold")
    ax2.grid(True, which="both", ls="--", alpha=0.5)

    plt.tight_layout()
    out_file = os.path.join(PLOTS_DIR, "fermat_vs_factor_distance.png")
    plt.savefig(out_file, dpi=200)
    plt.close()
    print(f"Generated: {out_file}")


def plot_rho_distributions(raw_df):
    """Plot 3: Pollard Rho iteration distributions across randomized seeds."""
    fig, ax = plt.subplots(figsize=(9, 5))

    df_rho = raw_df[(raw_df["algorithm"] == "rho") & (raw_df["category"].isin(["close_semiprimes", "distant_semiprimes"]))].copy()
    # Filter to cases with >1 trial
    trial_counts = df_rho.groupby("dataset_id")["trial_id"].count()
    multi_trial_ids = trial_counts[trial_counts >= 3].index
    df_rho_multi = df_rho[df_rho["dataset_id"].isin(multi_trial_ids)]

    # Prepare data for boxplot
    groups = []
    labels = []
    for d_id, grp in df_rho_multi.groupby("dataset_id"):
        groups.append(grp["iterations"].values)
        labels.append(d_id)

    bp = ax.boxplot(groups, tick_labels=labels, patch_artist=True)
    for box in bp['boxes']:
        box.set(facecolor='#ffbb78', color='#ff7f0e', linewidth=1.5)
    for median in bp['medians']:
        median.set(color='darkred', linewidth=2)

    ax.set_yscale("log")
    ax.set_xlabel("Dataset Case", fontsize=11, fontweight="bold")
    ax.set_ylabel("Iterations (log scale)", fontsize=11, fontweight="bold")
    ax.set_title("Pollard's Rho Iteration Variability Across Seeds", fontsize=12, fontweight="bold")
    plt.xticks(rotation=45, ha="right")
    ax.grid(True, which="both", ls="--", alpha=0.5)
    plt.tight_layout()

    out_file = os.path.join(PLOTS_DIR, "rho_iteration_distribution.png")
    plt.savefig(out_file, dpi=200)
    plt.close()
    print(f"Generated: {out_file}")


def plot_pminus1_smoothness():
    """Plot 4: Pollard p-1 success and runtime versus smoothness bound B."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    # Test cases: favorable (smooth p-1) vs safe prime (unfavorable)
    # p1_favorable_7681 (smooth: max factor 5), p1_favorable_65537 (smooth: 2)
    # p1_unfavorable_safe1 (safe prime: factor 500333)
    bounds = [10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000, 600000]

    bin_path = os.path.join(".", "build", "factorium_bench")
    if not os.path.exists(bin_path):
        bin_path = "factorium_bench"

    results = []
    tmp_raw = "tmp_p1_sweep.csv"
    tmp_sum = "tmp_p1_sum.csv"

    for b in bounds:
        cmd = [bin_path, "--smoke", "--bound", str(b), "--raw-out", tmp_raw, "--summary-out", tmp_sum, "--quiet"]
        try:
            subprocess.run(cmd, check=True)
            df = pd.read_csv(tmp_raw)
            p1_df = df[df["algorithm"] == "pminus1"]
            for _, row in p1_df.iterrows():
                results.append({
                    "bound": b,
                    "dataset_id": row["dataset_id"],
                    "status": row["status"],
                    "runtime_us": row["runtime_ns"] / 1000.0,
                    "modmul_ops": row["modmul_ops"]
                })
        except Exception as e:
            print(f"Warning during p-1 sweep at B={b}: {e}")

    # Clean temporary files
    for p in [tmp_raw, tmp_sum]:
        if os.path.exists(p):
            os.remove(p)

    if not results:
        print("Skipping p-1 plot: no sweep data collected")
        return

    res_df = pd.DataFrame(results)

    # Plot representative cases
    fav_case = "p1_favorable_65537"
    unfav_case = "p1_unfavorable_safe1"

    fav_data = res_df[res_df["dataset_id"] == fav_case]
    unfav_data = res_df[res_df["dataset_id"] == unfav_case]

    # Subplot 1: Runtime vs Bound B
    if not fav_data.empty:
        ax1.plot(fav_data["bound"], fav_data["runtime_us"], marker="o", color="#2ca02c", label="Favorable (65537, 2-smooth)", lw=2)
    if not unfav_data.empty:
        ax1.plot(unfav_data["bound"], unfav_data["runtime_us"], marker="x", color="#d62728", ls="--", label="Safe Prime (factor > 500k)", lw=2)

    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel("Smoothness Bound B (log scale)", fontsize=11, fontweight="bold")
    ax1.set_ylabel("Runtime (μs, log scale)", fontsize=11, fontweight="bold")
    ax1.set_title("Pollard p-1 Runtime vs. Bound B", fontsize=12, fontweight="bold")
    ax1.grid(True, which="both", ls="--", alpha=0.5)
    ax1.legend()

    # Subplot 2: ModMul Operations vs Bound B
    if not fav_data.empty:
        ax2.plot(fav_data["bound"], fav_data["modmul_ops"], marker="s", color="#2ca02c", label="Favorable", lw=2)
    if not unfav_data.empty:
        ax2.plot(unfav_data["bound"], unfav_data["modmul_ops"], marker="^", color="#d62728", ls="--", label="Safe Prime", lw=2)

    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel("Smoothness Bound B (log scale)", fontsize=11, fontweight="bold")
    ax2.set_ylabel("Modular Multiplications (log scale)", fontsize=11, fontweight="bold")
    ax2.set_title("Pollard p-1 Work vs. Bound B", fontsize=12, fontweight="bold")
    ax2.grid(True, which="both", ls="--", alpha=0.5)
    ax2.legend()

    plt.tight_layout()
    out_file = os.path.join(PLOTS_DIR, "pminus1_vs_smoothness_bound.png")
    plt.savefig(out_file, dpi=200)
    plt.close()
    print(f"Generated: {out_file}")


def plot_success_timeout_rates(summary_df):
    """Plot 5: Success and timeout rates by algorithm and dataset category."""
    fig, ax = plt.subplots(figsize=(10, 5))

    # Aggregate by category and algorithm
    cat_summary = summary_df.groupby(["category", "algorithm"])[["num_success", "num_trials"]].sum().reset_index()
    cat_summary["success_rate"] = (cat_summary["num_success"] / cat_summary["num_trials"]) * 100.0

    categories = cat_summary["category"].unique()
    x = np.arange(len(categories))
    width = 0.25

    colors = {"fermat": "#1f77b4", "rho": "#ff7f0e", "pminus1": "#2ca02c"}

    for idx, algo in enumerate(["fermat", "rho", "pminus1"]):
        algo_data = cat_summary[cat_summary["algorithm"] == algo]
        # Align with categories
        rates = []
        for cat in categories:
            match = algo_data[algo_data["category"] == cat]
            rates.append(match["success_rate"].values[0] if not match.empty else 0.0)

        offset = (idx - 1) * width
        ax.bar(x + offset, rates, width, label=algo.capitalize(), color=colors.get(algo, "gray"), alpha=0.85)

    ax.set_ylabel("Success Rate (%)", fontsize=11, fontweight="bold")
    ax.set_title("Success and Resolution Rates by Category and Algorithm", fontsize=12, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels([c.replace("_", "\n") for c in categories], fontsize=9)
    ax.set_ylim(0, 110)
    ax.grid(axis="y", ls="--", alpha=0.5)
    ax.legend(title="Algorithm", loc="lower right")

    plt.tight_layout()
    out_file = os.path.join(PLOTS_DIR, "success_timeout_rates.png")
    plt.savefig(out_file, dpi=200)
    plt.close()
    print(f"Generated: {out_file}")


def main():
    os.makedirs(PLOTS_DIR, exist_ok=True)
    ensure_data()

    raw_df = pd.read_csv(RAW_CSV_PATH)
    summary_df = pd.read_csv(SUMMARY_CSV_PATH)

    print(f"Loaded {len(raw_df)} raw benchmark records and {len(summary_df)} summary records.")

    plot_runtime_vs_bit_length(raw_df)
    plot_fermat_distance(raw_df)
    plot_rho_distributions(raw_df)
    plot_pminus1_smoothness()
    plot_success_timeout_rates(summary_df)

    print(f"\nAll benchmark plots successfully written to: {os.path.abspath(PLOTS_DIR)}")


if __name__ == "__main__":
    main()
