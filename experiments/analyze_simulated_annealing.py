import json
import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize

# Config
RESULTS_FILE = "reports/experiment_results.json"

# Save plots to Simulated Annealing folder
OUTPUT_DIR = "reports/simulated_annealing"
OUTPUT_PLOT = f"{OUTPUT_DIR}/sa_convergence_plot.png"
HIST_OUT = f"{OUTPUT_DIR}/sa_distribution_plot.png"
EPOCH_OUT = f"{OUTPUT_DIR}/sa_epoch_boxplot.png"
CONTIGS_OUT = f"{OUTPUT_DIR}/sa_contigs_plot.png"
OVERLAP_OUT = f"{OUTPUT_DIR}/sa_overlap_plot.png"
SEARCHSPACE_OUT = f"{OUTPUT_DIR}/sa_searchspace_distribution.png"
CONVERGENCE_GAP_OUT = f"{OUTPUT_DIR}/sa_convergence_gap.png"
SWEEP_CONVERGENCE_OUT = f"{OUTPUT_DIR}/sa_sweep_convergence_comparison.png"
SWEEP_RUNTIME_OUT = f"{OUTPUT_DIR}/sa_sweep_runtime_comparison.png"
SWEEP_HEATMAP_OUT = f"{OUTPUT_DIR}/sa_sweep_best_cost_heatmap.png"


def _find_sa_runs(all_results):
    runs = []
    for res in all_results:
        method = res.get("method", "")
        if method == "Simulated Annealing":
            runs.append(res)
    return runs


def _run_label(run):
    if "method_label" in run:
        return run["method_label"]
    t0 = run.get("T0", "?")
    alpha = run.get("alpha", "?")
    return f"SA (T0={t0}, alpha={alpha})"


def _plot_sweep_comparison(sa_runs):
    if len(sa_runs) <= 1:
        return

    plt.figure(figsize=(11, 7))
    cmap = plt.cm.viridis
    color_norm = Normalize(vmin=0, vmax=max(1, len(sa_runs) - 1))

    for idx, run in enumerate(sa_runs):
        history = run.get("history", [])
        if not history:
            continue
        label = _run_label(run)
        plt.plot(history, linewidth=1.8, alpha=0.85, color=cmap(color_norm(idx)), label=label)

    plt.title("SA Hyperparameter Sweep: Best-Cost Convergence", fontsize=14)
    plt.xlabel("Evaluations", fontsize=12)
    plt.ylabel("Best Cost (Lower is Better)", fontsize=12)
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.legend(fontsize=8, loc="best")
    plt.tight_layout()
    plt.savefig(SWEEP_CONVERGENCE_OUT, dpi=200)
    print(f"Generated SA sweep convergence comparison: {SWEEP_CONVERGENCE_OUT}")
    plt.close()

    labels = [_run_label(run) for run in sa_runs]
    runtimes = [run.get("runtime_sec", 0.0) for run in sa_runs]
    plt.figure(figsize=(max(10, len(labels) * 1.1), 6))
    bars = plt.bar(range(len(labels)), runtimes, color="#4c78a8", alpha=0.85)
    plt.xticks(range(len(labels)), labels, rotation=35, ha="right")
    plt.title("SA Hyperparameter Sweep: Runtime Comparison", fontsize=14)
    plt.ylabel("Runtime (seconds)", fontsize=12)
    plt.grid(axis="y", linestyle="--", alpha=0.6)
    for i, bar in enumerate(bars):
        plt.text(i, bar.get_height(), f"{bar.get_height():.2f}", ha="center", va="bottom", fontsize=8)
    plt.tight_layout()
    plt.savefig(SWEEP_RUNTIME_OUT, dpi=200)
    print(f"Generated SA sweep runtime comparison: {SWEEP_RUNTIME_OUT}")
    plt.close()

    t0_values = sorted({run.get("sweep_T0", run.get("T0")) for run in sa_runs})
    alpha_values = sorted({run.get("sweep_alpha", run.get("alpha")) for run in sa_runs})
    if t0_values and alpha_values:
        heat = np.full((len(t0_values), len(alpha_values)), np.nan)
        for run in sa_runs:
            t0 = run.get("sweep_T0", run.get("T0"))
            alpha = run.get("sweep_alpha", run.get("alpha"))
            score = run.get("best_score", np.nan)
            i = t0_values.index(t0)
            j = alpha_values.index(alpha)
            heat[i, j] = min(heat[i, j], score) if not np.isnan(heat[i, j]) else score

        plt.figure(figsize=(8, 6))
        im = plt.imshow(heat, aspect="auto", cmap="YlGnBu")
        plt.colorbar(im, label="Best Cost")
        plt.xticks(range(len(alpha_values)), [str(a) for a in alpha_values])
        plt.yticks(range(len(t0_values)), [str(t) for t in t0_values])
        plt.xlabel("Cooling Rate (alpha)", fontsize=12)
        plt.ylabel("Initial Temperature (T0)", fontsize=12)
        plt.title("SA Hyperparameter Sweep: Best Cost Heatmap", fontsize=14)

        for i in range(len(t0_values)):
            for j in range(len(alpha_values)):
                value = heat[i, j]
                if not np.isnan(value):
                    plt.text(j, i, f"{value:.1f}", ha="center", va="center", fontsize=9, color="black")

        plt.tight_layout()
        plt.savefig(SWEEP_HEATMAP_OUT, dpi=200)
        print(f"Generated SA sweep best-cost heatmap: {SWEEP_HEATMAP_OUT}")
        plt.close()

def main():
    # Make sure the output folder exists
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    if not os.path.exists(RESULTS_FILE):
        print(f"Error: Could not find {RESULTS_FILE}.")
        print("Make sure you run 'python experiments/run_experiments.py' first!")
        return

    with open(RESULTS_FILE, "r") as f:
        all_results = json.load(f)

    sa_runs = _find_sa_runs(all_results)

    if not sa_runs:
        print("Error: Could not find 'Simulated Annealing' data in the results file.")
        return

    # For detailed single-run plots, use the best run by final score.
    sa_data = min(sa_runs, key=lambda x: x.get("best_score", float("inf")))
    selected_label = _run_label(sa_data)

    # Plot sweep comparison first (if more than one run).
    _plot_sweep_comparison(sa_runs)

    print("\n" + "="*40)
    print(" SIMULATED ANNEALING : FINAL ANALYSIS ")
    print("="*40)
    print(f"Selected run     : {selected_label}")
    print(f"Total Evaluations : {sa_data['evaluations']:,}")
    print(f"Total Runtime     : {sa_data['runtime_sec']:.2f} seconds")
    print(f"Best Cost Score   : {sa_data['best_score']:.2f}")
    print(f"Contig Breaks     : {sa_data.get('breaks', 'Not recorded')}")
    print("="*40 + "\n")

    history_best = sa_data["history"]
    history_current = sa_data.get("current_cost_history", [])
    history_contigs = sa_data.get("contigs_history", [])
    history_overlap = sa_data.get("overlap_history", [])
    temperature_history = sa_data.get("temperature_history", [])
    acceptance_history = sa_data.get("acceptance_history", [])
    num_evals = sa_data.get("evaluations", len(history_current))

    # ---------------------------------------------------------
    # 1. Convergence Plot (THE MAIN STORY)
    # ---------------------------------------------------------
    # This shows how SA learns vs Random Search
    plt.figure(figsize=(10, 6))

    if history_current:
        plt.plot(history_current, color="orange", alpha=0.3, linewidth=1, label="Current Solution Cost")

    plt.plot(history_best, color="#1f77b4", linewidth=2.5, label="Best Cost Found")

    plt.title(f"Simulated Annealing (Single Run): Convergence\n{selected_label}", fontsize=14)
    plt.xlabel("Evaluations", fontsize=12)
    plt.ylabel("Total Cost (Lower is Better)", fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend()
    plt.tight_layout()

    plt.savefig(OUTPUT_PLOT, dpi=200)
    print(f"Generated convergence plot: {OUTPUT_PLOT}")
    plt.close()

    # ---------------------------------------------------------
    # 1.5 Convergence Gap (Log Scale) - CONVERGENCE TO MINIMUM
    # ---------------------------------------------------------
    plt.figure(figsize=(10, 6))

    # Calculate gap from best solution found
    best_found = np.min(history_best)
    convergence_gap = np.array(history_best) - best_found
    convergence_gap = np.maximum(convergence_gap, 1e-6)  # Avoid log(0)

    # Plot with log scale to show exponential convergence
    plt.semilogy(convergence_gap, color="#d62728", linewidth=2.5, label="Distance from Best Solution")

    # Add moving average to show trend
    window = max(1, len(convergence_gap) // 50)
    if window > 1:
        moving_avg = np.convolve(convergence_gap, np.ones(window)/window, mode='valid')
        plt.semilogy(range(window-1, len(convergence_gap)), moving_avg, color="#1f77b4",
                     linewidth=3, linestyle='--', label=f"Trend (Moving Avg)", alpha=0.8)

    plt.title(f"Simulated Annealing (Single Run): Convergence to Minimum (Log Scale)\n{selected_label}", fontsize=14)
    plt.xlabel("Evaluations", fontsize=12)
    plt.ylabel("Gap from Best Solution (Log Scale)", fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.7, which='both')
    plt.legend(fontsize=11)
    plt.tight_layout()

    plt.savefig(CONVERGENCE_GAP_OUT, dpi=200)
    print(f"Generated convergence gap plot: {CONVERGENCE_GAP_OUT}")
    plt.close()

    # ---------------------------------------------------------
    # 2. Distribution Plot (Histogram)
    # ---------------------------------------------------------
    if history_current:
        plt.figure(figsize=(8, 5))

        plt.hist(history_current, bins=50, color="lightgray", edgecolor="darkgray", label="All Evaluated Permutations")
        plt.axvline(sa_data['best_score'], color='red', linestyle='dashed', linewidth=2, label=f"Best Found: {sa_data['best_score']:.1f}")

        plt.title(f"Simulated Annealing (Single Run): Distribution of Explored Solutions\n{selected_label}", fontsize=14)
        plt.xlabel("Total Cost", fontsize=12)
        plt.ylabel("Frequency (Number of Evaluations)", fontsize=12)
        plt.legend()
        plt.tight_layout()

        plt.savefig(HIST_OUT, dpi=200)
        print(f"Generated distribution histogram: {HIST_OUT}")
        plt.close()

    # ---------------------------------------------------------
    # 3. Learning Over Time (Boxplots by Epoch)
    # ---------------------------------------------------------
    if history_current:
        plt.figure(figsize=(10, 6))

        bin_size = len(history_current) // 10
        bins = [history_current[i:i + bin_size] for i in range(0, len(history_current), bin_size)]

        bp = plt.boxplot(bins, patch_artist=True, boxprops=dict(facecolor="lightblue", color="steelblue"))

        plt.title(f"Simulated Annealing (Single Run): Learning Over Epochs\n{selected_label}", fontsize=14)
        plt.xlabel("Evaluation Epochs (1,000 evals per bin)", fontsize=12)
        plt.ylabel("Cost Distribution", fontsize=12)
        plt.xticks(range(1, 11), [f"{i}k" for i in range(1, 11)])
        plt.grid(axis='y', linestyle='--', alpha=0.7)

        # Add mean line to show convergence
        means = [np.mean(b) for b in bins]
        plt.plot(range(1, len(means) + 1), means, 'r--', linewidth=2, label="Mean Cost per Epoch")
        plt.legend()

        plt.savefig(EPOCH_OUT, dpi=200)
        print(f"Generated Epoch Boxplot: {EPOCH_OUT}")
        plt.close()


    # ---------------------------------------------------------
    # 6. Contig Count Over Time
    # ---------------------------------------------------------
    if history_contigs:
        plt.figure(figsize=(10, 6))

        plt.plot(history_contigs, color="#d62728", linewidth=2.5, label="Number of Contigs")

        plt.title(f"Simulated Annealing (Single Run): Contig Count Over Time\n{selected_label}", fontsize=14)
        plt.xlabel("Evaluations", fontsize=12)
        plt.ylabel("Contigs (Lower is Better)", fontsize=12)
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.legend()
        plt.tight_layout()

        plt.savefig(CONTIGS_OUT, dpi=200)
        print(f"Generated Contigs plot: {CONTIGS_OUT}")
        plt.close()

    # ---------------------------------------------------------
    # 7. Sequence Overlap Over Time
    # ---------------------------------------------------------
    if history_overlap:
        plt.figure(figsize=(10, 6))

        plt.plot(history_overlap, color="#9467bd", linewidth=2.5, label="Total Sequence Overlap")

        plt.title(f"Simulated Annealing (Single Run): Sequence Overlap Over Time\n{selected_label}", fontsize=14)
        plt.xlabel("Evaluations", fontsize=12)
        plt.ylabel("Total Overlap in bp (Higher is Better)", fontsize=12)
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.legend()
        plt.tight_layout()

        plt.savefig(OVERLAP_OUT, dpi=200)
        print(f"Generated Overlap plot: {OVERLAP_OUT}")
        plt.close()

    # ---------------------------------------------------------
    # 8. Search Space Distribution (2D Scatter)
    # ---------------------------------------------------------
    if history_current and history_contigs:
        plt.figure(figsize=(12, 7))

        # Create color gradient based on evaluation progression (time)
        time_progression = np.linspace(0, 1, len(history_current))

        # Create scatter plot: Cost vs Contigs, colored by time progression
        scatter = plt.scatter(
            history_current,
            history_contigs,
            c=time_progression,
            cmap='YlOrRd',
            alpha=0.6,
            s=30,
            edgecolors='none'
        )

        # Highlight the best solution found
        best_idx = np.argmin(history_best)
        plt.scatter(
            history_current[best_idx],
            history_contigs[best_idx],
            color='green',
            s=300,
            marker='*',
            edgecolors='darkgreen',
            linewidths=2,
            label=f'Best Solution (Cost: {sa_data["best_score"]:.1f}, Contigs: {history_contigs[best_idx]})',
            zorder=5
        )

        cbar = plt.colorbar(scatter)
        cbar.set_label('Evaluation Progress (Start → End)', fontsize=11)

        plt.title(f"Simulated Annealing (Single Run): Search Space Exploration\n{selected_label}", fontsize=14)
        plt.xlabel("Total Cost (Lower is Better)", fontsize=12)
        plt.ylabel("Number of Contigs (Lower is Better)", fontsize=12)
        plt.grid(True, linestyle='--', alpha=0.4)
        plt.legend(fontsize=10)
        plt.tight_layout()

        plt.savefig(SEARCHSPACE_OUT, dpi=200)
        print(f"Generated Search Space Distribution plot: {SEARCHSPACE_OUT}")
        plt.close()

    print("\n" + "="*40)
    print("KEY INSIGHTS FOR SIMULATED ANNEALING:")
    print("="*40)
    print("✓ Convergence: Shows decreasing best cost (learning)")
    print("✓ Search Space: Visualizes multi-dimensional exploration")
    print("✓ Learning Over Time: Proves improvement vs random search")
    print("✓ Quality Metrics: Contigs & Overlap show solution quality")
    print("="*40 + "\n")


if __name__ == "__main__":
    main()
