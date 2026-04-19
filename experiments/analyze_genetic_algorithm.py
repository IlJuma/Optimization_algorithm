import json
import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize


RESULTS_FILE = "reports/genetic_algorithm/ga_sweep_results.json"
OUTPUT_DIR = "reports/genetic_algorithm"
OUTPUT_PLOT = f"{OUTPUT_DIR}/ga_convergence_plot.png"
HIST_OUT = f"{OUTPUT_DIR}/ga_distribution_plot.png"
CONTIGS_OUT = f"{OUTPUT_DIR}/ga_contigs_plot.png"
OVERLAP_OUT = f"{OUTPUT_DIR}/ga_overlap_plot.png"
SWEEP_CONVERGENCE_OUT = f"{OUTPUT_DIR}/ga_sweep_convergence_comparison.png"
SWEEP_RUNTIME_OUT = f"{OUTPUT_DIR}/ga_sweep_runtime_comparison.png"
SWEEP_SCORE_OUT = f"{OUTPUT_DIR}/ga_sweep_best_cost_comparison.png"


def _run_label(run):
    if "method_label" in run:
        return run["method_label"]
    return (
        "GA "
        f"(pop={run.get('sweep_population_size', '?')}, "
        f"gen={run.get('sweep_num_generations', '?')}, "
        f"mut={run.get('sweep_mutation_rate', '?')}, "
        f"cross={run.get('sweep_crossover_rate', '?')})"
    )


def _plot_sweep_comparison(ga_runs):
    if len(ga_runs) <= 1:
        return

    plt.figure(figsize=(11, 7))
    cmap = plt.cm.viridis
    color_norm = Normalize(vmin=0, vmax=max(1, len(ga_runs) - 1))

    for idx, run in enumerate(ga_runs):
        history = run.get("history", [])
        if not history:
            continue
        plt.plot(
            history,
            linewidth=1.8,
            alpha=0.85,
            color=cmap(color_norm(idx)),
            label=_run_label(run),
        )

    plt.title("GA Hyperparameter Sweep: Best-Cost Convergence", fontsize=14)
    plt.xlabel("Generation", fontsize=12)
    plt.ylabel("Best Cost (Lower is Better)", fontsize=12)
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.legend(fontsize=8, loc="best")
    plt.tight_layout()
    plt.savefig(SWEEP_CONVERGENCE_OUT, dpi=200)
    print(f"Generated GA sweep convergence comparison: {SWEEP_CONVERGENCE_OUT}")
    plt.close()

    labels = [_run_label(run) for run in ga_runs]
    runtimes = [run.get("runtime_sec", 0.0) for run in ga_runs]
    plt.figure(figsize=(max(10, len(labels) * 1.1), 6))
    bars = plt.bar(range(len(labels)), runtimes, color="#4c78a8", alpha=0.85)
    plt.xticks(range(len(labels)), labels, rotation=35, ha="right")
    plt.title("GA Hyperparameter Sweep: Runtime Comparison", fontsize=14)
    plt.ylabel("Runtime (seconds)", fontsize=12)
    plt.grid(axis="y", linestyle="--", alpha=0.6)
    for i, bar in enumerate(bars):
        plt.text(i, bar.get_height(), f"{bar.get_height():.2f}", ha="center", va="bottom", fontsize=8)
    plt.tight_layout()
    plt.savefig(SWEEP_RUNTIME_OUT, dpi=200)
    print(f"Generated GA sweep runtime comparison: {SWEEP_RUNTIME_OUT}")
    plt.close()

    ordered_runs = sorted(ga_runs, key=lambda run: run.get("best_score", float("inf")))
    ordered_labels = [_run_label(run) for run in ordered_runs]
    ordered_scores = [run.get("best_score", np.nan) for run in ordered_runs]

    plt.figure(figsize=(max(10, len(ordered_labels) * 1.1), 6))
    bars = plt.bar(range(len(ordered_labels)), ordered_scores, color="#2ca02c", alpha=0.85)
    plt.xticks(range(len(ordered_labels)), ordered_labels, rotation=35, ha="right")
    plt.title("GA Hyperparameter Sweep: Best Cost Ranking", fontsize=14)
    plt.ylabel("Best Cost", fontsize=12)
    plt.grid(axis="y", linestyle="--", alpha=0.6)
    for i, bar in enumerate(bars):
        plt.text(i, bar.get_height(), f"{bar.get_height():.1f}", ha="center", va="bottom", fontsize=8)
    plt.tight_layout()
    plt.savefig(SWEEP_SCORE_OUT, dpi=200)
    print(f"Generated GA sweep best-cost comparison: {SWEEP_SCORE_OUT}")
    plt.close()


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    if not os.path.exists(RESULTS_FILE):
        print(f"Error: Could not find {RESULTS_FILE}.")
        print("Make sure you run 'python experiments/run_experiments.py' first!")
        return

    with open(RESULTS_FILE, "r") as f:
        ga_runs = json.load(f)

    if not ga_runs:
        print("Error: Could not find Genetic Algorithm sweep data in the results file.")
        return

    best_run = min(ga_runs, key=lambda run: run.get("best_score", float("inf")))
    selected_label = _run_label(best_run)

    _plot_sweep_comparison(ga_runs)

    print("\n" + "=" * 40)
    print(" GENETIC ALGORITHM : FINAL ANALYSIS ")
    print("=" * 40)
    print(f"Selected run      : {selected_label}")
    print(f"Total Evaluations : {best_run['evaluations']:,}")
    print(f"Total Runtime     : {best_run['runtime_sec']:.2f} seconds")
    print(f"Best Cost Score   : {best_run['best_score']:.2f}")
    print(f"Contig Breaks     : {best_run.get('breaks', 'Not recorded')}")
    print("=" * 40 + "\n")

    history_best = best_run.get("history", [])
    history_current = best_run.get("current_cost_history", [])
    history_contigs = best_run.get("contigs_history", [])
    history_overlap = best_run.get("overlap_history", [])

    if history_best:
        plt.figure(figsize=(10, 6))
        if history_current:
            plt.plot(history_current, color="lightgray", alpha=0.5, linewidth=1.2, label="Current Best per Generation")
        plt.plot(history_best, color="#2ca02c", linewidth=2.5, label="Best Cost Found")
        plt.title(f"Genetic Algorithm: Convergence\n{selected_label}", fontsize=14)
        plt.xlabel("Generation", fontsize=12)
        plt.ylabel("Total Cost (Lower is Better)", fontsize=12)
        plt.grid(True, linestyle="--", alpha=0.7)
        plt.legend()
        plt.tight_layout()
        plt.savefig(OUTPUT_PLOT, dpi=200)
        print(f"Generated convergence plot: {OUTPUT_PLOT}")
        plt.close()

    if history_current:
        plt.figure(figsize=(8, 5))
        plt.hist(history_current, bins=min(30, max(5, len(history_current))), color="lightgray", edgecolor="darkgray")
        plt.axvline(best_run["best_score"], color="red", linestyle="dashed", linewidth=2, label=f"Best Found: {best_run['best_score']:.1f}")
        plt.title(f"Genetic Algorithm: Distribution of Best Costs\n{selected_label}", fontsize=14)
        plt.xlabel("Cost", fontsize=12)
        plt.ylabel("Frequency", fontsize=12)
        plt.legend()
        plt.tight_layout()
        plt.savefig(HIST_OUT, dpi=200)
        print(f"Generated distribution histogram: {HIST_OUT}")
        plt.close()

    if history_contigs:
        plt.figure(figsize=(10, 6))
        plt.plot(history_contigs, color="#1f77b4", linewidth=2.5)
        plt.title(f"Genetic Algorithm: Contig Count Over Time\n{selected_label}", fontsize=14)
        plt.xlabel("Generation", fontsize=12)
        plt.ylabel("Contig Count", fontsize=12)
        plt.grid(True, linestyle="--", alpha=0.7)
        plt.tight_layout()
        plt.savefig(CONTIGS_OUT, dpi=200)
        print(f"Generated contigs plot: {CONTIGS_OUT}")
        plt.close()

    if history_overlap:
        plt.figure(figsize=(10, 6))
        plt.plot(history_overlap, color="#ff7f0e", linewidth=2.5)
        plt.title(f"Genetic Algorithm: Overlap Over Time\n{selected_label}", fontsize=14)
        plt.xlabel("Generation", fontsize=12)
        plt.ylabel("Total Overlap", fontsize=12)
        plt.grid(True, linestyle="--", alpha=0.7)
        plt.tight_layout()
        plt.savefig(OVERLAP_OUT, dpi=200)
        print(f"Generated overlap plot: {OVERLAP_OUT}")
        plt.close()


if __name__ == "__main__":
    main()
