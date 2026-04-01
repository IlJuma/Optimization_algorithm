import json
import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# config
RESULTS_FILE = "reports/experiment_results.json"
OUTPUT_DIR = "reports/comparisons/"

# Consistent color coding for the report
COLORS = {
    "Random Search": "#1f77b4",       # Blue
    "Simulated Annealing": "#ff7f0e", # Orange
    "Genetic Algorithm": "#2ca02c"    # Green
}

def ensure_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

# Plotting functions

def plot_convergence(results):
    """Plots the 'best cost' history over evaluations for all algorithms."""
    plt.figure(figsize=(10, 6))
    
    for res in results:
        method = res["method"]
        history = res.get("history", [])
        color = COLORS.get(method, "#333333")
        
        # Plot the staircase of best scores
        plt.plot(history, label=method, color=color, linewidth=2.5, alpha=0.9)
        
    plt.title("Optimization Convergence: Algorithm Comparison", fontsize=14, fontweight='bold')
    plt.xlabel("Evaluations", fontsize=12)
    plt.ylabel("Total Cost (Lower is Better)", fontsize=12)
    plt.legend(fontsize=11)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    
    out_path = os.path.join(OUTPUT_DIR, "1_convergence_comparison.png")
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Generated: {out_path}")

def plot_bar_metrics(results):
    """Creates a multi-panel bar chart comparing Final Cost, Breaks, and Runtime."""
    methods = [res["method"] for res in results]
    costs = [res["best_score"] for res in results]
    breaks = [res.get("breaks", 0) for res in results]
    runtimes = [res["runtime_sec"] for res in results]
    bar_colors = [COLORS.get(m, "#333333") for m in methods]
    
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
    
    # Panel 1: Best Cost
    bars1 = ax1.bar(methods, costs, color=bar_colors, alpha=0.8)
    ax1.set_title("Final Optimized Cost", fontweight='bold')
    ax1.set_ylabel("Cost")
    ax1.grid(axis='y', linestyle='--', alpha=0.7)
    for bar in bars1:
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height(), 
                 f'{bar.get_height():.1f}', ha='center', va='bottom')

    # Panel 2: Contig Breaks
    bars2 = ax2.bar(methods, breaks, color=bar_colors, alpha=0.8)
    ax2.set_title("Remaining Contig Breaks", fontweight='bold')
    ax2.set_ylabel("Number of Breaks")
    ax2.grid(axis='y', linestyle='--', alpha=0.7)
    for bar in bars2:
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height(), 
                 f'{int(bar.get_height())}', ha='center', va='bottom')

    # Panel 3: Runtime
    bars3 = ax3.bar(methods, runtimes, color=bar_colors, alpha=0.8)
    ax3.set_title("Execution Time", fontweight='bold')
    ax3.set_ylabel("Seconds")
    ax3.grid(axis='y', linestyle='--', alpha=0.7)
    for bar in bars3:
        ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height(), 
                 f'{bar.get_height():.2f}s', ha='center', va='bottom')

    plt.tight_layout()
    out_path = os.path.join(OUTPUT_DIR, "2_metrics_comparison.png")
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Generated: {out_path}")

def plot_search_space_distributions(results):
    """Plots Violin plots of the 'current_cost_history' to prove if an algorithm 'learns'."""
    plt.figure(figsize=(10, 6))
    
    data_to_plot = []
    labels = []
    colors_to_use = []
    
    for res in results:
        # Check if the algorithm exported its noisy exploration history
        if "current_cost_history" in res and len(res["current_cost_history"]) > 0:
            data_to_plot.append(res["current_cost_history"])
            labels.append(res["method"])
            colors_to_use.append(COLORS.get(res["method"], "#333333"))
            
    if not data_to_plot:
        return # Skip if no algorithms exported current_cost_history
        
    parts = plt.violinplot(data_to_plot, showmeans=True)
    
    # Color the violins
    for pc, color in zip(parts['bodies'], colors_to_use):
        pc.set_facecolor(color)
        pc.set_edgecolor('black')
        pc.set_alpha(0.6)
        
    parts['cmeans'].set_color('black')
    
    plt.xticks(np.arange(1, len(labels) + 1), labels, fontsize=11)
    plt.title("Search Space Exploration Distribution", fontsize=14, fontweight='bold')
    plt.ylabel("Evaluated Costs (Where did the algorithm spend its time?)", fontsize=12)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    
    out_path = os.path.join(OUTPUT_DIR, "3_exploration_distribution.png")
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Generated: {out_path}")

def generate_summary_table(results):
    """Creates a Markdown summary of the final metrics for your report."""
    data = []
    for res in results:
        data.append({
            "Method": res["method"],
            "Best Cost": round(res["best_score"], 2),
            "Breaks": res.get("breaks", "N/A"),
            "Evaluations": res["evaluations"],
            "Runtime (s)": round(res["runtime_sec"], 2)
        })
        
    df = pd.DataFrame(data)
    
    md_path = os.path.join(OUTPUT_DIR, "4_summary_metrics.md")
    with open(md_path, "w") as f:
        f.write("### Final Optimization Results\n\n")
        f.write(df.to_markdown(index=False))
        
    print(f"Generated: {md_path}")

# 
# MAIN EXECUTION ---------------

def main():
    ensure_dir(OUTPUT_DIR)
    
    if not os.path.exists(RESULTS_FILE):
        print(f"Error: {RESULTS_FILE} not found.")
        print("Run 'python experiments/run_experiments.py' first!")
        return
        
    print(f"Loading results from {RESULTS_FILE}...")
    with open(RESULTS_FILE, "r") as f:
        results = json.load(f)
        
    print("\nGenerating Comparison Plots...")
    plot_convergence(results)
    plot_bar_metrics(results)
    plot_search_space_distributions(results)
    generate_summary_table(results)
    
    print(f"\nPipeline complete! All plots saved to: {OUTPUT_DIR}")

if __name__ == "__main__":
    main()