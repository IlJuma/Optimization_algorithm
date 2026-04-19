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
    "Oracle Ground Truth": "#d62728", # Red
    "Random Search": "#1f77b4",       # Blue
    "Simulated Annealing": "#ff7f0e", # Orange
    "Genetic Algorithm": "#2ca02c"    # Green
}

def ensure_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

# Plotting functions --------------

def plot_metric_convergence(results, metric_key, title, ylabel, filename, oracle_value=None, lower_is_better=True):
    """Generic function to plot line graphs for Cost, Contigs, and Overlap."""
    plt.figure(figsize=(10, 6))
    has_data = False
    
    for res in results:
        method = res["method"]
        
        # skip oracle
        if method == "Oracle Ground Truth":
            continue
            
        history = res.get(metric_key, [])
        total_evals = res.get("evaluations", len(history)) # Get the true evaluation count
        
        if history:
            has_data = True
            color = COLORS.get(method, "#333333")
            
            # fix for GA
            # This generates a mathematically spaced X-axis. 
            # If GA has 200 history items but did 10,000 evals, this correctly 
            # spaces the points at 50, 100, 150... up to 10,000
            x_axis = np.linspace(1, total_evals, len(history))
            
            plt.plot(x_axis, history, label=method, color=color, linewidth=2.5, alpha=0.9)
            
    if not has_data:
        plt.close()
        return
    
    # Draw the horizontal dashed line
    if oracle_value is not None:
        plt.axhline(y=oracle_value, color='black', linestyle='--', linewidth=2, label=f'Oracle Baseline ({oracle_value:.1f})')

    plt.title(title, fontsize=14, fontweight='bold')
    plt.xlabel("Evaluations", fontsize=12)
    plt.ylabel(ylabel, fontsize=12)
    plt.legend(fontsize=11)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    
    out_path = os.path.join(OUTPUT_DIR, filename)
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Generated: {out_path}")

def plot_bar_metrics(results):
    """Creates a multi-panel bar chart comparing Final Cost, Breaks, and Runtime."""
    methods = [res["method"] for res in results]
    costs = [res["best_score"] for res in results]
    # Oracle uses 'best_breaks', algorithms use 'breaks'
    breaks = [res.get("breaks", res.get("best_breaks", 0)) for res in results]
    runtimes = [res["runtime_sec"] for res in results]
    bar_colors = [COLORS.get(m, "#333333") for m in methods]
    
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
    
    # Panel 1: Best Cost
    bars1 = ax1.bar(methods, costs, color=bar_colors, alpha=0.8)
    ax1.set_title("Final Optimized Cost", fontweight='bold')
    ax1.set_ylabel("Cost")
    ax1.grid(axis='y', linestyle='--', alpha=0.7)
    # Automatically rotate labels if they start overlapping (helpful now that we added Oracle)
    ax1.set_xticks(range(len(methods)))
    ax1.set_xticklabels(methods, rotation=15, ha="right") 
    for bar in bars1:
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height(), 
                 f'{bar.get_height():.1f}', ha='center', va='bottom')

    # Panel 2: Contig Breaks
    bars2 = ax2.bar(methods, breaks, color=bar_colors, alpha=0.8)
    ax2.set_title("Remaining Contig Breaks", fontweight='bold')
    ax2.set_ylabel("Number of Breaks")
    ax2.grid(axis='y', linestyle='--', alpha=0.7)
    ax2.set_xticks(range(len(methods)))
    ax2.set_xticklabels(methods, rotation=15, ha="right")
    for bar in bars2:
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height(), 
                 f'{int(bar.get_height())}', ha='center', va='bottom')

    # Panel 3: Runtime
    bars3 = ax3.bar(methods, runtimes, color=bar_colors, alpha=0.8)
    ax3.set_title("Execution Time", fontweight='bold')
    ax3.set_ylabel("Seconds")
    ax3.grid(axis='y', linestyle='--', alpha=0.7)
    ax3.set_xticks(range(len(methods)))
    ax3.set_xticklabels(methods, rotation=15, ha="right")
    for bar in bars3:
        ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height(), 
                 f'{bar.get_height():.2f}s', ha='center', va='bottom')

    plt.tight_layout()
    out_path = os.path.join(OUTPUT_DIR, "4_metrics_comparison.png")
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
        # The Oracle doesn't have a 'current_cost_history', so it naturally gets skipped here.
        if "current_cost_history" in res and len(res["current_cost_history"]) > 0:
            data_to_plot.append(res["current_cost_history"])
            labels.append(res["method"])
            colors_to_use.append(COLORS.get(res["method"], "#333333"))
            
    if not data_to_plot:
        plt.close()
        return 
        
    parts = plt.violinplot(data_to_plot, showmeans=True)
    
    for pc, color in zip(parts['bodies'], colors_to_use):
        pc.set_facecolor(color)
        pc.set_edgecolor('black')
        pc.set_alpha(0.6)
        
    parts['cmeans'].set_color('black')
    
    plt.xticks(np.arange(1, len(labels) + 1), labels, fontsize=11)
    plt.title("Search Space Exploration Distribution", fontsize=14, fontweight='bold')
    plt.ylabel("Evaluated Costs", fontsize=12)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    
    out_path = os.path.join(OUTPUT_DIR, "5_exploration_distribution.png")
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Generated: {out_path}")

def generate_summary_table(results):
    """Generates a clean markdown table summarizing the results."""
    data = []
    for res in results:
        # Handle the different dictionary keys between the Oracle and the Optimization Algorithms
        if res["method"] == "Oracle Ground Truth":
            final_contigs = res.get("best_contigs", "N/A")
            final_overlap = res.get("best_total_overlap", "N/A")
            breaks = res.get("best_breaks", "N/A")
        else:
            contigs_hist = res.get("contigs_history", [])
            overlap_hist = res.get("overlap_history", [])
            final_contigs = contigs_hist[-1] if contigs_hist else "N/A"
            final_overlap = overlap_hist[-1] if overlap_hist else "N/A"
            breaks = res.get("breaks", "N/A")
        
        data.append({
            "Method": res["method"],
            "Best Cost": round(res["best_score"], 2),
            "Breaks": breaks,
            "Final Contigs": final_contigs,
            "Final Overlap": final_overlap,
            "Evaluations": res["evaluations"],
            "Runtime (s)": round(res["runtime_sec"], 2)
        })
        
    df = pd.DataFrame(data)
    md_path = os.path.join(OUTPUT_DIR, "6_summary_metrics.md")
    with open(md_path, "w") as f:
        f.write("### Final Optimization Results\n\n")
        f.write(df.to_markdown(index=False))
        
    print(f"Generated: {md_path}")

# Main execution --------------

def main():
    ensure_dir(OUTPUT_DIR)
    
    if not os.path.exists(RESULTS_FILE):
        print(f"Error: {RESULTS_FILE} not found.")
        print("Run 'python experiments/run_experiments.py' first!")
        return
        
    print(f"Loading results from {RESULTS_FILE}...")
    with open(RESULTS_FILE, "r") as f:
        results = json.load(f)

    # extract Oracle metrics for horizontal lines in the plots
    oracle_data = {}
    for res in results:
        if res["method"] == "Oracle Ground Truth":
            oracle_data["best_score"] = res["best_score"]
            oracle_data["contigs"] = res.get("best_contigs")
            oracle_data["overlap"] = res.get("best_total_overlap")
            break # Once we find it, stop searching
        
    print("\nGenerating Comparison Plots...")
    
    # Generate the 3 convergence lines, passing our dynamically extracted Oracle targets
    plot_metric_convergence(results, "history", "Convergence: Overall Cost", "Cost (Lower is Better)", "1_cost_convergence.png", oracle_value=oracle_data.get("best_score"))
    plot_metric_convergence(results, "contigs_history", "Convergence: Contig Count", "Number of Contigs (Lower is Better)", "2_contigs_convergence.png", oracle_value=oracle_data.get("contigs"))
    plot_metric_convergence(results, "overlap_history", "Convergence: Total Sequence Overlap", "Total Overlap (Higher is Better)", "3_overlap_convergence.png", oracle_value=oracle_data.get("overlap"), lower_is_better=False)
   
    # Generate the bar chart, violin plot, and table
    plot_bar_metrics(results)
    plot_search_space_distributions(results)
    generate_summary_table(results)
    
    print(f"\nPipeline complete! All plots saved to: {OUTPUT_DIR}")

if __name__ == "__main__":
    main()