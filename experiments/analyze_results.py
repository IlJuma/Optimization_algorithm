import json
import os
import matplotlib.pyplot as plt
import pandas as pd

# -------------------------
# Config
# -------------------------
RESULTS_FILE = "reports/experiment_results.json"
OUTPUT_DIR = "reports/"

# Plot styling
COLORS = {
    "Random Search": "#1f77b4",       # Blue
    "Simulated Annealing": "#ff7f0e", # Orange
    "Genetic Algorithm": "#2ca02c"    # Green
}

def ensure_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

# -------------------------
# Plotting functions
# -------------------------

def plot_convergence(results):
    """Plots the cost history over evaluations for all algorithms."""
    plt.figure(figsize=(10, 6))
    
    for res in results:
        method = res["method"]
        history = res["history"]
        color = COLORS.get(method, "#333333") # Default to dark gray if method not in dict
        
        plt.plot(history, label=method, color=color, linewidth=2)
        
    plt.title("Optimization Convergence: Algorithm Comparison", fontsize=14)
    plt.xlabel("Evaluations", fontsize=12)
    plt.ylabel("Total Cost (Lower is Better)", fontsize=12)
    plt.legend(fontsize=11)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    
    out_path = os.path.join(OUTPUT_DIR, "convergence_plot.png")
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Saved convergence plot to {out_path}")

def plot_bar_chart(results, metric, title, ylabel, filename):
    """Generic function to plot bar charts for best scores or runtimes."""
    methods = [res["method"] for res in results]
    values = [res[metric] for res in results]
    bar_colors = [COLORS.get(m, "#333333") for m in methods]
    
    plt.figure(figsize=(8, 5))
    bars = plt.bar(methods, values, color=bar_colors, alpha=0.8)
    
    plt.title(title, fontsize=14)
    plt.ylabel(ylabel, fontsize=12)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    
    # Add value labels on top of the bars
    for bar in bars:
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2, yval, 
                 f'{yval:.2f}', ha='center', va='bottom', fontsize=10)
                 
    plt.tight_layout()
    out_path = os.path.join(OUTPUT_DIR, filename)
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Saved {metric} bar chart to {out_path}")

# ------------------------------
# Aggregation and Summary Table
# ------------------------------

def generate_summary_table(results):
    """Creates a clean CSV and Markdown summary of the final metrics."""
    data = []
    for res in results:
        data.append({
            "Method": res["method"],
            "Best Cost": res["best_score"],
            "Evaluations": res["evaluations"],
            "Runtime (sec)": round(res["runtime_sec"], 2)
        })
        
    df = pd.DataFrame(data)
    
    # save to CSV
    csv_path = os.path.join(OUTPUT_DIR, "summary_metrics.csv")
    df.to_csv(csv_path, index=False)
    print(f"Saved summary CSV to {csv_path}")
    
    # Save summary to md
    md_path = os.path.join(OUTPUT_DIR, "summary_metrics.md")
    with open(md_path, "w") as f:
        f.write(df.to_markdown(index=False))
    print(f"Saved summary Markdown to {md_path}")

# ------------------------
# Main execution
# ------------------------

def main():
    ensure_dir(OUTPUT_DIR)
    
    if not os.path.exists(RESULTS_FILE):
        print(f"Error: {RESULTS_FILE} not found.")
        print("Waiting for run_experiments.py to generate the data.")
        return
        
    print(f"Loading results from {RESULTS_FILE}...")
    with open(RESULTS_FILE, "r") as f:
        results = json.load(f)
        
    # Generate all outputs
    plot_convergence(results)
    plot_bar_chart(results, "best_score", "Final Optimized Cost by Algorithm", "Cost", "best_cost_comparison.png")
    plot_bar_chart(results, "runtime_sec", "Execution Time by Algorithm", "Seconds", "runtime_comparison.png")
    generate_summary_table(results)
    
    print("\nPipeline execution complete. Check the 'reports/' directory.")

if __name__ == "__main__":
    main()