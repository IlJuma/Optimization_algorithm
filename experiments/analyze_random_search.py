import json
import os
import matplotlib.pyplot as plt

#config
# 1. Read from the shared team folder
RESULTS_FILE = "reports/experiment_results.json" 

# 2. Save plots to your specific Random Search folder
OUTPUT_DIR = "reports/random_search"
OUTPUT_PLOT = f"{OUTPUT_DIR}/rs_convergence_plot.png"
HIST_OUT = f"{OUTPUT_DIR}/rs_distribution_plot.png"
EPOCH_OUT = f"{OUTPUT_DIR}/rs_epoch_boxplot.png"

def main():
    # Make sure your new folder actually exists before trying to save images into it!
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    if not os.path.exists(RESULTS_FILE):
        print(f"Error: Could not find {RESULTS_FILE}.")
        print("Make sure you run 'python experiments/run_experiments.py' first!")
        return

    with open(RESULTS_FILE, "r") as f:
        all_results = json.load(f)

    rs_data = next((res for res in all_results if res.get("method") == "Random Search"), None)

    if not rs_data:
        print("Error: Could not find 'Random Search' data in the results file.")
        return

    print("\n" + "="*40)
    print(" RANDOM SEARCH : FINAL ANALYSIS ")
    print("="*40)
    print(f"Total Evaluations : {rs_data['evaluations']:,}")
    print(f"Total Runtime     : {rs_data['runtime_sec']:.2f} seconds")
    print(f"Best Cost Score   : {rs_data['best_score']:.2f}")
    print(f"Contig Breaks     : {rs_data.get('breaks', 'Not recorded')}")
    print("="*40 + "\n")

    history_best = rs_data["history"]
    history_current = rs_data.get("current_cost_history", [])
    
    # ---------------------------------------------------------
    # 1. Convergence Plot
    # ---------------------------------------------------------
    plt.figure(figsize=(10, 6))
    
    if history_current:
        plt.plot(history_current, color="orange", alpha=0.3, linewidth=1, label="Evaluated Cost (Noise)")
    
    plt.plot(history_best, color="#1f77b4", linewidth=2.5, label="Best Cost Found")
    
    plt.title("Random Search Baseline: Exploration vs. Best Cost", fontsize=14)
    plt.xlabel("Evaluations", fontsize=12)
    plt.ylabel("Total Cost (Lower is Better)", fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend()
    plt.tight_layout()

    plt.savefig(OUTPUT_PLOT, dpi=200)
    print(f"Generated convergence plot: {OUTPUT_PLOT}")
    plt.close()

    # ---------------------------------------------------------
    # 2. Distribution Plot (Histogram)
    # ---------------------------------------------------------
    if history_current:
        plt.figure(figsize=(8, 5))
        
        plt.hist(history_current, bins=50, color="lightgray", edgecolor="darkgray", label="All Evaluated Permutations")
        plt.axvline(rs_data['best_score'], color='red', linestyle='dashed', linewidth=2, label=f"Best Found: {rs_data['best_score']:.1f}")
        
        plt.title("Random Search: Distribution of the Search Space", fontsize=14)
        plt.xlabel("Total Cost", fontsize=12)
        plt.ylabel("Frequency (Number of Guesses)", fontsize=12)
        plt.legend()
        plt.tight_layout()
        
        plt.savefig(HIST_OUT, dpi=200)
        print(f"Generated distribution histogram: {HIST_OUT}")
        plt.close()

    # ---------------------------------------------------------
    # 3. Lack of Learning Plot (Boxplots)
    # ---------------------------------------------------------
    if history_current:
        plt.figure(figsize=(10, 6))
        
        bin_size = len(history_current) // 10
        bins = [history_current[i:i + bin_size] for i in range(0, len(history_current), bin_size)]
        
        plt.boxplot(bins, patch_artist=True, boxprops=dict(facecolor="lightgray", color="gray"))
        
        plt.title("Random Search Epochs: Proving a Lack of 'Learning'", fontsize=14)
        plt.xlabel("Evaluation Epochs (1,000 evals per bin)", fontsize=12)
        plt.ylabel("Cost Distribution", fontsize=12)
        plt.xticks(range(1, 11), [f"{i}k" for i in range(1, 11)])
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        
        plt.savefig(EPOCH_OUT, dpi=200)
        print(f"Generated Epoch Boxplot: {EPOCH_OUT}")
        plt.close()

if __name__ == "__main__":
    main()