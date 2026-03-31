import matplotlib.pyplot as plt

# Simulate the output dictionaries from all three algorithms
mock_results = [
    {
        "method": "Random Search",
        "history": [1000, 950, 950, 900, 800, 800, 800, 750, 750, 700]
    },
    {
        "method": "Simulated Annealing",
        "history": [1000, 900, 850, 700, 650, 500, 450, 400, 380, 350]
    },
    {
        "method": "Genetic Algorithm",
        "history": [1000, 800, 600, 400, 300, 250, 200, 180, 150, 140]
    }
]

plt.figure(figsize=(10, 6))
for res in mock_results:
    plt.plot(res["history"], label=res["method"], marker='o')

plt.title("Optimization Convergence (Mock Data)")
plt.xlabel("Evaluations (x1000)")
plt.ylabel("Total Cost (Lower is Better)")
plt.legend()
plt.grid(True)
plt.savefig("mock_convergence_plot.png")
print("Mock plot generated...")