import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

output_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "Results"))

sns.set_theme(style="whitegrid", context="paper", font_scale=1.2)

def generate_figure2():
    import pickle
    results_path = os.path.join(os.path.dirname(__file__), "ibd_correct_results.pkl")
    with open(results_path, "rb") as f:
        result = pickle.load(f)
    consensus_counts = result['consensus_counts']
    filtered = sorted([(g, c) for g, c in consensus_counts.items() if c >= 5], key=lambda kv: -kv[1])
    features = [g for g, _ in filtered]
    counts = [c for _, c in filtered]

    features.reverse()
    counts.reverse()

    plt.figure(figsize=(10, 6))
    ax = sns.barplot(x=counts, y=features, hue=features, palette="Blues_d", legend=False)

    plt.axvline(x=36, color='red', linestyle='--', linewidth=2, label='Maximum Folds (n=36)')

    plt.title('Robustness of Multimodal Feature Selection Across LOO-CV Folds', fontsize=16, weight='bold')
    plt.xlabel('Number of Folds Selected (out of 36)', fontsize=14)
    plt.ylabel('Features', fontsize=14)
    plt.xlim(0, 40)
    plt.legend(loc='lower right')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'path1_figure2_stability.png'), dpi=300, bbox_inches='tight')
    plt.close()

def generate_figure4():
    import pickle
    results_path = os.path.join(os.path.dirname(__file__), "ibd_correct_results.pkl")
    with open(results_path, "rb") as f:
        result = pickle.load(f)
    null_dist = result['perm_accuracies']
    observed_accuracy = result['accuracy']
    p_value = result['p_value']

    plt.figure(figsize=(10, 6))

    sns.histplot(null_dist, bins=30, kde=True, color='lightgray', edgecolor='black')

    plt.axvline(x=observed_accuracy, color='red', linestyle='-', linewidth=4)

    plt.text(observed_accuracy - 0.015, plt.ylim()[1]*0.8,
             f'Observed Accuracy: {observed_accuracy*100:.1f}%\n$p = {p_value:.4f}$\n(Not Statistically Significant)',
             color='red', fontsize=12, weight='bold', ha='right',
             bbox=dict(facecolor='white', alpha=0.8, edgecolor='red', boxstyle='round,pad=0.5'))

    plt.title('Permutation Test: Null Distribution vs. Observed Accuracy', fontsize=16, weight='bold')
    plt.xlabel('Classification Accuracy', fontsize=14)
    plt.ylabel('Frequency (500 Permutations)', fontsize=14)
    plt.xlim(0.1, 0.9)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'path1_figure4_permutation.png'), dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    import sys
    try:
        print("Generating Figure 2: Feature Stability...")
        generate_figure2()
        print("Generating Figure 4: Permutation Test...")
        generate_figure4()
        print(f"Figures successfully saved to {output_dir}")
    except Exception as e:
        import traceback
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)
