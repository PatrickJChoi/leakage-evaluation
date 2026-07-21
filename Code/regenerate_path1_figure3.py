import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc

def main():
    # 1. Load the data
    data_path = os.path.join("Code", "roc_data.npz")
    print(f"Loading data from {data_path}...")
    data = np.load(data_path)
    y_true = data['y_true']
    y_probs = data['y_probs']
    
    # 2. Compute ROC and AUC directly using sklearn
    fpr, tpr, thresholds = roc_curve(y_true, y_probs)
    computed_auc = auc(fpr, tpr)
    
    print("\n" + "="*50)
    print("COMPUTED ROC METRICS:")
    print(f"  True labels shape : {y_true.shape}")
    print(f"  Probs shape       : {y_probs.shape}")
    print(f"  Computed AUC value: {computed_auc:.6f}")
    print(f"  Computed AUC (3dp): {computed_auc:.3f}")
    print("="*50)
    
    # 3. Regenerate path1_figure3.png from scratch
    plt.figure(figsize=(8, 6), dpi=300)
    
    # Clean, dark blue ROC line
    legend_label = f"Strict Nested LOO-CV (AUC = {computed_auc:.3f})"
    plt.plot(fpr, tpr, color='darkblue', lw=3, label=legend_label)
    
    # Diagonal gray dashed line representing random chance
    plt.plot([0, 1], [0, 1], color='gray', lw=2, linestyle='--', label='Random Chance (AUC = 0.50)')
    
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate', fontsize=14, fontweight='semibold')
    plt.ylabel('True Positive Rate', fontsize=14, fontweight='semibold')
    plt.title('Receiver Operating Characteristic (ROC)', fontsize=16, weight='bold', pad=15)
    plt.legend(loc="lower right", fontsize=12)
    plt.grid(True, alpha=0.3)
    
    output_path = "path1_figure3.png"
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"\nSuccessfully regenerated and saved figure to: {output_path}")
    print(f"Legend text written on image: '{legend_label}'")
    print("="*50)

if __name__ == "__main__":
    main()
