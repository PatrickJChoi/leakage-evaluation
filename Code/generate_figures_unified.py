import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import roc_curve, auc
import os

# Ensure output directory exists
os.makedirs("Results", exist_ok=True)

# Set global plotting style for professional look
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']

# ==============================================================================
# FIGURE 5: Melanoma ICI Feature Stability Bar Chart
# ==============================================================================
print("Generating Figure 5...")
fig5, ax5 = plt.subplots(figsize=(8, 5), dpi=300)

import pickle
_consensus_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "melanoma_consensus_counts.pkl")
with open(os.path.abspath(_consensus_path), "rb") as f:
    _consensus_counts = pickle.load(f)
genes = ['SELL', 'FOXP1', 'TNRC6B', 'TBC1D10B']
counts = [_consensus_counts[g] for g in genes]
# Threshold is 80% of 18 folds = 14.4
threshold = 14.4

colors = ['#1f77b4' if c >= threshold else '#d3d3d3' for c in counts]

bars = ax5.barh(genes, counts, color=colors, edgecolor='black', height=0.6)

# Add threshold line
ax5.axvline(x=threshold, color='red', linestyle='--', linewidth=1.5, label='80% Selection Threshold (14.4 folds)')

# Add counts on the bars
for bar in bars:
    width = bar.get_width()
    ax5.text(width + 0.3, bar.get_y() + bar.get_height()/2, f'{int(width)}/18', 
             va='center', ha='left', fontsize=10, fontweight='bold',
             color='black' if width < threshold else '#1f77b4')

ax5.set_xlim(0, 20)
ax5.set_xlabel('Selection Frequency (out of 18 folds)', fontsize=12, fontweight='bold')
ax5.set_title('Melanoma ICI Consensus Feature Selection Frequency', fontsize=13, fontweight='bold', pad=15)
ax5.spines['top'].set_visible(False)
ax5.spines['right'].set_visible(False)
ax5.tick_params(axis='both', which='major', labelsize=11)
ax5.legend(loc='lower right', fontsize=10)

plt.tight_layout()
output_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "Results"))
fig5.savefig(os.path.join(output_dir, 'path1_figure5_stability.png'), dpi=300, bbox_inches='tight')
plt.close(fig5)
print("Figure 5 saved as path1_figure5_stability.png")

# ==============================================================================
# FIGURE 6: Melanoma ROC Curve
# ==============================================================================
print("Generating Figure 6...")
if os.path.exists("melanoma_roc_data.npz"):
    data = np.load("melanoma_roc_data.npz")
    y_true = data['y_true']
    y_probs = data['y_probs']
    
    fpr, tpr, _ = roc_curve(y_true, y_probs)
    roc_auc = auc(fpr, tpr)
    
    fig6, ax6 = plt.subplots(figsize=(6, 6), dpi=300)
    ax6.plot(fpr, tpr, color='#1f77b4', lw=2.5, label=f'Random Forest LOO-CV (AUC = {roc_auc:.3f})')
    ax6.plot([0, 1], [0, 1], color='gray', lw=1.5, linestyle='--', label='Random Chance (AUC = 0.50)')
    
    ax6.set_xlim([-0.02, 1.02])
    ax6.set_ylim([-0.02, 1.02])
    ax6.set_xlabel('False Positive Rate (FPR)', fontsize=12, fontweight='bold', labelpad=8)
    ax6.set_ylabel('True Positive Rate (TPR)', fontsize=12, fontweight='bold', labelpad=8)
    ax6.set_title('Melanoma ICI LOO-CV ROC Curve', fontsize=13, fontweight='bold', pad=15)
    ax6.legend(loc='lower right', fontsize=10, frameon=True, shadow=False)
    ax6.grid(True, linestyle=':', alpha=0.6)
    
    # Square aspect ratio
    ax6.set_aspect('equal')
    ax6.spines['top'].set_visible(False)
    ax6.spines['right'].set_visible(False)
    
    plt.tight_layout()
    fig6.savefig(os.path.join(output_dir, 'figure6_melanoma_roc.png'), dpi=300, bbox_inches='tight')
    plt.close(fig6)
    print("Figure 6 saved as figure6_melanoma_roc.png")
else:
    print("WARNING: melanoma_roc_data.npz not found, skipping Figure 6.")

# ==============================================================================
# FIGURE 1: Pipeline Schematic
# ==============================================================================
print("Generating Figure 1...")
fig1, ax1 = plt.subplots(figsize=(14, 8), dpi=300)
ax1.axis('off')
ax1.set_xlim(0, 14)
ax1.set_ylim(0, 8)

# Draw Left Panel (Naive Pipeline - Leaky)
# Title Banner
ax1.text(3.5, 7.3, "Naive Pipeline (Leaky)", color='white', weight='bold', fontsize=14, 
         bbox=dict(facecolor='#d9534f', edgecolor='none', boxstyle='round,pad=0.5'), ha='center')

# Left Panel Boxes
box_width = 4.8
box_height = 0.7
left_x = 1.1

steps_left = [
    ("All Cells", 6.0),
    ("Global HVG Selection\n(Uses all cells)", 4.8),
    ("Global Scaling\n(Uses all cells)", 3.6),
    ("LOO-CV Split\n(Train / Test)", 2.4),
    ("Model Training & Evaluation", 1.2)
]

import matplotlib.patches as patches
for text, y in steps_left:
    # Draw box
    ax1.add_patch(patches.FancyBboxPatch((left_x, y), box_width, box_height, boxstyle='round,pad=0.1', facecolor='#fcf8f2', edgecolor='#d9534f', lw=1.5))
    # Draw text
    ax1.text(left_x + box_width/2, y + box_height/2, text, ha='center', va='center', color='black', fontsize=10, weight='semibold')
    
    # Draw arrow down (if not the last step)
    if y > 1.2:
        ax1.annotate('', xy=(left_x + box_width/2, y - 0.4), xytext=(left_x + box_width/2, y),
                    arrowprops=dict(arrowstyle="->", color='#d9534f', lw=1.5))

# Red dashed leak border around steps 2 and 3
ax1.add_patch(plt.Rectangle((left_x - 0.2, 3.4), box_width + 0.4, 2.3, facecolor='none', edgecolor='red', linestyle='--', lw=2))

# Leak callout box
ax1.text(0.1, 4.3, "TEST DATA\nLEAKS HERE!", color='white', weight='bold', fontsize=11, ha='center', va='center',
         bbox=dict(facecolor='red', edgecolor='none', boxstyle='round,pad=0.4'))
ax1.annotate('', xy=(left_x - 0.05, 4.5), xytext=(0.6, 4.5),
            arrowprops=dict(arrowstyle="->", color='red', lw=2))


# Draw Right Panel (Corrected Pipeline - Leakage-Free)
# Title Banner
ax1.text(10.5, 7.3, "Corrected Pipeline (Leakage-Free)", color='white', weight='bold', fontsize=14, 
         bbox=dict(facecolor='#5cb85c', edgecolor='none', boxstyle='round,pad=0.5'), ha='center')

right_x = 8.1

steps_right = [
    ("All Cells", 6.0),
    ("LOO-CV Split First\n(Hold out 1 patient)", 4.8),
    ("[Training Fold Only]\nHVG Selection", 3.6),
    ("[Training Fold Only]\nScaling", 2.4),
    ("Train on Train Fold\nPredict on Test Patient", 1.2)
]

for text, y in steps_right:
    # Draw box
    ax1.add_patch(patches.FancyBboxPatch((right_x, y), box_width, box_height, boxstyle='round,pad=0.1', facecolor='#f4fbf4', edgecolor='#5cb85c', lw=1.5))
    # Draw text
    ax1.text(right_x + box_width/2, y + box_height/2, text, ha='center', va='center', color='black', fontsize=10, weight='semibold')
    
    # Draw checkmark for corrected steps
    if y < 6.0:
        ax1.text(right_x + box_width - 0.35, y + box_height - 0.25, r'$\checkmark$', color='#5cb85c', fontsize=18)
        
    # Draw arrow down (if not the last step)
    if y > 1.2:
        ax1.annotate('', xy=(right_x + box_width/2, y - 0.4), xytext=(right_x + box_width/2, y),
                    arrowprops=dict(arrowstyle="->", color='#5cb85c', lw=1.5))

plt.tight_layout()
fig1.savefig(os.path.join(output_dir, 'figure1_pipeline_schematic.png'), dpi=300, bbox_inches='tight')
plt.close(fig1)
print("Figure 1 saved as figure1_pipeline_schematic.png")
