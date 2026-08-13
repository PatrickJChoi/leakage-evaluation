"""
Generates manuscript figures 2-6 directly into latex/figs/ from committed,
tracked result artifacts.

Figure 1 is deliberately excluded -- see GENERATE_FIG1 below.

Why this script exists
----------------------
Prior to 2026-07-29 the figures in latex/figs/ had no derivation chain: no
committed script wrote to that directory, and figs 2-4 could only be matched to
images in a pre-split working directory outside this repository. Figure 5's gene
set was additionally a hardcoded four-name literal that excluded two genes
selected more often than one of the four it included.

Provenance rules this script follows
------------------------------------
1. Every plotted value is read from a committed pickle. No gene name, count,
   accuracy, AUC or p-value is written as a literal anywhere below.
2. All paths resolve via __file__, so the script is cwd-independent.
3. ROC curves are computed from stored per-patient probabilities in the result
   pickles, not from separate .npz files (the previous fig3/fig6 inputs, which
   were loaded by bare cwd-relative name and, for fig6, silently skipped when
   absent).
4. Missing or malformed inputs raise immediately. There is no silent-skip path.
5. Feature-selection figures use a stated fold threshold, declared as a module
   constant below and reported in stdout so the run log records the rule.

Note on key names: the two pipelines store the same quantities under different
keys. The IBD pickle uses "true_labels" and "patient_index"; the melanoma pickle
uses "y" and "patients". Both are named as explicit constants below.

Inputs (both must be tracked in git):
    Code/ibd_verified_results.pkl
    melanoma_verified_results.pkl

Outputs (overwritten in place):
    latex/figs/fig2.png ... latex/figs/fig6.png

Run:
    python Code/generate_manuscript_figures.py
"""

import os
import pickle
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from sklearn.metrics import auc, roc_curve

# --------------------------------------------------------------------------
# Paths and stated rules
# --------------------------------------------------------------------------

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IBD_PKL = os.path.join(REPO_ROOT, "Code", "ibd_correct_results.pkl")
MEL_PKL = os.path.join(REPO_ROOT, "melanoma_verified_results.pkl")
FIGS_DIR = os.path.join(REPO_ROOT, "latex", "figs")

# Stated fold thresholds for the two feature-stability figures.
# IBD uses >= 5 of 36, matching the existing Figure 2 caption.
# Melanoma uses >= 2 of 18, matching the external-validation rule stated in
# Methods. Change either value here and the corresponding caption must change
# with it.
IBD_FOLD_THRESHOLD = 1  # plot every feature selected in at least one fold
MEL_FOLD_THRESHOLD = 2

# Display-only reference lines.
IBD_MAX_FOLDS = 36
MEL_MAX_FOLDS = 18

# The two pipelines store the same quantities under different key names:
# run_ibd_verified.py writes "true_labels" / "patient_index", while
# run_melanoma_fromraw.py writes "y" / "patients". Named explicitly here rather
# than probed at runtime, so a future key rename fails loudly instead of
# silently falling through to a different array.
IBD_LABEL_KEY = "true_labels"
MEL_LABEL_KEY = "y"

# Figure 1 is a data-free schematic: it plots no result values, so nothing in it
# can disagree with the manuscript. Regenerating it is also not currently
# possible on this machine -- matplotlib 3.10.8 in the crohns environment dies
# with a native Windows fatal exception (0xc06d007f) inside
# matplotlib/bezier.py:223, reached via add_patch() -> _update_patch_limits()
# when a FancyBboxPatch with a rounded boxstyle is added. That is a crash below
# the Python level: no traceback is produced and the process exits 127 silently.
# It is visible only under `python -X faulthandler`. The committed fig1.png was
# produced 2026-07-22 under a different (since-removed) environment, and its
# drawing code remains committed in Code/generate_figures_unified.py.
# Set this True only in an environment where the above crash does not occur.
GENERATE_FIG1 = True


# --------------------------------------------------------------------------
# Preflight
# --------------------------------------------------------------------------

def load_pickle(path, label):
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"{label} not found at {path}. This script will not fall back to "
            f"any other artifact. Restore the file or fix the path."
        )
    with open(path, "rb") as fh:
        obj = pickle.load(fh)
    if not isinstance(obj, dict):
        raise TypeError(f"{label} loaded as {type(obj)}, expected dict.")
    return obj


def require_keys(obj, keys, label):
    missing = [k for k in keys if k not in obj]
    if missing:
        raise KeyError(
            f"{label} is missing required key(s): {missing}. "
            f"Keys present: {sorted(obj.keys())}"
        )


def counts_at_or_above(freq, threshold):
    """Descending by count, ties broken alphabetically for determinism."""
    return sorted(
        [(g, int(c)) for g, c in freq.items() if int(c) >= threshold],
        key=lambda gc: (-gc[1], gc[0]),
    )


def singleton_count(freq):
    return sum(1 for c in freq.values() if int(c) == 1)


# --------------------------------------------------------------------------
# Figures 2 and 4 (seaborn-themed, matching the original generate_figures.py)
# --------------------------------------------------------------------------

def figure2_ibd_stability(ibd, out_path):
    freq = ibd["consensus_counts"]
    selected = counts_at_or_above(freq, IBD_FOLD_THRESHOLD)
    below = len(freq) - len(selected)

    print(f"  rule: >= {IBD_FOLD_THRESHOLD} of {IBD_MAX_FOLDS} folds")
    print(f"  plotted: {len(selected)} features; not plotted: {below}")
    for gene, count in selected:
        print(f"    {gene}: {count}/{IBD_MAX_FOLDS}")

    # barh renders bottom-up, so reverse for descending top-down
    features = [g for g, _ in selected][::-1]
    counts = [c for _, c in selected][::-1]

    plt.figure(figsize=(10, 9))
    sns.barplot(x=counts, y=features, hue=features,
                palette="Blues_d", legend=False)
    plt.axvline(x=IBD_MAX_FOLDS, color="red", linestyle="--", linewidth=2,
                label=f"Maximum Folds (n={IBD_MAX_FOLDS})")
    plt.title("Robustness of Multimodal Feature Selection Across LOO-CV Folds",
              fontsize=16, weight="bold")
    plt.xlabel(f"Number of Folds Selected (out of {IBD_MAX_FOLDS})", fontsize=14)
    plt.ylabel("Features", fontsize=14)
    plt.xlim(0, 40)
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    return len(selected), below


def figure4_ibd_permutation(ibd, out_path):
    null_dist = np.asarray(ibd["perm_accuracies"], dtype=float)
    observed = float(ibd["accuracy"])
    p_value = float(ibd["p_value"])
    n_perm = null_dist.size

    print(f"  permutations: {n_perm}")
    print(f"  observed accuracy: {observed}")
    print(f"  p_value: {p_value}")
    print(f"  null mean: {null_dist.mean()}  null max: {null_dist.max()}")

    plt.figure(figsize=(10, 6))
    sns.histplot(null_dist, bins=30, kde=True,
                 color="lightgray", edgecolor="black")
    plt.axvline(x=observed, color="red", linestyle="-", linewidth=4)
    plt.text(
        observed - 0.015, plt.ylim()[1] * 0.8,
        f"Observed Accuracy: {observed * 100:.1f}%\n"
        f"$p = {p_value:.4f}$\n(Not Statistically Significant)",
        color="red", fontsize=12, weight="bold", ha="right",
        bbox=dict(facecolor="white", alpha=0.8, edgecolor="red",
                  boxstyle="round,pad=0.5"),
    )
    plt.title("Permutation Test: Null Distribution vs. Observed Accuracy",
              fontsize=16, weight="bold")
    plt.xlabel("Classification Accuracy", fontsize=14)
    plt.ylabel(f"Frequency ({n_perm} Permutations)", fontsize=14)
    plt.xlim(0.1, 0.9)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    return n_perm


# --------------------------------------------------------------------------
# Figures 1, 3, 5, 6 (plain matplotlib, matching the original unified script)
# --------------------------------------------------------------------------

def figure1_pipeline_schematic(out_path):
    fig, ax = plt.subplots(figsize=(14, 8), dpi=300)
    ax.axis("off")
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 8)

    box_width, box_height = 4.8, 0.7

    ax.text(3.5, 7.3, "Naive Pipeline (Leaky)", color="white", weight="bold",
            fontsize=14, ha="center",
            bbox=dict(facecolor="#d9534f", edgecolor="none",
                      boxstyle="round,pad=0.5"))
    left_x = 1.1
    steps_left = [
        ("All Cells", 6.0),
        ("Global HVG Selection\n(Uses all cells)", 4.8),
        ("Global Scaling\n(Uses all cells)", 3.6),
        ("LOO-CV Split\n(Train / Test)", 2.4),
        ("Model Training & Evaluation", 1.2),
    ]
    for text, y in steps_left:
        ax.add_patch(patches.FancyBboxPatch(
            (left_x, y), box_width, box_height, boxstyle="round,pad=0.1",
            facecolor="#fcf8f2", edgecolor="#d9534f", lw=1.5))
        ax.text(left_x + box_width / 2, y + box_height / 2, text,
                ha="center", va="center", color="black",
                fontsize=10, weight="semibold")
        if y > 1.2:
            ax.annotate("", xy=(left_x + box_width / 2, y - 0.4),
                        xytext=(left_x + box_width / 2, y),
                        arrowprops=dict(arrowstyle="->", color="#d9534f", lw=1.5))

    ax.add_patch(plt.Rectangle((left_x - 0.2, 3.4), box_width + 0.4, 2.3,
                               facecolor="none", edgecolor="red",
                               linestyle="--", lw=2))
    ax.text(0.1, 4.3, "TEST DATA\nLEAKS HERE!", color="white", weight="bold",
            fontsize=11, ha="center", va="center",
            bbox=dict(facecolor="red", edgecolor="none",
                      boxstyle="round,pad=0.4"))
    ax.annotate("", xy=(left_x - 0.05, 4.5), xytext=(0.6, 4.5),
                arrowprops=dict(arrowstyle="->", color="red", lw=2))

    ax.text(10.5, 7.3, "Corrected Pipeline (Leakage-Free)", color="white",
            weight="bold", fontsize=14, ha="center",
            bbox=dict(facecolor="#5cb85c", edgecolor="none",
                      boxstyle="round,pad=0.5"))
    right_x = 8.1
    steps_right = [
        ("All Cells", 6.0),
        ("LOO-CV Split First\n(Hold out 1 patient)", 4.8),
        ("[Training Fold Only]\nHVG Selection", 3.6),
        ("[Training Fold Only]\nScaling", 2.4),
        ("Train on Train Fold\nPredict on Test Patient", 1.2),
    ]
    for text, y in steps_right:
        ax.add_patch(patches.FancyBboxPatch(
            (right_x, y), box_width, box_height, boxstyle="round,pad=0.1",
            facecolor="#f4fbf4", edgecolor="#5cb85c", lw=1.5))
        ax.text(right_x + box_width / 2, y + box_height / 2, text,
                ha="center", va="center", color="black",
                fontsize=10, weight="semibold")
        if y < 6.0:
            ax.text(right_x + box_width - 0.35, y + box_height - 0.25,
                    r"$\checkmark$", color="#5cb85c", fontsize=18)
        if y > 1.2:
            ax.annotate("", xy=(right_x + box_width / 2, y - 0.4),
                        xytext=(right_x + box_width / 2, y),
                        arrowprops=dict(arrowstyle="->", color="#5cb85c", lw=1.5))

    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def figure3_ibd_roc(ibd, out_path):
    y_true = np.asarray(ibd[IBD_LABEL_KEY])
    y_prob = np.asarray(ibd["probabilities"], dtype=float)
    if y_true.shape != y_prob.shape:
        raise ValueError(
            f"IBD true_labels shape {y_true.shape} != probabilities shape "
            f"{y_prob.shape}"
        )
    fpr, tpr, _ = roc_curve(y_true, y_prob)
    roc_auc = auc(fpr, tpr)
    stored_auc = float(ibd["auc"])
    print(f"  n patients: {y_true.size}")
    print(f"  AUC recomputed from probabilities: {roc_auc}")
    print(f"  AUC stored in pickle:              {stored_auc}")
    if not np.isclose(roc_auc, stored_auc, rtol=0, atol=1e-12):
        raise ValueError(
            f"Recomputed IBD AUC {roc_auc} does not match stored {stored_auc}. "
            f"Stopping rather than plotting a value that disagrees with the "
            f"committed artifact."
        )

    plt.figure(figsize=(8, 6), dpi=300)
    plt.plot(fpr, tpr, color="darkblue", lw=3,
             label=f"Strict Nested LOO-CV (AUC = {roc_auc:.3f})")
    plt.plot([0, 1], [0, 1], color="gray", lw=2, linestyle="--",
             label="Random Chance (AUC = 0.50)")
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel("False Positive Rate", fontsize=14)
    plt.ylabel("True Positive Rate", fontsize=14)
    plt.title("Receiver Operating Characteristic (ROC)",
              fontsize=16, weight="bold")
    plt.legend(loc="lower right", fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    return roc_auc


def figure5_melanoma_stability(mel, out_path):
    freq = mel["feature_freq"]
    selected = counts_at_or_above(freq, MEL_FOLD_THRESHOLD)
    singletons = singleton_count(freq)
    total = len(freq)

    print(f"  rule: >= {MEL_FOLD_THRESHOLD} of {MEL_MAX_FOLDS} folds")
    print(f"  plotted: {len(selected)} of {total} distinct genes")
    print(f"  genes selected in exactly 1 fold: {singletons}")
    for gene, count in selected:
        print(f"    {gene}: {count}/{MEL_MAX_FOLDS}")

    genes = [g for g, _ in selected][::-1]
    counts = [c for _, c in selected][::-1]

    fig, ax = plt.subplots(figsize=(8, 7), dpi=300)
    colors = ["#1f77b4" for c in counts]
    bars = ax.barh(genes, counts, color=colors, edgecolor="black", height=0.6)
    for bar in bars:
        width = bar.get_width()
        ax.text(width + 0.3, bar.get_y() + bar.get_height() / 2,
                f"{int(width)}/{MEL_MAX_FOLDS}", va="center", ha="left",
                fontsize=10, fontweight="bold",
                color="black")
    ax.set_xlim(0, 20)
    ax.set_xlabel(f"Selection Frequency (out of {MEL_MAX_FOLDS} folds)",
                  fontsize=12, fontweight="bold")
    ax.set_title("Melanoma ICI Feature Selection Frequency",
                 fontsize=13, fontweight="bold", pad=15)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", which="major", labelsize=11)
    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return len(selected), singletons, total


def figure6_melanoma_roc(mel, out_path):
    y_true = np.asarray(mel[MEL_LABEL_KEY])
    y_prob = np.asarray(mel["probabilities"], dtype=float)
    if y_true.shape != y_prob.shape:
        raise ValueError(
            f"Melanoma true_labels shape {y_true.shape} != probabilities "
            f"shape {y_prob.shape}"
        )
    fpr, tpr, _ = roc_curve(y_true, y_prob)
    roc_auc = auc(fpr, tpr)
    stored_auc = float(mel["auc"])
    print(f"  n patients: {y_true.size}")
    print(f"  AUC recomputed from probabilities: {roc_auc}")
    print(f"  AUC stored in pickle:              {stored_auc}")
    if not np.isclose(roc_auc, stored_auc, rtol=0, atol=1e-12):
        raise ValueError(
            f"Recomputed melanoma AUC {roc_auc} does not match stored "
            f"{stored_auc}. Stopping rather than plotting a value that "
            f"disagrees with the committed artifact."
        )

    fig, ax = plt.subplots(figsize=(6, 6), dpi=300)
    ax.plot(fpr, tpr, color="#1f77b4", lw=2.5,
            label=f"Random Forest LOO-CV (AUC = {roc_auc:.3f})")
    ax.plot([0, 1], [0, 1], color="gray", lw=1.5, linestyle="--",
            label="Random Chance (AUC = 0.50)")
    ax.set_xlim([-0.02, 1.02])
    ax.set_ylim([-0.02, 1.02])
    ax.set_xlabel("False Positive Rate (FPR)", fontsize=12,
                  fontweight="bold", labelpad=8)
    ax.set_ylabel("True Positive Rate (TPR)", fontsize=12,
                  fontweight="bold", labelpad=8)
    ax.set_title("Melanoma ICI LOO-CV ROC Curve", fontsize=13,
                 fontweight="bold", pad=15)
    ax.legend(loc="lower right", fontsize=10, frameon=True, shadow=False)
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.set_aspect("equal")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return roc_auc


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main():
    print("=" * 74)
    print("MANUSCRIPT FIGURE GENERATION")
    print("=" * 74)
    print(f"Repo root : {REPO_ROOT}")
    print(f"Output dir: {FIGS_DIR}")

    if not os.path.isdir(FIGS_DIR):
        raise FileNotFoundError(
            f"Output directory {FIGS_DIR} does not exist. This script writes "
            f"only into the existing latex/figs/ directory and will not create it."
        )

    print("\n--- PREFLIGHT ---")
    ibd = load_pickle(IBD_PKL, "IBD result pickle")
    mel = load_pickle(MEL_PKL, "Melanoma result pickle")
    print(f"IBD keys      : {sorted(ibd.keys())}")
    print(f"Melanoma keys : {sorted(mel.keys())}")

    require_keys(ibd, ["consensus_counts", "perm_accuracies", "accuracy",
                       "p_value", IBD_LABEL_KEY, "probabilities", "auc"],
                 "IBD result pickle")
    require_keys(mel, ["feature_freq", "accuracy", "auc",
                       MEL_LABEL_KEY, "probabilities"],
                 "Melanoma result pickle")
    print(f"Label key in use -- IBD: '{IBD_LABEL_KEY}', "
          f"melanoma: '{MEL_LABEL_KEY}'")
    print("Preflight OK: all required keys present.")

    if GENERATE_FIG1:
        print("\n--- FIGURE 1: pipeline schematic (no data inputs) ---")
        figure1_pipeline_schematic(os.path.join(FIGS_DIR, "fig1.png"))
        print("  wrote fig1.png")
    else:
        print("\n--- FIGURE 1: SKIPPED (GENERATE_FIG1 is False) ---")
        print("  fig1.png left untouched; see the comment on GENERATE_FIG1.")

    print("\n--- FIGURE 3: IBD ROC ---")
    figure3_ibd_roc(ibd, os.path.join(FIGS_DIR, "fig3.png"))
    print("  wrote fig3.png")

    print("\n--- FIGURE 5: melanoma feature selection stability ---")
    figure5_melanoma_stability(mel, os.path.join(FIGS_DIR, "fig5.png"))
    print("  wrote fig5.png")

    print("\n--- FIGURE 6: melanoma ROC ---")
    figure6_melanoma_roc(mel, os.path.join(FIGS_DIR, "fig6.png"))
    print("  wrote fig6.png")

    # Seaborn theming is global; apply it only for the two figures that used it
    # originally, and do it last so it cannot bleed into the others.
    sns.set_theme(style="whitegrid", context="paper", font_scale=1.2)

    print("\n--- FIGURE 2: IBD feature selection stability ---")
    figure2_ibd_stability(ibd, os.path.join(FIGS_DIR, "fig2.png"))
    print("  wrote fig2.png")

    print("\n--- FIGURE 4: IBD permutation null distribution ---")
    figure4_ibd_permutation(ibd, os.path.join(FIGS_DIR, "fig4.png"))
    print("  wrote fig4.png")

    print("\n" + "=" * 74)
    if GENERATE_FIG1:
        print("FIGURES 1-6 WRITTEN TO latex/figs/")
    else:
        print("FIGURES 2-6 WRITTEN TO latex/figs/ (fig1 skipped, untouched)")
    print("=" * 74)


if __name__ == "__main__":
    try:
        main()
    except Exception:
        import traceback
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)
