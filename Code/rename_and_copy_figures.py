import os
import shutil

def main():
    root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    results_dir = os.path.join(root_dir, "Results")

    # Mapping for existing single-cell leakage-evaluation figures in root
    # NOTE: the original "figure5_melanoma_feature_stability.png" -> "path1_figure5.png"
    # entry was REMOVED here — that source file does not exist anywhere on disk
    # (confirmed by direct search). The only real source for path1_figure5.png is
    # Results/path1_figure5_stability.png, handled by path1_results below.
    path1_existing = {
        "figure1_pipeline_schematic.png": "path1_figure1.png",
        "figure6_melanoma_roc.png": "path1_figure6.png"
    }

    # Mapping for single-cell leakage-evaluation figures still in Results folder
    path1_results = {
        "path1_figure2_stability.png": "path1_figure2.png",
        "path1_figure3_roc.png": "path1_figure3.png",
        "path1_figure4_permutation.png": "path1_figure4.png",
        "path1_figure5_stability.png": "path1_figure5.png"
    }

    print("\nRenaming leakage-evaluation figures in root...")
    for old_name, new_name in path1_existing.items():
        old_path = os.path.join(root_dir, old_name)
        new_path = os.path.join(root_dir, new_name)
        if os.path.exists(old_path):
            shutil.move(old_path, new_path)
            print(f"  Moved {old_name} -> {new_name}")
        elif os.path.exists(new_path):
            print(f"  {new_name} already exists.")
        else:
            print(f"  Warning: {old_name} not found.")

    print("\nCopying leakage-evaluation figures from Results to root...")
    for old_name, new_name in path1_results.items():
        src_path = os.path.join(results_dir, old_name)
        dst_path = os.path.join(root_dir, new_name)
        if not os.path.exists(src_path):
            print(f"  Warning: Results/{old_name} not found. We will need to generate it.")
            continue
        if os.path.exists(dst_path):
            print(f"  Skipped: {new_name} already exists at destination (not overwriting).")
            continue
        shutil.copy(src_path, dst_path)
        print(f"  Copied Results/{old_name} -> {new_name}")

    # Verify final set of figures in root
    print("\nVerification check:")
    all_path1 = [f"path1_figure{i}.png" for i in range(1, 7)]

    print("Leakage-evaluation figures status:")
    for f in all_path1:
        exists = os.path.exists(os.path.join(root_dir, f))
        print(f"  {f:<18}: {'FOUND' if exists else 'MISSING'}")

if __name__ == "__main__":
    main()
