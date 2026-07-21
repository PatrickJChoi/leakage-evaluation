# Leakage-Free Evaluation of Single-Cell Response Classifiers

Code and results backing the manuscript **"Leakage-Free Evaluation of Single-Cell
Classifiers Reveals Systematic Underpowering in Current Public Datasets"**
(`Leakage_Evaluation_Preprint.md`, in this repository root).

**Summary.** Two single-cell classifiers — anti-TNF response prediction in
inflammatory bowel disease, and immune-checkpoint-inhibitor response prediction in
melanoma — are evaluated under standard nested cross-validation with no data
leakage. Both collapse to near-chance performance. The paper argues that
leakage-free evaluation should be standard practice for this class of model, and
that some published results in this space may be inflated by evaluation artifacts
rather than genuine biological signal.

---

## 1. Headline results and where each number comes from

Every number below is reproduced from raw data by the named script in this
repository. Nothing is hardcoded.

| Claim | Value | Script | Saved output |
|---|---|---|---|
| IBD accuracy | 72.22% | `Code/run_ibd_fast.py` | `Code/ibd_correct_results.pkl` |
| IBD AUC | 0.7656 | same | same |
| IBD permutation p | 0.0878 (43/500) | same | same |
| Melanoma accuracy | 72.22% (13/18) | `Code/run_melanoma_fromraw.py` | `melanoma_verified_results.pkl` |
| Melanoma AUC | 0.7407 | same | same |
| Melanoma permutation p | 0.2040 (102/500) | same | same |
| Melanoma consensus counts | TBC1D10B 18/18, TNRC6B 17/18, FOXP1 10/18, SELL 2/18 | same | `melanoma_consensus_counts.pkl` |
| IBD consensus counts | TNK1 36, IGSF8 35, KCNQ1 31, CYP4F12 31, SSTR1 19, SIGMAR1 8 | `Code/run_ibd_fast.py` | `Code/ibd_correct_results.pkl` |
| GIMATS sensitivity | 90.9% (10 of 11) | `Code/run_gimats_correct.py` | `Code/gimats_correct_results.pkl` |
| GIMATS mean remission prob. | 0.1504 (range 0.00085–0.7505) | same | same |
| GIMATS cells | 32,458 total; smallest patient 831 | same | same |
| GSE91061 directional concordance | 10/13 (76.9%) | `Code/validate_gse91061.py` | `gse91061_validation_results.csv` |
| GSE91061 binomial p | 0.092, two-sided, not significant | `Code/binomial_directional_check.py` | — |

---

## 2. Figures

| Figure | File (in `Results/`) | Produced by | Source data |
|---|---|---|---|
| 1 — pipeline schematic | `figure1_pipeline_schematic.png` | `Code/generate_figures_unified.py` | diagram, no data |
| 2 — IBD feature stability | `figure2_stability.png` | `Code/generate_figures.py` | `ibd_correct_results.pkl['consensus_counts']`, filtered to ≥5 folds |
| 3 — IBD ROC | `figure3_roc.png` | `Code/generate_figure3.py` | `roc_data.npz` |
| 4 — IBD permutation null | `figure4_permutation.png` | `Code/generate_figures.py` | the actual 500-permutation null from `ibd_correct_results.pkl` |
| 5 — melanoma feature stability | `figure5_stability.png` | `Code/generate_figures_unified.py` | `melanoma_consensus_counts.pkl` |
| 6 — melanoma ROC | `figure6_melanoma_roc.png` | `Code/generate_figures_unified.py` | `melanoma_roc_data.npz` |

Figure 6's ROC and the cited melanoma AUC are computed via two independent code
paths (`run_melanoma_fromraw.py` directly, and
`ml_pipeline.py` → `run_melanoma_cv.py` for the figure). Both agree to full
float precision (0.7407407407407408).

---

## 3. Data — not included in this repository

Raw and derived data are excluded from version control (`.gitignore`). Obtain
from GEO and place as shown. Total footprint is roughly 50 GB, dominated by the
`.h5ad` files.

| Accession | File | Destination | Notes |
|---|---|---|---|
| GSE282122 | `taurus_lightweight.h5ad` | `Data/GSE282122/` | Preprocessed derivative, not the raw GEO download directly. See caveat below. |
| GSE134809 | `gimats_annotated.h5ad` | `Data/GSE134809/` | Annotated derivative of `GSE134809_RAW.tar` |
| GSE120575 | `GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz` | `Data/` | Raw TPM matrix (Sade-Feldman et al.) |
| GSE120575 | `GSE120575_patient_ID_single_cells.txt.gz` | `Data/` | Patient/response metadata — required |
| GSE91061 | `GSE91061_series_matrix.txt.gz` | `Data/` | Response metadata |
| GSE91061 | rld and FPKM expression tables | `Data/` | Riaz et al. 2017 external validation cohort |

> **Note on derived inputs.** `taurus_lightweight.h5ad` and `gimats_annotated.h5ad`
> are preprocessed derivatives, not files downloadable directly from GEO. The
> upstream preprocessing script that produced them is not currently included in
> this repository. All results are reproducible starting from these files, but a
> reader beginning from the raw GEO archives alone would need to reconstruct that
> preprocessing step. Contact the author for the preprocessing script if needed.

`Code/run_gimats_correct.py` additionally requires `ibd_feature_cache.npz`,
produced by `Code/run_ibd_fast.py`. Run the IBD pipeline first.

---

## 4. Running the analyses

Run from the repository root.

```
# IBD (slow — approx. 20 minutes; 500-permutation null)
python Code/run_ibd_fast.py

# GIMATS external validation (requires ibd_feature_cache.npz above)
python Code/run_gimats_correct.py

# Melanoma, from raw TPM (approx. 15 minutes)
python Code/run_melanoma_fromraw.py

# Melanoma ROC chain
python Code/filter_and_save_baseline_final.py
python Code/ml_pipeline.py
python Code/run_melanoma_cv.py

# GSE91061 external validation
python Code/validate_gse91061.py
python Code/binomial_directional_check.py

# Figures
python Code/generate_figures_unified.py
python Code/generate_figures.py
```

See `requirements.txt` for the Python environment.

---

## 5. Notes and limitations

**Permutation p-value convention.** The melanoma p-value is computed as plain
*k/N*; the IBD p-value uses *(k+1)/(N+1)*. Both are defensible conventions in the
literature; the inconsistency between the two analyses is noted here for
transparency.

**Pandas version fragility.** The derived pickles `baseline_expression_clean.pkl`
and `pseudobulk_matrix.pkl` were originally written under an older pandas version
and could not be unpickled under Python 3.12 / pandas 2.x (StringDtype
NotImplementedError). They were rebuilt from raw GSE120575 rather than shimmed.
If you encounter this error, regenerate from raw via
`Code/filter_and_save_baseline_final.py` rather than attempting a compatibility
workaround.

**Result artifacts kept on disk but not committed to git.** Three pickles are
gitignored because their evidentiary status is ambiguous, not because they are
unimportant:

- `Code/ibd_verified_results.pkl` — written by `run_ibd_verified.py`. Bit-identical
  to `ibd_correct_results.pkl`, which is the script the manuscript's IBD numbers
  are attributed to. Both are kept on disk; only the attributed one is committed.
- `Code/gimats_multimodal_results.pkl` — written by `run_gimats_multimodal.py`.
  Its sensitivity matches `gimats_correct_results.pkl` but its probability values
  differ, and no manuscript number is attributed to it.
- `melanoma_verified_results.pkl` — this filename was historically written by two
  different scripts. `run_melanoma_fromraw.py` is the valid from-raw source;
  `run_melanoma_verified.py` has been updated to write a distinctly named output
  to prevent future collision.

**Scope.** This analysis covers two datasets (IBD and melanoma) and one external
validation cohort (GSE91061). It is not a comprehensive survey of the single-cell
response-prediction literature.

---

## 6. Repository layout

```
Code/       analysis and figure scripts
Data/       raw and derived data (not included — see §3)
Results/    generated figures
Leakage_Evaluation_Preprint.md    the manuscript
```

Committed result artifacts — `Code/ibd_correct_results.pkl`,
`Code/gimats_correct_results.pkl`, `melanoma_consensus_counts.pkl` — are the
saved outputs each cited number in the manuscript is derived from.
