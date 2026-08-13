# Leakage-Free Evaluation of Single-Cell Response Classifiers

**License.** Code in this repository is licensed under the MIT License
(see `LICENSE`). The manuscript text and figures are licensed under
Creative Commons Attribution 4.0 International (CC BY 4.0).

Code and results backing the manuscript **"Leakage-Free Evaluation of Single-Cell
Classifiers Reveals Systematic Underpowering in Small Public Cohorts"**
(`Leakage_Evaluation_Preprint.md`, in this repository root).

**Summary.** Two single-cell classifiers, one for anti-TNF response prediction in
inflammatory bowel disease, and immune-checkpoint-inhibitor response prediction in
melanoma, are evaluated under standard nested cross-validation with no data
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
| Melanoma permutation p | 0.2056 (102/500) | same | same |
| Melanoma consensus counts | *TBC1D10B* 18/18, *TNRC6B* 17/18, *FOXP1* 10/18, *SELL* 2/18 | same | `melanoma_consensus_counts.pkl` |
| IBD consensus counts | *TNK1* 36, *IGSF8* 35, *KCNQ1* 31, *CYP4F12* 31, *SSTR1* 19, *SIGMAR1* 8 | `Code/run_ibd_fast.py` | `Code/ibd_correct_results.pkl` |
| GIMATS sensitivity | 90.9% (10 of 11) | `Code/run_gimats_correct.py` | `Code/gimats_correct_results.pkl` |
| GIMATS mean remission prob. | 0.1504 (range 0.00085–0.7505) | same | same |
| GIMATS cells | 32,458 total; smallest patient 831 | same | same |
| GSE91061 directional concordance | 5/6 (83.3%) | `Code/bulk_validation_threshold2.py` | `Results/bulk_validation_threshold2.json` |
| GSE91061 binomial p | 0.219, two-sided, not significant | `Code/bulk_validation_threshold2.py` | `Results/bulk_validation_threshold2.json` |

---

## 2. Figures

| Figure | File (in `latex/figs/`) | Produced by | Source data |
|---|---|---|---|
| 1: pipeline schematic | `fig1.png` | `Code/generate_manuscript_figures.py` | diagram, no data |
| 2: IBD feature stability | `fig2.png` | `Code/generate_manuscript_figures.py` | `ibd_correct_results.pkl['consensus_counts']`, filtered to ≥5 folds |
| 3: IBD ROC | `fig3.png` | `Code/generate_manuscript_figures.py` | `ibd_correct_results.pkl` (ROC recomputed from stored per-patient probabilities) |
| 4: IBD permutation null | `fig4.png` | `Code/generate_manuscript_figures.py` | the actual 500-permutation null from `ibd_correct_results.pkl` |
| 5: melanoma feature stability | `fig5.png` | `Code/generate_manuscript_figures.py` | `melanoma_consensus_counts.pkl` |
| 6: melanoma ROC | `fig6.png` | `Code/generate_manuscript_figures.py` | `melanoma_verified_results.pkl` (ROC recomputed from stored per-patient probabilities) |

Figure 6's ROC and the cited melanoma AUC are computed via two independent code
paths (`run_melanoma_fromraw.py` directly, and
`ml_pipeline.py` → `run_melanoma_cv.py` for the figure). Both agree to full
float precision (0.7407407407407408).

---

## 3. Data (not included in this repository)

Raw and derived data are excluded from version control (`.gitignore`). Obtain
from GEO and place as shown. Total footprint is roughly 50 GB, dominated by the
`.h5ad` files.

| Accession | File | Destination | Notes |
|---|---|---|---|
| GSE282122 | `taurus_lightweight.h5ad` | `Data/GSE282122/` | Preprocessed derivative, not the raw GEO download directly. See caveat below. |
| GSE134809 | `gimats_annotated.h5ad` | `Data/GSE134809/` | Annotated derivative of `GSE134809_RAW.tar` |
| GSE120575 | `GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz` | `Data/` | Raw TPM matrix (Sade-Feldman et al.) |
| GSE120575 | `GSE120575_patient_ID_single_cells.txt.gz` | `Data/` | Patient/response metadata (required) |
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

Working directory matters: these scripts do not all resolve paths the same way.
Some resolve inputs via `__file__` and run from anywhere ("any cwd"); others use
bare relative filenames and must be started from a specific directory. The
required cwd is noted per command; each `cd` below starts from the repository
root.

```
# IBD (slow, approx. 20 minutes; 500-permutation null)   [cwd: Code/, bare
# 'ibd_feature_cache.npz' and '../Data/...' resolve against Code/]
cd Code && python run_ibd_fast.py && cd ..

# GIMATS external validation (requires ibd_feature_cache.npz above)   [cwd: Code/]
cd Code && python run_gimats_correct.py && cd ..

# Melanoma, from raw TPM (approx. 15 minutes)   [cwd: Data/, reads the raw
# GSE120575 .txt.gz files, which live in Data/, by bare name]
cd Data && python ../Code/run_melanoma_fromraw.py && cd ..

# Melanoma ROC chain
#   filter_and_save_baseline_final.py  [cwd: Data/, bare raw GSE120575 inputs]
#   ml_pipeline.py                     [cwd: Data/, bare 'baseline_expression_clean.pkl']
#   run_melanoma_cv.py                 [any cwd, input via __file__; writes
#                                       'melanoma_roc_data.npz' to the cwd]
cd Data && python ../Code/filter_and_save_baseline_final.py && python ../Code/ml_pipeline.py && cd ..
python Code/run_melanoma_cv.py

# GSE91061 external validation
#   reproduce_directional_binomial.py  [cwd: Data/, bare raw GSE91061 inputs]
#   bulk_validation_threshold2.py      [any cwd, all paths via __file__/REPO_ROOT]
cd Data && python ../Code/reproduce_directional_binomial.py && cd ..
python Code/bulk_validation_threshold2.py

# Figures   [any cwd, all paths via __file__]
python Code/generate_manuscript_figures.py
```

Note on superseded scripts. `Code/validate_gse91061.py` and
`Code/binomial_directional_check.py` implement the earlier 13-gene bulk
validation, which was replaced by a stated selection rule (genes selected in
at least 2 of 18 folds and testable in GSE91061). They are retained for
historical reference and are not part of the current reproduction path.
`Code/generate_figures.py` and `Code/generate_figures_unified.py` are
likewise superseded by `Code/generate_manuscript_figures.py`.

See `requirements.txt` for the Python environment.

---

## 5. Notes and limitations

**Permutation p-value convention.** Permutation p-values in both case studies use
the exact estimator (b+1)/(m+1), where b is the number of permutations achieving
accuracy greater than or equal to the observed value and m is the total number of
permutations (Phipson & Smyth 2010). This estimator is conservative and cannot
return zero. The IBD p-value is (43+1)/(500+1) = 0.0878 and the melanoma p-value
is (102+1)/(500+1) = 0.2056.

**Pickle version fragility.** `Code/ibd_verified_results.pkl` cannot be
unpickled under Python 3.12 / pandas 2.x: it embeds a pandas `Categorical`
whose `__setstate__` raises `NotImplementedError` under current pandas.
`Code/ibd_correct_results.pkl` stores the same IBD results using only NumPy
objects and loads without error; it is the artifact the figure-generation and
GIMATS scripts read. If you encounter an unpickling error on a derived
artifact, regenerate it from raw rather than attempting a compatibility
workaround.

**Result artifacts kept on disk but not committed to git.** Three pickles are
gitignored because their evidentiary status is ambiguous, not because they are
unimportant:

- `Code/ibd_verified_results.pkl`: written by `run_ibd_verified.py`. Bit-identical
  to `ibd_correct_results.pkl`, which is the script the manuscript's IBD numbers
  are attributed to. Both are kept on disk; only the attributed one is committed.
- `Code/gimats_multimodal_results.pkl`: written by `run_gimats_multimodal.py`.
  Its sensitivity matches `gimats_correct_results.pkl` but its probability values
  differ, and no manuscript number is attributed to it.
- `melanoma_verified_results.pkl`: this filename was historically written by two
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
Data/       raw and derived data (not included, see §3)
Results/    generated figures
Leakage_Evaluation_Preprint.md    the manuscript
```

Committed result artifacts (`Code/ibd_correct_results.pkl`,
`Code/gimats_correct_results.pkl`, `melanoma_consensus_counts.pkl`) are the
saved outputs each cited number in the manuscript is derived from.
