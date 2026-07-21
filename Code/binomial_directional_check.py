"""
GSE91061 directional-concordance binomial check.

Purpose: provide a VERSIONED, reproducible source for the directional-consistency
p-value reported in the GSE91061 external-validation arm of Leakage_Evaluation_Preprint.md.

Discipline note: this script does NOT target or confirm any manuscript value.
It computes BOTH the one-sided and two-sided binomial p-values from the observed
directional counts and reports whatever it gets. Which value (if any) is reported
in the manuscript is a human decision that must be justified by an a-priori
directional hypothesis, not by which number is smaller.

Inputs (observed, from gse91061_validation_results.csv direction tally):
  n_genes  = 13 total protein-coding genes tested
  n_agree  = 10 genes trending in the predicted direction (higher in responders)
  p_null   = 0.5 (chance direction under the null)
"""

from scipy.stats import binomtest

n_genes = 13
n_agree = 10
p_null  = 0.5

# One-sided: P(>= n_agree successes in the predicted direction), requires a
# genuine a-priori directional hypothesis to be legitimate.
one_sided = binomtest(n_agree, n_genes, p_null, alternative='greater')

# Two-sided: P(a split at least this lopsided in EITHER direction), the default
# when there is no pre-registered directional hypothesis.
two_sided = binomtest(n_agree, n_genes, p_null, alternative='two-sided')

print(f"Observed: {n_agree} of {n_genes} genes trend the predicted direction "
      f"({n_agree/n_genes:.1%})")
print(f"One-sided binomial p (alternative='greater'):   {one_sided.pvalue:.6f}")
print(f"Two-sided binomial p (alternative='two-sided'): {two_sided.pvalue:.6f}")
print()
print("Interpretation is a human call:")
print(" - Report one-sided ONLY with a stated a-priori directional hypothesis,")
print("   and disclose the two-sided value alongside it.")
print(" - Otherwise report the two-sided value, or report the 10/13 trend")
print("   descriptively with no p-value.")

# Save a machine-readable record so the number traces to a re-runnable source.
import json
out = {
    "n_genes": n_genes,
    "n_agree": n_agree,
    "p_null": p_null,
    "one_sided_greater_p": one_sided.pvalue,
    "two_sided_p": two_sided.pvalue,
    "note": "Choice of test is a human decision requiring a-priori justification; "
            "this script reports both and targets neither.",
}
with open("binomial_directional_check_results.json", "w") as f:
    json.dump(out, f, indent=2)
print("\nSaved: binomial_directional_check_results.json")
