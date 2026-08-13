"""
Sample size / power calculation for detecting a single ROC AUC against the
null AUC = 0.5, using the Hanley & McNeil (1982) variance formula for a
single AUC estimate.

Hanley & McNeil (1982), "The Meaning and Use of the Area under a Receiver
Operating Characteristic (ROC) Curve," Radiology 143:29-36.

Variance formula (single AUC, positives n1, negatives n2):

    Q1 = A / (2 - A)
    Q2 = 2*A^2 / (1 + A)
    Var(A) = [A(1-A) + (n1-1)(Q1 - A^2) + (n2-1)(Q2 - A^2)] / (n1 * n2)

Sample size condition (two-sided alpha, target power), balanced design
n1 = n2 = N/2:

    (A - 0.5) >= z_(1-alpha/2) * SE0 + z_(1-beta) * SE1

where SE0 = sqrt(Var(A)) evaluated at A = 0.5 (null variance), and
SE1 = sqrt(Var(A)) evaluated at the alternative A (variance under the
effect being detected).

N is iterated upward from 4 in steps of 2, so that n1 = n2 = N/2 is always
an integer (a balanced design cannot be evaluated at odd N under this
n1 = n2 constraint). The smallest such N for which the condition holds is
reported, along with the corresponding achieved power.

Achieved power at a given N is the exact complement:

    power(N) = Phi( (A - 0.5 - z_(1-alpha/2) * SE0) / SE1 )

which is algebraically equivalent to the sample-size condition above
(power(N) >= 0.80 iff the condition holds), verified directly rather than
assumed.

No values are seeded or randomized; this is a fully deterministic
closed-form calculation.
"""

import sys
import platform
import json
from pathlib import Path

try:
    import numpy
except ImportError:
    numpy = None

try:
    import scipy
    from scipy.stats import norm
except ImportError:
    scipy = None
    norm = None

OUTPUT_PATH = Path(__file__).resolve().parent / "auc_power_output.txt"

# Committed machine-readable artifact, resolved via __file__ (cwd-independent),
# with Results/ as a sibling of Code/ -- same layout convention as
# generate_manuscript_figures.py.
REPO_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = REPO_ROOT / "Results"
JSON_PATH = RESULTS_DIR / "auc_power_calculation.json"

ALPHA = 0.05        # two-sided
TARGET_POWER = 0.80
AUC_VALUES = [0.75, 0.766, 0.741, 0.80]

N_START = 4
N_STEP = 2           # keeps n1 = n2 = N/2 an integer at every step
N_MAX = 1_000_000    # safety bound against a non-terminating search


def hanley_mcneil_variance(a, n1, n2):
    """Hanley & McNeil (1982) variance of a single AUC estimate A."""
    q1 = a / (2 - a)
    q2 = (2 * a ** 2) / (1 + a)
    numerator = a * (1 - a) + (n1 - 1) * (q1 - a ** 2) + (n2 - 1) * (q2 - a ** 2)
    return numerator / (n1 * n2)


def achieved_power(a, n_total, z_alpha, alpha_var_a=0.5):
    """Exact power at balanced total sample size n_total for detecting AUC a
    against the null alpha_var_a (0.5), given z_alpha = z_(1-alpha/2)."""
    n1 = n2 = n_total // 2
    se0 = hanley_mcneil_variance(alpha_var_a, n1, n2) ** 0.5
    se1 = hanley_mcneil_variance(a, n1, n2) ** 0.5
    z = (a - 0.5 - z_alpha * se0) / se1
    return norm.cdf(z), se0, se1, n1, n2


def find_required_n(a, z_alpha, z_beta, target_power):
    n_total = N_START
    while n_total <= N_MAX:
        n1 = n2 = n_total // 2
        se0 = hanley_mcneil_variance(0.5, n1, n2) ** 0.5
        se1 = hanley_mcneil_variance(a, n1, n2) ** 0.5
        condition_lhs = a - 0.5
        condition_rhs = z_alpha * se0 + z_beta * se1
        if condition_lhs >= condition_rhs:
            power, se0_check, se1_check, n1_check, n2_check = achieved_power(
                a, n_total, z_alpha
            )
            return {
                "auc": a,
                "n_per_group": n1,
                "n_total": n_total,
                "se0": se0,
                "se1": se1,
                "achieved_power": power,
            }
        n_total += N_STEP
    raise RuntimeError(
        f"No balanced N <= {N_MAX} satisfied the power condition for AUC={a}."
    )


def main():
    lines = []

    if scipy is None or norm is None:
        message = (
            "STOP: SciPy is not available in this Python environment "
            f"({sys.executable}, Python {platform.python_version()}). "
            "Per task instructions, the normal quantiles required for this "
            "calculation (z_(1-alpha/2), z_(1-beta), and the power CDF) were "
            "NOT approximated by any substitute. No sample-size results were "
            "computed. Install SciPy in this environment and re-run."
        )
        print(message)
        OUTPUT_PATH.write_text(message + "\n", encoding="utf-8")
        sys.exit(1)

    z_alpha = norm.ppf(1 - ALPHA / 2)   # z_(1-alpha/2), two-sided alpha=0.05
    z_beta = norm.ppf(TARGET_POWER)     # z_(1-beta), power=0.80

    lines.append("AUC POWER CALCULATION — Hanley & McNeil (1982) single-AUC variance")
    lines.append("=" * 78)
    lines.append("")
    lines.append(f"alpha (two-sided): {ALPHA}")
    lines.append(f"target power: {TARGET_POWER}")
    lines.append(f"z_(1-alpha/2) = {z_alpha!r}")
    lines.append(f"z_(1-beta)   = {z_beta!r}")
    lines.append(f"Null AUC: 0.5")
    lines.append(f"N search: start={N_START}, step={N_STEP} (balanced n1=n2=N/2), max={N_MAX}")
    lines.append("")

    results = []
    for a in AUC_VALUES:
        result = find_required_n(a, z_alpha, z_beta, TARGET_POWER)
        results.append(result)

    lines.append("RESULTS")
    lines.append("-" * 78)
    for r in results:
        lines.append(f"AUC = {r['auc']}")
        lines.append(f"  required n per group : {r['n_per_group']}")
        lines.append(f"  required total N      : {r['n_total']}")
        lines.append(f"  SE0 (Var at A=0.5)^0.5 : {r['se0']!r}")
        lines.append(f"  SE1 (Var at A)^0.5     : {r['se1']!r}")
        lines.append(f"  achieved power at N    : {r['achieved_power']!r}")
        lines.append("")

    lines.append("ENVIRONMENT")
    lines.append("-" * 78)
    lines.append(f"Python executable : {sys.executable}")
    lines.append(f"Python version    : {platform.python_version()} ({sys.version})")
    lines.append(f"SciPy version     : {scipy.__version__}")
    if numpy is not None:
        lines.append(f"NumPy version     : {numpy.__version__}")
    else:
        lines.append("NumPy version     : not available (not required by this script's computation path)")

    output_text = "\n".join(lines) + "\n"
    OUTPUT_PATH.write_text(output_text, encoding="utf-8")
    print(output_text)

    # Full-precision machine-readable result, committed under Results/.
    # Records every input assumption and the computed n for each AUC.
    json_payload = {
        "method": "Hanley & McNeil (1982) single-AUC variance, balanced design",
        "assumptions": {
            "null_auc": 0.5,
            "alpha_two_sided": ALPHA,
            "target_power": TARGET_POWER,
            "class_balance": "balanced (n1 = n2 = N/2)",
            "auc_values": AUC_VALUES,
            "z_alpha_half": z_alpha,
            "z_beta": z_beta,
            "n_search": {"start": N_START, "step": N_STEP, "max": N_MAX},
        },
        "results": results,
        "environment": {
            "python_executable": sys.executable,
            "python_version": platform.python_version(),
            "scipy_version": scipy.__version__,
            "numpy_version": (numpy.__version__ if numpy is not None else None),
        },
    }
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    JSON_PATH.write_text(json.dumps(json_payload, indent=2), encoding="utf-8")
    print("JSON written to:", JSON_PATH)
    print(json.dumps(json_payload, indent=2))


if __name__ == "__main__":
    main()
