# Power-weighted tree leaves — reproducibility code

Python code for numerical checks, two explanatory figures, and a table associated
with *Explicit centering for leaf counts in power-weighted recursive trees*.

## Run

Python 3.11 or newer is required by the pinned dependencies in `requirements.txt`.
The revised checks were run with Python 3.12.14. The preceding version was also
reported to run successfully with Python 3.11 in the referee's independent check.

```bash
python -m venv .venv
# Linux/macOS:
source .venv/bin/activate
# Windows PowerShell:
# .venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python checks/verify_and_plot.py
```

The script performs:

- exhaustive enumeration for eight small-tree parameter cases;
- 1,000 seeded checks of coupling identities and the deletion bound;
- twelve large-n evaluations of the leading mean correction;
- twelve evaluations of the third-term ratio, reproducing the manuscript table;
- three evaluations of the logarithmic correction at `s = 2`;
- four independent numerical integrations of the third coefficient;
- generation of both explanatory figures.

Run ordinary Python without `-O`, so assertions remain enabled. A failed assertion
stops execution. The large-n ratios are reported as diagnostics; their finite-n
closeness to an asymptotic limit is not an assertion or a proof.

## Outputs

Running the script creates local files ignored by Git:

- `checks/verification_results.json` — all numerical results;
- `figures/growth_example.pdf` and `.png` — the fixed coupling example;
- `figures/centering_scales.pdf` and `.png` — analytic exponent curves.

## Numerical scope

The code evaluates the number of internal vertices, the complement of the leaf
count. The coefficient `third_coefficient(s)` enters the internal-vertex mean
with a positive sign and the leaf mean with a negative sign. The table reports
`Q_n = (E[I_n] - A_s*n**(1/s) + D_s*n**(2-s)) / n**(2/s-1)`.
Its convergence can be slow: at `s = 1.3`, values near `0.57` for the tested sizes
do not estimate the limiting coefficient `0.333543...` accurately.

The logarithmic-series error estimate bounds truncation in exact arithmetic,
not floating-point roundoff. The random trials exercise deterministic identities;
they are not a distributional accuracy study. Figure 1 is a fixed realization,
and Figure 2 plots analytic exponents. Numerical checks supplement the proofs.

## Citation

Software citation metadata is provided in `CITATION.cff`. No archival DOI has
been assigned in this repository version.

Marek Gałązka and Hanna Wdowicka.
