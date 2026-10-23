# Validation and reproducibility

Umediation has two distinct validation tasks: checking the numerical implementation and assessing the statistical performance of its inference procedures. Numerical agreement alone does not establish confidence-interval coverage or control of false reversal decisions.

## Implementation checks

From the source root, run in Stata 17 or later:

```stata
do validation/check.do
```

The script uses synthetic data. It checks coefficients against explicit dummy-variable least squares; fixed effects and factor controls; reproducibility; separation of primary and wild inference; exposure translation; numeric and string identifiers; root-set edge cases; finite-support enumeration; handling of failed deletions; output-path collisions; preservation of data, frames, and random-number state; help rendering; and the standalone example. It writes results to `validation/output/` and clears the session at the start. Run it in a separate Stata session.

All numerical matrices, saved draw datasets, and classifications were also compared with the reference calculation source in seven representative designs. All 17,730 compared numeric cells agreed exactly. This comparison covers observation bootstrap, clustered fixed effects, jackknife projection, WCR-S, WCR-C, and Webb and Rademacher enumeration.

Actual Stata execution has been checked on Windows with Stata/MP 17. Other operating systems and Stata versions remain untested. Installation and ancillary-file retrieval have been checked from a local folder and a local HTTP server, including the downloaded example and license file. Public GitHub download URLs can be checked after the repository is published.

## Monte Carlo results

`simulation_summary.csv` records the completed simulation study. Rates are proportions. Columns ending in `_lo` and `_hi` give 95% Wilson Monte Carlo intervals for binary performance measures. Intervals for mediated curvature and the indirect-effect curve have nominal level 95%; reversal decisions use nominal level 5% and `direction(either)`.

Selected results for the balanced G=5, N=400 design, with 5,000 outer replications and 1,999 resampling draws:

| Performance measure | Pairs | Jackknife projection | WCR-S |
|---|---:|---:|---:|
| Mediated-curvature coverage, regular alternative | 84.46% | 99.96% | No product interval |
| Indirect-effect band coverage, regular alternative | 71.90% | 99.90% | No band |
| False reversal rate, b=0 | 16.08% | 0.62% | 1.08% |
| Correct-direction power, regular alternative | 93.28% | 23.02% | 40.02% |
| False reversal rate, root at interval boundary | 8.32% | 0.04% | 0.96% |

At G=100, N=400, false reversal rates under b=0 were 6.30%, 1.32%, and 2.76%; regular-alternative power was 98.68%, 85.20%, and 97.26%, respectively (5,000 outer replications). More clusters did not yield identical calibration across methods.

For the highly imbalanced G=10, N=400 design, with cluster sizes 220 and nine times 20, the nominal 95% four-component region covered 91.15% of the time (95% Monte Carlo interval 89.83% to 92.32%; 2,000 replications). High coverage for a projected product can coexist with undercoverage of the component region. In the high-leverage design, WCR-S had zero observed power under the regular alternative; projection power was 0.75%. Low false-positive rates in such a setting provide little evidence about useful power.

The full design comprises 74 structural scenarios across 16 data strata. It includes weak individual paths, a double-zero null, roots inside, outside, and at the interval boundary, different cluster counts, unequal cluster sizes, heavy-tailed errors, and concentrated exposure information. Parameters and seeds are recorded in `design.json`.

The simulation phases share samples. The `precision` phase uses repetitions 1 to 5,000, including the first 2,000 used in overlapping core scenarios. The `screen` phase uses subsets of the same samples. Different structural paths are evaluated on common generated data through exact least-squares parameter translations. There are 41,000 distinct base samples and 614,000 method-by-scenario records; these counts do not represent 614,000 independent samples or native Stata calls. Phase comparisons must account for this dependence.

These results support reporting the methods' design-specific tradeoffs. They do not establish uniform calibration for few clusters. Jackknife projection and the wild branch test remain experimental inference procedures. The default pairs method also has documented calibration limitations.

## Reproducing the simulations

The optional Python reference implementation is supplied in `ols.py`, `inference.py`, and `simulate.py`. These files are used only for validation; Python is not required to install or run the Stata commands. They reproduce the specified least-squares and resampling calculations with common random numbers. The Stata implementation has separately been compared with observation-level regression calculations.

Use the versions in `requirements.txt` for numerical replication. From the source root:

```text
python -m pip install -r validation/requirements.txt
python validation/simulate.py --smoke --workers 1
python validation/simulate.py --phase all --workers 4
```

The smoke run checks execution with three samples and 199 draws; it is not an accuracy study. The full schedule is computationally intensive. Individual phases can be run with `--phase core`, `--phase screen`, `--phase extensions`, or `--phase precision`. Generated files go to `validation/results/`. Completed tasks are reused only when their task and output hashes match. The core phase recomputes both pairs and jackknife results on the recorded samples. Full reruns may differ in their last floating-point digits across numerical-library versions and hardware.

Results include confidence-interval coverage, false reversal rates, correct-direction power, interval widths, and component coverage. Full simulation output and internal preparation records are excluded from the repository; the scripts generate the synthetic data required for reproduction.
