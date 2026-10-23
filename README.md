# Umediation

Mediation analysis with a quadratic exposure-to-mediator relationship in Stata.

Version 1.0.0.

Umediation fits a quadratic X-to-M path and a linear M-to-Y path on a common sample. It estimates the local indirect effect `IE(x) = b(a1 + 2a2x)`, decomposes outcome curvature into mediated and direct components, and assesses sign reversals over a specified exposure interval.

## Features and requirements

- Common controls, factor variables, and categorical fixed effects in both equations.
- Observation or whole-cluster pairs bootstrap, with percentile intervals and a simultaneous grid band.
- Cluster jackknife projection intervals, root confidence sets, and deletion diagnostics.
- An optional wild cluster bootstrap test of indirect-effect reversal.
- A plotting command and a self-contained example using synthetic data.

Stata 17 or later is required. The commands use built-in Stata and Mata routines, with no external package or Python dependency. The model uses continuous X, M, and Y, one mediator, unweighted least squares, and at most one clustering dimension.

## Installation

Install version 1.0.0:

```stata
net install umediation, from("https://raw.githubusercontent.com/shj0512/Umediation/v1.0.0")
help umediation
help umediationplot
```

For installation from an extracted source folder, use the absolute path to the folder containing `umediation.pkg`:

```stata
net install umediation, from("C:/path/to/Umediation")
```

Add `replace` when deliberately replacing an installed copy. To download the example, run the following in an empty working directory, using the same installation source:

```stata
net get umediation, from("https://raw.githubusercontent.com/shj0512/Umediation/v1.0.0")
do umediation_example.do
```

The example preserves the existing data and random-number state. Its default 199 resampling draws demonstrate operation; `do umediation_example.do 1999` requests more draws. No SSC installation address is available until the package is separately accepted by SSC.

## Basic use

For an outcome `y`, exposure `x`, mediator `m`, control `z`, firm identifier `firm`, and year identifier `year`:

```stata
umediation y, x(x) mediator(m) controls(z) absorb(firm year) cluster(firm) ///
    reps(1999) seed(42)
umediationplot
```

To use cluster jackknife projection and request the additional wild test:

```stata
umediation y, x(x) mediator(m) controls(z) absorb(firm year) cluster(firm) ///
    inference(jackknife) wildreps(1999) seed(42)
matrix list e(wild_reversal)
matrix list e(root_set)
```

The primary inference method determines `e(ci)` and `e(classification)`. The optional wild test reports separate p-values and leaves those results unchanged. The default evaluation positions are the 5th, 25th, 50th, 75th, and 95th exposure percentiles; the default research interval spans the 5th to 95th percentiles. These positions are fixed across replications. See `help umediation` for all options, stored results, and methods.

## Interpretation and validation

Coverage and reversal error rates depend on the design. Pairs inference can under-cover and produce excess false reversals with few or imbalanced clusters. Jackknife projection and the wild branch test are experimental inference procedures. Projection intervals can be highly conservative. The four-component confidence region may under-cover when clusters are imbalanced, and high leverage can substantially reduce power. Weak paths and roots near the interval boundary are also difficult. Increasing bootstrap repetitions improves simulation precision without resolving these issues.

Choose the research interval, direction, and primary method before inspecting results. A root estimate inside the interval alone does not establish a reversal, and failure to detect a reversal does not establish its absence. Fixed effects do not by themselves identify a causal mediation effect.

The [validation documentation](validation/README.md) gives simulation results, reproducibility instructions, and implementation checks. Tests were run on Windows with Stata/MP 17. Other operating systems and Stata versions have not been tested. The statistical limitations apply independently of the software version number.

## Citation

SUN, Huijie. Umediation: Mediation with a quadratic exposure-to-mediator relationship in Stata. Version 1.0.0. University of Science and Technology of China.

Machine-readable citation information is provided in [CITATION.cff](CITATION.cff). Cite the inference methods used alongside the software; references are listed in the help file.

## Author

Huijie SUN  
University of Science and Technology of China  
[shj0512@mail.ustc.edu.cn](mailto:shj0512@mail.ustc.edu.cn)  
[github.com/shj0512](https://github.com/shj0512)

## License

Copyright (c) 2026 Huijie SUN. Distributed under the [MIT License](LICENSE).
