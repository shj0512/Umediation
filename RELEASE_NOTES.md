# Umediation 1.0.0

First public release.

Umediation estimates local indirect effects for a quadratic exposure-to-mediator path and a linear mediator-to-outcome path. It supports common controls, fixed effects, pairs bootstrap inference, cluster jackknife projection, and a separate wild cluster bootstrap reversal test. `umediationplot` displays local effects and confidence limits.

Requires Stata 17 or later; no external packages are needed. Installation:

```stata
net install umediation, from("https://raw.githubusercontent.com/shj0512/Umediation/v1.0.0")
help umediation
```

Finite-sample performance depends on the cluster design and path strengths. The projection and wild procedures are experimental; validation results and limitations are documented in the help file and `validation/README.md`.
