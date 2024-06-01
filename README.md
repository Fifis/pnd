# pnd

R package for computing fast and accurate numerical derivatives

This package is accompanied by a vignette:

* Kostyrka, A. V. Using parallel numerical derivatives for fast optimisation and reliable inference in R. *In progress.*

The following article provide the theory behind the methods implemented in this package:

* Kostyrka, A. V. (2024). What are you doing, step size: fast computation of accurate numerical derivatives with finite precision. *In progress.*

To install and load this package, run:
```
if (!"devtools" %in% rownames(installed.packages())) install.packages("devtools")
devtools::install_github("Fifis/pnd")
library(pnd)
```

There are no dependencies (the `parallel` package belongs to the `base` group).

## TODO

* BUG: pnd::Grad(func = sin, x = 1:4, deriv.order = 1) fails
* BUG: zero-order derivatives do not work with the BP algorithm
* MISC: Check with `lintr::lint_package()`!
* FEATURE: Extend the step selection routines to gradients
* FEATURE: Auto-detect parallel type, create a cluster in `Grad` for Windows machines inside the function
* SYNTAX: Add side -1, 0, 1 for compatibility with `numDeriv`
* SYNTAX: Align with the syntax of `optimParallel`
* FEATURE: Add absolute or relative step size
* FEATURE: Step selection in Curtis--Reid: parallelise the evaluation
* FEATURE: add `diagnostics` to return an attribute of grad containing original f evaluations
* FEATURE: Return an attribute for the remainder term in `fdCoef`
* FEATURE: Add the arguments f0 and precomputed list(stencil, f) to reuse existing values
* FEATURE: Add a vector of step sizes for different arguments
* FEATURE: Pass arguments from Grad to `fdCoef`, e.g. allow stencils in `Grad` or stencil matrices like in the presentation
* FEATURE: If `f(x)` is present in multiple stencils, do not compute multiple times
* FEATURE: Add safety checks: func(x) must be numeric of length 1; if NA, fail; if vector, warn (Jacobian)
* FEATURE: Create `method.args` for `Grad` with automatic step selection
* FEATURE: Hessian via direct 4-term difference for a faster evaluation
* FEATURE: Return attribute of the estimated absolute error
* FEATURE: Add print-out showing the order of h and derivative from the Taylor expansion
* FEATURE: Arbitrary mixed orders
* MISC: Write unit tests

