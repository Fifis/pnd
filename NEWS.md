# pnd 0.dev (2024-XX-XX)

- BUG: Derivatives of vectorised functions do not work. Check compatibility with `numDeriv` for functions that take vector `x` are return `f(x)` of the same length. Solution: is f vectorised? Is the output the same length as x? Example: `Grad(sin, 1:4)`
- BUG: Check the example with neural networks where `gr` does not accept the arguments of `f`
- BUG: Matching in the Hessian is too slow -- de-duplicate first
- BUG: 1x1 Hessians?
- SYNTAX: Split the gradient into 1D vectorised input and multi-d non-vectorised? Think about it.
- SYNTAX: Align with the syntax of `optimParallel`, borrow ideas from `forecast::auto.arima`
- FEATURE: Replace the long formula in the default step with zero tolerance: `abs(x) < zero.tol`
- FEATURE: disable parallelisation if `f(x)` takes less than 0.001 s
- FEATURE: homogenise handling of missing values (FUN1, FUN2, ...)
- FEATURE: plug-in step size with an estimated `f'''` with a rule of thumb (add this as a one-step DV, e.g. via `maxit = 1`)
- FEATURE: SW algorithm for arbitrary derivative and accuracy orders
- FEATURE: update the rounding error as the estimated sum of differences and evaluations of f
- FEATURE: Handle NA in step size selection
- FEATURE: Auto-shrink the step size at the beginning of all procedures if FUN(x) is not finite
- FEATURE: Extend the step selection routines to gradients
- FEATURE: Auto-detect parallel type, create a cluster in `Grad` for Windows machines inside the function
- FEATURE: Add absolute or relative step size
- FEATURE: Step selection in Curtis--Reid: parallelise the evaluation
- FEATURE: add `diagnostics` to return an attribute of grad containing original f evaluations
- FEATURE: Add the arguments f0 and precomputed list(stencil, f) to reuse existing values
- FEATURE: Add a vector of step sizes for different arguments
- FEATURE: Pass arguments from Grad to `fdCoef`, e.g. allow stencils in `Grad` or stencil matrices like in the presentation
- FEATURE: If `f(x)` is present in multiple stencils, do not compute multiple times
- FEATURE: Add safety checks: func(x) must be numeric of length 1; if NA, fail; if vector, warn (Jacobian)
- FEATURE: auto-detecting the number of cores available on multi-core machines to speed up computations
- FEATURE: Create `control` or `method.args` for `Grad` with automatic step selection
- FEATURE: Hessian via direct 4-term difference for a faster evaluation
- FEATURE: Functions for fast and reliable Hessian computation based on parallel gradient evaluation
- FEATURE: Return attribute of the estimated absolute error
- FEATURE: Add print-out showing the order of h and derivative from the Taylor expansion
- FEATURE: Arbitrary mixed orders
- DOCUMENTATION: Compare with `stats::numericDeriv`
- MISC: Write the list of controls on the help page of `gradstep()` explicitly!
- MISC: Check which packages depend on `numDeriv` and check compatibility with 10 top!
- MISC: Add links to documentation and tutorials onto the GitHub page.
- MISC: Detailed vignette explaining the mathematics behind the functions with full transparency about the choice of parameters
- DEV: add examples for large matrices (200 x 200)
- DEV: ensure that `Grad` takes all the arguments of `GenD` and `Jacobian`, and vice versa
- DEV: Ensure unit-test coverage >90%
- DEV: Check the compatibility between the function and its documentation
- DEV: Check the release with `todor::todor_package()`, `lintr::lint_package()`, `R CMD check --as-cran`, and `goodpractice::gp()`

# pnd 0.next (2023-06-XX)
- Feature: parallelisation of steps in all selection algorithms
- Feature: Mathur's AutoDX algorithm for step size selection `step.M()`
- Feature: added `Hessian()` that supports central differences (for the moment) and arbitrary accuracy
- Feature: separate `Grad()` and `Jacobian()` that call the workhors, `GenD()`, for compatibility with `numDeriv`

# pnd 0.0.4 (2023-06-10)
- Feature: Stepleman--Winarsky algorithm for step size selection `step.SW()`
- Feature: automated wrapper for step size selection `gradstep()`
- Improvement: Safe handling of function errors and non-finite returns in step-size selection procedures
- Improvement: Finite-difference coefficients gained attributes: Taylor expansion, coefficient on the largest truncated term, and effective accuracy (useful for custom stencils)
- Improvement: added unit tests for core features

# pnd 0.0.3 (2023-06-01)
- Feature: `solveVandermonde()` to solve ill-conditioned problems that arise in weight calculation
- Feature: Dumontet--Vignes algorithm for step size selection `step.DV()`
- Feature: Curtis--Reid algorithm for step size selection `step.CR()` and its modification
- Feature: different step sizes for the gradient
- Fix: if the user supplies a short custom stencil and requests a high accuracy order, it will provide the best available order and produce a warning
- Fix: the output of `Grad()` preserves the names of `x` and `FUN(x)`, which prevents errors in cases where names are required

# pnd 0.0.2 (2023-12-06)
- Fix: bug in stencil calculation

# pnd 0.0.1 (2023-09-01)
- Initial release
- Computing finite-difference coefficients on arbitrary stencils
- Computing numerical gradients with reasonable default step sizes
- Numerical Jacobians
- Support for `mclapply()` on *nix systems only
