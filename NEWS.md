# pnd 0.next (2024-XX-XX)
- SYNTAX: should `x` bo before `FUN` in `Grad`? Check for 2 first unnamed argument types for compatibility.
- SYNTAX: `Grad()` fails if `h` is text and a Jacobian is needed.
- BUG: pnd::Grad(func = sin, x = 1:4, deriv.order = 1) fails
- BUG: zero-order derivatives do not work with the BP algorithm
- BUG: the attributes on long stencils, e.g. `fdCoef(deriv.order = 3, stencil = -4:4)`, are full of digits -- `zero.tol` must be attached to the coefficient magnitude
- MISC: Write the list of controls on the help page of `gradstep()` explicitly!
- MISC: Check which packages depend on `numDeriv` and check compatibility with 10 top!
- MISC: Add links to documentation and tutorials onto the GitHub page.
- Evaluating the accuracy of numerical approximations
- FEATURE: SW algorithm for arbitrary derivative and accuracy orders
- FEATURE: Handle NA in step size selection
- FEATURE: Auto-shrink the step size at the beginning of all procedures if FUN(x) is not finite
- FEATURE: Extend the step selection routines to gradients
- FEATURE: Auto-detect parallel type, create a cluster in `Grad` for Windows machines inside the function
- SYNTAX: Add side -1, 0, 1 for compatibility with `numDeriv`
- SYNTAX: Align with the syntax of `optimParallel`
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
- MISC: Write unit tests
- MISC: Detailed vignette explaining the mathematics behind the functions with full transparency about the choice of parameters
- MISC: Check the release with `todor`, `lintr::lint_package()`, `R CMD check --as-cran`, and `library(goodpractice)`

# pnd 0.0.4 (2023-06-10)
- Feature: Stepleman--Winarsky algorithm for step size selection `step.SW()`
- Feature: automated wrapper for step size selection `gradstep()`
- Improvement: Safe handling of function errors and non-finite returns in step-size selection procedures.
- Improvement: Finite-difference coefficients gained attributes: Taylor expansion, coefficient on the largest truncated term, and effective accuracy (useful for custom stencils).
- Improvement: added unit tests for core features

# pnd 0.0.3 (2023-06-01)
- Feature: `solveVandermonde()` to solve ill-conditioned problems that arise in weight calculation
- Feature: Dumontet--Vignes algorithm for step size selection `step.DV()`
- Feature: Curtis--Reid algorithm for step size selection `step.CR()` and its modification
- Feature: different step sizes for the gradient
- Fix: if the user supplies a short custom stencil and requests a high accuracy order, it will provide the best available order and produce a warning
- Fix: the output of `Grad()` preserves the names of `x` and `FUN(x)`, which prevents errors in cases where names are required.

# pnd 0.0.2 (2023-12-06)
- Fix: bug in stencil calculation

# pnd 0.0.1 (2023-09-01)
- Initial release
- Computing finite-difference coefficients on arbitrary stencils
- Computing numerical gradients with reasonable default step sizes
- Numerical Jacobians
- Support for `mclapply()` on *nix systems only
