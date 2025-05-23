# pnd 0.dev roadmap (2025-XX-XX)

- UX: make the warnings once-per-session; print the first error in `runParallel` in Grad
- FEATURE: add options for default accuracy (maybe 4?)
- FEATURE: If `h` is a character in `Grad`, extract the gradient directly if the order is 2
- FEATURE: disable parallelisation if `f(x)` takes less than 0.002 s
- FEATURE: CR, DW, SW, M, K algorithm for arbitrary derivative and accuracy orders
- FEATURE: if any algorithm returns a step size larger than `|x| / max(stencil)`, throw a warning
- FEATURE: Create `control` or `method.args` for `Grad` with automatic step selection
- FEATURE: Arbitrary mixed orders
- MISC: Write the list of controls on the help page of `gradstep()` explicitly!
- MISC: Check which packages depend on `numDeriv` and check compatibility with 10 top!
- MISC: Add links to documentation and tutorials onto the GitHub page.
- MISC: Detailed vignette explaining the mathematics behind the functions with full transparency about the choice of parameters
- DEV: ensure that `Grad` takes all the arguments of `GenD` and `Jacobian`, and vice versa
- DEV: Ensure unit-test coverage >90%
- DEV: Check the compatibility between the function and its documentation
- DEV: Check the release with `todor::todor_package()`, `lintr::lint_package()`, `R CMD check --as-cran`, and `goodpractice::gp(checks = all_checks()[!grepl("^lintr", all_checks())])`

# pnd 0.1.1 (2025-07-XX)
- Fix: added a simpler and more reliable fall-back option for `step.M`
- Feature: added plotting methods for `step...` functions
- Feature: added arbitrary derivation and accuracy order to Curtis--Reid method

# pnd 0.1.0 (2025-05-20)
- Feature: original kink-based algorithm for step size selection `step.K()`
- Feature: added safety shrinking if `FUN(x)` is finite but `FUN(x+h)` is not in all SSS routines
- Feature: added S3 printing methods for derivatives and step sizes
- Feature: removed the `diagnostics` and `report` arguments; the iteration information is always saved, but not printed
- Feature: added support for `max.rel.error` for all step-selection methods
- Feature: all step-search methods now return both the truncation and the rounding-error estimate
- Fix: corrected the wrong formula for the plug-in step size
- Fix: now 1x1 Hessians can be computed (why, though, if second derivatives exist?)
- Fix: added the `v` argument for `numDeriv` compatibility

# pnd 0.0.10 (2025-04-02)
- Fix: GitHub issue #2 -- `checkDimensions` could not handle character `h` passed for auto-selection
- Fix: GitHub issue #1 -- function arguments in `...` did non propagate properly to `step...` functions

# pnd 0.0.9 (2025-03-10)
- Fix: fixed a regression with the default step size
- Fix: parallelised Hessians in the same manner as gradients
- Feature: compatibility of `Hessian()` with the arguments for methods `"Richardson"` and `"simple"` from `numDeriv`

# pnd 0.0.8 (2025-03-05)
- Fix: sped up CPU core request diagnostics for 1-core operations
- Fix: Using full paths on Macs

# pnd 0.0.7 (2025-03-01)
- Fix: removed obsolete environment creation for cluster export
- Fix: changed physical core detection on Macs
- Misc: the package has been added to CRAN, fewer syntax changes are expected

# pnd 0.0.6 (2025-01-25)
- Fix: Derivatives of vectorised functions are working. Example: `Grad(sin, 1:4)`
- Feature: Auto-detecting the number of cores available on multi-core machines to speed up computations
- Feature: Added plug-in step size selection with an estimated `f'''` with a rule of thumb
- Feature: Auto-detection of parallel type
- Feature: Added zero tolerance to the default step for a fixed step

# pnd 0.0.5 (2025-01-14)
- Feature: Extended the step-selection routines to gradients (vector input `x`)
- Feature: Parallelisation of step selection in all algorithms
- Feature: Mathur's AutoDX algorithm for step size selection `step.M()`
- Feature: Added `Hessian()` that supports central differences (for the moment) and arbitrary accuracy
- Feature: Separate `Grad()` and `Jacobian()` that call the workhorse, `GenD()`, for compatibility with `numDeriv`

# pnd 0.0.4 (2024-06-10)
- Feature: Stepleman--Winarsky algorithm for step size selection `step.SW()`
- Feature: Automated wrapper for step size selection `gradstep()`
- Improvement: Safe handling of function errors and non-finite returns in step-size selection procedures
- Improvement: Finite-difference coefficients gained attributes: Taylor expansion, coefficient on the largest truncated term, and effective accuracy (useful for custom stencils)
- Improvement: added unit tests for core features

# pnd 0.0.3 (2024-06-01)
- Feature: `solveVandermonde()` to solve ill-conditioned problems that arise in weight calculation
- Feature: Dumontet--Vignes algorithm for step size selection `step.DV()`
- Feature: Curtis--Reid algorithm for step size selection `step.CR()` and its modification
- Feature: Different step sizes for the gradient
- Fix: If the user supplies a short custom stencil and requests a high accuracy order, it will provide the best available order and produce a warning
- Fix: The output of `Grad()` preserves the names of `x` and `FUN(x)`, which prevents errors in cases where names are required

# pnd 0.0.2 (2023-12-06)
- Fix: bug in stencil calculation

# pnd 0.0.1 (2023-09-01)
- Initial release
- Computing finite-difference coefficients on arbitrary stencils
- Computing numerical gradients with reasonable default step sizes
- Numerical Jacobians
- Support for `mclapply()` on *nix systems only
