# pnd next
- Detailed vignette explaining the mathematics behind the functions with full transparency about the choice of parameters
- Evaluating the accuracy of numerical approximations
- Functions for fast and reliable Hessian computation based on parallel gradient evaluation
- Feature: auto-detecting the number of cores available on multi-core machines to speed up computations

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
