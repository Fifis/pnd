<!-- badges: start -->
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![R-CMD-check](https://github.com/Fifis/pnd/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Fifis/pnd/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/Fifis/pnd/graph/badge.svg?token=2ZTHBCRLBR)](https://app.codecov.io/gh/Fifis/pnd)
<!-- badges: end -->

# pnd

An R package for computing fast and accurate numerical derivatives.

<img src="https://kostyrka.lu/user/pages/05.programming/05.pnd.package/parallel-numerical-derivatives-R-package.png" alt="Parallel numerical derivatives in R" width="640"/>

In the past, I was using [numDeriv](https://CRAN.R-project.org/package=numDeriv) to compute numerical gradients.
However, the results were not stable for some function, and I could not investigate the source of this instability.
Different step sizes yielded different results. Small step sizes were sometimes better, sometimes worse.

The `pnd` package was designed to offer a comprehensive tool-kit containing popular algorithms for finite differences, numerical gradients, Jacobians, and Hessians.

Optimal step sizes and parallel evaluation of numerical derivatives translate directly to faster numerical optimisation and statistical inference.


## Features
- **Robust numerical differentiation:** effortlessly compute derivatives while controlling the accuracy-speed trade-off.
- **Gradient and Hessian calculations:** obtain the direction and curvature required by most quasi-Newton optimisation algorithms.
- **Parallel capabilities:** evaluate multiple values under the best parallelisation scheme that reduces overhead. For example, on a 12-core machine, a 4th-order accurate Jacobian of a 3-dimensional function takes almost the same amount of time as one function evaluation.
- **Optimal step size selection:** obtain adaptive step size to ensure the best trade-off between mathematical truncation error and computer floating-point rounding error for the best overall accuracy.
- **Five optimal step selection algorithms:** choose between Curtis–Reid (1974) and its modern (2025) modification, Dumontet–Vignes (1977), Stepleman–Winarsky (1979), and Mathur (2012) algorithms. Future versions will feature parallelised algorithms.

## Getting started

This package has `numDeriv`-compatible syntax.
Simply replace the first letter of `numDeriv` commands with a capital one to get the improved commands: `Grad`, `Jacobian`, and `Hessian`.

Here is how to compute the gradient of `f(x) = sum(sin(x))` at the point `x = (1, 2, 3, 4)`.

```r
f <- function(x) sum(sin(x))
x <- 1:4
names(x) <- c("Jan", "Feb", "Mar", "Apr")

numDeriv::grad(f, x)
#> 0.5403023 -0.4161468 -0.9899925 -0.6536436

pnd::Grad(f, x)
#>       Jan        Feb        Mar        Apr
#> 0.5403023 -0.4161468 -0.9899925 -0.6536436
#> attr(,"step.size")
#>          Jan          Feb          Mar          Apr
#> 6.055454e-06 1.211091e-05 1.816636e-05 2.422182e-05
#> attr(,"step.size.method")
#> "default"
```

The output contains diagnostic information about the chosen step size. Our function
preserved the names of the input argument, unlike `grad`.

The default step size in many implementations is proportional to the argument value, and this is reflected in the default output.
Should the user desire a fixed step size, this can be easily achieved with an extra argument named `h`:

```r
pnd::Grad(f, x, h = c(1e-5, 1e-5, 1e-5, 2e-5))
#>       Jan        Feb        Mar        Apr 
#> 0.5403023 -0.4161468 -0.9899925 -0.6536436 
#> attr(,"step.size")
#>   Jan   Feb   Mar   Apr 
#> 1e-05 1e-05 1e-05 2e-05 
attr(,"step.size.method")
#> "user-supplied"
```

Finally, it is easy to request an algorithmically chosen optimal step size -- here is how to do it with the Stepleman--Winarsky (1979) rule, named `"SW"`, that works well in practice:

```r
pnd::Grad(f, x, h = "SW")
#>       Jan        Feb        Mar        Apr 
#> 0.5403023 -0.4161468 -0.9899925 -0.6536436 
#> attr(,"step.size")
#>          Jan          Feb          Mar          Apr 
#> 5.048535e-06 1.000000e-05 7.500000e-06 1.000000e-05 
#> attr(,"step.size.method")
#> "SW"
```

Extensive diagnostics and error estimates can be requested at any time:
`pnd::Grad(f, x, h = "SW", report = 2)` will contain the step-search path for each coordinate of `x`.
Use `report = 0` to produce just the numerical gradient without any attributes, like `numDeriv::grad` would.

## Learning resources

- [PDF of a 2025 presentation at the University of Luxembourg.](https://kostyrka.lu/en/education/presentations/2025-dem-internal-seminar.pdf)
- [PDF of an early 2024 presentation at the University of Luxembourg.](https://kostyrka.lu/en/education/presentations/2024-brown-bag-seminar.pdf) *(Obsolete – check the one above or the vignettes for up-to-date examples!)*

## Literature

This package is supported by 3 vignettes:

* Kostyrka, A. V. Fast and accurate parallel numerical derivatives in R. *In progress.*
* Kostyrka, A. V. Compatilibility of pnd with the syntax of numDeriv. *In progress.*
* Kostyrka, A. V. Step-size-selection algorithm benchmark. *In progress.*

The following articles provide the theory behind the methods implemented in this package:

* Kostyrka, A. V. (2025). What are you doing, step size: fast computation of accurate numerical derivatives with finite precision. *In progress.*

## Installation

This package currently exists only on GitHub. To install it, run the following two commands:
```r
install.packages("devtools")
devtools::install_github("Fifis/pnd")
```

To load this package, include this line in the code:
```r
library(pnd)
```

This package is almost dependency-free; the `parallel` library belongs to the `base`
group and is included in most R distributions.

## Licence

This software is released under the free/open-source [EUPL 1.2 licence](https://interoperable-europe.ec.europa.eu/collection/eupl/eupl-text-eupl-12).
