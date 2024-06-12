<!-- badges: start -->
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
- **Four optimal step selection algorithms:** choose between Curtis–Reid (1974) and its modern (2024) modification, Dumontet–Vignes (1977), and Stepleman–Winarsky (1979) algorithms. Future versions will feature parallelised algorithms.

## Learning resources

- [PDF of an early 2024 presentation at the University of Luxembourg.](https://kostyrka.lu/en/education/presentations/2024-brown-bag-seminar.pdf) *(Some functions might be outdated – check the package vignette for up-to-date examples!)*

## Literature

This package is supported by a vignette:

* Kostyrka, A. V. Using parallel numerical derivatives for fast optimisation and reliable inference in R. *In progress.*

The following articles provide the theory behind the methods implemented in this package:

* Kostyrka, A. V. (2024). What are you doing, step size: fast computation of accurate numerical derivatives with finite precision. *In progress.*

## Installation

This package currently exists only on GitHub. To install it, run the following two commands:
```{r}
install.packages("devtools")
devtools::install_github("Fifis/pnd")
```

To load this package, include this line in the code:
```{r}
library(pnd)
```

This package is almost dependency-free; the `parallel` library belongs to the `base`
group and is included in most R distributions.

## Licence

This software is released under the free/open-source [EUPL 1.2 licence](https://joinup.ec.europa.eu/collection/eupl/eupl-text-eupl-12).
