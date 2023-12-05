# pnd

R package for computing fast and accurate numerical derivatives

This package accompanies the following article (currently provided as a vignette):

* Kostyrka, A. V. (2023). What are you doing, step size: computing accurate numerical derivatives with finite precision. *In progress.*

## TODO

### High priority

* Align with the syntax of `optimParallel`
* Auto-detect parallel type, create a cluster in `Grad` for Windows machines inside the function
* Ignore method an `method.args`
* Add stencil to `Grad`
* `fdCoef` for arbitrary orders
* Hessian via symmetrisation of the gradient or directly -- or a faster evaluation
* Return attribute of the estimated absolute error

### Medium priority

* Add print-out showing the order of h and derivative from the Taylor expansion
* Arbitrary mixed orders
* Add side -1, 0, 1 for compatibility.

### Low priority

* Mention that PND does not support 0th-order differences, i.e. interpolation
