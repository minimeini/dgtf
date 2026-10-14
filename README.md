# dgtf

> Variational and MCMC inference for Bayesian state-space models of count data.

`dgtf` provides a unified R interface to variational Bayes and MCMC for a general class of dynamic generalized transfer function models for count data.

## Installation

1. Clone the repository from a terminal.
```
git clone git@github.com:minimeini/dgtf.git
```
2. Change into the repository with `cd dgtf`.
3. Install the R package.
    - In terminal: `R CMD build .`
    - Or in R: `devtools::install(".")`


### System requirements

- C++17 compiler
- LAPACK and BLAS, both provided by R

The current release builds without OpenMP, so the parallel sections of the
sequential Monte Carlo filter are evaluated serially.

### R package dependencies

- `LinkingTo`: `Rcpp`, `RcppArmadillo`, `BH`, `RcppProgress`, `pg`.
- `Imports`: `Rcpp`, `ggplot2` (plus base `stats`, `utils`, `graphics`, and
  `grDevices`).

## License

MIT (see [`LICENSE`](LICENSE)).
