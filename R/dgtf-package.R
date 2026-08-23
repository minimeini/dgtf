#' dgtf: Variational and Monte Carlo Inference for Bayesian State-Space Models of Count Data
#'
#' The `dgtf` package provides a unified interface to several inference
#' algorithms (variational Bayes, MCMC, sequential Monte Carlo, and linear
#' Bayes) for a general class of dynamic generalized transfer function models
#' for count data.
#'
#' @keywords internal
#' @useDynLib dgtf, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom graphics plot
#' @importFrom stats coef confint fitted logLik nobs predict residuals vcov
"_PACKAGE"
