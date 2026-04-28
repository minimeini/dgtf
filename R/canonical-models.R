#' Canonical DGTF model templates
#'
#' Convenience constructors for the canonical model variants discussed in
#' the manuscript: Poisson autoregressions, discretised Hawkes, Koyck's
#' first-order distributed lag, and the higher-order distributed-lag model.
#' Each is a thin wrapper around [`dgtf_model()`] with the appropriate
#' components pre-filled, so the same components can be swapped out
#' downstream to study close variants.
#'
#' @name dgtf-canonical
NULL

#' Poisson autoregression / integer ARCH model
#'
#' \deqn{y_t \mid \lambda_t \sim \mathrm{Pois}(\lambda_t),\quad
#'       \lambda_t = a + \sum_{i=1}^p b_i y_{t-i}.}
#'
#' Implemented via `obs_poisson() + link_identity() + sys_identity() +
#' gain_identity() + lag_uniform(window = p) + err_constant()`.
#'
#' @param p Order of the autoregression.
#' @param intercept Baseline `a`.
#' @param ... Forwarded to [`dgtf_model()`].
#'
#' @return A `dgtf_model`.
#' @export
#' @examples
#' dgtf_poisson_ar(p = 3)
dgtf_poisson_ar <- function(p = 1L, intercept = 0, ...) {
    dgtf_model(
        obs       = obs_poisson(),
        link      = link_identity(),
        sys       = sys_identity(),
        gain      = gain_identity(),
        lag       = lag_uniform(window = as.integer(p)),
        err       = err_constant(),
        intercept = intercept,
        ...
    )
}

#' Discretised Hawkes model (Koyama-style)
#'
#' \deqn{y_t \mid \lambda_t \sim \mathrm{NB}(\lambda_t, \rho),\quad
#'       \lambda_t = a + \sum_{l=1}^L \phi_l\, h(\psi_{t+1-l})\, y_{t-l},\quad
#'       \psi_t = \psi_{t-1} + w_t.}
#'
#' @param lag A lag-distribution component (default
#'   `lag_lognormal(meanlog = 1.386, sdlog = 0.323, window = 30)`).
#' @param gain Gain function (default `gain_softplus()`).
#' @param obs Observation distribution (default `obs_nbinom()`).
#' @param err System-error distribution (default `err_gaussian(W = 0.01)`).
#' @param intercept Baseline `a`.
#' @param ... Forwarded to [`dgtf_model()`].
#'
#' @return A `dgtf_model`.
#' @export
#' @examples
#' dgtf_hawkes()
#' dgtf_hawkes(lag = lag_nbinom(r = 2, kappa = 0.5, window = 30))
dgtf_hawkes <- function(lag       = lag_lognormal(),
                        gain      = gain_softplus(),
                        obs       = obs_nbinom(),
                        err       = err_gaussian(W = 0.01),
                        intercept = 0,
                        ...) {
    dgtf_model(
        obs       = obs,
        link      = link_identity(),
        sys       = sys_shift(),
        gain      = gain,
        lag       = lag,
        err       = err,
        intercept = intercept,
        ...
    )
}

#' Koyck's first-order distributed lag model
#'
#' \deqn{\beta_t = \kappa\beta_{t-1} + (1-\kappa) h(\psi_t) y_{t-1},\quad
#'       \lambda_t = a + \beta_t.}
#'
#' Equivalent to a discretised Hawkes with a geometric lag distribution.
#'
#' @param kappa Memory parameter.
#' @param gain Gain function (default `gain_softplus()`).
#' @param obs Observation distribution (default `obs_nbinom()`).
#' @param err System-error distribution (default `err_gaussian(W = 0.01)`).
#' @param intercept Baseline `a`.
#' @param ... Forwarded to [`dgtf_model()`].
#'
#' @return A `dgtf_model`.
#' @export
#' @examples
#' dgtf_koyck(kappa = 0.5)
dgtf_koyck <- function(kappa     = 0.5,
                       gain      = gain_softplus(),
                       obs       = obs_nbinom(),
                       err       = err_gaussian(W = 0.01),
                       intercept = 0,
                       ...) {
    dgtf_model(
        obs       = obs,
        link      = link_identity(),
        sys       = sys_nbinom(kappa = kappa),
        gain      = gain,
        lag       = lag_nbinom(r = 1, kappa = kappa, window = 1L),
        err       = err,
        intercept = intercept,
        ...
    )
}

#' Generalised distributed lag model (Solow form)
#'
#' Negative-binomial lag distribution with `r` successes and failure
#' probability `kappa`; reduces to Koyck when `r = 1`. Uses the recursive
#' `sys_nbinom()` evolution for efficiency on long histories.
#'
#' @param r Number of successes (positive integer).
#' @param kappa Failure probability in (0, 1).
#' @param gain Gain function (default `gain_softplus()`).
#' @param obs Observation distribution (default `obs_poisson()`).
#' @param err System-error distribution (default `err_gaussian(W = 0.01)`).
#' @param intercept Baseline `a`.
#' @param ... Forwarded to [`dgtf_model()`].
#'
#' @return A `dgtf_model`.
#' @export
#' @examples
#' dgtf_distributed_lag(r = 2, kappa = 0.5)
dgtf_distributed_lag <- function(r         = 2L,
                                 kappa     = 0.5,
                                 gain      = gain_softplus(),
                                 obs       = obs_poisson(),
                                 err       = err_gaussian(W = 0.01),
                                 intercept = 0,
                                 ...) {
    dgtf_model(
        obs       = obs,
        link      = link_identity(),
        sys       = sys_nbinom(kappa = kappa),
        gain      = gain,
        lag       = lag_nbinom(r = r, kappa = kappa, window = 1L),
        err       = err,
        intercept = intercept,
        ...
    )
}
