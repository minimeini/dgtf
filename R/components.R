## Component constructors --------------------------------------------------
##
## Each constructor returns a small list with class
##   c("dgtf_<comp>_<type>", "dgtf_<comp>", "dgtf_component")
##
## where `type` is the C++ engine's string identifier. The `params` slot
## carries any continuous-valued component parameters (lag distribution
## hyperparameters, error variance, etc.). The C++ side is fed via
## `as_settings()` (see `R/model.R`).

new_component <- function(kind, type, params = list(), ...) {
    structure(
        c(list(type = type, params = params), list(...)),
        class = c(paste0("dgtf_", kind, "_", type),
                  paste0("dgtf_", kind),
                  "dgtf_component")
    )
}

#' @export
print.dgtf_component <- function(x, ...) {
    cls <- class(x)[1]
    nm <- sub("^dgtf_", "", cls)
    cat(sprintf("<%s>", nm))
    if (length(x$params)) {
        cat(":")
        for (p in names(x$params))
            cat(sprintf(" %s=%s", p, format(x$params[[p]])))
    }
    cat("\n")
    invisible(x)
}


# ---- Observation distribution -------------------------------------------

#' Observation-distribution components
#'
#' Constructors for the observational level of a `dgtf_model()`.
#' Counts: `obs_poisson()`, `obs_nbinom()` (mean/dispersion form),
#' `obs_nbinom_p()` (number-of-successes / probability form).
#' Continuous: `obs_normal()`.
#'
#' @return A `dgtf_obs` component object.
#' @name dgtf-obs
NULL

#' @rdname dgtf-obs
#' @export
obs_poisson <- function() new_component("obs", "poisson")

#' @rdname dgtf-obs
#' @export
obs_nbinom <- function() new_component("obs", "nbinom")

#' @rdname dgtf-obs
#' @export
obs_nbinom_p <- function() new_component("obs", "nbinomp")

#' @rdname dgtf-obs
#' @export
obs_normal <- function() new_component("obs", "gaussian")


# ---- Link function ------------------------------------------------------

#' Link-function components
#'
#' Maps the conditional mean to a linear predictor.
#'
#' @return A `dgtf_link` component object.
#' @name dgtf-link
NULL

#' @rdname dgtf-link
#' @export
link_identity <- function() new_component("link", "identity")

#' @rdname dgtf-link
#' @export
link_exponential <- function() new_component("link", "exponential")

#' @rdname dgtf-link
#' @export
link_logistic <- function() new_component("link", "logistic")


# ---- Gain function ------------------------------------------------------

#' Gain-function components
#'
#' Differentiable transforms applied to the first state element inside the
#' transfer function.
#'
#' @return A `dgtf_gain` component object.
#' @name dgtf-gain
NULL

#' @rdname dgtf-gain
#' @export
gain_identity <- function() new_component("gain", "identity")

#' @rdname dgtf-gain
#' @export
gain_softplus <- function() new_component("gain", "softplus")

#' @rdname dgtf-gain
#' @export
gain_ramp <- function() new_component("gain", "ramp")

#' @rdname dgtf-gain
#' @export
gain_exponential <- function() new_component("gain", "exponential")

#' @rdname dgtf-gain
#' @export
gain_logistic <- function() new_component("gain", "logistic")


# ---- Lag distribution ---------------------------------------------------

#' Lag-distribution components
#'
#' The discretised PMF that weights past observations in a sliding-window
#' transfer function.
#'
#' \describe{
#'   \item{\code{lag_lognormal(meanlog, sigma2, window)}}{Discretised
#'     log-normal distribution. \strong{Note:} \code{sigma2} is the
#'     \emph{variance} on the log scale, i.e. \eqn{\sigma^2} where
#'     \eqn{\log X \sim N(\mu, \sigma^2)}. This matches the C++ engine
#'     and the CSDA paper. R's \code{dlnorm()} takes the standard
#'     deviation \code{sdlog = sqrt(sigma2)}, so do not pass \code{sdlog}
#'     directly.}
#'   \item{\code{lag_nbinom(r, kappa, window)}}{Discretised
#'     negative-binomial (Solow form), with \code{r} successes and
#'     \code{kappa} failure probability. Reduces to geometric when
#'     \code{r = 1}.}
#'   \item{\code{lag_uniform(window)}}{Discrete uniform over the window
#'     (used for AR models where each lag gets its own time-varying
#'     coefficient).}
#' }
#'
#' @param meanlog Numeric. Log-normal mean parameter \eqn{\mu}.
#' @param sigma2 Numeric. Log-normal \strong{variance} parameter
#'   \eqn{\sigma^2} on the log scale. Default 0.32.
#' @param r,kappa Negative-binomial hyperparameters.
#' @param window Integer truncation length \eqn{L}. For lognormal and
#'   nbinom lags the C++ engine auto-computes the truncation from the
#'   99.5\% quantile; this parameter is stored as R-side metadata.
#'   For uniform lags it sets the AR order \eqn{p}.
#'
#' @return A `dgtf_lag` component object.
#' @name dgtf-lag
NULL

#' @rdname dgtf-lag
#' @export
lag_lognormal <- function(meanlog = 1.386, sigma2 = 0.32, window = 30L) {
    if (sigma2 <= 0) stop("`sigma2` must be positive.", call. = FALSE)
    if (window < 1L) stop("`window` must be a positive integer.", call. = FALSE)
    new_component("lag", "lognorm",
                  params = list(par1 = meanlog, par2 = sigma2),
                  window = as.integer(window))
}

#' @rdname dgtf-lag
#' @export
lag_nbinom <- function(r = 1, kappa = 0.5, window = 30L) {
    if (r <= 0) stop("`r` must be positive.", call. = FALSE)
    if (kappa <= 0 || kappa >= 1)
        stop("`kappa` must lie strictly in (0, 1).", call. = FALSE)
    if (window < 1L) stop("`window` must be a positive integer.", call. = FALSE)
    new_component("lag", "nbinomp",
                  params = list(par1 = kappa, par2 = r),
                  window = as.integer(window))
}

#' @rdname dgtf-lag
#' @export
lag_uniform <- function(window = 1L) {
    if (window < 1L) stop("`window` must be a positive integer.", call. = FALSE)
    new_component("lag", "uniform",
                  params = list(),
                  window = as.integer(window))
}


# ---- System equation ----------------------------------------------------

#' System-equation components
#'
#' Determines the form of the latent state evolution \eqn{\theta_t = g_t(\cdot) + w_t}.
#'
#' - `sys_identity()`: \eqn{G = I}, no state evolution (e.g. Poisson AR).
#' - `sys_shift()`: shift matrix appropriate for windowed Hawkes / DL models.
#' - `sys_nbinom()`: nonlinear evolution for the Solow distributed-lag model
#'   (`kappa` is the geometric / NB decay parameter).
#'
#' @param kappa Decay parameter for the nonlinear `sys_nbinom()` evolution.
#'
#' @return A `dgtf_sys` component object.
#' @name dgtf-sys
NULL

#' @rdname dgtf-sys
#' @export
sys_identity <- function() new_component("sys", "identity")

#' @rdname dgtf-sys
#' @export
sys_shift <- function() new_component("sys", "shift")

#' @rdname dgtf-sys
#' @export
sys_nbinom <- function(kappa = 0.5) {
    if (kappa <= 0 || kappa >= 1)
        stop("`kappa` must lie strictly in (0, 1).", call. = FALSE)
    new_component("sys", "nbinom", params = list(kappa = kappa))
}


# ---- System-error distribution ------------------------------------------

#' System-error distribution components
#'
#' - `err_gaussian(W, full_rank, w0)`: Gaussian system error with variance `W`
#'   (scalar or matrix). `full_rank = TRUE` allows correlated state errors.
#' - `err_constant()`: degenerate (zero) system error.
#'
#' @param W Variance (scalar or square matrix).
#' @param full_rank If `TRUE`, treat `W` as a full covariance matrix.
#' @param w0 Initial state variance offset.
#'
#' @return A `dgtf_err` component object.
#' @name dgtf-err
NULL

#' @rdname dgtf-err
#' @export
err_gaussian <- function(W = 0.01, full_rank = FALSE, w0 = 0) {
    Wm <- if (is.matrix(W)) W else matrix(W)
    if (any(diag(Wm) < 0))
        stop("Diagonal of `W` must be non-negative.", call. = FALSE)
    new_component("err", "gaussian",
                  params = list(var = Wm, w0 = w0, full_rank = isTRUE(full_rank)))
}

#' @rdname dgtf-err
#' @export
err_constant <- function() new_component("err", "constant")


# ---- Coercion of string shortcuts ---------------------------------------

#' @keywords internal
as_component <- function(x, kind) {
    if (inherits(x, paste0("dgtf_", kind))) return(x)
    if (is.character(x) && length(x) == 1L) {
        ctor <- switch(
            kind,
            obs  = list(poisson = obs_poisson, nbinom = obs_nbinom,
                        nbinomp = obs_nbinom_p, normal = obs_normal,
                        gaussian = obs_normal),
            link = list(identity = link_identity, exponential = link_exponential,
                        logistic = link_logistic),
            gain = list(identity = gain_identity, softplus = gain_softplus,
                        ramp = gain_ramp, exponential = gain_exponential,
                        logistic = gain_logistic),
            lag  = list(lognorm = lag_lognormal, lognormal = lag_lognormal,
                        nbinom = lag_nbinom, nbinomp = lag_nbinom,
                        uniform = lag_uniform),
            sys  = list(identity = sys_identity, shift = sys_shift,
                        nbinom = sys_nbinom),
            err  = list(gaussian = err_gaussian, normal = err_gaussian,
                        constant = err_constant)
        )
        f <- ctor[[tolower(x)]]
        if (is.null(f))
            stop(sprintf("Unknown %s component: \"%s\".", kind, x),
                 call. = FALSE)
        return(f())
    }
    stop(sprintf("`%s` must be a dgtf_%s object or a recognised string.",
                 kind, kind), call. = FALSE)
}
