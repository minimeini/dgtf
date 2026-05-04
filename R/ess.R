#' Effective sample size of an MCMC chain
#'
#' Computes the effective sample size of a 1-D vector of posterior draws
#' using Geyer's initial positive sequence applied to the FFT-based
#' autocorrelation function. Non-finite entries are dropped before the
#' computation.
#'
#' @param x A numeric vector of posterior draws (a single chain).
#' @return A scalar in `[1, length(x)]`.
#' @examples
#' set.seed(1); ess(rnorm(2000))           # ~ 2000 (i.i.d.)
#' x <- as.numeric(arima.sim(list(ar = 0.9), 2000))
#' ess(x)                                  # ~ 2000 * (1-0.9)/(1+0.9) ≈ 105
#' @export
ess <- function(x) {
    if (!is.numeric(x) || !is.null(dim(x)))
        stop("`x` must be a numeric vector.", call. = FALSE)
    ess_cpp(as.numeric(x))
}
