#' Posterior predictive draws
#'
#' Generates posterior predictive samples for an in-sample series under
#' the fitted model. Use this for posterior-predictive checks and
#' goodness-of-fit assessment.
#'
#' @param object A `dgtf_fit`.
#' @param nrep Number of replicate draws **per posterior sample**. The
#'   total number of posterior-predictive draws is
#'   `nrep * <number of posterior samples>` (number of samples is the
#'   second dimension of `object$fit$psi_stored` / `object$fit$Theta`).
#' @param Rt Optional length-`length(y)` vector of "true" intensity
#'   values used as a reference (e.g. `softplus(sim$psi)` in a
#'   simulation study). Pass `NULL` (default) to skip the reference
#'   comparison.
#' @param ... Unused.
#'
#' @return A list as returned by [`dgtf_posterior_predictive()`]:
#'   posterior predictive intensity samples, replicate `y` draws, and
#'   summary statistics.
#'
#' @export
#' @examples
#' \dontrun{
#' fit <- dgtf(sim$y, mod, prior, method = "hva", control = vb_control())
#' pp  <- posterior_predict(fit, nrep = 100)
#' }
posterior_predict <- function(object, ...) UseMethod("posterior_predict")

#' @rdname posterior_predict
#' @export
posterior_predict.dgtf_fit <- function(object,
                                        nrep = 100L,
                                        Rt   = NULL,
                                        ...) {
    dgtf_posterior_predictive(
        output     = object$fit,
        model_opts = as_settings(object$model),
        y          = as.numeric(object$y),
        nrep       = as.integer(nrep),
        Rt         = if (is.null(Rt)) NULL else as.numeric(Rt)
    )
}


#' Out-of-sample forecasts from a DGTF fit
#'
#' \eqn{k}-step-ahead posterior predictive forecasts. Can also score
#' forecasts against held-out true values when supplied.
#'
#' @param object A `dgtf_fit`.
#' @param h Forecast horizon (positive integer).
#' @param nrep Number of replicate draws per posterior sample.
#' @param ypred_true Optional numeric vector of held-out future
#'   observations of length `h`. When supplied, the returned list
#'   includes per-step coverage and error metrics.
#' @param only_last If `TRUE`, only the final `h`-step-ahead
#'   prediction is evaluated (skips intermediate horizons).
#' @param return_samples If `TRUE`, returns the full
#'   `(ntime + h + 1) x nrep x nsample` cube of posterior predictive
#'   draws. Otherwise just summaries (median + 95% CI).
#' @param ... Unused.
#'
#' @return A list with components produced by [`dgtf_forecast()`]:
#'   typically `ymean`, `yvar`, posterior predictive bands, and,
#'   when `ypred_true` is supplied, scoring tables.
#'
#' @export
#' @examples
#' \dontrun{
#' fit <- dgtf(sim$y[1:180], mod, prior, method = "hva")
#' fc  <- dgtf_forecast_fit(fit, h = 20, nrep = 200,
#'                          ypred_true = sim$y[181:200])
#' }
dgtf_forecast_fit <- function(object,
                              h               = 1L,
                              nrep            = 100L,
                              ypred_true      = NULL,
                              only_last       = FALSE,
                              return_samples  = FALSE,
                              ...) {
    if (!inherits(object, "dgtf_fit"))
        stop("`object` must be a `dgtf_fit`.", call. = FALSE)
    if (h < 1L)
        stop("`h` must be a positive integer.", call. = FALSE)

    dgtf_forecast(
        output                = object$fit,
        model_settings        = as_settings(object$model),
        y                     = as.numeric(object$y),
        nrep                  = as.integer(nrep),
        k                     = as.integer(h),
        ypred_true            = if (is.null(ypred_true)) NULL
                                else as.numeric(ypred_true),
        only_evaluate_last_one = isTRUE(only_last),
        return_all_samples    = isTRUE(return_samples)
    )
}
