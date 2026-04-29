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
#' @param level Credible-interval level for coverage / width / interval-score
#'   reporting (default `0.95`).
#'
#' @return An object of class `dgtf_ppc` (a list). Always present:
#'   * `chi` — averaged chi-square discrepancy of `y`
#'   * `Rt` — `(ntime + 1) x 3` matrix of `R_t` posterior quantiles at
#'     `(alpha/2, 0.5, 1 - alpha/2)` with `alpha = 1 - level`
#'   * `width_Rt` — mean width of the `R_t` credible bands
#'   * `level` — the credible-interval level used
#'
#'   When `nrep > 0`:
#'   * `crps` — averaged continuous ranked probability score for `y`
#'   * `yhat` — `(ntime + 1) x 3` quantile matrix of posterior-predictive `y`
#'   * `coverage_yhat`, `width_yhat`, `interval_score_yhat` — Winkler interval
#'     score and its components for `y`
#'
#'   When a reference `Rt` is supplied (simulation studies):
#'   * `mae_Rt`, `rmse_Rt`, `coverage_Rt`, `interval_score_Rt`
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
                                       Rt = NULL,
                                       level = 0.95,
                                       ...) {
    if (!is.numeric(level) || length(level) != 1L ||
        level <= 0 || level >= 1) {
        stop("`level` must be a single number in (0, 1).", call. = FALSE)
    }

    out <- dgtf_posterior_predictive(
        output     = object$fit,
        model_opts = as_settings(object$model),
        y          = as.numeric(object$y),
        nrep       = as.integer(nrep),
        Rt         = if (is.null(Rt)) NULL else as.numeric(Rt),
        level      = level
    )
    class(out) <- c("dgtf_ppc", "list")
    out
}


#' @export
print.dgtf_ppc <- function(x, digits = 3L, ...) {
    lvl <- x$level %||% 0.95
    cat(sprintf("<dgtf posterior predictive check, level = %.0f%%>\n",
                100 * lvl))
    rows <- list(
        c("chi-sq discrepancy",   x$chi),
        c("CRPS (y)",              x$crps),
        c("coverage (y)",          x$coverage_yhat),
        c("width (y)",             x$width_yhat),
        c("interval score (y)",    x$interval_score_yhat),
        c("MAE (R_t)",             x$mae_Rt),
        c("RMSE (R_t)",            x$rmse_Rt),
        c("coverage (R_t)",        x$coverage_Rt),
        c("width (R_t)",           x$width_Rt),
        c("interval score (R_t)",  x$interval_score_Rt)
    )
    rows <- Filter(function(r) !is.null(r[[2]]), rows)
    tab <- data.frame(
        metric = vapply(rows, `[[`, character(1), 1L),
        value  = vapply(rows, function(r)
            formatC(as.numeric(r[[2]]), digits = digits, format = "f"),
            character(1))
    )
    print(tab, row.names = FALSE)
    invisible(x)
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
