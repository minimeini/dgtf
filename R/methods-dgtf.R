## S3 methods for `dgtf_fit` objects --------------------------------------
##
## The C++ engines write into `fit$fit` slightly different structures
## depending on the method. The methods below extract what they need
## defensively.


#' Latent state posterior draws
#'
#' Returns a matrix of posterior draws of psi (one column per draw),
#' or, for autoregressive HVA fits, a 3-D array of Theta draws
#' (state dimension x length(y) x draws). See Value.
#'
#' @param object A `dgtf_fit`.
#' @param what One of `"psi"` (the scalar `Rt`-driver) or `"Theta"`
#'   (full state matrix).
#' @return For `what = "psi"`, a numeric matrix of posterior draws of
#'   the latent state \eqn{\psi_t}, with one row per time point (that is,
#'   per element of the fitted series) and one column per posterior draw.
#'   For autoregressive models fitted by HVA, the rows instead hold the
#'   lag coefficients at the final time point. For `what = "Theta"`, a
#'   numeric 3-D array of draws of the full state vector, with dimensions
#'   state dimension x time points x posterior draws; it is stored only
#'   for autoregressive models fitted by HVA. An error is raised when the
#'   requested component is not stored for the fit.
#' @export
extract_state_draws <- function(object, what = c("psi", "Theta")) {
    what <- match.arg(what)
    out <- switch(what,
        psi   = object$fit$psi_stored, # n x nsample (nP x nsample for AR fits from HVA)
        Theta = object$fit$Theta_stored %||% object$fit$Theta
    )
    if (is.null(out)) {
        stop(sprintf("State component `%s` not found in fit output.", what),
            call. = FALSE
        )
    }
    out
}


# Returns a list:
#   inferred     : character, names of static params with posterior draws
#   draws        : named list, each element `nsample x p_dim` with
#                  column names like "W", "seas[1]", "seas[2]", "rho", ...
#   draws_matrix : `nsample x total_dim` matrix (concatenation of the above)
#   nsample      : integer
.dgtf_static_draws <- function(object) {
    f <- object$fit
    method <- object$method

    # 1. Unified `inferred` regardless of engine.
    inferred <- f$inferred
    if (is.null(inferred)) {
        flags <- grep("^infer_", names(f), value = TRUE)
        inferred <- sub(
            "^infer_", "",
            flags[vapply(f[flags], isTRUE, logical(1))]
        )
    }
    inferred <- inferred[inferred %in% names(f)] # only keep those with draws

    # 2. Determine canonical nsample from the control object (most reliable
    #    signal that survives both methods).
    nsample <- object$control$nsample %||% object$control$n_sample

    # 3. Normalize each draw matrix to nsample x p_dim.
    draws <- lapply(inferred, function(nm) {
        x <- f[[nm]]
        if (is.null(dim(x))) x <- matrix(x, ncol = 1)
        if (nrow(x) != nsample && ncol(x) == nsample) x <- t(x)
        if (ncol(x) == 1L) {
            colnames(x) <- nm
        } else {
            colnames(x) <- paste0(nm, "[", seq_len(ncol(x)), "]")
        }
        x
    })
    names(draws) <- inferred

    list(
        inferred = inferred,
        draws = draws,
        draws_matrix = do.call(cbind, draws),
        nsample = nsample
    )
}

# Internal: scale-free convergence diagnostics for an HVA marglik trace.
#
# `ml` is the per-iteration SMC log marginal likelihood estimate
# (`fit$marglik`). It is heavy-tailed, so every summary uses medians
# and MADs rather than means and SDs.
#
# Returns NULL if `ml` has fewer than 20 finite values; otherwise a list:
#   iterations      : number of finite values in the trace
#   tail_median     : median over the last `tail_window` iterations
#   tail_mad        : MAD over the same window (1 if zero or non-finite)
#   early_median    : median over an early window (after a 1% warm-up,
#                     capped at the first half of the trace)
#   early_mad       : MAD over the early window
#   mad_ratio       : tail_mad / early_mad (< 0.5 substantial, 0.5-0.8
#                     moderate, ~1 no stabilization)
#   tail_drift      : OLS slope over the last 25% of iterations times the
#                     window length, in units of tail_mad (|.| < 1 = flat)
#   stabilized_iter : first window start after which the rolling MAD stays
#                     below 2 x tail_mad for the rest of the trace; NA if never
#   last_spike_iter : last iteration more than `spike_k` x tail_mad away
#                     from tail_median; NA if none
#   spike_k         : spike threshold used
#   tail_window     : max(20, ceiling(tail_frac * n))
#   roll_win        : rolling-MAD window length, max(50, ceiling(0.05 * n))
#                     by default
.dgtf_marglik_convergence <- function(ml,
                                      tail_frac = 0.05,
                                      spike_k = 5,
                                      roll_win = NULL) {
    ml <- as.numeric(ml)
    ml <- ml[is.finite(ml)]
    n <- length(ml)
    if (n < 20L) {
        return(NULL)
    }

    win_tail <- max(20L, ceiling(tail_frac * n))
    tail_ix <- (n - win_tail + 1L):n

    # Robust late-window summary. The SMC log marginal estimator has
    # heavy tails (a few iters can sit hundreds of units away from the
    # bulk), so the mean and sd are dominated by outliers while the
    # median and MAD are not.
    tail_median <- stats::median(ml[tail_ix])
    tail_mad <- stats::mad(ml[tail_ix])
    if (!is.finite(tail_mad) || tail_mad == 0) tail_mad <- 1

    # Early window: skip a 1% warmup, take next `win_tail` iters but
    # cap at the first half of the trace.
    head_skip <- min(
        max(as.integer(ceiling(0.01 * n)), 5L),
        n %/% 4L
    )
    head_end <- min(head_skip + win_tail, n %/% 2L)
    head_ix <- (head_skip + 1L):head_end
    early_median <- stats::median(ml[head_ix])
    early_mad <- stats::mad(ml[head_ix])

    # MAD ratio is a "did variance shrink?" diagnostic.
    # < 0.5 = substantial; 0.5-0.8 = moderate; ~1 = no stabilization.
    mad_ratio <- if (is.finite(early_mad) && early_mad > 0) {
        tail_mad / early_mad
    } else {
        NA_real_
    }

    # Tail drift over last 25%, expressed as a multiple of tail MAD.
    # |tail_drift| < 1 means the trace drifted by less than one
    # tail-noise-unit over the entire last quarter, that is, flat.
    q_ix <- max(2L, ceiling(0.75 * n)):n
    tslope <- if (length(q_ix) >= 2L) {
        xx <- seq_along(q_ix)
        unname(stats::coef(stats::lm(ml[q_ix] ~ xx))[2L])
    } else {
        NA_real_
    }
    tail_drift <- if (is.finite(tslope)) {
        tslope * length(q_ix) / tail_mad
    } else {
        NA_real_
    }

    # Stabilization detector based on rolling MAD.
    # Fire only when rolling-MAD drops below 2x the final MAD AND
    # stays below for the rest of the trace (avoids matching a quiet
    # window between two outlier clusters).
    rwin <- if (is.null(roll_win)) {
        max(50L, ceiling(0.05 * n))
    } else {
        as.integer(roll_win)
    }
    stabilized_iter <- NA_integer_
    if (n >= rwin + 10L) {
        thr <- 2 * tail_mad
        roll_mad <- vapply(
            seq_len(n - rwin + 1L),
            function(i) stats::mad(ml[i:(i + rwin - 1L)]),
            numeric(1)
        )
        below <- roll_mad < thr
        all_below <- rev(cumprod(rev(below))) # 1 iff all later < thr
        min_ok <- max(rwin, as.integer(ceiling(0.10 * n)))
        cands <- which(all_below == 1L & seq_along(below) >= min_ok)
        if (length(cands)) stabilized_iter <- as.integer(cands[1L])
    }

    # Tail-aware companion to `stabilized_iter`. MAD is robust enough
    # that a single large spike inside a rolling window does not move it
    # much, so `stabilized_iter` reports when the *bulk* settled. This
    # variable instead reports the *latest* iter whose value is more
    # than `spike_k` final-MADs away from `tail_median`. The two
    # answer different questions:
    #   stabilized_iter  -> variational parameters have settled
    #   last_spike_iter  -> SMC estimator has stopped going degenerate
    spike_thr <- spike_k * tail_mad
    big_dev <- which(abs(ml - tail_median) > spike_thr)
    last_spike_iter <- if (length(big_dev)) {
        as.integer(max(big_dev))
    } else {
        NA_integer_
    }

    list(
        iterations      = n,
        tail_median     = tail_median,
        tail_mad        = tail_mad,
        early_median    = early_median,
        early_mad       = early_mad,
        mad_ratio       = mad_ratio,
        tail_drift      = tail_drift,
        stabilized_iter = stabilized_iter,
        last_spike_iter = last_spike_iter,
        spike_k         = spike_k,
        tail_window     = as.integer(win_tail),
        roll_win        = as.integer(rwin)
    )
}


#' @export
print.dgtf_fit <- function(x, ...) {
    cat("<dgtf_fit>\n")
    cat(sprintf("  method   : %s\n", x$method))
    cat(sprintf("  obs eq   : %s + link %s\n",
                x$model$obs$type, x$model$link$type))
    cat(sprintf("  state eq : %s\n", x$model$sys$type))
    cat(sprintf("  n        : %d\n", length(x$y)))
    cat(sprintf("  elapsed  : %.2f s\n", as.numeric(x$elapsed)))
    invisible(x)
}

# Counts every element of `y`, including the fixed initial value, to match
# `print()`, `summary()`, and the length of `fitted()` and `residuals()`.
#' @export
nobs.dgtf_fit <- function(object, ...) length(object$y)

# Internal: posterior fitted moments of y[t], from the per-draw intensity,
# mean and variance that `dgtf_posterior_predictive()` computes. Every
# element is aligned with `object$y`, with NA at the initial value.
.dgtf_fitted_moments <- function(object, level = 0.95) {
    if (isTRUE(object$model$zi))
        stop("Fitted values are not yet available for zero-inflated models.",
             call. = FALSE)
    if (identical(object$model$obs$type, "gaussian"))
        stop("Fitted values are not available for Gaussian observations.",
             call. = FALSE)

    pp <- posterior_predict(object, nrep = 0L, level = level)
    first_na <- function(x) {
        if (is.matrix(x)) x[1L, ] <- NA else x[1L] <- NA
        x
    }
    list(
        mean        = first_na(as.numeric(pp$fitted_mean)),
        var         = first_na(as.numeric(pp$fitted_var)),
        mean_qt     = first_na(as.matrix(pp$fitted_qt)),
        lambda      = first_na(as.numeric(pp$lambda_mean)),
        lambda_qt   = first_na(as.matrix(pp$lambda_qt))
    )
}

#' Fitted values, residuals and predictions of a fitted DGTF model
#'
#' `fitted()` returns the posterior mean of the expected count
#' \eqn{E[y_t]}, `residuals()` compares it with the observed series, and
#' `predict()` returns fitted values or forecasts, optionally with
#' intervals. In-sample results are aligned with `object$y`; the first
#' element is `NA` because the initial value of the series is not
#' modelled.
#'
#' @param object A `dgtf_fit` object.
#' @param type For `fitted()` and `predict()`, `"response"` (default)
#'   for the expected count \eqn{E[y_t]}, or `"intensity"` for the
#'   intensity \eqn{\lambda_t}. The two coincide except for
#'   [`obs_nbinom_p()`] observations. For `residuals()`, `"response"`
#'   (default) for \eqn{y_t - E[y_t]}, or `"pearson"` for the response
#'   residual divided by the posterior predictive standard deviation of
#'   \eqn{y_t}.
#' @param ... For `predict()` with `horizon > 0`, passed to
#'   [`dgtf_forecast_fit()`]; otherwise unused.
#'
#' @details The posterior predictive variance used by Pearson residuals
#'   is the posterior mean of the per-draw variance of \eqn{y_t} plus the
#'   posterior variance of the per-draw mean. Zero-inflated models and
#'   Gaussian observations are not yet supported.
#'
#' @return `fitted()` and `residuals()` return a numeric vector of the
#'   same length as `object$y`. See `predict()` below for its return
#'   value.
#' @seealso [`posterior_predict()`], [`dgtf_forecast_fit()`]
#' @export
fitted.dgtf_fit <- function(object, type = c("response", "intensity"), ...) {
    type <- match.arg(type)
    m <- .dgtf_fitted_moments(object)
    if (type == "response") m$mean else m$lambda
}

#' @rdname fitted.dgtf_fit
#' @export
residuals.dgtf_fit <- function(object,
                                type = c("response", "pearson"),
                                ...) {
    type <- match.arg(type)
    m <- .dgtf_fitted_moments(object)
    res <- as.numeric(object$y) - m$mean
    if (type == "pearson")
        res <- res / sqrt(pmax(m$var, .Machine$double.eps))
    res
}

#' @export
coef.dgtf_fit <- function(object, ...) {
    d <- .dgtf_static_draws(object)
    apply(d$draws_matrix, 2, stats::median)
}

#' @export
vcov.dgtf_fit <- function(object, ...) {
    d <- .dgtf_static_draws(object)
    stats::cov(d$draws_matrix)
}

#' @export
confint.dgtf_fit <- function(object, parm = NULL, level = 0.95, ...) {
    d <- .dgtf_static_draws(object)
    alpha <- 1 - level
    t(apply(d$draws_matrix, 2, stats::quantile,
        probs = c(alpha / 2, 1 - alpha / 2)
    ))
}

#' @param horizon Number of steps ahead. `0` (default) returns in-sample
#'   predictions; a positive value returns forecasts of the next
#'   `horizon` observations, computed by [`dgtf_forecast_fit()`].
#' @param interval `"none"` (default) for point predictions only,
#'   `"credible"` for a credible interval of the expected count (in
#'   sample only), or `"prediction"` for a posterior predictive interval
#'   of a new observation.
#' @param level Interval level. Forecasts (`horizon > 0`) always use
#'   `0.95`.
#' @param nrep Number of posterior predictive replicates per posterior
#'   draw, used for prediction intervals and forecasts.
#'
#' @return `predict()` returns a numeric vector when
#'   `interval = "none"`, and otherwise a matrix with columns `fit`,
#'   `lwr`, and `upr`. In sample, it has one entry or row per element of
#'   `object$y`. For forecasts, it has `horizon` entries or rows, and
#'   `fit` is the median of the forecast draws.
#' @rdname fitted.dgtf_fit
#' @export
predict.dgtf_fit <- function(object,
                              horizon  = 0L,
                              type     = c("response", "intensity"),
                              interval = c("none", "credible", "prediction"),
                              level    = 0.95,
                              nrep     = 100L,
                              ...) {
    type <- match.arg(type)
    interval <- match.arg(interval)
    if (!is.numeric(horizon) || length(horizon) != 1L || is.na(horizon) ||
        horizon < 0 || horizon != round(horizon))
        stop("`horizon` must be a single non-negative integer.", call. = FALSE)

    if (horizon > 0) {
        if (type != "response")
            stop("Forecasts are only available for `type = \"response\"`.",
                 call. = FALSE)
        if (interval == "credible")
            stop("Forecast intervals are predictive; use ",
                 "`interval = \"prediction\"`.", call. = FALSE)
        if (interval != "none" && !isTRUE(all.equal(level, 0.95)))
            warning("Forecast intervals are computed at the 95% level; ",
                    "`level` is ignored.", call. = FALSE)
        fc <- dgtf_forecast_fit(object,
                                h    = as.integer(horizon),
                                nrep = as.integer(nrep),
                                ...)
        q <- fc$quantiles
        if (interval == "none")
            return(unname(q[, "median"]))
        return(cbind(fit = q[, "median"], lwr = q[, "lower"],
                     upr = q[, "upper"]))
    }

    if (interval == "prediction" && type != "response")
        stop("Prediction intervals are only available for ",
             "`type = \"response\"`.", call. = FALSE)

    m <- .dgtf_fitted_moments(object, level = level)
    fit <- if (type == "response") m$mean else m$lambda
    if (interval == "none")
        return(fit)
    q <- switch(interval,
        credible   = if (type == "response") m$mean_qt else m$lambda_qt,
        prediction = {
            yq <- as.matrix(posterior_predict(object, nrep = nrep,
                                              level = level)$yhat)
            yq[1L, ] <- NA
            yq
        })
    cbind(fit = fit, lwr = q[, 1L], upr = q[, 3L])
}

#' Summarize a fitted DGTF model
#'
#' Summarizes inferred static parameters, latent states, runtime information,
#' and method-specific convergence diagnostics.
#'
#' @param object A `dgtf_fit` object.
#' @param truth Optional named numeric vector or list of true values,
#' named as in [param_recovery()] (e.g. `"W"`, `"rho"`, `"par1"`, `"seas[1]"`).
#' @param ... Additional arguments reserved for future methods.
#'
#' @details For MCMC fits, the table also has ess and rhat columns;
#'   rhat is single-chain split-Rhat (Gelman--Rubin
#'   on the two halves of the chain). It detects within-chain drift but
#'   cannot diagnose multi-modality. For a stronger check, run multiple
#'   independent fits and call `rhat(list(fit1, fit2, ...))`.
#'   HVA fits instead get a convergence summary of the SMC log marginal likelihood trace.
#'
#' @return An object of class `summary.dgtf_fit`.
#' @export
summary.dgtf_fit <- function(object, truth = NULL, ...) {
    f       <- object$fit
    method  <- object$method
    nsample <- object$control$nsample %||% object$control$n_sample

    # Static parameter draws, handling the orientation differences
    # between the HVA and MCMC engines.
    d <- tryCatch(.dgtf_static_draws(object), error = function(e) NULL)

    # Per-parameter table.
    param_table <- NULL
    if (!is.null(d) && !is.null(d$draws_matrix) && ncol(d$draws_matrix) > 0L) {
        dm <- d$draws_matrix
        param_table <- data.frame(
            parameter = colnames(dm),
            mean = apply(dm, 2, mean),
            median = apply(dm, 2, stats::median),
            sd = apply(dm, 2, stats::sd),
            q025 = apply(dm, 2, stats::quantile, probs = 0.025),
            q975 = apply(dm, 2, stats::quantile, probs = 0.975),
            row.names = NULL,
            stringsAsFactors = FALSE
        )

        if (identical(method, "mcmc")) {
            param_table$ess <- apply(dm, 2, ess)
            param_table$rhat <- apply(dm, 2, function(col) {
                tryCatch(rhat(matrix(col, ncol = 1L)),
                    error = function(e) NA_real_
                )
            })
        }

        # Truth-based recovery columns (simulation studies).
        if (!is.null(truth)) {
            rec <- tryCatch(
                param_recovery(object, truth, level = 0.95),
                error = function(e) {
                    warning("param_recovery() failed: ", conditionMessage(e),
                            call. = FALSE)
                    NULL
                })
            if (!is.null(rec) && nrow(rec) > 0L) {
                ix <- match(param_table$parameter, rec$parameter)
                param_table$truth    <- rec$truth[ix]
                param_table$bias     <- rec$bias[ix]
                param_table$mae      <- rec$mae[ix]
                param_table$rmse     <- rec$rmse[ix]
                param_table$crps     <- rec$crps[ix]
                param_table$coverage <- rec$coverage[ix]
                # Move `truth` next to `parameter` for readability.
                param_table <- param_table[, c(
                    "parameter", "truth", "mean", "median", "sd",
                    "q025", "q975",
                    setdiff(names(param_table),
                            c("parameter","truth","mean","median","sd",
                              "q025","q975"))
                )]
            }
        }
    } # param_table

    # Latent state summary: median band of psi (before the gain function).
    state_summary <- NULL
    if (!is.null(f$psi) && is.matrix(f$psi) && ncol(f$psi) >= 2L) {
        med <- as.numeric(f$psi[, 2])
        state_summary <- list(
            median_range = range(med, na.rm = TRUE),
            iqr_median   = stats::IQR(med, na.rm = TRUE)
        )
    }

    # Model components (best-effort; tolerant of partial specs).
    m <- object$model
    model_components <- list(
        obs         = m$obs$type      %||% NA_character_,
        link        = m$link$type     %||% NA_character_,
        sys         = m$sys$type      %||% NA_character_,
        gain        = m$gain$type     %||% NA_character_,
        lag         = m$lag$type      %||% NA_character_,
        lag_window  = m$lag$window    %||% NA_integer_,
        lag_residual_prob = m$lag$residual_prob %||% NA_real_,
        lag_truncation    = lag_truncation_length(m$lag),
        seasonality = if (!is.null(m$seasonality))
                          sprintf("period %s",
                                  m$seasonality$period %||% NA)
                      else NA_character_,
        zi          = isTRUE(m$zi)
    )

    # HVA-only convergence block (uses fit$marglik trace).
    convergence <- NULL
    if (identical(method, "hva") && !is.null(f$marglik)) {
        convergence <- .dgtf_marglik_convergence(f$marglik)
    }

    # MCMC-only HMC block (uses fit$hmc).
    hmc <- NULL
    if (identical(method, "mcmc") && !is.null(f$hmc)) {
        h  <- f$hmc
        ed <- h$diagnostics$energy_diff
        gn <- h$diagnostics$grad_norm
        hmc <- list(
            acceptance      = f$hmc_accept %||% h$acceptance_rate,
            leapfrog_step   = h$leapfrog_step_size,
            n_leapfrog      = h$n_leapfrog,
            energy_diff_max = if (length(ed)) max(abs(ed), na.rm = TRUE)         else NA_real_,
            energy_diff_q99 = if (length(ed)) stats::quantile(abs(ed), 0.99,
                                                              na.rm = TRUE)      else NA_real_,
            grad_norm_mean  = if (length(gn)) mean(gn, na.rm = TRUE)             else NA_real_,
            mass_diag       = h$mass_diag
        )
    }

    # MCMC-only disturbance MH block (uses fit$wt_accept).
    disturbance_mh <- NULL
    if (identical(method, "mcmc") && !is.null(f$wt_accept)) {
        wa <- as.numeric(f$wt_accept)
        disturbance_mh <- list(
            mean  = mean(wa, na.rm = TRUE),
            range = range(wa, na.rm = TRUE),
            n_low = sum(wa < 0.1, na.rm = TRUE)
        )
    }

    out <- list(
        method           = method,
        n                = length(object$y),
        nsample          = nsample,
        elapsed          = object$elapsed,
        elapsed_optimization = object$elapsed_optimization,
        elapsed_sampling     = object$elapsed_sampling,
        model_components = model_components,
        param_table      = param_table,
        state_summary    = state_summary,
        convergence      = convergence,
        hmc              = hmc,
        disturbance_mh   = disturbance_mh,
        ppc              = attr(object, "ppc"),
        error            = object$error
    )
    class(out) <- c("summary.dgtf_fit", "list")
    out
}

#' @param x A `summary.dgtf_fit` object.
#' @param digits Significant digits for most numeric output.
#' @param ... Additional arguments reserved for future methods.
#' @return `x`, invisibly.
#' @rdname summary.dgtf_fit
#' @export
print.summary.dgtf_fit <- function(x, digits = 3L, ...) {
    cat(sprintf("DGTF model fit (method = %s)\n",
                toupper(x$method  %||% "?")))
    cat(strrep("-", 50), "\n", sep = "")

    mc <- x$model_components
    if (!is.null(mc)) {
        if (!is.na(mc$obs))
            cat(sprintf("Observation : %s + link %s\n", mc$obs, mc$link))
        if (!is.na(mc$sys)) {
            lag_tag <- if (identical(mc$lag, "uniform")) {
                if (!is.na(mc$lag_window))
                    sprintf(" [window = %d]", as.integer(mc$lag_window))
                else ""
            } else if (!is.na(mc$lag_residual_prob)) {
                if (!is.na(mc$lag_truncation))
                    sprintf(" [residual = %g, L = %d]",
                            mc$lag_residual_prob,
                            as.integer(mc$lag_truncation))
                else
                    sprintf(" [residual = %g]", mc$lag_residual_prob)
            } else if (!is.na(mc$lag_window)) {
                sprintf(" [window = %d]", as.integer(mc$lag_window))
            } else ""
            cat(sprintf("System      : %s  (gain = %s, lag = %s%s)\n",
                        mc$sys, mc$gain, mc$lag, lag_tag))
        }
        if (!is.na(mc$seasonality))
            cat(sprintf("Seasonality : %s\n", mc$seasonality))
        if (isTRUE(mc$zi))
            cat("Zero-infl.  : enabled\n")
    }
    cat(sprintf("Data        : n = %d observations\n", x$n))
    if (!is.null(x$nsample))
        cat(sprintf("Posterior   : %s draws\n",
                    format(x$nsample, big.mark = ",")))
    if (!is.null(x$elapsed))
        cat(sprintf("Elapsed     : %.2f s\n", as.numeric(x$elapsed)))
    if (!is.null(x$elapsed_optimization)) {
        cat(sprintf(
            "  optim     : %.2f s\n",
            as.numeric(x$elapsed_optimization)
        ))
    }
    if (!is.null(x$elapsed_sampling)) {
        cat(sprintf(
            "  sampling  : %.2f s\n",
            as.numeric(x$elapsed_sampling)
        ))
    }

    if (!is.null(x$param_table) && nrow(x$param_table) > 0L) {
        cat("\nStatic parameters\n")
        pt <- x$param_table
        for (nm in setdiff(names(pt), "parameter")) {
            pt[[nm]] <- if (nm == "ess") {
                formatC(round(pt[[nm]]), format = "d", big.mark = ",")
            } else if (nm == "rhat") {
                formatC(pt[[nm]], digits = 3L, format = "f")
            } else if (nm == "coverage") {
                ifelse(is.na(pt[[nm]]), "",
                    ifelse(pt[[nm]], "yes", "no")
                )
            } else {
                formatC(pt[[nm]], digits = digits, format = "g")
            }
        }
        print(pt, row.names = FALSE)
    }

    if (!is.null(x$state_summary)) {
        cat("\nLatent state psi (median band)\n")
        cat(sprintf("  median range : [%s, %s]   IQR median : %s\n",
                    formatC(x$state_summary$median_range[1], digits = digits, format = "g"),
                    formatC(x$state_summary$median_range[2], digits = digits, format = "g"),
                    formatC(x$state_summary$iqr_median,      digits = digits, format = "g")))
    }

    if (!is.null(x$convergence)) {
        c0 <- x$convergence
        cat("\nConvergence (HVA marglik trace)\n")
        cat(sprintf("  iterations             : %d\n", c0$iterations))
        cat(sprintf("  late window (last %-3d) : median %s, MAD %s\n",
                    c0$tail_window,
                    formatC(c0$tail_median, digits = digits, format = "g"),
                    formatC(c0$tail_mad,    digits = digits, format = "g")))
        cat(sprintf("  early window           : median %s, MAD %s\n",
                    formatC(c0$early_median, digits = digits, format = "g"),
                    formatC(c0$early_mad,    digits = digits, format = "g")))
        mr <- c0$mad_ratio
        mr_tag <- if (!is.finite(mr)) ""
                  else if (mr < 0.5) "  [substantial stabilization]"
                  else if (mr < 0.8) "  [moderate stabilization]"
                  else               "  [no clear stabilization]"
        cat(sprintf("  MAD ratio (late/early) : %s%s\n",
                    formatC(mr, digits = 3L, format = "g"), mr_tag))
        flat <- is.finite(c0$tail_drift) && abs(c0$tail_drift) < 1
        cat(sprintf("  tail drift (in MADs)   : %s%s\n",
                    formatC(c0$tail_drift, digits = 2L, format = "f"),
                    if (flat) "  [flat]" else ""))
        cat(sprintf("  stabilized             : %s\n",
                    if (is.na(c0$stabilized_iter)) "no"
                    else sprintf("iter ~ %d  (rolling MAD < 2x final MAD)",
                                 c0$stabilized_iter)))
        cat(sprintf("  last large spike       : %s\n",
                    if (is.na(c0$last_spike_iter)) "none"
                    else sprintf("iter ~ %d  (|deviation| > %dx final MAD)",
                                 c0$last_spike_iter, c0$spike_k)))
    }

    if (!is.null(x$hmc)) {
        h <- x$hmc
        cat("\nHMC sampler\n")
        cat(sprintf("  acceptance rate      : %.3f\n", as.numeric(h$acceptance)))
        cat(sprintf("  leapfrog step        : %.4g\n", as.numeric(h$leapfrog_step)))
        cat(sprintf("  leapfrog steps       : %d\n",   as.integer(h$n_leapfrog)))
        cat(sprintf("  max |dH|             : %.4g\n", h$energy_diff_max))
        cat(sprintf("  99%% |dH|             : %.4g\n", h$energy_diff_q99))
        cat(sprintf("  mean ||grad log pi|| : %.4g\n", h$grad_norm_mean))
    }

    if (!is.null(x$disturbance_mh)) {
        d <- x$disturbance_mh
        cat("\nDisturbance MH\n")
        cat(sprintf("  mean acceptance      : %.3f\n", d$mean))
        cat(sprintf("  range acceptance     : [%.3f, %.3f]\n",
                    d$range[1], d$range[2]))
        cat(sprintf("  t with accept < 0.1  : %d\n", as.integer(d$n_low)))
    }

    if (!is.null(x$ppc) && inherits(x$ppc, "dgtf_ppc")) {
        cat("\nPosterior predictive check\n")
        print(summary(x$ppc))
    }

    invisible(x)
}
