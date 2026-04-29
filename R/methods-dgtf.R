## S3 methods for `dgtf_fit` objects --------------------------------------
##
## The C++ engines write into `fit$fit` slightly different structures
## depending on the method. The methods below extract what they need
## defensively.


# Returns a list:
#   inferred     : character — names of static params with posterior draws
#   draws        : named list — each element is `nsample x p_dim` with
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

#' @export
nobs.dgtf_fit <- function(object, ...) length(object$y)

#' @export
fitted.dgtf_fit <- function(object, ...) {
    lam <- object$fit$lambda
    if (is.null(lam))
        stop("This fit does not expose `lambda` (fitted intensity).",
             call. = FALSE)
    if (is.matrix(lam) || (is.array(lam) && length(dim(lam)) > 1L)) {
        # Posterior summary: take the median row / column
        if (nrow(lam) == 3L) return(lam[2L, ])
        if (ncol(lam) == 3L) return(lam[, 2L])
        return(apply(lam, 2L, stats::median))
    }
    as.numeric(lam)
}

#' @export
residuals.dgtf_fit <- function(object,
                                type = c("response", "pearson"),
                                ...) {
    type <- match.arg(type)
    yhat <- fitted(object)
    res <- object$y - yhat
    if (type == "pearson") {
        # Approximate variance ~ lambda for Poisson; for NB use lambda(1+lambda/rho)
        v <- yhat
        if (identical(object$model$obs$type, "nbinom")) {
            rho <- coef(object)["rho"]
            if (!is.na(rho))
                v <- yhat * (1 + yhat / rho)
        }
        res <- res / sqrt(pmax(v, .Machine$double.eps))
    }
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
    alpha <- 1 - level
    t(apply(d$draws_matrix, 2, stats::quantile,
        probs = c(alpha / 2, 1 - alpha / 2)
    ))
}

#' @export
confint.dgtf_fit <- function(object, parm = NULL, level = 0.95, ...) {
    d <- .dgtf_static_draws(object)
    alpha <- 1 - level
    t(apply(d$draws_matrix, 2, stats::quantile,
        probs = c(alpha / 2, 1 - alpha / 2)
    ))
}

#' @export
predict.dgtf_fit <- function(object,
                              horizon = 0L,
                              nrep    = 100L,
                              type    = c("response", "intensity", "state"),
                              ...) {
    type <- match.arg(type)
    if (horizon > 0L) {
        return(dgtf_forecast_fit(object,
                                 h    = as.integer(horizon),
                                 nrep = as.integer(nrep),
                                 ...))
    }
    switch(type,
           response  = fitted(object),
           intensity = fitted(object),
           state     = object$fit$psi)
}

#' @export
logLik.dgtf_fit <- function(object, ...) {
    ll <- object$fit$loglik %||% object$fit$log_lik
    if (is.null(ll))
        stop("`logLik` not stored for this method.", call. = FALSE)
    val <- if (length(ll) == 1L) ll else sum(ll)
    structure(val,
              df = length(coef(object)),
              nobs = nobs(object),
              class = "logLik")
}

#' @export
summary.dgtf_fit <- function(object, ...) {
    out <- list(
  method            = "hva" | "mcmc" | ...,
  engine            = control$method,        # "vb" | "mcmc"
  n                 = length(y),
  nsample           = nsample,
  elapsed           = elapsed,
  model_components  = list(obs = "...", link = "...", sys = "...",
                           gain = "...", lag = "...", lag_window = ...,
                           seasonality = "period N (in_state = T/F)",
                           zi = TRUE/FALSE),
  param_table       = data.frame(parameter, mean, median, sd, q025, q975,
                                 ess_bulk, ess_tail, rhat),
  state_summary     = list(
    median_range = c(min, max),     # of post-gain R_t = h(psi)
    iqr_median   = ...
  ),
  convergence       = NULL | list(   # HVA only
    final_log_marglik   = ...,
    last_delta          = ...,
    iterations          = niter,
    plateau_iter        = first iter where rolling SD of marglik drops below
                          a heuristic threshold; NA if not detected
  ),
  hmc               = NULL | list(   # MCMC only
    acceptance        = hmc_accept,
    leapfrog_step     = hmc$leapfrog_step_size,
    n_leapfrog        = hmc$n_leapfrog,
    energy_diff_max   = max(abs(diagnostics$energy_diff)),
    energy_diff_q99   = quantile(abs(diagnostics$energy_diff), 0.99),
    grad_norm_mean    = mean(diagnostics$grad_norm)
  ),
  disturbance_mh    = NULL | list(   # MCMC only
    mean      = mean(wt_accept),
    range     = range(wt_accept),
    n_low     = sum(wt_accept < 0.1)
  ),
  ppc               = attr(object, "ppc")    # NULL unless attached
)
    class(out) <- "summary.dgtf_fit"
    out
}

#' @export
print.summary.dgtf_fit <- function(x, ...) {
    cat(sprintf("<dgtf_fit summary>\n  method = %s, n = %d, elapsed = %.2fs\n\n",
                x$method, x$n, as.numeric(x$elapsed)))
    if (!is.null(x$coef)) {
        cat("Static parameter point estimates:\n")
        print(x$coef)
        cat("\n")
    }
    if (!is.null(x$ci)) {
        cat("Posterior credible intervals:\n")
        print(x$ci)
        cat("\n")
    }
    if (!is.null(x$error)) {
        if (!is.null(x$error$fitted))
            cat("In-sample error:    ", format(x$error$fitted), "\n")
        if (!is.null(x$error$forecast))
            cat("Out-of-sample error:", format(x$error$forecast), "\n")
    }
    invisible(x)
}

#' @export
plot.dgtf_fit <- function(x, ...) {
    lam <- x$fit$lambda
    if (is.null(lam))
        stop("Nothing to plot: `lambda` missing from this fit.",
             call. = FALSE)
    if (is.matrix(lam)) {
        # Try to infer (lower, median, upper) layout
        if (nrow(lam) == 3L) {
            lo <- lam[1L, ]; me <- lam[2L, ]; hi <- lam[3L, ]
        } else if (ncol(lam) == 3L) {
            lo <- lam[, 1L]; me <- lam[, 2L]; hi <- lam[, 3L]
        } else {
            me <- apply(lam, 2L, stats::median)
            lo <- apply(lam, 2L, stats::quantile, probs = 0.025)
            hi <- apply(lam, 2L, stats::quantile, probs = 0.975)
        }
    } else {
        me <- as.numeric(lam); lo <- NULL; hi <- NULL
    }
    t  <- seq_along(me)
    yl <- range(c(x$y, lo, me, hi), na.rm = TRUE)
    plot(t, x$y, type = "h",
         xlab = "time", ylab = "y / intensity", ylim = yl,
         col  = "grey50",
         main = sprintf("dgtf fit (%s)", x$method))
    if (!is.null(lo)) {
        graphics::polygon(c(t, rev(t)), c(lo, rev(hi)),
                          col = grDevices::rgb(0.2, 0.4, 0.8, 0.3),
                          border = NA)
    }
    graphics::lines(t, me, col = "blue", lwd = 2)
    invisible(NULL)
}
