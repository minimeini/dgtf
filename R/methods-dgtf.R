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
#' @export
summary.dgtf_fit <- function(object, ...) {
    f       <- object$fit
    method  <- object$method
    engine  <- object$control$method
    nsample <- object$control$nsample %||% object$control$n_sample

    # Static-parameter draws (handles HVA / MCMC orientation differences).
    d <- tryCatch(.dgtf_static_draws(object), error = function(e) NULL)

    # Per-parameter table.
    param_table <- NULL
    if (!is.null(d) && !is.null(d$draws_matrix) && ncol(d$draws_matrix) > 0L) {
        dm <- d$draws_matrix
        param_table <- data.frame(
            parameter = colnames(dm),
            mean      = apply(dm, 2, mean),
            median    = apply(dm, 2, stats::median),
            sd        = apply(dm, 2, stats::sd),
            q025      = apply(dm, 2, stats::quantile, probs = 0.025),
            q975      = apply(dm, 2, stats::quantile, probs = 0.975),
            ess_bulk  = NA_real_,
            ess_tail  = NA_real_,
            rhat      = NA_real_,
            row.names = NULL,
            stringsAsFactors = FALSE
        )
        if (requireNamespace("posterior", quietly = TRUE)) {
            ds <- tryCatch(
                posterior::summarise_draws(
                    posterior::as_draws_matrix(dm),
                    "ess_bulk", "ess_tail", "rhat"),
                error = function(e) NULL)
            if (!is.null(ds)) {
                idx <- match(param_table$parameter, ds$variable)
                param_table$ess_bulk <- ds$ess_bulk[idx]
                param_table$ess_tail <- ds$ess_tail[idx]
                param_table$rhat     <- ds$rhat[idx]
            }
        }
    }

    # Latent-state summary (post-link R_t median band).
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
        seasonality = if (!is.null(m$seasonality))
                          sprintf("period %s (in_state = %s)",
                                  m$seasonality$period %||% NA,
                                  isTRUE(m$seasonality$in_state))
                      else NA_character_,
        zi          = isTRUE(m$zi$enabled %||% FALSE)
    )

    # HVA-only convergence block (uses fit$marglik trace).
    convergence <- NULL
    if (identical(method, "hva") && !is.null(f$marglik)) {
        ml <- as.numeric(f$marglik)
        last_delta <- if (length(ml) >= 2L) abs(diff(utils::tail(ml, 2L))) else NA_real_
        # Plateau: first index where 50-iter rolling SD < 1% * |median(ml)|.
        plateau_iter <- NA_integer_
        if (length(ml) >= 100L) {
            win <- 50L
            thr <- 0.01 * abs(stats::median(ml))
            sds <- vapply(seq_len(length(ml) - win + 1L),
                          function(i) stats::sd(ml[i:(i + win - 1L)]),
                          numeric(1))
            hit <- which(sds < thr)
            if (length(hit)) plateau_iter <- hit[1L]
        }
        convergence <- list(
            final_log_marglik = utils::tail(ml, 1L),
            last_delta        = last_delta,
            iterations        = length(ml),
            plateau_iter      = plateau_iter
        )
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
        engine           = engine,
        n                = length(object$y),
        nsample          = nsample,
        elapsed          = object$elapsed,
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

#' @export
#' @export
print.summary.dgtf_fit <- function(x, digits = 3L, ...) {
    cat(sprintf("DGTF model fit (method = %s, engine = %s)\n",
                toupper(x$method  %||% "?"),
                toupper(x$engine  %||% "?")))
    cat(strrep("-", 50), "\n", sep = "")

    mc <- x$model_components
    if (!is.null(mc)) {
        if (!is.na(mc$obs))
            cat(sprintf("Observation : %s + link %s\n", mc$obs, mc$link))
        if (!is.na(mc$sys))
            cat(sprintf("System      : %s  (gain = %s, lag = %s%s)\n",
                        mc$sys, mc$gain, mc$lag,
                        if (!is.na(mc$lag_window))
                            sprintf(" [window = %d]", as.integer(mc$lag_window))
                        else ""))
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

    if (!is.null(x$param_table) && nrow(x$param_table) > 0L) {
        cat("\nStatic parameters\n")
        pt <- x$param_table
        for (nm in setdiff(names(pt), "parameter"))
            pt[[nm]] <- formatC(pt[[nm]], digits = digits, format = "g")
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
        cat("\nConvergence (HVA)\n")
        cat(sprintf("  final log p(y | gamma) : %s\n",
                    formatC(c0$final_log_marglik, digits = digits, format = "g")))
        cat(sprintf("  last delta             : %s\n",
                    formatC(c0$last_delta,        digits = digits, format = "g")))
        cat(sprintf("  iterations             : %d\n", c0$iterations))
        cat(sprintf("  plateau detected       : %s\n",
                    if (is.na(c0$plateau_iter)) "no"
                    else sprintf("iter ~ %d", c0$plateau_iter)))
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
        cat("\nPosterior predictive check (attached)\n")
        print(x$ppc)
    } else {
        cat("\n(no posterior predictive check attached - call ",
            "posterior_predict(fit) and assign with ",
            "`attr(fit, \"ppc\") <- ppc`)\n", sep = "")
    }

    if (!is.null(x$error)) {
        if (!is.null(x$error$fitted))
            cat("In-sample error    :", format(x$error$fitted), "\n")
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
