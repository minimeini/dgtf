# Multi-fit comparison table for `dgtf_fit` objects.
#
# Each row is one fit; columns are goodness-of-fit and runtime metrics
# extracted from the fit object directly or from a posterior_predict()
# result (cached on `attr(fit, "ppc")` if available, otherwise run on
# the fly).

#' Side-by-side comparison of fitted DGTF models
#'
#' Compares an arbitrary named list of fits on a common set of
#' goodness-of-fit and runtime metrics. Useful for selecting between
#' candidate models or comparing inference engines (HVA vs MCMC).
#'
#' For each fit, an attached PPC (`attr(fit, "ppc")`) is reused when
#' available, otherwise [`posterior_predict()`] is run on the fly
#' with `nrep`. To control the level / Rt reference, attach the PPC
#' yourself first or pass `force_recompute = TRUE` / `Rt_truth = ...`.
#'
#' @param fits Named list of `dgtf_fit` objects. Names appear in the
#'   `model` column of the returned data frame.
#' @param level Credible-interval level used for coverage / interval
#'   scoring (default 0.95). Used only for fresh PPC computations.
#' @param nrep Posterior-predictive replicate draws per posterior
#'   sample, forwarded to [`posterior_predict()`] when a fit has no
#'   cached PPC.
#' @param Rt_truth Optional reference `R_t` vector (length matches
#'   the input series, e.g. `softplus(sim$psi)`). When supplied, the
#'   table gains `mae_Rt`, `rmse_Rt`, `coverage_Rt`, and
#'   `interval_score_Rt` columns; PPCs are recomputed to bind this
#'   reference.
#' @param force_recompute If `TRUE`, ignore any attached PPC and run
#'   [`posterior_predict()`] fresh for every fit.
#'
#' @return A data frame of class `dgtf_comparison`, one row per fit.
#'   Columns:
#'   * `model` -- name from the input list.
#'   * `n_param` -- count of inferred static-parameter scalars
#'     (cross-engine, via [`.dgtf_static_draws()`]).
#'   * `log_marglik` -- HVA-only model-evidence proxy: stable tail
#'     median of the SMC log-marginal trace (estimate of
#'     `log p(y | gamma_final)`). `NA` for MCMC fits, which don't
#'     currently store this. **Note**: this is the marginal
#'     likelihood, not the pointwise log-likelihood -- treat it as a
#'     ranking signal for HVA fits, not as a frequentist
#'     log-likelihood.
#'   * `crps_y`, `chi_y`, `coverage_y`, `width_y`,
#'     `interval_score_y` -- in-sample posterior-predictive scoring
#'     for the observed counts `y` (see [`posterior_predict()`]).
#'   * `mae_Rt`, `rmse_Rt`, `coverage_Rt`, `interval_score_Rt` --
#'     present only when `Rt_truth` is supplied.
#'   * `elapsed_total`, `elapsed_optim`, `elapsed_sample` -- runtimes
#'     in seconds (NA when the fit doesn't expose the split).
#'
#' @examples
#' \dontrun{
#' fit_hva  <- dgtf(y, mod, prior, method = "hva")
#' fit_mcmc <- dgtf(y, mod, prior, method = "mcmc")
#'
#' # Cache PPCs once -- reused on every comparison call.
#' attr(fit_hva,  "ppc") <- posterior_predict(fit_hva,  nrep = 100)
#' attr(fit_mcmc, "ppc") <- posterior_predict(fit_mcmc, nrep = 100)
#'
#' dgtf_compare(list(HVA = fit_hva, MCMC = fit_mcmc))
#'
#' # In a sim study, pass a reference R_t to also score Rt recovery:
#' dgtf_compare(list(HVA = fit_hva, MCMC = fit_mcmc),
#'              Rt_truth = log1p(exp(sim$psi)))
#' }
#'
#' @export
# Multi-fit comparison table for `dgtf_fit` objects.
#
# Each row is one fit; columns are goodness-of-fit and runtime metrics
# extracted from the fit object directly or from a posterior_predict()
# result (cached on `attr(fit, "ppc")` if available, otherwise run on
# the fly).

#' Side-by-side comparison of fitted DGTF models
#'
#' Compares an arbitrary named list of fits on a common set of
#' goodness-of-fit and runtime metrics. Useful for selecting between
#' candidate models or comparing inference engines (HVA vs MCMC).
#'
#' For each fit, an attached PPC (`attr(fit, "ppc")`) is reused when
#' available, otherwise [`posterior_predict()`] is run on the fly
#' with `nrep`. To control the level / Rt reference, attach the PPC
#' yourself first or pass `force_recompute = TRUE` / `Rt_truth = ...`.
#'
#' @param fits Named list of `dgtf_fit` objects. Names appear in the
#'   `model` column of the returned data frame.
#' @param level Credible-interval level used for coverage / interval
#'   scoring (default 0.95). Used only for fresh PPC computations.
#' @param nrep Posterior-predictive replicate draws per posterior
#'   sample, forwarded to [`posterior_predict()`] when a fit has no
#'   cached PPC.
#' @param Rt_truth Optional reference `R_t` vector (length matches
#'   the input series, e.g. `softplus(sim$psi)`). When supplied, the
#'   table gains `mae_Rt`, `rmse_Rt`, `coverage_Rt`, and
#'   `interval_score_Rt` columns; PPCs are recomputed to bind this
#'   reference.
#' @param force_recompute If `TRUE`, ignore any attached PPC and run
#'   [`posterior_predict()`] fresh for every fit.
#'
#' @return A data frame of class `dgtf_comparison`, one row per fit.
#'   Columns:
#'   * `model` -- name from the input list.
#'   * `n_param` -- count of inferred static-parameter scalars
#'     (cross-engine, via [`.dgtf_static_draws()`]).
#'   * `log_marglik` -- HVA-only model-evidence proxy: stable tail
#'     median of the SMC log-marginal trace (estimate of
#'     `log p(y | gamma_final)`). `NA` for MCMC fits, which don't
#'     currently store this. **Note**: this is the marginal
#'     likelihood, not the pointwise log-likelihood -- treat it as a
#'     ranking signal for HVA fits, not as a frequentist
#'     log-likelihood.
#'   * `crps_y`, `chi_y`, `coverage_y`, `width_y`,
#'     `interval_score_y` -- in-sample posterior-predictive scoring
#'     for the observed counts `y` (see [`posterior_predict()`]).
#'   * `mae_Rt`, `rmse_Rt`, `coverage_Rt`, `interval_score_Rt` --
#'     present only when `Rt_truth` is supplied.
#'   * `elapsed_total`, `elapsed_optim`, `elapsed_sample` -- runtimes
#'     in seconds (NA when the fit doesn't expose the split).
#'
#' @examples
#' \dontrun{
#' fit_hva  <- dgtf(y, mod, prior, method = "hva")
#' fit_mcmc <- dgtf(y, mod, prior, method = "mcmc")
#'
#' # Cache PPCs once -- reused on every comparison call.
#' attr(fit_hva,  "ppc") <- posterior_predict(fit_hva,  nrep = 100)
#' attr(fit_mcmc, "ppc") <- posterior_predict(fit_mcmc, nrep = 100)
#'
#' dgtf_compare(list(HVA = fit_hva, MCMC = fit_mcmc))
#'
#' # In a sim study, pass a reference R_t to also score Rt recovery:
#' dgtf_compare(list(HVA = fit_hva, MCMC = fit_mcmc),
#'              Rt_truth = log1p(exp(sim$psi)))
#' }
#'
#' @export
dgtf_compare <- function(fits,
                         level = 0.95,
                         nrep = 100L,
                         Rt_truth = NULL,
                         truth = NULL,
                         force_recompute = FALSE) {
    if (!is.list(fits) || length(fits) == 0L ||
        is.null(names(fits)) || any(!nzchar(names(fits)))) {
        stop("`fits` must be a named, non-empty list of `dgtf_fit` ",
            "objects.",
            call. = FALSE
        )
    }
    if (!is.numeric(level) || length(level) != 1L ||
        level <= 0 || level >= 1) {
        stop("`level` must be a single number in (0, 1).", call. = FALSE)
    }

    rows <- lapply(seq_along(fits), function(i) {
        fit <- fits[[i]]
        name <- names(fits)[[i]]
        if (!inherits(fit, "dgtf_fit")) {
            stop(
                sprintf(
                    "Element `%s` is not a `dgtf_fit` (got class: %s).",
                    name, paste(class(fit), collapse = "/")
                ),
                call. = FALSE
            )
        }
        .dgtf_compare_row(fit,
            name = name, level = level,
            nrep = nrep, Rt_truth = Rt_truth,
            truth = truth,
            force_recompute = force_recompute
        )
    })

    df <- do.call(rbind, rows)
    rownames(df) <- NULL
    class(df) <- c("dgtf_comparison", "data.frame")
    attr(df, "level") <- level
    attr(df, "has_Rt_truth") <- !is.null(Rt_truth)
    attr(df, "has_truth") <- !is.null(truth)
    df
}


# Internal: extract one row of the comparison table.
.dgtf_compare_row <- function(fit, name, level, nrep,
                              Rt_truth, truth, force_recompute) {

    # 1. Static-parameter count.
    d <- tryCatch(.dgtf_static_draws(fit), error = function(e) NULL)
    n_param <- if (!is.null(d) && !is.null(d$draws_matrix))
                   as.integer(ncol(d$draws_matrix))
               else NA_integer_

    # 2. HVA-only model-evidence proxy.
    log_marglik <- NA_real_
    if (identical(fit$method, "hva") && !is.null(fit$fit$marglik)) {
        cv <- .dgtf_marglik_convergence(fit$fit$marglik)
        if (!is.null(cv)) log_marglik <- cv$tail_median
    }

    # 3. PPC: cached when safe, fresh otherwise. We always recompute
    #    when the caller passes `Rt_truth` (the cached PPC was likely
    #    built without it) or when force_recompute is set.
    use_cache <- !isTRUE(force_recompute) && is.null(Rt_truth)
    ppc <- if (use_cache) attr(fit, "ppc") else NULL
    if (is.null(ppc) || !inherits(ppc, "dgtf_ppc")) {
        ppc <- tryCatch(
            posterior_predict(fit, nrep = nrep, level = level,
                              Rt = Rt_truth),
            error = function(e) {
                warning(sprintf(
                    "posterior_predict() failed for `%s`: %s",
                    name, conditionMessage(e)), call. = FALSE)
                NULL
            })
    }

    pull <- function(x, default = NA_real_)
        if (is.null(x)) default else as.numeric(x)
    secs <- function(x) if (is.null(x)) NA_real_ else as.numeric(x)

    out <- data.frame(
        model              = name,
        n_param            = n_param,
        log_marglik        = log_marglik,
        crps_y             = pull(ppc$crps),
        chi_y              = pull(ppc$chi),
        coverage_y         = pull(ppc$coverage_yhat),
        width_y            = pull(ppc$width_yhat),
        interval_score_y   = pull(ppc$interval_score_yhat),
        elapsed_total      = secs(fit$elapsed),
        elapsed_optim      = secs(fit$elapsed_optimization),
        elapsed_sample     = secs(fit$elapsed_sampling),
        stringsAsFactors   = FALSE
    )
    if (!is.null(Rt_truth)) {
        out$mae_Rt            <- pull(ppc$mae_Rt)
        out$rmse_Rt           <- pull(ppc$rmse_Rt)
        out$coverage_Rt       <- pull(ppc$coverage_Rt)
        out$interval_score_Rt <- pull(ppc$interval_score_Rt)
    }

    # Static-parameter recovery (aggregate). The detailed per-parameter
    # view comes from `dgtf_compare_params()`. CRPS and bias are not
    # aggregated here because they mix scales across parameters
    # (W ~ 1e-3 vs rho ~ 30); coverage is a scale-free fraction and
    # reads well as an at-a-glance calibration signal.
    if (!is.null(truth)) {
        rec <- tryCatch(
            param_recovery(fit, truth, level = level),
            error = function(e) {
                warning(sprintf(
                    "param_recovery() failed for `%s`: %s",
                    name, conditionMessage(e)), call. = FALSE)
                NULL
                })
        out$n_param_truth <- if (is.null(rec)) NA_integer_ else as.integer(nrow(rec))
        out$param_coverage_rate <- if (is.null(rec) || !nrow(rec)) NA_real_ else mean(rec$coverage)
    }
    out
}


#' @export
print.dgtf_comparison <- function(x, digits = 3L, ...) {
    cat(sprintf("<dgtf_comparison: %d model%s, level = %.0f%%>\n",
                nrow(x), if (nrow(x) == 1L) "" else "s",
                100 * (attr(x, "level") %||% 0.95)))

    pt <- as.data.frame(unclass(x), stringsAsFactors = FALSE)

    for (nm in setdiff(names(pt), "model")) {
        v <- pt[[nm]]
        if (!is.numeric(v)) next
        pt[[nm]] <- if (nm == "n_param") {
            ifelse(is.na(v), "", formatC(v, format = "d"))
        } else if (grepl("^elapsed", nm)) {
            ifelse(is.na(v), "", sprintf("%.1fs", v))
        } else {
            formatC(v, digits = digits, format = "g")
        }
    }
    print(pt, row.names = FALSE)

    if (any(is.na(x$log_marglik)))
        cat("\nNote: `log_marglik` is the HVA SMC tail estimate of\n",
            "log p(y | gamma_final); MCMC fits have no equivalent stored\n",
            "and show NA. Pointwise log-likelihood is not currently\n",
            "exposed by either engine.\n", sep = "")
    invisible(x)
}

#' @export
print.dgtf_comparison <- function(x, digits = 3L, ...) {
    cat(sprintf("<dgtf_comparison: %d model%s, level = %.0f%%>\n",
                nrow(x), if (nrow(x) == 1L) "" else "s",
                100 * (attr(x, "level") %||% 0.95)))

    pt <- as.data.frame(unclass(x), stringsAsFactors = FALSE)

    for (nm in setdiff(names(pt), "model")) {
        v <- pt[[nm]]
        if (!is.numeric(v)) next
        pt[[nm]] <- if (nm == "n_param") {
            ifelse(is.na(v), "", formatC(v, format = "d"))
        } else if (nm == "n_param_truth") {
            ifelse(is.na(v), "", formatC(v, format = "d"))
        } else if (grepl("^elapsed", nm)) {
            ifelse(is.na(v), "", sprintf("%.1fs", v))
        } else {
            formatC(v, digits = digits, format = "g")
        }
    }
    print(pt, row.names = FALSE)

    if (any(is.na(x$log_marglik)))
        cat("\nNote: `log_marglik` is the HVA SMC tail estimate of\n",
            "log p(y | gamma_final); MCMC fits have no equivalent stored\n",
            "and show NA. Pointwise log-likelihood is not currently\n",
            "exposed by either engine.\n", sep = "")
    if (isTRUE(attr(x, "has_truth")))
        cat("\nUse `dgtf_compare_params(fits, truth)` for the\n", 
            "per-parameter recovery table.\n", sep = "")
        
    invisible(x)
}

#' Per-parameter posterior summaries across multiple fits
#'
#' Long-form table stacking per-parameter posterior summaries for every
#' fit in `fits`, with an added `model` column. One row per
#' (parameter, fit) pair. Use this for the detailed cross-fit view of
#' static-parameter estimation.
#'
#' When `truth` is supplied, the table additionally reports recovery
#' metrics (`bias`, `mae`, `rmse`, `crps`, `coverage`) from
#' [`param_recovery()`]. When `truth` is `NULL`, only the posterior
#' summaries (`mean`, `median`, `sd`, `q_lo`, `q_hi`) are reported,
#' giving a side-by-side view of the marginal posteriors without
#' assuming a known data-generating process.
#'
#' For an at-a-glance summary across fits, use [`dgtf_compare()`]
#' (with `truth` supplied for the recovery columns).
#'
#' @param fits Named list of `dgtf_fit` objects.
#' @param truth Optional named numeric vector or named list of true
#'   static parameter values (same convention as
#'   [`param_recovery()`]). When `NULL` (default), only posterior
#'   summaries are reported.
#' @param level Credible-interval level for the `q_lo`, `q_hi`, and
#'   (when `truth` is supplied) `coverage` columns (default 0.95).
#'
#' @return A data frame of class `dgtf_param_comparison` with columns
#'   `parameter`, `model`, `mean`, `median`, `sd`, `q_lo`, `q_hi`,
#'   and, when `truth` is supplied, `truth`, `bias`, `mae`, `rmse`,
#'   `crps`, `coverage`. Rows are ordered by parameter, then by the
#'   input order of `fits`.
#'
#' @export
dgtf_compare_params <- function(fits, truth = NULL, level = 0.95) {
    if (!is.list(fits) || length(fits) == 0L ||
        is.null(names(fits)) || any(!nzchar(names(fits))))
        stop("`fits` must be a named, non-empty list of `dgtf_fit` ",
             "objects.", call. = FALSE)
    if (!is.numeric(level) || length(level) != 1L ||
        level <= 0 || level >= 1)
        stop("`level` must be a single number in (0, 1).", call. = FALSE)

    rows <- lapply(seq_along(fits), function(i) {
        fit  <- fits[[i]]
        name <- names(fits)[[i]]
        if (!inherits(fit, "dgtf_fit"))
            stop(sprintf("Element `%s` is not a `dgtf_fit`.", name),
                 call. = FALSE)
        rec <- tryCatch(
            .compare_params_row(fit, truth, level),
            error = function(e) {
                warning(sprintf(
                    "compare_params row failed for `%s`: %s",
                    name, conditionMessage(e)), call. = FALSE)
                NULL
            })
        if (is.null(rec) || nrow(rec) == 0L) return(NULL)
        rec$model <- name
        rec
    })
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (!length(rows)) return(NULL)

    df <- do.call(rbind, rows)
    other <- setdiff(names(df), c("parameter", "model"))
    df <- df[, c("parameter", "model", other)]
    df$model <- factor(df$model, levels = names(fits))
    df <- df[order(df$parameter, df$model), , drop = FALSE]
    df$model <- as.character(df$model)
    rownames(df) <- NULL
    class(df) <- c("dgtf_param_comparison", "data.frame")
    attr(df, "level")     <- level
    attr(df, "has_truth") <- !is.null(truth)
    df
}

# Internal: produce one fit's worth of summary rows. Dispatches to
# `param_recovery()` when truth is provided, otherwise computes the
# posterior summaries directly from `.dgtf_static_draws()`.
.compare_params_row <- function(fit, truth, level) {
    if (!is.null(truth))
        return(param_recovery(fit, truth, level = level))

    d <- tryCatch(.dgtf_static_draws(fit), error = function(e) NULL)
    if (is.null(d) || is.null(d$draws_matrix) ||
        ncol(d$draws_matrix) == 0L)
        return(NULL)

    dm    <- d$draws_matrix
    alpha <- 1 - level
    qs    <- apply(dm, 2, stats::quantile,
                   probs = c(alpha / 2, 1 - alpha / 2))

    data.frame(
        parameter = colnames(dm),
        mean      = round(apply(dm, 2, mean),         3),
        median    = round(apply(dm, 2, stats::median), 3),
        sd        = round(apply(dm, 2, stats::sd),     3),
        q_lo      = round(qs[1, ],                     3),
        q_hi      = round(qs[2, ],                     3),
        row.names = NULL,
        stringsAsFactors = FALSE
    )
}

#' @export
print.dgtf_param_comparison <- function(x, digits = 3L, ...) {
    np  <- length(unique(x$parameter))
    nm  <- length(unique(x$model))
    lvl <- attr(x, "level") %||% 0.95
    has_truth <- isTRUE(attr(x, "has_truth"))

    cat(sprintf("<dgtf_param_comparison: %d parameter%s x %d model%s, level = %.0f%%>\n",
                np, if (np == 1L) "" else "s",
                nm, if (nm == 1L) "" else "s",
                100 * lvl))

    pt <- as.data.frame(unclass(x), stringsAsFactors = FALSE)

    # Convert the logical coverage column into a blank-named asterisk
    # column before the numeric formatting loop runs. Only present
    # when truth was supplied.
    if ("coverage" %in% names(pt)) {
        sig <- ifelse(is.na(pt$coverage), " ",
                      ifelse(pt$coverage, "*", " "))
        pt$coverage <- NULL
        pt[[" "]]   <- sig
    }

    for (col in setdiff(names(pt), c("parameter", "model", " "))) {
        pt[[col]] <- formatC(pt[[col]], digits = digits, format = "g")
    }

    print(pt, row.names = FALSE)
    if (has_truth)
        cat(sprintf("---\n'*' indicates the truth lies within the %g%% credible interval\n",
                    100 * lvl))
    invisible(x)
}
