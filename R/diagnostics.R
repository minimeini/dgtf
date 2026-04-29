## Diagnostic helpers -----------------------------------------------------

#' Diagnostics for a DGTF fit
#'
#' Returns method-appropriate convergence and goodness-of-fit
#' diagnostics:
#'
#' - For `method = "vb"`: ELBO trajectory, final ELBO, change-per-iter.
#' - For `method = "mcmc"`: per-parameter \eqn{\hat{R}} and effective
#'   sample size (computed via the **posterior** package if installed).
#' - In-sample / out-of-sample loss tables (always).
#'
#' @param object A `dgtf_fit`.
#'
#' @return A named list of diagnostics; printed by its own method.
#' @export
dgtf_diagnostics <- function(object) {
    if (!inherits(object, "dgtf_fit"))
        stop("`object` must be a `dgtf_fit`.", call. = FALSE)
    out <- list(method = object$method, error = object$error)

    if (object$method == "vb") {
        elbo <- object$fit$elbo %||% object$fit$elbo_trace
        if (!is.null(elbo)) {
            out$elbo_final <- utils::tail(elbo, 1L)
            out$elbo_diff  <- diff(utils::tail(elbo, 2L))
            out$elbo_trace <- elbo
        }
    }

    if (object$method == "mcmc" &&
        requireNamespace("posterior", quietly = TRUE)) {
        s <- object$fit$samples %||% object$fit$static_samples
        if (!is.null(s)) {
            d <- posterior::as_draws_array(s)
            out$summary <- posterior::summarise_draws(
                d, "mean", "median", "sd",
                ess = posterior::ess_bulk,
                rhat = posterior::rhat
            )
        }
    }

    ppc <- attr(object, "ppc") %||% NULL
    if (!is.null(ppc) && inherits(ppc, "dgtf_ppc")) {
        out$ppc <- ppc
    }

    structure(out, class = "dgtf_diagnostics")
}

#' @export
print.dgtf_diagnostics <- function(x, ...) {
    cat(sprintf("<dgtf diagnostics: method = %s>\n", x$method))
    if (!is.null(x$elbo_final)) {
        cat(sprintf("  final ELBO     : %.4f\n", x$elbo_final))
        cat(sprintf("  last delta ELBO: %.2e\n", x$elbo_diff))
        cat(sprintf("  iterations     : %d\n", length(x$elbo_trace)))
    }
    if (!is.null(x$summary)) {
        cat("  posterior summary:\n")
        print(x$summary)
    }
    if (!is.null(x$error)) {
        if (!is.null(x$error$fitted))
            cat("  fitted error  :", format(x$error$fitted), "\n")
        if (!is.null(x$error$forecast))
            cat("  forecast error:", format(x$error$forecast), "\n")
    }

    if (!is.null(x$ppc)) {
        cat("  posterior predictive check:\n")
        print(x$ppc)
    }
    invisible(x)
}
