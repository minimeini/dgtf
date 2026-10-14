## Bridges to the wider Bayesian R ecosystem ------------------------------
##
## These methods let users hand a `dgtf_fit` to packages like `posterior`
## (`as_draws_array`) and `bayesplot`. They are thin extractors over the
## posterior draws the C++ engine deposited in `fit$fit`.

#' Coerce a fitted DGTF model to a posterior::draws object
#'
#' Returns posterior samples of the static parameters, as a
#' `draws_array`, or as a `draws_matrix` when only one chain is
#' available. Latent state trajectories are held separately and are
#' retrieved with [`extract_state_draws()`], which keeps the returned
#' object small.
#'
#' @param x A `dgtf_fit`.
#' @param ... Unused.
#'
#' @return A `posterior::draws_array` if the **posterior** package is
#'   installed and samples are available; otherwise a plain matrix.
#' @export
as_draws <- function(x, ...) UseMethod("as_draws")

#' @rdname as_draws
#' @export
as_draws.dgtf_fit <- function(x, ...) {
    d <- .dgtf_static_draws(x)
    if (is.null(d$draws_matrix) || ncol(d$draws_matrix) == 0L) {
        stop("This fit has no posterior samples for static parameters.",
            call. = FALSE
        )
    }
    if (requireNamespace("posterior", quietly = TRUE)) {
        return(posterior::as_draws_matrix(d$draws_matrix))
    }
    d$draws_matrix
}

#' Latent state posterior draws
#'
#' Returns an `n x nsample` matrix for `psi`, or a 3-D array (`nP x (n+1) x nsample`) for `Theta`.
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
        psi   = object$fit$psi_stored, # n x nsample
        Theta = object$fit$Theta_stored %||% object$fit$Theta
    )
    if (is.null(out)) {
        stop(sprintf("State component `%s` not found in fit output.", what),
            call. = FALSE
        )
    }
    out
}