## S3 methods for `dgtf_fit` objects --------------------------------------
##
## The C++ engines write into `fit$fit` slightly different structures
## depending on the method. The methods below extract what they need
## defensively.

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
    s <- object$fit$static
    if (is.null(s)) {
        # Fall back to whatever scalar parameters live at the top level
        out <- vapply(object$fit, function(x)
            if (is.numeric(x) && length(x) == 1L) x else NA_real_,
            numeric(1))
        return(out[!is.na(out)])
    }
    if (is.list(s))
        return(unlist(lapply(s, function(z)
            if (is.numeric(z) && length(z) == 1L) z
            else if (is.numeric(z)) stats::median(z)
            else NA_real_)))
    s
}

#' @export
vcov.dgtf_fit <- function(object, ...) {
    s <- object$fit$static_vcov
    if (is.null(s))
        stop("Variance-covariance matrix not available for this method.",
             call. = FALSE)
    s
}

#' @export
confint.dgtf_fit <- function(object, parm = NULL, level = 0.95, ...) {
    s <- object$fit$static_quantiles
    if (is.null(s))
        stop("Posterior quantiles for static parameters not available.",
             call. = FALSE)
    if (is.null(parm)) parm <- rownames(s) %||% seq_len(nrow(s))
    s[parm, , drop = FALSE]
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
        method  = object$method,
        n       = length(object$y),
        elapsed = object$elapsed,
        coef    = tryCatch(coef(object), error = function(e) NULL),
        ci      = tryCatch(confint(object), error = function(e) NULL),
        error   = object$error
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
