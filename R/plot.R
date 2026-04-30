.okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                "#0072B2", "#D55E00", "#CC79A7", "#000000")
.truth_pattern <- "^(true|truth)$"
.default_ylab  <- function(what) switch(what,
    psi  = expression(psi[t]),
    Rt   = expression(R[t]),
    yhat = expression(hat(y)[t]))


.ts_targets <- c("Rt", "psi", "yhat")

# Soft replacement for match.arg(): time-series targets are an enum, but
# any other single string is treated as a static-parameter name and
# validated downstream against `.dgtf_static_draws(fit)$inferred`.
.what_kind <- function(what) {
    if (length(what) != 1L || !is.character(what))
        stop("`what` must be a single string.", call. = FALSE)
    if (what %in% .ts_targets) "ts" else "param"
}


# Mirror of src/GainFunc.hpp:43 (psi2hpsi). Logistic is registered C++-side
# but falls through to identity in the switch, so we omit it here.
.gain_fun <- function(name) {
    switch(name,
        identity    = identity,
        softplus    = function(z) log1p(exp(z)),
        exponential = exp,
        ramp        = function(z) pmax(z, .Machine$double.eps),
        stop(sprintf("Unknown gain '%s'; extend .gain_fun().", name),
             call. = FALSE))
}

# Row-wise quantiles of an n x nsample draws matrix at the central
# (1-level)/2, 0.5, 1-(1-level)/2 probabilities. Returns an n x 3 matrix.
.row_quantiles <- function(draws, level = 0.95) {
    probs <- c((1 - level) / 2, 0.5, 1 - (1 - level) / 2)
    out   <- t(apply(draws, 1L, stats::quantile,
                     probs = probs, na.rm = TRUE, names = FALSE))
    colnames(out) <- c("lower", "median", "upper")
    out
}

# Lazy PPC fallback (only triggered when Rt/yhat is requested but no PPC is
# attached and -- for Rt -- psi draws are unavailable).
.require_ppc <- function(fit, level) {
    warning("No PPC attached to fit; running posterior_predict() with ",
            "default settings (nrep = 100, no Rt reference). Attach one ",
            "with `attr(fit, 'ppc') <- posterior_predict(fit, ...)` to ",
            "avoid this on every plot call.", call. = FALSE)
    posterior_predict(fit, level = level)
}

# Normalize any input to a long-form band data.frame with columns:
#   time (numeric), lower, central, upper, name (character), has_ribbon (logical)
# Dispatches on class(x):
#   dgtf_fit  + what == "psi"             -> x$fit$psi          (T x 3)
#   dgtf_fit  + what %in% c("Rt","yhat")  -> attr(x, "ppc"), or
#                                            posterior_predict(x, level = level)
#                                            on the fly (with a warning)
#   dgtf_ppc  + what %in% c("Rt","yhat")  -> x$Rt or x$yhat     (T+1 x 3)
#   matrix (T x 3)                        -> use as-is (escape hatch for
#                                            external methods like EpiEstim)
#   numeric vector                        -> central only,
#                                            has_ribbon = FALSE  (truth)
.dgtf_band <- function(x, what, level = 0.95, time = NULL, name = NULL) {
    # Returns data.frame:
    #   time       numeric
    #   lower      numeric (NA when has_ribbon = FALSE)
    #   central    numeric (median)
    #   upper      numeric (NA when has_ribbon = FALSE)
    #   name       character (recycled)
    #   has_ribbon logical (recycled)

    # Branch 1: dgtf_fit
    #   psi  -> use the precomputed quantile matrix `fit$fit$psi`.
    #   Rt   -> apply the model's gain function elementwise to the
    #           `n x nsample` posterior draws of psi (via the public
    #           `extract_state_draws()` from R/draws.R:65) and take row-wise
    #           quantiles. No C++ posterior_predict round-trip needed.
    #   yhat -> requires posterior_predict (observation noise sampling), so
    #           use any attached PPC, otherwise warn + run on the fly.
    if (inherits(x, "dgtf_fit")) {
        if (what == "psi") {
            qm <- x$fit$psi                         # n x 3
        } else if (what == "Rt") {
            psi_draws <- tryCatch(extract_state_draws(x, "psi"),
                                  error = function(e) NULL)
            if (!is.null(psi_draws)) {
                gain_name <- x$model$gain$type %||% "identity"
                h         <- .gain_fun(gain_name)
                Rt_draws  <- h(psi_draws)           # n x nsample
                qm        <- .row_quantiles(Rt_draws, level = level)
            } else {
                # Fallback: fits without stored psi draws fall back to PPC.
                ppc <- attr(x, "ppc") %||% .require_ppc(x, level = level)
                qm  <- ppc$Rt
                if (is.null(qm))
                    stop("PPC has no `Rt` field. Re-run posterior_predict() ",
                         "with a non-NULL `Rt` argument or attach a PPC ",
                         "that contains it.", call. = FALSE)
            }
        } else if (what == "yhat") {
            ppc <- attr(x, "ppc") %||% .require_ppc(x, level = level)
            qm  <- ppc$yhat
            if (is.null(qm))
                stop("PPC has no `yhat` field. Re-run posterior_predict() ",
                     "with `nrep > 0`.", call. = FALSE)
        } else {
            stop(sprintf("Unknown `what = %s` for dgtf_fit.", what),
                 call. = FALSE)
        }
        return(.dgtf_band_from_matrix(qm, time = time, name = name,
                                      has_ribbon = TRUE))
    }

    # Branch 2: dgtf_ppc
    if (inherits(x, "dgtf_ppc")) {
        if (!what %in% c("Rt", "yhat"))
            stop(sprintf("Unknown `what = %s` for dgtf_ppc.", what),
                 call. = FALSE)
        qm <- x[[what]]
        if (is.null(qm))
            stop(sprintf("PPC object has no `%s` field.", what),
                 call. = FALSE)
        return(.dgtf_band_from_matrix(qm, time = time, name = name,
                                      has_ribbon = TRUE))
    }

    # Branch 3: bare matrix (T x 3)
    if (is.matrix(x) || (is.data.frame(x) && ncol(x) == 3L)) {
        return(.dgtf_band_from_matrix(as.matrix(x), time = time, name = name,
                                      has_ribbon = TRUE))
    }

    # Branch 4: numeric vector -> truth (no ribbon)
    if (is.numeric(x) && is.null(dim(x))) {
        n <- length(x)
        if (is.null(time)) {
            time <- seq(0, n - 1)
        } else {
            stopifnot(length(time) == n)
        }
        return(data.frame(
            time       = as.numeric(time),
            lower      = NA_real_,
            central    = as.numeric(x),
            upper      = NA_real_,
            name       = name %||% "Truth",
            has_ribbon = FALSE,
            stringsAsFactors = FALSE
        ))
    }

    stop("Unsupported input to .dgtf_band(): ",
         paste(class(x), collapse = "/"), call. = FALSE)
}

# Helper: turn a T x 3 quantile matrix into a long band data.frame.
.dgtf_band_from_matrix <- function(qm, time, name, has_ribbon) {
    stopifnot(is.matrix(qm), ncol(qm) == 3L)
    n <- nrow(qm)
    if (is.null(time)) {
        time <- seq(0, n - 1)        # matches plot_ts() convention in
                                     # R/diagnostics.R:109 (time = 0:(T-1))
    } else {
        # Allow callers to pass `time = band$time` for re-indexing truth onto
        # an existing band's axis. Recycle / truncate to length n if needed.
        if (length(time) != n) {
            if (length(time) > n) time <- time[seq_len(n)]
            else time <- c(time, seq(max(time) + 1L,
                                      length.out = n - length(time)))
        }
    }
    data.frame(
        time       = as.numeric(time),
        lower      = as.numeric(qm[, 1]),
        central    = as.numeric(qm[, 2]),
        upper      = as.numeric(qm[, 3]),
        name       = name %||% "Estimate",
        has_ribbon = has_ribbon,
        stringsAsFactors = FALSE
    )
}

# Deterministic palette assignment.
# names: character vector of method names (already filtered to drop entries
#        matching .truth_pattern -- truth is rendered separately).
# palette: NULL | character vector | named character vector. Named entries
#        pin specific names to specific hexes; unnamed entries fill from
#        the default generator in alphabetical order of remaining names.
# Returns a named character vector keyed by sort(unique(names)).
.dgtf_palette <- function(names, palette = NULL) {
    names <- unique(names)
    if (!length(names))
        return(stats::setNames(character(0), character(0)))
    n    <- length(names)
    auto <- if (n <= 8L) .okabe_ito[seq_len(n)]
            else grDevices::hcl.colors(n, palette = "Set 2")
    # Default assignment: alphabetical order over `names` -> Okabe-Ito sequence.
    out <- stats::setNames(auto, sort(names))[names]

    if (is.null(palette)) return(out)

    if (is.null(names(palette))) {
        # Unnamed vector: pin by position over `names`; if shorter than n,
        # the remaining entries keep the auto colors.
        m <- min(length(palette), n)
        out[seq_len(m)] <- palette[seq_len(m)]
        return(out)
    }

    # Named vector: pin the named entries; auto fills the rest. Unrecognized
    # names in `palette` are silently ignored.
    pinned <- intersect(names(palette), names)
    if (length(pinned)) out[pinned] <- palette[pinned]
    out
}

# Shared ggplot renderer. `bands` is the long data.frame from .dgtf_band
# (possibly rbind'd from many sources).
.plot_dgtf_band <- function(bands, palette = NULL, alpha = 0.5,
                            xlab = "Time", ylab = NULL, main = NULL,
                            theme = ggplot2::theme_minimal(),
                            legend.position = "right",
                            legend.name = "Method",
                            legend.nrow = NULL) {
    is_truth     <- grepl(.truth_pattern, bands$name, ignore.case = TRUE)
    truth        <- bands[ is_truth, , drop = FALSE]
    methods      <- bands[!is_truth, , drop = FALSE]
    method_names <- unique(methods$name)
    cols         <- .dgtf_palette(method_names, palette)

    p <- ggplot2::ggplot(methods, ggplot2::aes(x = time, y = central,
                                               group = name)) +
        ggplot2::geom_ribbon(
            data = methods[methods$has_ribbon, , drop = FALSE],
            ggplot2::aes(ymin = lower, ymax = upper, fill = name),
            alpha = alpha, na.rm = TRUE) +
        ggplot2::geom_line(ggplot2::aes(color = name), na.rm = TRUE) +
        ggplot2::scale_color_manual(name = legend.name, values = cols) +
        ggplot2::scale_fill_manual( name = legend.name, values = cols) +
        theme +
        ggplot2::labs(x = xlab, y = ylab, title = main) +
        ggplot2::theme(legend.position = legend.position)

    if (!is.null(legend.nrow))
        p <- p + ggplot2::guides(
            color = ggplot2::guide_legend(nrow = legend.nrow),
            fill  = ggplot2::guide_legend(nrow = legend.nrow))

    if (nrow(truth)) {
        truth_names <- unique(truth$name)
        p <- p + ggplot2::geom_line(
            data = truth,
            ggplot2::aes(x = time, y = central,
                         linetype = name, group = name),
            color = "black", inherit.aes = FALSE) +
            ggplot2::scale_linetype_manual(
                name   = NULL,
                values = stats::setNames(rep("dashed", length(truth_names)),
                                         truth_names))
    }
    p
}

# Returns a data.frame with columns:
#   name   character  -- method/source label (recycled)
#   index  integer    -- 1..p_dim; identifies which entry of a vector param
#   value  numeric    -- a single draw (or the truth value)
#   kind   character  -- "draws" or "truth"
#
# Dispatches on class(x):
#   dgtf_fit   -> pulls .dgtf_static_draws(x)$draws[[param]] (nsample x p_dim)
#                 and reshapes to long.
#   numeric    -> length 1 OR length p_dim; treated as truth.
#   matrix     -> nsample x p_dim, treated as draws (escape hatch for
#                 externally-fitted methods).
#
# `p_dim_hint` is consulted only for the numeric truth branch (we need to
# know how many indices to broadcast a length-1 truth across, OR validate
# that a length-k truth matches the expected p_dim). It is filled in by the
# compare wrapper after at least one fit has been resolved.
.dgtf_param_long <- function(x, param, name = NULL, p_dim_hint = NULL) {
    if (inherits(x, "dgtf_fit")) {
        sd <- .dgtf_static_draws(x)
        if (!param %in% sd$inferred)
            stop(sprintf(
                "Parameter '%s' is not inferred for this fit. Inferred: %s",
                param, paste(sd$inferred, collapse = ", ")), call. = FALSE)
        m <- sd$draws[[param]]                 # nsample x p_dim
        p <- ncol(m)
        return(data.frame(
            name  = name %||% "Estimate",
            index = rep(seq_len(p), each = nrow(m)),
            value = as.numeric(m),
            kind  = "draws",
            stringsAsFactors = FALSE
        ))
    }
    if (is.matrix(x)) {
        return(data.frame(
            name  = name %||% "Estimate",
            index = rep(seq_len(ncol(x)), each = nrow(x)),
            value = as.numeric(x),
            kind  = "draws",
            stringsAsFactors = FALSE
        ))
    }
    if (is.numeric(x) && is.null(dim(x))) {
        n <- length(x)
        if (!is.null(p_dim_hint) && n != 1L && n != p_dim_hint)
            stop(sprintf(
                "Truth for '%s' has length %d, expected 1 or %d.",
                param, n, p_dim_hint), call. = FALSE)
        idx <- if (n == 1L && !is.null(p_dim_hint)) seq_len(p_dim_hint)
               else seq_len(n)
        val <- if (n == 1L && !is.null(p_dim_hint)) rep(x, p_dim_hint) else x
        return(data.frame(
            name  = name %||% "Truth",
            index = idx,
            value = as.numeric(val),
            kind  = "truth",
            stringsAsFactors = FALSE
        ))
    }
    stop("Unsupported input to .dgtf_param_long(): ",
         paste(class(x), collapse = "/"), call. = FALSE)
}

.plot_dgtf_param_hist <- function(long, palette = NULL, alpha = 0.5,
                                  bins = 50, xlab = NULL, ylab = "Count",
                                  main = NULL,
                                  theme = ggplot2::theme_minimal(),
                                  legend.position = "right",
                                  legend.name = "Method") {
    is_truth     <- long$kind == "truth"
    truth        <- long[ is_truth, , drop = FALSE]
    methods      <- long[!is_truth, , drop = FALSE]
    method_names <- unique(methods$name)
    cols         <- .dgtf_palette(method_names, palette)

    p <- ggplot2::ggplot(methods, ggplot2::aes(x = value, fill = name)) +
        ggplot2::geom_histogram(position = "identity", alpha = alpha,
                                bins = bins, na.rm = TRUE) +
        ggplot2::scale_fill_manual(name = legend.name, values = cols) +
        theme +
        ggplot2::labs(x = xlab, y = ylab, title = main) +
        ggplot2::theme(legend.position = legend.position)

    if (nrow(truth))
        p <- p + ggplot2::geom_vline(
            xintercept = truth$value,
            color = "black", linetype = "dashed")
    p
}

.plot_dgtf_param_box <- function(long, palette = NULL, alpha = 0.5,
                                 xlab = "Index", ylab = "Value",
                                 main = NULL,
                                 theme = ggplot2::theme_minimal(),
                                 legend.position = "right",
                                 legend.name = "Method") {
    is_truth     <- long$kind == "truth"
    truth        <- long[ is_truth, , drop = FALSE]
    methods      <- long[!is_truth, , drop = FALSE]
    methods$index_f <- factor(methods$index)
    method_names <- unique(methods$name)
    cols         <- .dgtf_palette(method_names, palette)

    p <- ggplot2::ggplot(methods, ggplot2::aes(x = index_f, y = value,
                                               fill = name)) +
        ggplot2::geom_boxplot(alpha = alpha, outlier.size = 0.6) +
        ggplot2::scale_fill_manual(name = legend.name, values = cols) +
        theme +
        ggplot2::labs(x = xlab, y = ylab, title = main) +
        ggplot2::theme(legend.position = legend.position)

    if (nrow(truth)) {
        truth$index_f <- factor(truth$index,
                                levels = levels(methods$index_f))
        # Drop method/group from inheritance so the truth marker is
        # rendered once per index regardless of how many methods are
        # plotted. Use shape = 23 (filled diamond) to read clearly on
        # both light and dark themes.
        p <- p + ggplot2::geom_point(
            data = truth,
            ggplot2::aes(x = index_f, y = value),
            inherit.aes = FALSE,
            color = "black", fill = "black",
            shape = 23, size = 3)
    }
    p
}

.plot_dgtf_param_single <- function(x, param, truth = NULL,
                                    color = NULL, alpha = 0.5,
                                    bins = 50, ...) {
    sd <- .dgtf_static_draws(x)
    if (!param %in% sd$inferred)
        stop(sprintf(
            "Parameter '%s' is not inferred for this fit. Inferred: %s",
            param, paste(sd$inferred, collapse = ", ")), call. = FALSE)
    p_dim <- ncol(sd$draws[[param]])
    long  <- .dgtf_param_long(x, param, name = "Estimate")
    if (!is.null(truth))
        long <- rbind(long,
                      .dgtf_param_long(truth, param,
                                       name = "Truth",
                                       p_dim_hint = p_dim))
    palette <- if (!is.null(color)) c(Estimate = color) else NULL
    legend  <- if (is.null(truth)) "none" else "right"
    if (p_dim == 1L) {
        .plot_dgtf_param_hist(long, palette = palette, alpha = alpha,
                              bins = bins, xlab = param,
                              legend.position = legend, ...)
    } else {
        .plot_dgtf_param_box(long, palette = palette, alpha = alpha,
                             ylab = param, legend.position = legend, ...)
    }
}

.plot_dgtf_param_compare <- function(fits, param,
                                     palette = NULL, alpha = 0.5,
                                     bins = 50,
                                     legend.position = "top", ...) {
    stopifnot(is.list(fits), length(fits) >= 1L,
              !is.null(names(fits)), all(nzchar(names(fits))))

    # Resolve p_dim from the first dgtf_fit found; needed to broadcast /
    # validate any numeric-truth entries in the list.
    p_dim_hint <- NULL
    for (e in fits) {
        if (inherits(e, "dgtf_fit")) {
            sd <- .dgtf_static_draws(e)
            if (param %in% sd$inferred) {
                p_dim_hint <- ncol(sd$draws[[param]])
                break
            }
        } else if (is.matrix(e)) {
            p_dim_hint <- ncol(e); break
        }
    }
    if (is.null(p_dim_hint))
        stop(sprintf(
            "No fit in the input has '%s' inferred (and no draws matrix ",
            "supplied to seed p_dim).", param), call. = FALSE)

    long <- do.call(rbind, Map(
        function(x, nm) .dgtf_param_long(x, param, name = nm,
                                         p_dim_hint = p_dim_hint),
        fits, names(fits)))

    if (p_dim_hint == 1L)
        .plot_dgtf_param_hist(long, palette = palette, alpha = alpha,
                              bins = bins, xlab = param,
                              legend.position = legend.position, ...)
    else
        .plot_dgtf_param_box(long, palette = palette, alpha = alpha,
                             ylab = param,
                             legend.position = legend.position, ...)
}

#' @export
plot.dgtf_fit <- function(x,
                          what  = "Rt",
                          level = 0.95,
                          truth = NULL,
                          time  = NULL,
                          color = NULL,
                          alpha = 0.5,
                          xlab  = "Time",
                          ylab  = NULL,
                          main  = NULL,
                          theme = ggplot2::theme_minimal(),
                          ...) {
    kind <- .what_kind(what)
    if (kind == "param") {
        return(.plot_dgtf_param_single(x, what, truth = truth, ...))
    }

    band <- .dgtf_band(x, what = what, level = level, time = time,
                       name = "Estimate")
    if (!is.null(truth))
        band <- rbind(band,
                      .dgtf_band(truth, what = what,
                                 time = band$time, name = "Truth"))
    .plot_dgtf_band(
        band,
        palette = if (!is.null(color)) c(Estimate = color) else NULL,
        alpha   = alpha,
        xlab    = xlab,
        ylab    = ylab %||% .default_ylab(what),
        main    = main,
        theme   = theme,
        legend.position = if (is.null(truth)) "none" else "right"
    )
}

#' @export
dgtf_compare_plot <- function(fits,
                              what            = "Rt",
                              level           = 0.95,
                              time            = NULL,
                              palette         = NULL,
                              alpha           = 0.5,
                              legend.position = "right",
                              legend.name     = "Method",
                              legend.nrow     = NULL,
                              xlab            = "Time",
                              ylab            = NULL,
                              main            = NULL,
                              theme           = ggplot2::theme_minimal(),
                              ...) {
    stopifnot(is.list(fits), length(fits) >= 1L,
              !is.null(names(fits)), all(nzchar(names(fits))))
    kind <- .what_kind(what)
    if (kind == "param") {
        return(.plot_dgtf_param_compare(fits, what,
            palette = palette, alpha = alpha,
            legend.position = legend.position,
            legend.name = legend.name,
            theme = theme, ...
        ))
    }

    bands <- Map(
        function(x, nm) .dgtf_band(x, what = what, level = level,
                                   time = time, name = nm),
        fits, names(fits)
    )
    bands <- do.call(rbind, bands)
    .plot_dgtf_band(bands,
                    palette         = palette,
                    alpha           = alpha,
                    xlab            = xlab,
                    ylab            = ylab %||% .default_ylab(what),
                    main            = main,
                    theme           = theme,
                    legend.position = legend.position,
                    legend.name     = legend.name,
                    legend.nrow     = legend.nrow)
}

#' @export
plot.dgtf_ppc <- function(x, what = c("Rt", "yhat"),
                          truth = NULL, ...) {
    what <- match.arg(what)
    band <- .dgtf_band(x, what = what, name = "Estimate")
    if (!is.null(truth))
        band <- rbind(band,
                      .dgtf_band(truth, what = what,
                                 time = band$time, name = "Truth"))
    .plot_dgtf_band(band, ...,
                    legend.position = if (is.null(truth)) "none" else "right",
                    ylab = .default_ylab(what))
}