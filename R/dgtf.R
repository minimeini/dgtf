#' Fit a DGTF model
#'
#' The user-facing entry point for inference on a DGTF model. Wraps the
#' C++ engines exposed by [`dgtf_infer()`] and returns an S3 `dgtf_fit`
#' object with familiar methods (`print`, `summary`, `coef`, `predict`,
#' ...).
#'
#' @section Choosing a method:
#'
#' - `"vb"`: mean-field variational Bayes. The simplest variational
#'   family. Use [`vb_control()`] to tune.
#' - `"hva"`: hybrid variational approximation that combines SMC and
#'   variational updates (Tang et al., the engine featured in the
#'   manuscript). Fast; good for real-time applications. Also tuned via
#'   [`vb_control()`].
#' - `"mcmc"`: disturbance Metropolis-Hastings within Gibbs.
#'   Gold-standard uncertainty quantification; slower. Use
#'   [`mcmc_control()`] to tune.
#' - `"smc"`, `"ffbs"`, `"tfs"`, `"mcs"`, `"pl"`: sequential Monte Carlo
#'   smoothers. Use [`smc_control()`] to tune.
#' - `"linear_bayes"`: linear-Bayes approximation (closed-form, fastest).
#'   Use [`lba_control()`] to tune.
#'
#' @section Inference rule:
#' Anything supplied a prior in `prior` is treated as unknown and
#' inferred. Anything left out of `prior` is treated as fixed at the
#' value supplied to [`dgtf_model()`].
#'
#' @param y Numeric vector of observed counts (or, in future versions,
#'   a matrix for multivariate observations).
#' @param model A [`dgtf_model`] object.
#' @param prior A [`dgtf_prior`] object (default: empty -> all parameters fixed).
#' @param method One of `"vb"`, `"hva"`, `"mcmc"`, `"smc"`, `"ffbs"`,
#'   `"tfs"`, `"mcs"`, `"pl"`, `"linear_bayes"`.
#' @param control A method-specific control object. If `NULL`, sensible
#'   defaults are used.
#' @param loss_func Loss for the in-sample / forecast error tables.
#'   `"quadratic"` (RMSE) or `"absolute"` (MAE).
#' @param horizon Forecast horizon `k` for the forecast-error table.
#' @param seed Optional integer seed.
#' @param verbose If `TRUE`, returns the lowered model and method
#'   settings *without* calling the C++ engine. Useful for debugging
#'   what gets handed to the engines.
#'
#' @return An object of class `dgtf_fit`.
#' @export
#' @examples
#' \dontrun{
#' mod <- dgtf_hawkes()
#' set.seed(1)
#' sim <- dgtf_simulate(mod, ntime = 200)
#' fit <- dgtf(sim$y, mod,
#'             prior = dgtf_prior(W = inv_gamma(1, 1)),
#'             method = "vb",
#'             control = vb_control(iter = 200))
#' }
dgtf <- function(y,
                 model,
                 prior   = dgtf_prior(),
                 method  = c("vb", "hva", "mcmc", "smc", "ffbs", "tfs",
                             "mcs", "pl", "linear_bayes"),
                 control = NULL,
                 loss_func = c("quadratic", "absolute"),
                 horizon   = 1L,
                 seed      = NULL,
                 verbose   = FALSE) {

    if (!inherits(model, "dgtf_model"))
        stop("`model` must be a `dgtf_model` (see `?dgtf_model`).",
             call. = FALSE)
    if (!inherits(prior, "dgtf_prior"))
        stop("`prior` must be a `dgtf_prior` (see `?dgtf_prior`).",
             call. = FALSE)
    method    <- match.arg(method)
    loss_func <- match.arg(loss_func)
    if (!is.numeric(y))
        stop("`y` must be numeric.", call. = FALSE)
    if (!is.null(seed)) set.seed(seed)

    if (is.null(control)) {
        control <- switch(
            method,
            vb            = vb_control(),
            mcmc          = mcmc_control(),
            linear_bayes  = lba_control(),
            smc_control()
        )
    }
    if (!inherits(control, "dgtf_control"))
        stop("`control` must be a `dgtf_control` object.", call. = FALSE)

    model_settings  <- as_settings(model)
    method_settings <- control_to_method_settings(control, prior, method, model)

    if (isTRUE(verbose))
        return(list(model_settings = model_settings,
                    method_settings = method_settings,
                    method = method))

    t0 <- Sys.time()
    raw <- dgtf_infer(
        model_settings = model_settings,
        y_in           = as.numeric(y),
        method         = method,
        method_settings = method_settings,
        loss_func      = loss_func,
        k              = as.integer(horizon)
    )
    elapsed <- difftime(Sys.time(), t0, units = "secs")

    structure(
        list(
            fit     = raw$fit,
            error   = raw$error,
            y       = as.numeric(y),
            model   = model,
            prior   = prior,
            control = control,
            method  = method,
            elapsed = elapsed
        ),
        class = "dgtf_fit"
    )
}
