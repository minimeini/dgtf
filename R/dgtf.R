#' Fit a DGTF model
#'
#' The user-facing entry point for inference on a DGTF model. The
#' function wraps the C++ inference engines and returns an S3
#' `dgtf_fit` object supporting the standard methods (`print`,
#' `summary`, `coef`, `predict`, and others).
#'
#' @section Choosing a method:
#'
#' - `"hva"` (or its synonym `"vb"`): the hybrid variational
#'   approximation, which combines sequential Monte Carlo smoothing with
#'   variational updates. It is the faster of the two methods and is well
#'   suited to exploratory analysis and model screening. Tuned through
#'   [`vb_control()`].
#' - `"mcmc"`: a blocked sampler that updates the latent state
#'   disturbances by Metropolis--Hastings and the static parameters by
#'   Hamiltonian Monte Carlo. It targets the exact posterior and is the
#'   reference method for uncertainty quantification, at a higher
#'   computational cost. Tuned through [`mcmc_control()`].
#'
#' @section Inference rule:
#' Anything supplied a prior in `prior` is treated as unknown and
#' inferred. Anything left out of `prior` is treated as fixed at the
#' value supplied to [`dgtf_model()`].
#'
#' @param y Numeric vector of observed counts (or, in future versions,
#'   a matrix for multivariate observations).
#' @param model A [`dgtf_model`] object.
#' @param prior A [`dgtf_prior`] object. The default is an empty
#'   specification, under which all static parameters are held fixed.
#' @param method `"hva"` or `"mcmc"`. `"vb"` is accepted as a synonym
#'   for `"hva"`, and the fitted object records the method as `"hva"`
#'   either way.
#' @param control A method-specific control object. When `NULL`, the
#'   defaults of the selected method are used.
#' @param seed Optional integer seed.
#' @param verbose If `TRUE`, returns the lowered model and method
#'   settings *without* calling the C++ engine, which is useful for
#'   inspecting the arguments passed to the engines.
#'
#' @return An object of class `dgtf_fit`.
#' @export
#' @examples
#' \dontrun{
#' mod <- dgtf_hawkes()
#' sim <- dgtf_simulate_model(mod, ntime = 200, seed = 1)
#' fit <- dgtf(sim$y, mod,
#'             prior = dgtf_prior(W = inv_gamma(1, 1)),
#'             method = "hva",
#'             control = vb_control(iter = 200))
#' }
dgtf <- function(y,
                 model,
                 prior   = dgtf_prior(),
                 method  = c("hva", "mcmc", "vb"),
                 control = NULL,
                 seed      = NULL,
                 verbose   = FALSE) {

    if (!inherits(model, "dgtf_model"))
        stop("`model` must be a `dgtf_model` (see `?dgtf_model`).",
             call. = FALSE)
    if (!inherits(prior, "dgtf_prior"))
        stop("`prior` must be a `dgtf_prior` (see `?dgtf_prior`).",
             call. = FALSE)
    method <- match.arg(method)
    # "vb" is a synonym for the hybrid variational engine. Normalizing it
    # here keeps the rest of the function, and the `method` field recorded
    # on the fit, to the two canonical names.
    if (identical(method, "vb")) method <- "hva"
    if (!is.numeric(y))
        stop("`y` must be numeric.", call. = FALSE)
    if (!is.null(seed)) set.seed(seed)

    if (is.null(control)) {
        control <- switch(
            method,
            hva  = vb_control(),
            mcmc = mcmc_control()
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
        method_settings = method_settings
    )
    elapsed <- difftime(Sys.time(), t0, units = "secs")

    elapsed_optimization <- NULL
    elapsed_sampling <- NULL
    if (!is.null(raw$fit$elapsed_opt_us)) {
        elapsed_optimization <- as.difftime(
            raw$fit$elapsed_opt_us / 1e6,
            units = "secs"
        )
        raw$fit$elapsed_opt_us <- NULL
    }
    if (!is.null(raw$fit$elapsed_sample_us)) {
        elapsed_sampling <- as.difftime(
            raw$fit$elapsed_sample_us / 1e6,
            units = "secs"
        )
        raw$fit$elapsed_sample_us <- NULL
    }

    structure(
        list(
            fit     = raw$fit,
            y       = as.numeric(y),
            model   = model,
            prior   = prior,
            control = control,
            method  = method,
            elapsed = elapsed,
            elapsed_optimization = elapsed_optimization,
            elapsed_sampling     = elapsed_sampling
        ),
        class = "dgtf_fit"
    )
}
