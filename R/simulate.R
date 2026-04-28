#' Simulate from a DGTF model
#'
#' Draws an observation series and the corresponding latent trajectories
#' from a [`dgtf_model()`].
#'
#' @param model A `dgtf_model` object.
#' @param ntime Length of the simulated series.
#' @param y0 Initial observation (default 0).
#' @param z Optional zero-inflation indicator vector.
#' @param seed Optional RNG seed.
#'
#' @return A list with components `y`, `lambda`, `psi`, `Theta`, `ft`,
#'   `nlag`, and (if zero-inflated) `prob`, `z`.
#' @export
#' @examples
#' mod <- dgtf_hawkes()
#' sim <- dgtf_simulate_model(mod, ntime = 200, seed = 1)
#' str(sim, max.level = 1)
dgtf_simulate_model <- function(model, ntime, y0 = 0, z = NULL, seed = NULL) {
    if (!inherits(model, "dgtf_model"))
        stop("`model` must be a `dgtf_model`.", call. = FALSE)
    if (!is.null(seed)) set.seed(seed)
    settings <- as_settings(model)
    dgtf_simulate(
        settings = settings,
        ntime    = as.integer(ntime),
        y0       = as.numeric(y0),
        z        = if (is.null(z)) NULL else as.numeric(z)
    )
}
