#' Low-level C++ interface functions
#'
#' These functions are generated from the package's Rcpp interface and are
#' exported for advanced use and for the higher-level R wrappers. Most users
#' should call \code{dgtf()}, \code{dgtf_model()}, \code{posterior_predict()},
#' or \code{dgtf_forecast_fit()} instead.
#'
#' @param mu,sd2 Log-scale mean \eqn{\mu} and log-scale variance
#'   \eqn{\sigma^2} of a lognormal distribution.
#' @param m,s2 Mean and variance of the serial interval on the original
#'   scale.
#' @param nL Number of lags, which is the length of the returned PMF.
#' @param kappa,r Failure probability and number of successes of the
#'   negative binomial lag distribution.
#' @param method Name of the inference algorithm, either \code{"hva"}
#'   (or \code{"vb"}) or \code{"mcmc"}.
#' @param settings,model_settings,model_opts A model settings list, as
#'   returned by \code{as_settings()}.
#' @param ntime Length of the simulated series.
#' @param y0 Initial observation.
#' @param z Optional numeric 0/1 vector of structural zeros, used only by
#'   zero-inflated models.
#' @param y_in,y Numeric vector of observations.
#' @param method_settings Algorithm settings list, typically the result
#'   of \code{dgtf_default_algo_settings()} with some fields overridden.
#' @param X Optional numeric covariate matrix for the zero-inflation
#'   component, with one column per time point.
#' @param add_y0 If \code{TRUE}, a zero is prepended to \code{y_in} as the
#'   initial observation.
#' @param output The raw engine output, that is, the \code{fit} element of
#'   a \code{dgtf_infer()} result or of a \code{dgtf_fit} object.
#' @param nrep Number of replicate draws per posterior sample.
#' @param Rt Optional reference \eqn{R_t} trajectory of the same length
#'   as \code{y}.
#' @param level Credible interval level.
#' @param k Forecast horizon.
#' @param ypred_true Optional numeric vector of held-out observations
#'   used to score the forecast.
#' @param only_evaluate_last_one If \code{TRUE}, only the
#'   \code{k}-step-ahead forecast is evaluated.
#' @param return_all_samples If \code{TRUE}, the forecast draws are
#'   returned as well.
#' @param verbose If \code{TRUE}, a progress bar is shown.
#'
#' @return
#' \itemize{
#'   \item \code{lognorm2serial()}: a numeric vector of length 2 holding
#'     the mean and variance of the serial interval implied by a
#'     lognormal distribution with parameters \code{mu} and \code{sd2}.
#'   \item \code{serial2lognorm()}: a numeric vector of length 2 holding
#'     the lognormal parameters \eqn{(\mu, \sigma^2)} that match a serial
#'     interval with mean \code{m} and variance \code{s2}.
#'   \item \code{dlognorm0()}, \code{dnbinom0()}: a one-column numeric
#'     matrix with \code{nL} rows, the discretized lag PMF at lags
#'     \eqn{1, \dots, nL}.
#'   \item \code{dgtf_default_algo_settings()}: a named list of the default
#'     tuning options and prior blocks of the selected algorithm.
#'   \item \code{dgtf_default_model()}: a default model settings list with
#'     elements \code{model}, \code{param}, \code{season}, and \code{zero},
#'     in the format returned by \code{as_settings()}.
#'   \item \code{dgtf_simulate()}: a list with the simulated series
#'     \code{y} and the latent quantities \code{psi}, \code{Rt},
#'     \code{lambda}, and \code{ft} (each a \eqn{1 \times (ntime + 1)}
#'     matrix), the state matrix \code{Theta}, the lag truncation length
#'     \code{nlag}, and, for zero-inflated models, \code{prob} and
#'     \code{z}. \code{Rt} is \eqn{h(\psi_t)} and is omitted for
#'     autoregressive models.
#'   \item \code{dgtf_infer()}: a list with the single element \code{fit},
#'     the raw engine output (posterior draws and summaries) that
#'     \code{dgtf()} wraps into a \code{dgtf_fit} object.
#'   \item \code{dgtf_posterior_predictive()}: a list of posterior
#'     predictive summaries and scores, described in the Value section of
#'     \code{posterior_predict()}.
#'   \item \code{dgtf_forecast()}: a list with the forecast quantile
#'     matrix \code{yqt}; the scores \code{coverage}, \code{chi}, and
#'     \code{crps} together with \code{ypred_true} when held-out
#'     observations are supplied; and the draws \code{ypred_samples} when
#'     \code{return_all_samples = TRUE}.
#' }
#'
#' @name dgtf-cpp-lowlevel
#' @order 1
#' @keywords internal
NULL
