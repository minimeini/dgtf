# Posterior summaries are sample quantiles taken across draws. With 20
# draws, the 97.5% quantile falls exactly on the largest order statistic;
# arma::quantile() reads one element past the end of its work vector
# there. These tests fix the summaries to R's definition 5 quantiles at
# that boundary.

# A small HVA fit that stores exactly 20 draws. It is built once, on first
# use.
boundary_fit <- local({
    cache <- NULL
    function() {
        if (is.null(cache)) {
            mod <- dgtf_hawkes(
                meanlog = 1.35,
                sigma2 = 0.32,
                W = 0.01,
                seasonality = seas_period(1, init = 2)
            )
            sim <- suppressMessages(dgtf_simulate_model(mod, ntime = 40, seed = 7))
            prior <- dgtf_prior(W = inv_gamma(1, 1), rho = inv_gamma(1, 1))
            control <- vb_control(iter = 5, n_sample = 20, n_particle = 30)
            cache <<- dgtf(sim$y, mod, prior, method = "hva",
                           control = control, seed = 123)
        }
        cache
    }
})

# Row-wise definition 5 quantiles (Hyndman and Fan, 1996).
row_quantile_r <- function(draws, probs) {
    t(apply(draws, 1, stats::quantile, probs = probs, type = 5,
            names = FALSE))
}

interval_probs <- function(level) {
    alpha <- 1 - level
    c(alpha / 2, 0.5, 1 - alpha / 2)
}

test_that("state quantiles from 20 draws match definition 5 quantiles", {
    fit <- boundary_fit()
    draws <- fit$fit$psi_stored

    expect_identical(ncol(draws), 20L)
    expect_equal(unname(as.matrix(fit$fit$psi)),
                 row_quantile_r(draws, c(0.025, 0.5, 0.975)))
})

test_that("posterior_predict() quantiles match definition 5 at every level", {
    fit <- boundary_fit()
    # The model uses the softplus gain, so R[t] = log(1 + exp(psi[t])).
    Rt_draws <- log1p(exp(fit$fit$psi_stored))

    # With N draws, the upper probability reaches the largest order
    # statistic when 1 - level = 1 / N: level 0.95 for N = 20. The other
    # levels fall inside and outside the range of the order statistics.
    for (level in c(0.5, 0.8, 0.9, 0.95, 0.99)) {
        pp <- posterior_predict(fit, nrep = 0, level = level)
        expect_equal(unname(as.matrix(pp$Rt)),
                     row_quantile_r(Rt_draws, interval_probs(level)),
                     info = paste("level =", level))
    }
})

test_that("predictive quantiles from 20 draws are finite and ordered", {
    fit <- boundary_fit()

    # One replicate per stored draw gives 20 predictive draws.
    set.seed(1)
    pp <- posterior_predict(fit, nrep = 1)
    qt <- as.matrix(pp$yhat)[-1, ]
    expect_true(all(is.finite(qt)))
    expect_true(all(qt[, 1] <= qt[, 2] & qt[, 2] <= qt[, 3]))

    set.seed(1)
    fc <- dgtf_forecast_fit(fit, h = 2, nrep = 1)
    expect_identical(fc$nsample, 20L)
    expect_true(all(is.finite(fc$quantiles)))
    expect_true(all(fc$quantiles[, "lower"] <= fc$quantiles[, "median"] &
                    fc$quantiles[, "median"] <= fc$quantiles[, "upper"]))
})
