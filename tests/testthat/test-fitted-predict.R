# A small HVA fit shared by the tests in this file. It is built once,
# on first use.
small_fit <- local({
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
            control <- vb_control(iter = 5, n_sample = 30, n_particle = 30)
            cache <<- dgtf(sim$y, mod, prior, method = "hva",
                           control = control, seed = 123)
        }
        cache
    }
})

test_that("fitted() returns the posterior mean of E[y_t] aligned with y", {
    fit <- small_fit()
    f <- fitted(fit)

    expect_type(f, "double")
    expect_length(f, length(fit$y))
    expect_true(is.na(f[1]))
    expect_true(all(is.finite(f[-1])))
    expect_true(all(f[-1] >= 0))

    pp <- posterior_predict(fit, nrep = 0)
    expect_equal(f[-1], as.numeric(pp$fitted_mean)[-1])

    # For mean-dispersion negative binomial observations the expected
    # count is the intensity itself.
    expect_equal(fitted(fit, type = "intensity"), f)
})

test_that("posterior_predict() returns the fitted moments", {
    fit <- small_fit()
    pp <- posterior_predict(fit, nrep = 0)
    n <- length(fit$y)

    expect_length(pp$fitted_mean, n)
    expect_length(pp$fitted_var, n)
    expect_length(pp$lambda_mean, n)
    expect_identical(dim(as.matrix(pp$fitted_qt)), c(n, 3L))
    expect_identical(dim(as.matrix(pp$lambda_qt)), c(n, 3L))
    expect_true(all(as.numeric(pp$fitted_var)[-1] > 0))

    qt <- as.matrix(pp$fitted_qt)[-1, ]
    expect_true(all(qt[, 1] <= qt[, 2] & qt[, 2] <= qt[, 3]))
})

test_that("residuals() are y minus fitted, and Pearson residuals are scaled", {
    fit <- small_fit()
    r <- residuals(fit)
    expect_equal(r, as.numeric(fit$y) - fitted(fit))
    expect_true(is.na(r[1]))

    pp <- posterior_predict(fit, nrep = 0)
    rp <- residuals(fit, type = "pearson")
    expect_true(is.na(rp[1]))
    expect_equal(rp[-1], r[-1] / sqrt(as.numeric(pp$fitted_var)[-1]))
})

test_that("nobs() matches the length of fitted() and residuals()", {
    fit <- small_fit()
    expect_identical(nobs(fit), length(fit$y))
    expect_length(fitted(fit), nobs(fit))
    expect_length(residuals(fit), nobs(fit))
})

test_that("predict() returns in-sample fits and intervals", {
    fit <- small_fit()
    n <- length(fit$y)

    expect_equal(predict(fit), fitted(fit))
    expect_equal(predict(fit, type = "intensity"),
                 fitted(fit, type = "intensity"))

    cr <- predict(fit, interval = "credible")
    expect_identical(colnames(cr), c("fit", "lwr", "upr"))
    expect_identical(nrow(cr), n)
    expect_equal(unname(cr[, "fit"]), fitted(fit))
    expect_true(all(is.na(cr[1, ])))
    expect_true(all(cr[-1, "lwr"] <= cr[-1, "upr"]))

    set.seed(1)
    pr <- predict(fit, interval = "prediction", nrep = 20)
    expect_identical(colnames(pr), c("fit", "lwr", "upr"))
    expect_identical(nrow(pr), n)
    expect_true(all(is.na(pr[1, c("lwr", "upr")])))
    expect_true(all(pr[-1, "lwr"] <= pr[-1, "upr"]))
    # Prediction intervals include observation noise, so they are wider
    # than credible intervals for the mean.
    expect_gt(mean(pr[-1, "upr"] - pr[-1, "lwr"]),
              mean(cr[-1, "upr"] - cr[-1, "lwr"]))
})

test_that("predict() with horizon > 0 returns forecasts", {
    fit <- small_fit()

    set.seed(1)
    fc <- predict(fit, horizon = 3, interval = "prediction", nrep = 10)
    expect_identical(dim(fc), c(3L, 3L))
    expect_identical(colnames(fc), c("fit", "lwr", "upr"))
    expect_true(all(fc[, "lwr"] <= fc[, "fit"] & fc[, "fit"] <= fc[, "upr"]))

    set.seed(1)
    point <- predict(fit, horizon = 3, nrep = 10)
    expect_type(point, "double")
    expect_length(point, 3)
})

test_that("predict() validates its arguments", {
    fit <- small_fit()

    expect_error(predict(fit, horizon = -1), "non-negative integer")
    expect_error(predict(fit, horizon = 1.5), "non-negative integer")
    expect_error(predict(fit, horizon = c(1, 2)), "non-negative integer")
    expect_error(predict(fit, horizon = 2, type = "intensity"),
                 "only available for `type = \"response\"`")
    expect_error(predict(fit, horizon = 2, interval = "credible"),
                 "predictive")
    expect_error(predict(fit, type = "intensity", interval = "prediction"),
                 "only available for `type = \"response\"`")
    expect_warning(
        predict(fit, horizon = 2, interval = "prediction", level = 0.8,
                nrep = 5),
        "95%"
    )
})

test_that("unsupported observation models give a clear error", {
    fit <- small_fit()

    zi <- fit
    zi$model$zi <- TRUE
    expect_error(fitted(zi), "zero-inflated")
    expect_error(residuals(zi), "zero-inflated")
    expect_error(predict(zi), "zero-inflated")

    gauss <- fit
    gauss$model$obs$type <- "gaussian"
    expect_error(fitted(gauss), "Gaussian")
})

test_that("fitted() works for MCMC fits and autoregressive models", {
    skip_on_cran()
    fit <- small_fit()

    mfit <- dgtf(fit$y, fit$model, fit$prior, method = "mcmc",
                 control = mcmc_control(iter = 50, warmup = 20), seed = 1)
    f <- fitted(mfit)
    expect_length(f, length(mfit$y))
    expect_true(all(is.finite(f[-1])))
    expect_equal(residuals(mfit), as.numeric(mfit$y) - f)

    # Autoregressive models store Theta draws instead of psi draws.
    ar <- dgtf_nb_ar(p = 2L)
    afit <- dgtf(fit$y, ar, dgtf_prior(W = inv_gamma(1, 1), rho = inv_gamma(1, 1)),
                 method = "hva",
                 control = vb_control(iter = 5, n_sample = 30, n_particle = 30),
                 seed = 1)
    fa <- fitted(afit)
    expect_length(fa, length(afit$y))
    expect_true(all(is.finite(fa[-1])))
})
