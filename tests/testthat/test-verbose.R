test_that("the C++ engines print nothing unless verbose = TRUE", {
    mod <- dgtf_hawkes(seasonality = seas_period(1, init = 2))
    sim <- dgtf_simulate_model(mod, ntime = 25, seed = 7)
    prior <- dgtf_prior(W = inv_gamma(1, 1))

    # Collect everything written to stdout (Rprintf) and stderr (the
    # progress bars) while fitting, predicting, and forecasting.
    console_output <- function(verbose) {
        out <- capture.output(
            msg <- capture.output(
                {
                    fit <- dgtf(
                        sim$y, mod, prior,
                        control = vb_control(iter = 3, n_sample = 10,
                                             n_particle = 25,
                                             verbose = verbose),
                        seed = 1
                    )
                    posterior_predict(fit, nrep = 2, verbose = verbose)
                    dgtf_forecast_fit(fit, h = 2, nrep = 2,
                                      verbose = verbose)
                    invisible(NULL)
                },
                type = "message"
            ),
            type = "output"
        )
        c(out, msg)
    }

    expect_length(console_output(FALSE), 0L)
    expect_true(length(console_output(TRUE)) > 0L)
})
