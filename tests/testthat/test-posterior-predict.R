test_that("dgtf_forecast_fit validates inputs", {
    expect_error(dgtf_forecast_fit("not a fit"),    "must be a `dgtf_fit`")
    fake_fit <- structure(list(fit = list(), model = NULL,
                               y = numeric(0)),
                          class = "dgtf_fit")
    expect_error(dgtf_forecast_fit(fake_fit, h = 0), "positive integer")
})

test_that("print.dgtf_forecast shows the scoring metrics", {
    fc <- structure(
        list(
            quantiles  = cbind(lower = c(1, 2), median = c(3, 4),
                               upper = c(5, 6)),
            level      = 0.95,
            horizon    = 2L,
            only_last  = FALSE,
            nsample    = 100L,
            ypred_true = c(3, 7),
            coverage   = 0.5,
            chi        = 0.25,
            crps       = 1.5
        ),
        class = c("dgtf_forecast", "list")
    )
    out <- capture.output(print(fc))
    expect_true(any(grepl("coverage\\s+:\\s+0.5", out)))
    expect_true(any(grepl("chi-sq\\s+:\\s+0.25", out)))
    expect_true(any(grepl("CRPS\\s+:\\s+1.5", out)))
})

test_that("print.dgtf_ppc works when optional metrics are absent", {
    # Without `Rt` and with `nrep = 0`, only these elements are present.
    pp <- structure(
        list(chi = 1.2, Rt = matrix(1, 3, 3), width_Rt = 0.3,
             level = 0.95),
        class = c("dgtf_ppc", "list")
    )
    out <- capture.output(print(pp))
    expect_true(any(grepl("chi-sq discrepancy", out)))
    expect_true(any(grepl("width \\(R_t\\)", out)))
    expect_false(any(grepl("MAE", out)))
})
