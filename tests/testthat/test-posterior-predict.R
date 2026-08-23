test_that("dgtf_forecast_fit validates inputs", {
    expect_error(dgtf_forecast_fit("not a fit"),    "must be a `dgtf_fit`")
    fake_fit <- structure(list(fit = list(), model = NULL,
                               y = numeric(0)),
                          class = "dgtf_fit")
    expect_error(dgtf_forecast_fit(fake_fit, h = 0), "positive integer")
})
