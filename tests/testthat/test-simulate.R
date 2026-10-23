test_that("dgtf_simulate_model() returns Rt = h(psi)", {
    mod <- dgtf_hawkes(seasonality = seas_period(1, init = 2))
    sim <- dgtf_simulate_model(mod, ntime = 30, seed = 1)

    expect_false(is.null(sim$Rt))
    expect_identical(dim(sim$Rt), dim(sim$psi))
    expect_equal(as.numeric(sim$Rt), log1p(exp(as.numeric(sim$psi))))

    mod_exp <- dgtf_hawkes(gain = "exponential",
                           seasonality = seas_period(1, init = 2))
    sim_exp <- dgtf_simulate_model(mod_exp, ntime = 30, seed = 1)
    expect_equal(as.numeric(sim_exp$Rt), exp(as.numeric(sim_exp$psi)))
})

test_that("autoregressive simulations have no Rt", {
    sim <- dgtf_simulate_model(dgtf_nb_ar(p = 2L), ntime = 30, seed = 1)
    expect_null(sim$Rt)
    expect_false(is.null(sim$psi))
})

test_that("print() and plot() use the stored Rt", {
    mod <- dgtf_hawkes(seasonality = seas_period(1, init = 2))
    sim <- dgtf_simulate_model(mod, ntime = 30, seed = 1)
    n <- length(sim$y)

    sim$Rt[] <- 5
    out <- capture.output(print(sim))
    expect_true(any(grepl("R_t\\s+:\\s+range \\[\\s*5,\\s*5\\]", out)))
    expect_equal(plot(sim, what = "Rt")$data$value, rep(5, n))

    # Without a stored Rt, both fall back to applying the gain to psi.
    sim$Rt <- NULL
    expect_equal(plot(sim, what = "Rt")$data$value,
                 log1p(exp(as.numeric(sim$psi))))
})

test_that("as_dgtf_sim() keeps a supplied Rt and checks its length", {
    mod <- dgtf_hawkes(seasonality = seas_period(1, init = 2))
    draw <- dgtf_simulate_model(mod, ntime = 30, seed = 1)
    raw <- list(y = as.numeric(draw$y), psi = as.numeric(draw$psi),
                Rt = as.numeric(draw$Rt))

    sim <- suppressMessages(as_dgtf_sim(raw, model = mod))
    expect_equal(sim$Rt, raw$Rt)

    raw$Rt <- raw$Rt[-1]
    expect_error(suppressMessages(as_dgtf_sim(raw, model = mod)),
                 "`x\\$Rt` has length")
})
