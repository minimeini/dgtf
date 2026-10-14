test_that("control constructors return tagged objects", {
    expect_s3_class(vb_control(),    c("dgtf_vb_control", "dgtf_control"))
    expect_s3_class(mcmc_control(),  c("dgtf_mcmc_control", "dgtf_control"))
    expect_s3_class(smc_control(),   c("dgtf_smc_control", "dgtf_control"))
    expect_s3_class(lba_control(),   c("dgtf_lba_control", "dgtf_control"))
})

test_that("lba_control validates discount factor", {
    expect_error(lba_control(0),   "must lie")
    expect_error(lba_control(1.1), "must lie")
})

test_that("smc_control matches resample arg", {
    expect_error(smc_control(resample = "weird"), "should be one of")
})

test_that("an unrecognized method name is rejected", {
    # The engine used to resolve method names through `algo_list[name]`,
    # which inserts a value-initialized entry for a missing key and so
    # returned the first enumerator, linear Bayes. A misspelled method
    # passed to the low-level interface then ran the wrong engine
    # silently rather than failing.
    expect_error(dgtf_default_algo_settings("not_a_method"),
                 "Unknown inference method")
})
