test_that("control constructors return tagged objects", {
    expect_s3_class(vb_control(),    c("dgtf_vb_control", "dgtf_control"))
    expect_s3_class(mcmc_control(),  c("dgtf_mcmc_control", "dgtf_control"))
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
