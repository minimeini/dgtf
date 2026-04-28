#' Validate parameter names for inference
#'
#' Checks whether the supplied parameter names are recognised inference
#' targets. Used internally by the fitting functions.
#'
#' @param param_infer Character vector of parameter names to inspect.
#' @param strict Logical. If `TRUE` (default) an unknown name raises an error;
#'   otherwise it issues a warning.
#' @return `TRUE` if every name is recognised, `FALSE` otherwise (and a
#'   warning is signalled when `strict = FALSE`).
#' @keywords internal
check_params <- function(param_infer, strict = TRUE) {
    param_list <- c("seas", "rho", "W", "par1", "par2", "zintercept", "zzcoef")
    out <- vapply(param_infer, function(x) x %in% param_list, logical(1))
    flag <- all(out)
    if (!flag) {
        msg <- "Unknown parameters."
        if (strict) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
    }
    flag
}
