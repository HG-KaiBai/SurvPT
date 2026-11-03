pt_scalar <- function(value, name, lower, upper = Inf, integer = FALSE) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < lower || value > upper || (integer && value != floor(value)))
        stop(name, ' is outside its allowed range.')
}

pt_model <- function(object) {
    if (!inherits(object,'survPT') || is.null(object$Model) ||
        !identical(object$Model$version,1L))
        stop('Use a model returned by the current survPT(); older fitted lists must be refitted.')
    invisible(object)
}
