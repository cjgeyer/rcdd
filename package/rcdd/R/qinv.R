qinv <- function(x) {

    if (! (is.numeric(x) || is.character(x)))
        stop("'x' not numeric or character")

    if (is.numeric(x))
        x <- d2q(x)

    .Call("qo", x, as.integer(3), PACKAGE = "rcdd")
}