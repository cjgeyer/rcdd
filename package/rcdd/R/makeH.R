makeH <- function(a1, b1, a2, b2, x = NULL) {

    if (missing(a1) != missing(b1))
        stop("must have either both or neither of 'a1' and 'b1' missing")
    if (missing(a2) != missing(b2))
        stop("must have either both or neither of 'a2' and 'b2' missing")
    if (missing(a1) && missing(a2))
        stop("all arguments missing")
    if ((! missing(a1)) && (! is.matrix(a1)))
        a1 <- matrix(a1, nrow = 1)
    if ((! missing(a2)) && (! is.matrix(a2)))
        a2 <- matrix(a2, nrow = 1)
    if ((! missing(a1)) && (! missing(a2)))
        if (ncol(a1) != ncol(a2))
            stop("ncol(a1) != ncol(a2)")

    fred <- attr(x, "representation")
    if ((! is.null(fred)) && (fred != "H"))
        stop("\"representation\" attribute of argument 'x' not \"H\"")
    if ((! is.null(x)) && (! is.matrix(x)))
        stop("argument 'x' must be NULL or matrix")

    foo <- NULL
    if (! missing(a1))
        foo <- cbind(0, as.vector(b1), - a1)

    bar <- NULL
    if (! missing(a2))
        bar <- cbind(1, as.vector(b2), - a2)

    mcol <- NULL
    if (! is.null(x))
        mcol <- ncol(x)
    if (! is.null(foo)) {
        if (is.null(mcol)) {
            mcol <- ncol(foo)
        } else {
            if (mcol != ncol(foo))
                stop("ncol(a1) + 2 != ncol(x)")
        }
    }
    if (! is.null(bar)) {
        if (is.null(mcol)) {
            mcol <- ncol(bar)
        } else {
            if (mcol != ncol(bar))
                stop("ncol(a2) + 2 != ncol(x) || ncol(a2) != ncol(a1)")
        }
    }
    baz <- rbind(x, bar, foo)
    dimnames(baz) <- NULL
    attr(baz, "representation") <- "H"
    validcdd(baz)
    return(baz)
}

addHeq <- function(a, b, x) makeH(a2 = a, b2 = b, x = x)

addHin <- function(a, b, x) makeH(a1 = a, b1 = b, x = x)

