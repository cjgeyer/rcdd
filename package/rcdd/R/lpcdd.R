lpcdd <- function(hrep, objgrd, objcon, minimize = TRUE,
    solver = c("DualSimplex", "CrissCross")) {

    solver <- match.arg(solver)
    stopifnot(is.character(hrep) || is.numeric(hrep))
    stopifnot(is.character(objgrd) || is.numeric(objgrd))
    stopifnot(missing(objcon) || is.character(objcon) || is.numeric(objcon))
    stopifnot(is.character(hrep) == is.character(objgrd))
    stopifnot(is.logical(minimize))

    stopifnot(ncol(hrep) - 2 == length(objgrd))
    if (missing(objcon))
        objcon <- as(0, class(objgrd))
    stopifnot(length(objcon) == 1)
    stopifnot(is.character(objcon) == is.character(objgrd))
    stopifnot(length(minimize) == 1)

    validcdd(hrep, representation = "H")

    if (is.character(hrep)) {
        .Call(C_lpcdd, hrep, c(objcon, objgrd), minimize, solver)
    } else {
        storage.mode(hrep) <- "double"
        objgrd <- as.double(objgrd)
        objcon <- as.double(objcon)
        .Call(C_lpcdd_f, hrep, c(objcon, objgrd), minimize, solver)
    }
}
