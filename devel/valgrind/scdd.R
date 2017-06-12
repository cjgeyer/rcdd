
 library("rcdd", lib.loc = "../../package/rcdd.Rcheck")

 set.seed(42)

 d <- 15
 ng <- 25 

 x <- rnorm(d * ng)
 x <- matrix(x, ncol = d)
 vrep <- cbind(as.numeric(runif(ng) < 0.5), x)
 vrep <- cbind(as.numeric(runif(ng) < 0.5), vrep)
 vrep <- d2q(vrep)

 hout <- scdd(vrep, representation = "V")
 vout <- scdd(hout$output)
 hout <- scdd(vout$output, adjacency = TRUE)
 vout <- scdd(hout$output, inputadjacency = TRUE)
 hout <- scdd(vout$output, inputincidence = TRUE)
 vout <- scdd(hout$output, incidence = TRUE)

