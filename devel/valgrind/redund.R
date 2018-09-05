
 library("rcdd", lib.loc = "../../package/rcdd.Rcheck")

 set.seed(42)

 d <- 15
 ng <- 50

 for (i in 1:3) {

     x <- rnorm(d * ng)
     x <- matrix(x, ncol = d)
     vrep <- cbind(as.numeric(runif(ng) < 0.5), x)
     vrep <- cbind(as.numeric(runif(ng) < 0.5), vrep)
     vrep <- d2q(vrep)

     vout <- redundant(vrep, representation = "V")
     print(nrow(vout$output))

 }

 for (i in 1:3) {

     x <- rnorm(d * ng)
     x <- matrix(x, ncol = d)
     hrep <- cbind(rnorm(ng), x)
     hrep <- cbind(as.numeric(runif(ng) < 0.5), hrep)
     hrep <- d2q(hrep)

     hout <- redundant(hrep, representation = "H")
     print(nrow(hout$output))

 }

