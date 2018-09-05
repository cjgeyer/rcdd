
 library("rcdd", lib.loc = "../../package/rcdd.Rcheck")

 set.seed(42)

 d <- 15
 ng <- 15

 for (i in 1:3) {

     x <- rnorm(d * ng)
     x <- matrix(x, ncol = d)
     hrep <- cbind(rnorm(ng), x)
     hrep <- cbind(as.numeric(runif(ng) < 0.5), hrep)
     hrep <- d2q(hrep)

     lout <- lpcdd(hrep, rep("1",d))
     print(lout$solution.type)

 }

