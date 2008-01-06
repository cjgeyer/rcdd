
 library(rcdd)

 ##### H-representation #####

 hrep <- scan()
 0  0  1  1  0
 0  0 -1  0  0
 0  0  0 -1  0
 0  0  0  0 -1
 0  0 -1 -1 -1

 hrep <- matrix(hrep, ncol = 5, byrow = TRUE)

 redundant(d2q(hrep), representation = "H")

 redundant(hrep, representation = "H")

 ##### V-representation #####

 foo <- c(1, 0, -1)
 hrep <- cbind(0, 1, rep(foo, each = 9), rep(foo, each = 3), foo)
 print(hrep)

 redundant(d2q(hrep), representation = "V")

 redundant(hrep, representation = "V")

 ##### another V-representation #####

 hrep <- scan()
 0  0  1  0  0
 0  0  0  1  0
 0  0  0  0  1
 0  0 -1 -1 -1

 hrep <- matrix(hrep, ncol = 5, byrow = TRUE)

 redundant(d2q(hrep), representation = "V")

 redundant(hrep, representation = "V")

 ##### negative new position #####

 hrep <- scan()
 1  0  1  0  0
 1  0  0  1  0
 1  0  0  0  1
 1  0  0  1  0
 1  0  0  0  1
 1  0  0  0  1

 hrep <- matrix(hrep, ncol = 5, byrow = TRUE)

 redundant(hrep, representation = "V")

