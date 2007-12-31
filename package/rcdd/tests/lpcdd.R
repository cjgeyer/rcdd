
 library(rcdd)

 ### optimal solution exists -- file samplelp1.ine in cddlib ###

 hrep <- scan()
 0  1  1  0  0
 0  1  0  1  0
 0  1  0  0  1
 0  1 -1  0  0
 0  1  0 -1  0
 0  1  0  0 -1

 hrep <- matrix(hrep, nrow = 6, byrow = TRUE)
 a <- c(1, 1, 1)

 lpcdd(hrep, a, minimize = FALSE)

 ### optimal solution exists -- file samplelp2.ine in cddlib ###

 hrep <- scan(what = character(0))
 0  0  1  1  0   0
 0  0  0  2  0   0
 1  3  0 -1  0   0
 1  9/2  0  0 -1  -1

 hrep <- matrix(hrep, nrow = 4, byrow = TRUE)
 a <- c("2", "3/5", "0", "0")

 lpcdd(hrep, a)

 ### primal inconsistent problem ###

 hrep <- scan(what = character(0))
 0  0  1  0
 0  0  0  1
 0  -2  -1  -1

 hrep <- matrix(hrep, nrow = 3, byrow = TRUE)
 a <- c("1", "1")

 lpcdd(hrep, a)

 lpcdd(q2d(hrep), q2d(a))

 ### dual inconsistent problem ###

 hrep <- scan(what = character(0))
 0  0  1  0
 0  0  0  1

 hrep <- matrix(hrep, nrow = 2, byrow = TRUE)
 a <- c("1", "1")

 lpcdd(hrep, a, minimize = FALSE)

 lpcdd(q2d(hrep), q2d(a), minimize = FALSE)
