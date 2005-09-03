
 library(rcdd)

 d <- 5
 nface <- 20
 set.seed(42)
 a <- matrix(rnorm(nface * d), nface)

 foo <- makeH(a, 1)
 foo[1:3, ]

 ## for simplicity, round to 4 decimal places and convert to GMP rational

 bar <- round(1e4 * as.vector(foo))
 baz <- z2q(bar, rep(1e4, length(bar)))
 attributes(baz) <- attributes(foo)
 baz[1:3, ]
 all.equal(round(foo, 4), q2d(baz))

 out <- scdd(baz, incidence = TRUE, inputincidence = TRUE)
 names(out)
 attributes(out$input)
 attributes(out$output)
 class(out$incidence)
 length(out$incidence)
 class(out$inputincidence)
 length(out$inputincidence)

 unique(out$output[ , 1])
 unique(out$output[ , 2])

 v <- out$output[ , - c(1, 2)]
 v[1:3, ]

 v <- q2d(v)
 v[1:3, ]

 sally <- out$inputincidence
 sally[[length(sally)]] <- NULL
 length(sally)

 ##### 4-D faces (facets)

 inies <- maximal(sally)
 all(inies)
 face4d <- sally

 ##### 3-D faces

 sally <- face4d
 fred <- all.intersect(sally)
 length(fred)
 inies <- maximal(fred)
 length(inies)
 sum(inies)
 face3d <- fred[inies]

 ##### 3-D faces check

 outies <- fred[!inies]
 lout <- sapply(outies, length)
 all(lout == 0)

 ##### 2-D faces

 sally <- face3d
 fred <- all.intersect(sally)
 length(fred)
 inies <- maximal(fred)
 length(inies)
 sum(inies)
 face2d <- fred[inies]

 ##### 2-D faces check

 outies <- fred[!inies]
 lout <- sapply(outies, length)
 all(lout == 0)

 outies <- outies[lout > 0]
 length(outies)

 # each outie is redundant ???

 for (i in seq(along = outies)) {
     set1 <- outies[[i]]
     foo <- FALSE
     for (j in seq(along = face2d)) {
         set2 <- face2d[[j]]
         if (all(is.element(set1, set2))) {
             foo <- TRUE
         }
     }
     if (! foo)
         cat("outie set, number", i, "not redundant\n")
 }

 # each face2d is non-redundant ???

 for (i in seq(along = face2d)) {
     set1 <- face2d[[i]]
     for (j in seq(along = face2d)) {
         set2 <- face2d[[j]]
         if (i != j && all(is.element(set1, set2)))
             cat("face2d set number", i, "subset of face2d set number", j, "\n")
     }
 }

 ##### 1-D faces

 sally <- face2d
 fred <- all.intersect(sally)
 length(fred)
 inies <- maximal(fred)
 length(inies)
 sum(inies)
 face1d <- fred[inies]

 ##### 1-D faces check

 outies <- fred[!inies]
 lout <- sapply(outies, length)
 all(lout == 0)

 outies <- outies[lout > 0]
 length(outies)

 # each outie is redundant ???

 for (i in seq(along = outies)) {
     set1 <- outies[[i]]
     foo <- FALSE
     for (j in seq(along = face1d)) {
         set2 <- face1d[[j]]
         if (all(is.element(set1, set2))) {
             foo <- TRUE
         }
     }
     if (! foo)
         cat("outie set, number", i, "not redundant\n")
 }

 # each face1d is non-redundant ???

 for (i in seq(along = face2d)) {
     set1 <- face1d[[i]]
     for (j in seq(along = face2d)) {
         set2 <- face1d[[j]]
         if (i != j && all(is.element(set1, set2)))
             cat("face1d set number", i, "subset of face1d set number", j, "\n")
     }
 }

 ##### 0-D faces

 sally <- face1d
 fred <- all.intersect(sally)
 length(fred)
 inies <- maximal(fred)
 length(inies)
 sum(inies)
 face0d <- fred[inies]

 ##### 0-D faces check

 unique(sapply(face0d, length))

 all(sort(unlist(face0d)) == 1:nrow(out$output))

 ##### summary #####

 length(face4d)
 length(face3d)
 length(face2d)
 length(face1d)
 length(face0d)

 sort(unique(sapply(face4d, length)))
 sort(unique(sapply(face3d, length)))
 sort(unique(sapply(face2d, length)))
 sort(unique(sapply(face1d, length)))
 sort(unique(sapply(face0d, length)))


