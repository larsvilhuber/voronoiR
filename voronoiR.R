
install.packages(c("deldir","gpclib"),repos = getOption("repos")[[1]])
library(deldir)
x <- 1000*runif(10)
y <- 1000*runif(10)
vt <- deldir(x, y)
plot(x, y, axes=FALSE, ann=FALSE, pch=16)
plot(vt, wlines="tess", lty="solid", add=TRUE)
w <- runif(10, 1, 100)
library("gpclib")
unitSquare <- as(list(x=c(0, 0, 1000, 1000, 0),
                      y=c(0, 1000, 1000, 0, 0)),
                 "gpc.poly")
awvt <- awv(list(x=x, y=y), w, unitSquare)
