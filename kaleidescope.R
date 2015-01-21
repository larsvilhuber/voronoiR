
#  Copyright (C) 2012 Paul Murrell
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.gnu.org/licenses/gpl.txt

library(gpclib) # Non-commercial only!

debugPlotGen <- function() {
    plotStarted <- FALSE
    x <- 1
    prevYoverall <- NA
    prevYmax <- NA
    
    function(yoverall, ymax) {
        if (!plotStarted) {
            plot.new()
            plot.window(c(0, 200), c(0, .2))
            abline(h=c(0, .001), col="grey", lty=c("solid", "dashed"))
            axis(1)
            axis(2)
            box()
            plotStarted <<- TRUE
        }
        segments(x, prevYoverall, x + 1, yoverall, col="grey")
        segments(x, prevYmax, x + 1, ymax, col="black")
        x <<- x + 1
        prevYoverall <<- yoverall
        prevYmax <<- ymax
    }
}

debugPlot <- debugPlotGen()

prevError <- 0

cellError <- function(a, target) {
    normA <- a/sum(a)
    diff <- abs(normA - target)
    max(diff)
}

breaking <- function(a, target, max=TRUE, debug=FALSE, dplot=FALSE) {
    if (max) {
        # Stop when largest individual cell error is less than .1%
        # (the default)
        err <- cellError(a, target)
        if (debug)
            cat(paste("Difference: ", err,
                      " (", abs(err - prevError), ")",
                      "\n", sep=""))
        stopping <- err < .001
        prevError <<- err
        if (!debug && dplot) {
            debugPlot(err, err)
        }
    } else {
        normA <- a/sum(a)
        diff <- abs(normA - target)
        # Stop when *change* in *total* cell error is tiny
        # (i.e., we are not improving the solution)
        if (debug)
            cat(paste("Difference: ", round(sum(diff), 2),
                      " (", round(abs(sum(diff) - prevError), 3), ")",
                      "\n", sep=""))
        stopping <- abs(sum(diff) - prevError) < .0001
        prevError <<- sum(diff)
        if (!debug && dplot) {
            debugPlot(sum(diff)/sum(target), max(diff/target))
        }
    }
    stopping
}

OLDadjustWeights <- function(w, a, target) {
    # The choice of 'eps' is important
    # Either too large (e.g., 10) or too small (e.g., .001)
    # can result in failure to converge
    eps <- scale/10000
    # Watch out for zero-area cells (set to target value)
    # a <- ifelse(a == 0, target/sum(target)*sum(a), a)
    # Watch out for zero-area cells (set to small value)
    a <- ifelse(a == 0, .01*sum(a), a)
    normA <- a/sum(a)
    w <- ifelse(abs(w) < eps,
                ifelse(w < 0, -eps, eps),
                w)
    w + abs(w)*((target - normA)/target)
}

# Instead of adjusting by a multiple of current weight
# adjust by multiple of average absolute weights
# This avoids problem of getting stuck at a tiny weight
# (and stabilizes the algorithm generally)
adjustWeights <- function(w, a, target) {
    # Watch out for zero-area cells (set to small value)
    a <- ifelse(a == 0, .01*sum(a), a)
    normA <- a/sum(a)
    w + mean(abs(w))*((target - normA)/target)
}

# Use centroid finding code in 'soiltexture' for now
# May be able to find a better source later (?)
library(soiltexture)
shiftSites <- function(s, k) {
    newSites <- mapply(function(poly, x, y) {
                           # Handle empty polygons
                           if (length(poly@pts)) {
                               pts <- getpts(poly)
                               TT.polygon.centroids(pts$x,
                                                    pts$y)
                           } else {
                               c(x=x, y=y)
                           }
                       },
                       k, as.list(s$x), as.list(s$y),
                       SIMPLIFY=FALSE)
    list(x=sapply(newSites, "[", "x"),
         y=sapply(newSites, "[", "y"))
}

OLDshiftWeights <- function(s, w) {
    n <- length(s$x)
    factor <- Inf
    for (i in 1:n) {
        for (j in 1:n) {
            if (i != j) {
                # Deviation from published algorithm here
                # to use abs(w) so that ensure non-overlapping
                # circles even when weights are negative
                if (FALSE) {
                    f = sqrt((s$x[i] - s$x[j])^2 +
                        (s$y[i] - s$y[j])^2)/(abs(w[i]) + abs(w[j]))
                }
                # Original published algorithm
                if (FALSE) {
                    f = sqrt((s$x[i] - s$x[j])^2 +
                        (s$y[i] - s$y[j])^2)/(w[i] + w[j])
                }
                # Another variation where max(abs(w)) used
                # to avoid circle overlapping other site
                if (FALSE) {
                    f = sqrt((s$x[i] - s$x[j])^2 +
                        (s$y[i] - s$y[j])^2)/
                            max(abs(w[i]), abs(w[j]))
                }
                
                if (f > 0 && f < factor)
                    factor <- f
            }
        }
    }
    if (factor < 1) {
        w <- w*factor
    }
    w
}

# Variation from published algorithm
# Allow factor adjustment per pair of sites
shiftWeights <- function(s, w) {
    n <- length(s$x)
    for (i in 1:n) {
        for (j in 1:n) {
            if (i != j) {
                # Deviation from published algorithm here
                # to use abs(w) so that ensure non-overlapping
                # circles even when weights are negative
                f = sqrt((s$x[i] - s$x[j])^2 +
                    (s$y[i] - s$y[j])^2)/(abs(w[i]) + abs(w[j]))
                
                if (f > 0 && f < 1) {
                    w[i] <- w[i]*f
                    w[j] <- w[j]*f
                }
            }
        }
    }
    w
}

# The algorithm can fail to converge sometimes so
# just give up after 'maxIteration's
allocate <- function(names, s, w, outer, target, maxIteration=200,
                     debug=FALSE, dplot=FALSE, debugCell=FALSE) {
    count <- 1
    debugPlot <<- debugPlotGen()
    repeat {
        k <- awv(s, w, outer, debug, debugCell)
        areas <- lapply(k, area.poly)
        if (debug) {
            drawRegions(list(names=names, k=k,
                             s=s, w=w, a=areas, t=target),
                        debug)
            info <- rbind(area=round(unlist(areas)/sum(unlist(areas)), 4),
                          target=round(target, 4),
                          weight=round(w, 1))
            colnames(info) <- names
            print(info)
        }
        if (count > maxIteration ||
            breaking(unlist(areas), target, debug=debug, dplot=dplot)) {
            return(list(names=names, k=k, s=s, w=w, a=areas, t=target))
        } else {
            w <- adjustWeights(w, unlist(areas), target)
            s <- shiftSites(s, k)
            w <- shiftWeights(s, w)
        }
        count <- count + 1
    }
}


