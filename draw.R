
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

drawPoly <- function(p, name, gp=gpar()) {
    if (length(p@pts)) {
        pts <- getpts(p)
        grid.polygon(pts$x, pts$y, default="native",
                     gp=gp, name=name)
    }
}

library(grid)

polyRangeX <- function(p) {
    if (length(p@pts)) {
        pts <- getpts(p)
        range(pts$x)
    } else {
        NA
    }
}

polyRangeY <- function(p) {
    if (length(p@pts)) {
        pts <- getpts(p)
        range(pts$y)
    } else {
        NA
    }
}

drawRegions <- function(result, label=FALSE, top=TRUE,
                        gp=gpar(), newpage=TRUE, debug=FALSE) {
    names <- result$names
    k <- result$k
    sites <- result$s
    weights <- result$w
    if (newpage) {
        grid.newpage()
        xrange <- range(unlist(lapply(k, polyRangeX)), na.rm=TRUE)
        yrange <- range(unlist(lapply(k, polyRangeY)), na.rm=TRUE)
        pushViewport(viewport(width=.8, height=.8,
                              xscale=xrange, yscale=yrange))
    }
    if (top) {
        invisible(mapply(drawPoly, k, names, MoreArgs=list(gp=gpar(lwd=8))))
        invisible(mapply(drawPoly, k, names,
                         MoreArgs=list(gp=gpar(lwd=6, col="grey90"))))
    } else {
        invisible(mapply(drawPoly, k, names, MoreArgs=list(gp=gp)))
    }
    if (label) {
        if (top) {
            cex=1
        } else {
            cex=.5
        }
        grid.text(result$names, sites$x, sites$y, default="native",
                  gp=gpar(cex=cex))
        if (debug) {
            col <- ifelse(weights < 0, "red", "black")
            grid.circle(sites$x, sites$y, default="native",
                        r=unit(weights, "native"),
                        gp=gpar(col=col, fill=NA))
        }
    }
}
