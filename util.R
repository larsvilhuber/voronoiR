
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

# Extracting coordinates from "gpc.poly"

# It is possible for us to end up with a cell containing a hole
# (because when the cell boundary is either zero area or
#  contains a zero-area indent, gpclib::intersect() can
#  produce isolated islands)

# To handle these cases, we just ignore the holes (which should
# be zero-area anyway) and take the outer boundary.

# There is defensive code in there to warn if the unexpected
# happens and there is more than one outer boundary.
getpts <- function(x) {
    pts <- get.pts(x)
    nonholes <- sapply(pts, function(y) !y$hole)
    if (sum(nonholes) < 1 || sum(nonholes) > 1) {
        warning("Cell or region with more than one boundary")
    }
    border <- which(nonholes)[1]
    list(x=pts[[border]]$x,
         y=pts[[border]]$y)
}
