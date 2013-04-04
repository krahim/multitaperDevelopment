##     The multitaper R package (development code)
##     Multitaper and spectral analysis package for R
##     Copyright (C) 2013 Karim Rahim 
##
##     Written by Karim Rahim.
##
##     This file is part of the multitaper package for R.
##     http://cran.r-project.org/web/packages/multitaper/index.html
## 
##     The multitaper package is free software: you can redistribute it and 
##     or modify it under the terms of the GNU General Public License as 
##     published by the Free Software Foundation, either version 2 of the 
##     License, or any later version.
##
##     The multitaper package is distributed in the hope that it will be 
##     useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
##     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##     GNU General Public License for more details.
##
##     You should have received a copy of the GNU General Public License
##     along with multitaper.  If not, see <http://www.gnu.org/licenses/>.
##
##     If you wish to report bugs please contact the author:
## 
##     Karim Rahim
##     karim.rahim@gmail.com
##     Jeffery Hall, Queen's University, Kingston Ontario
##     Canada, K7L 3N6

library("multitaper")

## Utilities for and qi code

sizeplotForRightAxis <- function() {
    par(mar=c(5.1,4.1,4.1,4.1))
}

scaleFnToFit <- function(fnToFit, fnPlotted) {
    const1 <- min(fnPlotted) - min(fnToFit)
    scale1 <- diff(range(fnPlotted)) / diff(range(fnToFit))
    const1+scale1*fnToFit
}
