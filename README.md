# kcopula

## Overview
kcopula provides the bivariate K-copula by Wollschläger and Schäfer (2016).
It provides two functions:
* `pkcopula()` gives the distribution function of the bivariate K-copula.
* `dkcopula()` gives the density of the bivariate K-copula.

## Installation
Install release version from CRAN:
```r
install.packages("kcopula")
```
Install development version from GitHub:
```r
# install.packages("devtools")
devtools::install_github("mlkremer/kcopula")
```

## Usage
This example plots the bivariate K-copula density and distribution function.
```r
library(kcopula)

## Parameters
b <- .05   # step size
u <- seq(b, 1 - b, b)
v <- u
rho <- .2
N <- 4

## K-copula CDF
pkcopula(.5, .5, rho, N)

## Plot full K-copula CDF
kcopula_cdf <- pkcopula(u, v, rho, N, output = "matrix")
persp(u, v, kcopula_cdf, zlim = c(0, 1), xlab = "\n\n u", ylab = "\n\n v",
      zlab = "\n\n CDF", theta = 30, phi = 30, ticktype = "detailed")

## K-copula PDF
dkcopula(.5, .5, rho, N)

## Plot full K-copula PDF
kcopula_pdf <- dkcopula(u, v, rho, N, output = "matrix")
persp(u, v, kcopula_pdf, zlim = c(0, max(kcopula_pdf)), xlab = "\n\n u", ylab = "\n\n v",
      zlab = "\n\n PDF", theta = 30, phi = 30, ticktype = "detailed")
```

## References
Wollschläger, M. and Schäfer, R. (2016). Impact of nonstationarity on estimating 
and modeling empirical copulas of daily stock returns. 
*Journal of Risk*, 19(1):1--23. https://doi.org/10.21314/JOR.2016.342. 
SSRN version: https://ssrn.com/abstract=3533903.

Chetalova, D., Wollschläger, M., and Schäfer, R. (2015). 
Dependence structure of market states.
*Journal of Statistical Mechanics: Theory and Experiment*, 2015(8):P08012.
https://doi.org/10.1088/1742-5468/2015/08/P08012.
SSRN version: https://ssrn.com/abstract=3533951.
