#' The Bivariate K-Copula
#'
#' Density and distribution function of the bivariate K-copula
#' by Wollschläger and Schäfer (2016).
#'
#' @param u,v Numeric vectors with values in \eqn{[0, 1]}.
#' @param c Numeric; Pearson correlation coefficient in \eqn{[-1, 1]}.
#' @param N Numeric; inverse fluctuation strength of correlations around
#'     their average \code{c}, \eqn{N>0}. The larger \code{N} the smaller the
#'     fluctuations around \code{c}, and vice versa.
#' @param output Character; output as "vector" (default) for single values of
#'     the K-copula, or "matrix" for the full K-copula.
#' @param method Character; method to be used for
#'     \code{pkcopula(..., output = "vector")}.
#'     If \code{method = "interpolate"} (default), values are computed by
#'     interpolating the bivariate K-copula \emph{distribution function}
#'     (computationally faster); returns \code{NA}, if \code{u, v} are out of
#'     range (here: outside of \eqn{[.025, .975]}).
#'     If \code{method = "integrate"}, values are computed by integrating
#'     the bivariate K-copula \emph{density} (computationally slower).
#'
#' @return \code{dkcopula} gives the density (PDF), \code{pkcopula}
#'     gives the distribution function (CDF) of the bivariate K-copula.
#'
#' @inherit kcopula author references
#'
#' @keywords distribution multivariate
#'
#' @name kcop
#'
#' @examples
#' ## Parameters
#' u <- seq(.05, .95, .05)
#' v <- u
#' rho <- .2
#' N <- 4
#'
#' ## K-copula PDF
#' dkcopula(.5, .5, rho, N)
#'
#' \donttest{
#' ## Plot full K-copula PDF
#' kcopula_pdf <- dkcopula(u, v, rho, N, output = "matrix")
#' persp(u, v, kcopula_pdf)
#'
#' ## K-copula CDF
#' pkcopula(.5, .5, rho, N)
#'
#' ## Plot full K-copula CDF
#' kcopula_cdf <- pkcopula(u, v, rho, N, output = "matrix")
#' persp(u, v, kcopula_cdf)
#' }
#'
NULL

#' @rdname kcop
#' @importFrom pracma interp2
#' @export
#'
pkcopula <- function(u, v, c, N, output = "vector", method = "interpolate") {

  cop_cdf <- Vectorize(function(u, v, c, N) {
    Fj(Fmi(u, N), Fmi(v, N), c, N)
  })

  if (output == "vector") {

    if (method == "interpolate") {
      b <- .05   # step size
      x <- seq(b/2, 1 - b/2, b)
      y <- x
      Z <- (apply(apply(dkcopula(x, y, c, N, "matrix"), 2, cumsum), 1, cumsum)
            * diff(x)[1] * diff(y)[1])
      pracma::interp2(x, y, Z, u, v)

    } else if (method == "integrate") {
      pracma::integral2(dkcopula, 1e-3, u, 1e-3, v, c = c, N = N)$Q
      # cubature::hcubature(dkcopula, c(1e-3, 1e-3), c(u, v), c = c, N = N)$integral
      # cop_cdf(u, v, c, N)   # does not work for all parameter sets

    }

  } else if (output == "matrix") {
    (apply(apply(dkcopula(u, v, c, N, "matrix"), 2, cumsum), 1, cumsum)
     * diff(u)[1] * diff(v)[1])
  }
}

#' @rdname kcop
#' @export
#'
dkcopula <- function(u, v, c, N, output = "vector") {

  cop_pdf <- Vectorize(function(u, v, c, N) {
    fj(Fmi(u, N), Fmi(v, N), c, N) / (fm(Fmi(u, N), N) * fm(Fmi(v, N), N))
  })

  if (output == "vector") {
    cop_pdf(u, v, c, N)

  } else if (output == "matrix") {
    outer(u, v, cop_pdf, c, N)
  }
}
