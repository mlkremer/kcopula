## Marginal CDF of the K-distribution
Fm <- Vectorize(function(x, N) {
  integrand <- function(z, x, N) {
    z^(N/2-1)*exp(-z)*pnorm(x*sqrt(N/(2*z)))
    }
  1/(gamma(N/2))*integrate(integrand, 0, Inf, x = x, N = N)$value
})

## Marginal PDF of the K-distribution
fm <- Vectorize(function(x, N) {
  integrand <- function(z, x, N) {
    z^(N/2-1)*exp(-z)*exp(-N/(4*z)*x^2)/sqrt(4*pi*z/N)
    }
  1/(gamma(N/2))*integrate(integrand, 0, Inf, x = x, N = N)$value
})

## Inverse marginal CDF of the K-distribution
Fmi <- Vectorize(function(u, N) {
  x <- seq(-4, 4, .1)  # sufficient for m \in [0.1, 500] (i.e. always)
  y <- Fm(x, N)  # sapply(x, function(a) {Fm(a, N)})
  approx(y, x, xout = u)$y
})

## Joint CDF of the K-distribution
#' @importFrom pracma integral3
Fj <- Vectorize(function(x, y, c, N) {
  integrand <- function(z, x, y, c, N) {
    1/(gamma(N/2))*z^(N/2-1)*exp(-z)/sqrt(1-c^2)*N/(4*pi*z)*exp(-N/(4*z)*(x^2-2*c*x*y+y^2)/(1-c^2))
  }
  a <- x
  b <- y
  integrand_tf <- function(z, x, y, c, N) {
    x_t <- a - (1-x)/x
    y_t <- b - (1-y)/y
    z_t <- z/(1-z)
    integrand(z_t, x_t, y_t, c = c, N = N)*(1/(1-z)^2)*(1/x^2)*(1/y^2)
  }
  pracma::integral3(integrand_tf, 0, 1, 0, 1, 0, 1, c = c, N = N)#$Q
  # cubature::hcubature(integrand_tf, c(0, 0, 0), c(1, 1, 1), c = c, N = N)$integral
})

## Joint PDF of the K-distribution
fj <- Vectorize(function(x, y, c, N) {
  integrand <- function(z, x, y, c, N) {
    z^(N/2-1)*exp(-z)/sqrt(1-c^2)*N/(4*pi*z)*exp(-N/(4*z)*(x^2-2*c*x*y+y^2)/(1-c^2))
    }
  1/(gamma(N/2))*integrate(integrand, 0, Inf, x = x, y = y, c = c, N = N)$value
  # 1/(gamma(N/2))*pracma::integral(integrand, 0, Inf, x = x, y = y, c = c, N = N)$Q
  # 1/(gamma(N/2))*cubature::hcubature(integrand, 0, Inf, x = x, y = y, c = c, N = N)$integral
})
