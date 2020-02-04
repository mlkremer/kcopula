library(kcopula)
invisible(lapply(c("plot3D", "rgl", "plotly"), library, character.only = T))

## Parameters
b <- .05   # step size
u <- seq(b, 1 - b, b)   # seq(b/2, 1 - b/2, b)
v <- u
rho <- .2
N <- 4

## K-copula PDF
dkcopula(.5, .5, rho, N)
kcopula_pdf <- dkcopula(u, v, rho, N, "matrix")   # takes 7 seconds (system.time())

# 3D plots using different packages
# graphics
persp(u, v, kcopula_pdf, zlim = c(0, max(kcopula_pdf)), xlab = "\n\n u", ylab = "\n\n v",
      zlab = "\n\n PDF", theta = 30, phi = 30, ticktype = "detailed")
# plot3D
persp3D(u, v, kcopula_pdf, zlim = c(0, max(kcopula_pdf)), xlab = "\n\n u", ylab = "\n\n v",
        zlab = "\n\n PDF", ticktype = "detailed")
# rgl
persp3d(u, v, kcopula_pdf)
# plotly
plot_ly(x = u, y = u, z = kcopula_pdf) %>% add_surface()


## K-copula CDF
pkcopula(.5, .5, rho, N)  # takes 7 or 94 seconds
kcopula_cdf <- pkcopula(u, v, rho, N, "matrix")   # takes 7 seconds

persp(u, v, kcopula_cdf, zlim = c(0, 1), xlab = "\n\n u", ylab = "\n\n v",
      zlab = "\n\n CDF", theta = 30, phi = 30, ticktype = "detailed")
