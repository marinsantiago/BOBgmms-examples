y.norm <- mvnfast::rmvn(500, 1:5, diag(5))
x.norm <- mvnfast::rmvn(500, 1:5, diag(5))
z.gamm <- matrix(rgamma(500 * 5, 1, 1), nrow = 500)

plot(x.norm[,1], x.norm[,2])
points(y.norm[,1], y.norm[,2], col = 2)
points(z.gamm[,1], z.gamm[,2], col = 4)

library(T4transport)

swdist(y.norm, x.norm)$distance
swdist(y.norm, z.gamm)$distance

