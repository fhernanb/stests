# Example 5.2.2 from Rencher & Christensen (2012) page 127
# Test H0: mu = (70, 170) versus H1: mu != (70, 170)
# with known Sigma

Sigma <- matrix(c(20, 100, 100, 1000), ncol=2, nrow=2)

res1 <- one_mean_vector_test(mu0=c(70, 170), xbar=c(71.45, 164.7),
                             n=20, Sigma=Sigma)
res1
plot(res1, from=4, to=10, shade.col='dodgerblue2')

# Repeating the last example with raw data
x1 <- c(69, 74, 68, 70, 72, 67, 66, 70, 76, 68,
        72, 79, 74, 67, 66, 71, 74, 75, 75, 76)
x2 <- c(153, 175, 155, 135, 172, 150, 115, 137, 200, 130,
        140, 265, 185, 112, 140, 150, 165, 185, 210, 220)
dt <- data.frame(x1, x2)
mu0 <- c(70, 170)
Sigma <- matrix(c(20, 100, 100, 1000), ncol=2, nrow=2)

res2 <- one_mean_vector_test(mu0=mu0, xbar=colMeans(dt),
                             n=nrow(dt), Sigma=Sigma)
res2
plot(res2, from=4, to=10, shade.col='lightpink1')

# Example 5.2 from Johnson and Wichern (2012) page 214
# Test H0: mu = (4, 50, 10) versus H1: mu != (4, 50, 10)
# with unknown Sigma

S <- matrix(c(2.879, 10.010, -1.810,
              10.010, 199.788, -5.640,
              -1.810, -5.640, 3.628), ncol=3, nrow=3)

res3 <- one_mean_vector_test(mu0=c(4, 50, 10),
                             xbar=c(4.640, 45.400, 9.965),
                             n=20, S=S)
res3
plot(res3, from=0, to=5, shade.col='aquamarine3')

\dontrun{
library(rrcov)
data(delivery)
delivery.x <- delivery[, 1:2]

# Using T2.test from rrcov package
T2.test(delivery.x)

one_mean_vector_test(mu0=c(0, 0), xbar=colMeans(delivery.x),
                     n=25, S=var(delivery.x))
}
