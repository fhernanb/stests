# Example 5.2.3 from Diaz and Morales (2015) page 200
S1 <- matrix(c(12.65, -16.45,
               -16.45, 73.04), ncol=2, nrow=2)
S2 <- matrix(c(11.44, -27.77,
               -27.77, 100.64), ncol=2, nrow=2)
S3 <- matrix(c(14.46, -31.26,
               -31.26, 101.03), ncol=2, nrow=2)
N1 <- 26
N2 <- 23
N3 <- 25
S <- list(S1, S2, S3)
N <- list(N1, N2, N3)

res <- mult_var_matrices_test(S=S, N=N, method="box")
res
plot(res, shade.col="tomato")

# Example 5.3.4 from Mardia (1979) page 141
S1 <- matrix(c(132.99, 75.85, 35.82,
               75.85, 47.96, 20.75,
               35.82, 20.75, 10.79), ncol=3, nrow=3)
S2 <- matrix(c(432.58, 259.87, 161.67,
               259.87, 164.57, 98.99,
               161.67, 98.99, 63.87), ncol=3, nrow=3)
N1 <- 24
N2 <- 24
S <- list(S1, S2)
N <- list(N1, N2)

res <- mult_var_matrices_test(S=S, N=N, method="modified_LRT")
res
plot(res, from=20, to=30, shade.col="pink")

res <- mult_var_matrices_test(S=S, N=N, method="wald_schott")
res
plot(res, from=100, to=110, shade.col="lightblue")

