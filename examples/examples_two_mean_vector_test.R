# Example 5.4.2 from Rencher & Christensen (2012) page 137,
# using Hotelling's test
n1 <- 32
xbar1 <- c(15.97, 15.91, 27.19, 22.75)
s1 <- matrix(c(5.192, 4.545, 6.522, 5.25,
               4.545, 13.18, 6.76, 6.266,
               6.522, 6.76, 28.67, 14.47,
               5.25, 6.266, 14.47, 16.65), ncol = 4)

n2 <- 32
xbar2 <- c(12.34, 13.91, 16.66, 21.94)
s2 <- matrix(c(9.136, 7.549, 4.864, 4.151,
               7.549, 18.6, 10.22, 5.446,
               4.864, 10.22, 30.04, 13.49,
               4.151, 5.446, 13.49, 28), ncol = 4)

res1 <- two_mean_vector_test(xbar1 = xbar1, s1 = s1, n1 = n1,
                             xbar2 = xbar2, s2 = s2, n2 = n2,
                             method = "T2")
res1
plot(res1, from=21, to=25, shade.col='tomato')

# Example 3.7 from Seber (1984) page 116.
# using the James first order test (1954).
n1 <- 16
xbar1 <- c(9.82, 15.06)
s1 <- matrix(c(120, -16.3,
               -16.3, 17.8), ncol = 2)

n2 <- 11
xbar2 <- c(13.05, 22.57)
s2 <- matrix(c(81.8, 32.1,
               32.1, 53.8), ncol = 2)

res2 <- two_mean_vector_test(xbar1 = xbar1, s1 = s1, n1 = n1,
                             xbar2 = xbar2, s2 = s2, n2 = n2,
                             method = 'james', alpha=0.05)
res2
plot(res2, from=5, to=10, shade.col="lightgreen")

# Example from page 141 from Yao (1965),
# using Yao's test

n1 <- 16
xbar1 <- c(9.82, 15.06)
s1 <- matrix(c(120, -16.3,
               -16.3, 17.8), ncol = 2)

n2 <- 11
xbar2 <- c(13.05, 22.57)
s2 <- matrix(c(81.8, 32.1,
               32.1, 53.8), ncol = 2)

res3 <- two_mean_vector_test(xbar1 = xbar1, s1 = s1, n1 = n1,
                             xbar2 = xbar2, s2 = s2, n2 = n2,
                             method = 'yao')
res3
plot(res3, from=2, to=6, shade.col="pink")


# Example for Johansen's test using the data from
# Example from page 141 from Yao (1965)

res4 <- two_mean_vector_test(xbar1 = xbar1, s1 = s1, n1 = n1,
                             xbar2 = xbar2, s2 = s2, n2 = n2,
                             method = 'johansen')
res4
plot(res4, from=2, to=6, shade.col="aquamarine1")

# Example 4.1 from Nel and Van de Merwe (1986) page 3729
# Test H0: mu1 = mu2 versus H1: mu1 != mu2
n1 <- 45
xbar1 <- c(204.4, 556.6)
s1 <- matrix(c(13825.3, 23823.4,
               23823.4, 73107.4), ncol=2)

n2 <- 55
xbar2 <- c(130.0, 355.0)
s2 <- matrix(c(8632.0, 19616.7,
               19616.7, 55964.5), ncol=2)

res5 <- two_mean_vector_test(xbar1 = xbar1, s1 = s1, n1 = n1,
                             xbar2 = xbar2, s2 = s2, n2 = n2,
                             method = 'nvm')
res5
plot(res5, from=6, to=10, shade.col='cyan2')

# using mnvm method for same data
res6 <- two_mean_vector_test(xbar1 = xbar1, s1 = s1, n1 = n1,
                             xbar2 = xbar2, s2 = s2, n2 = n2,
                             method = 'mnvm')
res6
plot(res6, from=6, to=10, shade.col='lightgoldenrodyellow')

# using gamage method for same data
res7 <- two_mean_vector_test(xbar1 = xbar1, s1 = s1, n1 = n1,
                             xbar2 = xbar2, s2 = s2, n2 = n2,
                             method = 'gamage')
res7
plot(res7, from=0, to=20, shade.col='mediumpurple4')
text(x=10, y=0.30, "The curve corresponds to an empirical density",
     col='orange')

# using Yanagihara and Yuan method for same data
res8 <- two_mean_vector_test(xbar1 = xbar1, s1 = s1, n1 = n1,
                             xbar2 = xbar2, s2 = s2, n2 = n2,
                             method = 'yy')
res8
plot(res8, from=6, to=10, shade.col='gold2')

# using Bartlett Correction method for same data
res9 <- two_mean_vector_test(xbar1 = xbar1, s1 = s1, n1 = n1,
                             xbar2 = xbar2, s2 = s2, n2 = n2,
                             method = 'byy')
res9

# using modified Bartlett Correction method for same data
res10 <- two_mean_vector_test(xbar1 = xbar1, s1 = s1, n1 = n1,
                              xbar2 = xbar2, s2 = s2, n2 = n2,
                              method = 'mbyy')
res10

# using Second Order Procedure (Kawasaki and Seo (2015)) method for same data
res11 <- two_mean_vector_test(xbar1 = xbar1, s1 = s1, n1 = n1,
                              xbar2 = xbar2, s2 = s2, n2 = n2,
                              method = 'ks1')
res11
plot(res11, from=6, to=10, shade.col='hotpink1')

# using Bias Correction Procedure (Kawasaki and Seo (2015)) method for same data
res12 <- two_mean_vector_test(xbar1 = xbar1, s1 = s1, n1 = n1,
                              xbar2 = xbar2, s2 = s2, n2 = n2,
                              method = 'ks2')
res12
plot(res12, from=6, to=10, shade.col='seagreen2')
