% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/one_covar_matrix_test.R
\name{one_covar_matrix_test}
\alias{one_covar_matrix_test}
\title{Test for \eqn{\Sigma} in a \eqn{Np(\mu, \Sigma)}}
\usage{
one_covar_matrix_test(Sigma0, S, n, method = "lrt")
}
\arguments{
\item{Sigma0}{a matrix indicating the hypothesized value of the covariance matrix \eqn{\Sigma}.}

\item{S}{a matrix with sample variances and covariances.}

\item{n}{sample size.}

\item{method}{a character string specifying the method, it must be one of \code{"lrt"} (default), \code{"modlrt1"} (modified LRT test) or \code{"modlrt2"} (modified LRT test for moderate \code{"n"}). You can specify just the initial letter. See details.}
}
\value{
A list with class \code{"htest"} containing the following components:
\item{statistic}{the value of the statistic.}
\item{parameter}{the degrees of freedom for the test.}
\item{p.value}{the p-value for the test.}
\item{estimate}{the estimated covariance matrix S.}
\item{method}{a character string indicating the type of test performed.}
}
\description{
This function can be used to test \eqn{H_0: \Sigma = \Sigma_0} versus \eqn{H_1: \Sigma} not = \eqn{\Sigma_0}.
}
\details{
When \code{method="lrt"} (default) the function performs the LRT test given in Mardia et. al (1979), page 126, expression 5.2.7. For \code{method="modlrt1"} or \code{method="modlrt2"} the function performs the LRT test given in Rencher and Christensen (2012), page 260, expressions 7.2 and 7.4.
}
\examples{
# Example 5.3.2 from Mardia (1979) page 127
# Test H0: Sigma = diag(100, 100) versus H1: Sigma != diag(100, 100)

Sigma0 <- matrix(c(100, 0, 0, 100), ncol=2)
S <- matrix(c(91.481, 66.875, 66.875, 96.775), ncol=2)

res1 <- one_covar_matrix_test(Sigma0=Sigma0, S=S, n=25, method='lrt')
res1
plot(res1, from=12, to=20, shade.col='dodgerblue2')

# Example from Morrison (1990) page 293
# Test H0: Sigma = Sigma0 versus H1: Sigma != Sigma0
# using the modified LRT test versions

n <- 20

Sigma0 <- matrix(c(4, 3, 2,
                   3, 6, 5,
                   2, 5, 10), ncol=3)

S <- matrix(c(3.42, 2.60, 1.89,
              2.60, 8.00, 6.51,
              1.89, 6.51, 9.62), ncol=3)

res2 <- one_covar_matrix_test(Sigma0=Sigma0, S=S, n=n, method='modlrt1')
res2
plot(res2, from=0, to=20, shade.col='indianred1')

res3 <- one_covar_matrix_test(Sigma0=Sigma0, S=S, n=n, method='modlrt2')
res3
plot(res3, from=0, to=20, shade.col='aquamarine3')

}
\seealso{
\link{one_mean_vector_test} for test \eqn{\mu} in a \eqn{Np(\mu, \Sigma)}.
}
\author{
Freddy Hernandez
}
