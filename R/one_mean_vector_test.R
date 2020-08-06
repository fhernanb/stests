#' Test for \eqn{\mu} in a \eqn{Np(\mu, \Sigma)}
#'
#' This function can be used to test \eqn{H_0: \mu = \mu_0} versus \eqn{H_1: \mu} not = \eqn{\mu_0} under \eqn{\Sigma} known or unkwon.
#'
#' @param mu0 a vector indicating the hypothesized value of \eqn{\mu}.
#' @param xbar a vector with the sample mean.
#' @param n sample size.
#' @param S a matrix with sample variances and covariances.
#' @param Sigma the matrix \eqn{\Sigma} if known.
#'
#' @details The user must provide only one matrix, S to perform the T2 test or \eqn{\Sigma} to perform the X2 test. When \eqn{\Sigma} is unkwon, T2 is perform and two values are provided in the print output, the T2 and F value.
#'
#' @seealso \link{one_covar_matrix_test} for test \eqn{\Sigma} in a \eqn{Np(\mu, \Sigma)}.
#'
#' @return A list with class \code{"htest"} containing the following components:
#' \item{statistic}{the value of the statistic.}
#' \item{parameter}{the degrees of freedom for the test.}
#' \item{p.value}{the p-value for the test.}
#' \item{estimate}{the estimated covariance matrix S.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{method}{a character string indicating what type of test was performed.}
#'
#' @author Freddy Hernandez
#' @example examples/examples_one_mean_vector_test.R
#'
#' @importFrom stats pf pchisq
#' @export
#'
one_mean_vector_test <- function (mu0, xbar, n, S=NULL, Sigma=NULL) {

  if (length(mu0) != length(xbar))
    stop("The vectors mu0 and xbar do not have same dimension")

  type <- sum(c(is.null(S), is.null(Sigma)))
  if (type == 2)
    stop("Please provide S or Sigma !!!")
  if (type == 0)
    stop("Please provide only 1 covariance matrix !!!")

  p <- length(mu0)

  mu0  <- matrix(mu0, ncol=1)
  xbar <- matrix(xbar, ncol=1)

  # T2 test
  if (! is.null(S)) {
    T2 <- n * t(xbar - mu0) %*% solve(S) %*% (xbar - mu0)
    T2 <- as.numeric(T2)
    p.value <- pf(q=(n - p) * T2/(p * (n - 1)), df1=p,
                  df2=n - p, lower.tail=FALSE)

    parameter <- c(p, n-p)
    names(parameter) <- c('df1', 'df2')
    method <- 'T2 test for mean vector'
    statistic <- c(T2, (n - p) * T2/(p * (n - 1)))
    names(statistic) <- c('T2', 'F')

  }

  # X2 test
  if (! is.null(Sigma)) {
    X2 <- n * t(xbar - mu0) %*% solve(Sigma) %*% (xbar - mu0)
    X2 <- as.numeric(X2)
    p.value <- pchisq(q=X2, df=p, lower.tail=FALSE)

    parameter <- p
    names(parameter) <- 'df'
    method <- 'X2 test for mean vector'
    statistic <- X2
    names(statistic) <- 'X2'
  }

  alternative <- paste0("true mean vector is not equal to (", paste0(mu0, collapse=", "), ") \n", sep="")
  estimate <- c(xbar)
  names(estimate) <- paste('xbar', 1:p, sep='_')
  data.name <- 'this test uses summarized data'

  res <- list(statistic=statistic,
              parameter=parameter,
              p.value=p.value,
              estimate=estimate,
              #null.value=null.message,
              alternative=alternative,
              method=method,
              data.name=data.name
  )

  class(res) <- "htest"
  res
}

