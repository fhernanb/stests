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
#' @return A list with class \code{"htest"}.
#' @author Freddy Hernandez
#' @examples
#' # Example 5.2.2 from Rencher & Christensen (2012) page 127
#' # Test H0: mu = (70, 170) versus H1: mu != (70, 170)
#' # with known Sigma
#'
#' Sigma <- matrix(c(20, 100, 100, 1000), ncol=2, nrow=2)
#'
#' res1 <- one_mean_vector_test(mu0=c(70, 170), xbar=c(71.45, 164.7),
#'                              n=20, Sigma=Sigma)
#' res1
#' plot(res1, from=4, to=10, shade.col='dodgerblue2')
#'
#' # Repeating the last example with raw data
#' x1 <- c(69, 74, 68, 70, 72, 67, 66, 70, 76, 68,
#'         72, 79, 74, 67, 66, 71, 74, 75, 75, 76)
#' x2 <- c(153, 175, 155, 135, 172, 150, 115, 137, 200, 130,
#'         140, 265, 185, 112, 140, 150, 165, 185, 210, 220)
#' dt <- data.frame(x1, x2)
#' mu0 <- c(70, 170)
#' Sigma <- matrix(c(20, 100, 100, 1000), ncol=2, nrow=2)
#'
#' res2 <- one_mean_vector_test(mu0=mu0, xbar=colMeans(dt),
#'                              n=nrow(dt), Sigma=Sigma)
#' res2
#' plot(res2, from=4, to=10, shade.col='lightpink1')
#'
#' # Example 5.2 from Johnson and Wichern (2012) page 214
#' # Test H0: mu = (4, 50, 10) versus H1: mu != (4, 50, 10)
#' # with unknown Sigma
#'
#' S <- matrix(c(2.879, 10.010, -1.810,
#'               10.010, 199.788, -5.640,
#'               -1.810, -5.640, 3.628), ncol=3, nrow=3)
#'
#' res3 <- one_mean_vector_test(mu0=c(4, 50, 10),
#'                              xbar=c(4.640, 45.400, 9.965),
#'                              n=20, S=S)
#' res3
#' plot(res3, from=0, to=5, shade.col='aquamarine3')
#'
#' \dontrun{
#' library(rrcov)
#' data(delivery)
#' delivery.x <- delivery[, 1:2]
#'
#' # Using T2.test from rrcov package
#' T2.test(delivery.x)
#'
#' one_mean_vector_test(mu0=c(0, 0), xbar=colMeans(delivery.x),
#'                      n=25, S=var(delivery.x))
#' }
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

