#' Test for \eqn{\mu} in a \eqn{Np(\mu, \Sigma)}
#'
#' This function can be used to test \eqn{H0: \mu = \mu0} versus \eqn{H1: \mu not = \mu0}.
#'
#' @param mu0 a vector indicating the hypothesized value of the mean.
#' @param xbar a vector with the sample mean.
#' @param n sample size.
#' @param S a matrix with sample variances and covariances.
#' @param Sigma the matrix \eqn{\Sigma} if known.
#' @param alpha the significance level.
#'
#' @return a list.
#'
#' @examples
#' # Example 5.2.2 from Rencher & Christensen (2012) page 127
#' # Test H0: mu = (70, 170) versus H1: mu != (70, 170)
#' # with known Sigma
#'
#' Sigma <- matrix(c(20, 100, 100, 1000), ncol=2, nrow=2)
#'
#' one_mean_vector_test(mu0=c(70, 170), xbar=c(71.45, 164.7),
#'                      n=20, Sigma=Sigma)
#'
#' # Example 5.2 from Johnson and Wichern (2012) page 214
#' # Test H0: mu = (4, 50, 10) versus H1: mu != (4, 50, 10)
#' # with unknown Sigma
#'
#' S <- matrix(c(2.879, 10.010, -1.810,
#'               10.010, 199.788, -5.640,
#'               -1.810, -5.640, 3.628), ncol=3, nrow=3)
#'
#' one_mean_vector_test(mu0=c(4, 50, 10), xbar=c(4.640, 45.400, 9.965),
#'                      n=20, S=S)
#'
#' @importFrom stats qf pf qchisq pchisq
#' @export
#'
one_mean_vector_test <- function (mu0, xbar, n,
                                  S=NULL, Sigma=NULL,
                                  alpha=0.05) {
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

  if (! is.null(S)) {
    T2 <- n * t(xbar - mu0) %*% solve(S) %*% (xbar - mu0)
    T2 <- as.numeric(T2)
    critic_value <- (n - 1) * p * qf(p=alpha, df1=p,
                                     df2=n-p, lower.tail=FALSE)/(n-p)
    p.value <- pf(q=(n - p) * T2/(p * (n - 1)), df1=p,
                  df2=n - p, lower.tail=FALSE)

    method <- 'T2 test for mean vector'
    statistic <- T2
    names(statistic) <- 'T2'

  }
  if (! is.null(Sigma)) {
    X2 <- n * t(xbar - mu0) %*% solve(Sigma) %*% (xbar - mu0)
    X2 <- as.numeric(X2)
    critic_value <- qchisq(p=alpha, df=p, lower.tail=FALSE)
    p.value <- pchisq(q=X2, df=p, lower.tail=FALSE)

    method <- 'X2 test for mean vector'
    statistic <- X2
    names(statistic) <- 'X2'
  }

  alternative <- "two.sided"
  null.message <- paste0("(", paste0(mu0, collapse=", "), ")", sep="")
  estimate <- c(xbar)
  names(estimate) <- paste('xbar', 1:p, sep='_')
  data.name <- 'this test used summarized data'

  res <- list(statistic=statistic,
              p.value=p.value,
              #conf.int=conf.int,
              estimate=estimate,
              null.value=null.message,
              alternative=alternative,
              method=method,
              data.name=data.name
  )

  class(res) <- "htest"
  res
}

