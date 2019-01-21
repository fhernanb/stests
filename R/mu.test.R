#' Test for mean vector with covariance matrix known
#' 
#' This function performs the test for mean vector of a multivariate normal population with covariance matrix known.
#' 
#' @param data a data frame.
#' @param mu0 the reference vector, it could be an usual vector obtain with `c()` or a row/column matrix.
#' @param Sigma the covariance matrix.
#' @param alpha by default its value is 0.05.
#'  
#' @return A list with the following components:
#' \item{x2}{the value of the statistic.}
#' \item{critic_value}{the critic value to make a decision.}
#' \item{pvalue}{the p-value for the test.}
#' 
#' @examples
#' # Example 3.4.1 from Diaz & Morales (2016)
#' # H0: mu = (70, 170)
#' # H1: mu != (70, 170)
#' x1 <- c(69, 74, 68, 70, 72, 67, 66, 70, 76, 68,
#'         72, 79, 74, 67, 66, 71, 74, 75, 75, 76)
#' x2 <- c(153, 175, 155, 135, 172, 150, 115, 137, 200, 130,
#'         140, 265, 185, 112, 140, 150, 165, 185, 210, 220)
#' dt <- data.frame(x1, x2)
#' 
#' mu0 <- c(70, 170)
#' Sigma <- matrix(c(20, 100,
#'                   100, 1000), byrow=TRUE, ncol=2)
#' mu.test(data=dt, mu0=mu0, Sigma=Sigma, alpha=0.05)
#' 
#' @importFrom stats qchisq pchisq
#' @export
mu.test <- function(data, mu0, Sigma, alpha=0.05) {
  if (is.null(Sigma)) stop('Sigma matrix is needed.')
  p <- ncol(data)
  n <- nrow(data)
  mu0 <- matrix(mu0, ncol=1)
  mu <- matrix(colMeans(data), ncol=1)
  x2 <- n * t(mu-mu0) %*% solve(Sigma) %*% (mu-mu0)
  x2 <-  as.numeric(x2)
  critic_value <-  qchisq(p=alpha, df=p, lower.tail=FALSE)
  pvalue <- pchisq(q=x2, df=p, lower.tail=FALSE)
  list(x2=x2, critic_value=critic_value, pvalue=pvalue)
}
