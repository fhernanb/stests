#' Tests for Equality of Two Normal Mean Vectors
#'
#' Implements the test for \eqn{H_0: \mu_1 = \mu_2} versus \eqn{H_1: \mu_1} not = \eqn{\mu_2} when both random samples are from two p-variate normal populations \eqn{Np(\mu_1, \Sigma_1)} and \eqn{Np(\mu_2, \Sigma_2)}. By default, this function performs the Hotelling test for two normal mean vectors assumming equality in the covariance matrices. Also testing the multivariate Behrens-Fisher problem, when the assumption of equal covariance matrices is violated, including James, Yao, Johansen, ... approaches.
#'
#' @param xbar1 a vector with the sample mean from population 1.
#' @param Sigma1 a matrix with sample variances and covariances from population 1.
#' @param n1 sample size 1.
#' @param xbar2 a vector with the sample mean from population 2.
#' @param Sigma2 a matrix with sample variances and covariances from population 2.
#' @param n2 sample size 2.
#' @param delta0 a number indicating the true value of the difference in means.
#' @param alpha the significance level, by default its value is 0.05.
#' @param method a character string specifying the method, it must be one of \code{"T2"} (default), \code{"james"} (James first order test), \code{"yao"} (Yao test), \code{"johansen"} (Johansen test), ...
#'
#' @return A list with class \code{"htest"} containing the following components:
#' \item{statistic}{the value of the statistic.}
#' \item{parameter}{the degrees of freedom for the test.}
#' \item{p.value}{the p-value for the test.}
#' \item{estimate}{the estimated mean vectors.}
#' \item{method}{a character string indicating the type of test performed.}
#'
#' @author Freddy Hernandez
#' @examples
#' # Example 5.4.2 from Rencher & Christensen (2012) page 137.
#' n1 <- 32
#' xbar1 <- c(15.97, 15.91, 27.19, 22.75)
#' s1 <- matrix(c(5.192, 4.545, 6.522, 5.25, 4.545, 13.18, 6.76, 6.266, 6.522,
#'                6.76, 28.67, 14.47, 5.25, 6.266, 14.47, 16.65), ncol = 4)
#'
#' n2 <- 32
#' xbar2 <- c(12.34, 13.91, 16.66, 21.94)
#' s2 <- matrix(c(9.136, 7.549, 4.864, 4.151, 7.549, 18.6, 10.22, 5.446, 4.864,
#'                10.22, 30.04, 13.49, 4.151, 5.446, 13.49, 28), ncol = 4)
#'
#' two_mean_vector_test(xbar1 = xbar1, Sigma1 = s1, n1 = n1,
#'                      xbar2 = xbar2, Sigma2 = s2, n2 = n2, alpha = 0.01)
#'
#' # Example 3.7 from Seber (1984) page 116.
#' # using the James first order test (1954).
#' n1 <- 16
#' xbar1 <- c(9.82, 15.06)
#' s1 <- matrix(c(120, -16.3, -16.3, 17.8), ncol = 2)
#'
#' n2 <- 11
#' xbar2 <- c(13.05, 22.57)
#' s2 <- matrix(c(81.8, 32.1, 32.1, 53.8), ncol = 2)
#'
#' two_mean_vector_test(xbar1 = xbar1, Sigma1 = s1, n1 = n1,
#'                      xbar2 = xbar2, Sigma2 = s2, n2 = n2,
#'                      method = 'james')
#'
#' @importFrom stats pf
#' @export
two_mean_vector_test <- function(xbar1, Sigma1, n1,
                                 xbar2, Sigma2, n2,
                                 delta0=NULL, alpha=0.05,
                                 method="T2") {

  if (! identical(dim(Sigma1), dim(Sigma2)))
    stop("The matrices Sigma1 and Sigma2 do not have same dimension")

  method <- match.arg(arg=method,
                      choices=c("T2", "james"))

  # To generate the code for evaluating, without using cases
  my_code <- paste0("two_mean_vector_test_", method,
                    "(xbar1, Sigma1, n1, xbar2, Sigma2, n2, delta0=delta0, alpha=alpha)")

  # To obtain the result
  result <- eval(parse(text=my_code))
  class(result) <- "htest"
  return(result)
}
#' @importFrom stats pf qf
two_mean_vector_test_T2 <- function(xbar1, Sigma1, n1,
                                    xbar2, Sigma2, n2,
                                    delta0=NULL, alpha=0.05) {
  p <- ncol(Sigma1)
  v <- n1 + n2 - 2
  if (is.null(delta0)) delta0 <- matrix(0, ncol=1, nrow=p)
  xbar1 <- matrix(xbar1, ncol=1)
  xbar2 <- matrix(xbar2, ncol=1)

  sp <- ((n1 - 1) * Sigma1 + (n2 - 1) * Sigma2) / (n1 + n2 - 2)
  T2 <- t(xbar1-xbar2-delta0) %*% solve(sp * (1/n1+1/n2)) %*% (xbar1-xbar2-delta0)
  T2 <- as.numeric(T2)
  p.value <- pf(q=T2 * (v-p+1) / (v*p), df1=p, df2=v-p+1, lower.tail=FALSE)
  critic.value <- qf(p=alpha, df1=p, df2=v-p+1, lower.tail=F) * (v*p) / (v-p+1)

  method <- 'T2 test for two mean vectors'
  statistic <- c(T2, T2 * (v-p+1) / (v*p))
  names(statistic) <- c('T2', 'F')
  parameter <- c(p, v-p+1)
  names(parameter) <- c('df1', 'df2')
  alternative <- "mu1 is not equal to mu2 \n"
  estimate <- cbind(xbar1, xbar2)
  colnames(estimate) <- c('Sample 1', 'Sample 2')
  rownames(estimate) <- paste('xbar', 1:p, sep = '_')
  data.name <- 'this test uses summarized data'

  return(list(statistic = statistic,
              parameter = parameter,
              p.value = p.value,
              estimate = estimate,
              alternative = alternative,
              method = method,
              data.name = data.name,
              sp = sp,
              critic.value = critic.value))
}
#' @importFrom stats pf qf
two_mean_vector_test_james <- function(xbar1, Sigma1, n1,
                                       xbar2, Sigma2, n2,
                                       delta0=NULL, alpha=0.05) {
  p <- ncol(Sigma1)
  xbar1 <- matrix(xbar1, ncol=1)
  xbar2 <- matrix(xbar2, ncol=1)

  S1 <- Sigma1/n1 # Represents S1 tilde
  S2 <- Sigma2/n2 # Represents S2 tilde
  S <- S1 + S2    # Represents S tilde

  T2 <- t(xbar1-xbar2) %*% solve(S) %*% (xbar1-xbar2)
  T2 <- as.numeric(T2)

  # To obtain delta
  Sinv <- solve(S)
  b1 <- Sinv %*% S1
  b2 <- Sinv %*% S2
  trb1 <- sum(diag(b1))
  trb2 <- sum(diag(b2))
  A <- 1 + (trb1^2/(n1 - 1) + trb2^2/(n2 - 1))/(2 * p)
  B <- (sum(b1^2)/(n1 - 1) + sum(b2^2)/(n2 - 1) + 0.5 *
          trb1^2/(n1 - 1) + 0.5 * trb2^2/(n2 - 1))/(p * (p + 2))
  x2 <- qchisq(1 - alpha, p)
  delta <- (A + B * x2)
  twoha <- x2 * delta
  p.value <- pchisq(q=T2/delta, df=p, lower.tail=FALSE)

  method <- 'James test for two mean vectors'
  statistic <- c(T2,T2/delta)
  names(statistic) <- c('T2', 'X-squared')
  parameter <- p
  names(parameter) <- 'df'
  alternative <- "mu1 is not equal to mu2 \n"
  estimate <- cbind(xbar1, xbar2)
  colnames(estimate) <- c('Sample 1', 'Sample 2')
  rownames(estimate) <- paste('xbar', 1:p, sep = '_')
  data.name <- 'this test uses summarized data'
  return(list(statistic = statistic,
              parameter = parameter,
              p.value = p.value,
              estimate = estimate,
              alternative = alternative,
              method = method,
              data.name = data.name))
}
