#' Tests for Equality of Two Normal Mean Vectors
#'
#' Implements the test for \eqn{H_0: \mu_1 = \mu_2} versus \eqn{H_1: \mu_1} not = \eqn{\mu_2} when both random samples are from two p-variate normal populations \eqn{Np(\mu_1, \Sigma_1)} and \eqn{Np(\mu_2, \Sigma_2)}. By default, this function performs the Hotelling test for two normal mean vectors assumming equality in the covariance matrices. Tests for the multivariate Behrens-Fisher problem are also implemented, the tests are James, Yao, Johansen, Nel and Van der Merwe.
#'
#' @param xbar1 a vector with the sample mean from population 1.
#' @param s1 a matrix with sample variances and covariances from population 1.
#' @param n1 sample size 1.
#' @param xbar2 a vector with the sample mean from population 2.
#' @param s2 a matrix with sample variances and covariances from population 2.
#' @param n2 sample size 2.
#' @param delta0 a number indicating the true value of the difference in means.
#' @param method a character string specifying the method, it must be one of \code{"T2"} (default), \code{"james"} (James' first order test), \code{"yao"} (Yao's test), \code{"johansen"} (Johansen's test), \code{"nvm"} (Nel and Van Der Merwe test).
#' @param alpha the significance level for method \code{"james"}, by default its value is 0.05.
#'
#' @details For James test the critic value is reported, if T2 > critic_value we reject H0.
#'
#' @return A list with class \code{"htest"} containing the following components:
#' \item{statistic}{the value of the statistic.}
#' \item{parameter}{the degrees of freedom for the test.}
#' \item{p.value}{the p-value for the test.}
#' \item{estimate}{the estimated mean vectors.}
#' \item{method}{a character string indicating the type of test performed.}
#'
#' @author Freddy Hernandez, Jean Paul Piedrahita, Valentina Garcia.
#' @examples
#' # Example 5.4.2 from Rencher & Christensen (2012) page 137,
#' # using Hotelling's test
#' n1 <- 32
#' xbar1 <- c(15.97, 15.91, 27.19, 22.75)
#' s1 <- matrix(c(5.192, 4.545, 6.522, 5.25,
#'                4.545, 13.18, 6.76, 6.266,
#'                6.522, 6.76, 28.67, 14.47,
#'                5.25, 6.266, 14.47, 16.65), ncol = 4)
#'
#' n2 <- 32
#' xbar2 <- c(12.34, 13.91, 16.66, 21.94)
#' s2 <- matrix(c(9.136, 7.549, 4.864, 4.151,
#'                7.549, 18.6, 10.22, 5.446,
#'                4.864, 10.22, 30.04, 13.49,
#'                4.151, 5.446, 13.49, 28), ncol = 4)
#'
#' res1 <- two_mean_vector_test(xbar1 = xbar1, s1 = s1, n1 = n1,
#'                              xbar2 = xbar2, s2 = s2, n2 = n2,
#'                              method = "T2")
#' res1
#' plot(res1, from=21, to=25, shade.col='tomato')
#'
#'
#' # Example 3.7 from Seber (1984) page 116.
#' # using the James first order test (1954).
#' n1 <- 16
#' xbar1 <- c(9.82, 15.06)
#' s1 <- matrix(c(120, -16.3,
#'                -16.3, 17.8), ncol = 2)
#'
#' n2 <- 11
#' xbar2 <- c(13.05, 22.57)
#' s2 <- matrix(c(81.8, 32.1,
#'                32.1, 53.8), ncol = 2)
#'
#' res2 <- two_mean_vector_test(xbar1 = xbar1, s1 = s1, n1 = n1,
#'                              xbar2 = xbar2, s2 = s2, n2 = n2,
#'                              method = 'james')
#' res2
#' plot(res2, from=5, to=10, shade.col="lightgreen")
#'
#'
#' # Example from page 141 from Yao (1965),
#' # using Yao's test
#'
#' n1 <- 16
#' xbar1 <- c(9.82, 15.06)
#' s1 <- matrix(c(120, -16.3,
#'                -16.3, 17.8), ncol = 2)
#'
#' n2 <- 11
#' xbar2 <- c(13.05, 22.57)
#' s2 <- matrix(c(81.8, 32.1,
#'                32.1, 53.8), ncol = 2)
#'
#' res3 <- two_mean_vector_test(xbar1 = xbar1, s1 = s1, n1 = n1,
#'                              xbar2 = xbar2, s2 = s2, n2 = n2,
#'                              method = 'yao')
#' res3
#' plot(res3, from=2, to=6, shade.col="pink")
#'
#'
#' # Example for Johansen's test using the data from
#' # Example from page 141 from Yao (1965)
#'
#' res4 <- two_mean_vector_test(xbar1 = xbar1, s1 = s1, n1 = n1,
#'                              xbar2 = xbar2, s2 = s2, n2 = n2,
#'                              method = 'johansen')
#' res4
#' plot(res4, from=2, to=6, shade.col="aquamarine1")
#'
#' # Example 4.1 from Nel and Van de Merwe (1986) page 3729
#' # Test H0: mu1 = mu2 versus H1: mu1 != mu2
#' n1 <- 45
#' xbar1 <- c(204.4, 556.6)
#' s1 <- matrix(c(13825.3, 23823.4,
#'                23823.4, 73107.4), ncol=2)
#'
#' n2 <- 55
#' xbar2 <- c(130.0, 355.0)
#' s2 <- matrix(c(8632.0, 19616.7,
#'                19616.7, 55964.5), ncol=2)
#'
#' res5 <- two_mean_vector_test(xbar1 = xbar1, s1 = s1, n1 = n1,
#'                              xbar2 = xbar2, s2 = s2, n2 = n2,
#'                              method = 'nvm')
#' res5
#' plot(res5, from=6, to=10, shade.col='pink')
#'
#' @importFrom stats pf
#' @export
two_mean_vector_test <- function(xbar1, s1, n1, xbar2, s2, n2,
                                 delta0=NULL, method="T2", alpha=0.05) {

  if (! identical(dim(s1), dim(s2)))
    stop("The matrices s1 and s2 do not have same dimension")

  method <- match.arg(arg=method,
                      choices=c("T2", "james", "yao", "johansen",
                                "nvm"))

  # To generate the code for evaluating, without using cases
  my_code <- paste0("two_mean_vector_test_", method,
                    "(xbar1, s1, n1, xbar2, s2, n2, delta0=delta0, alpha=alpha)")

  # To obtain the result
  result <- eval(parse(text=my_code))
  class(result) <- "htest"
  return(result)
}
#' @importFrom stats pf qf
two_mean_vector_test_T2 <- function(xbar1, s1, n1,
                                    xbar2, s2, n2,
                                    delta0=NULL, alpha=0.05) {
  p <- ncol(s1)
  v <- n1 + n2 - 2
  if (is.null(delta0)) delta0 <- matrix(0, ncol=1, nrow=p)
  xbar1 <- matrix(xbar1, ncol=1)
  xbar2 <- matrix(xbar2, ncol=1)

  sp <- ((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2)
  T2 <- t(xbar1-xbar2-delta0) %*% solve(sp * (1/n1+1/n2)) %*% (xbar1-xbar2-delta0)
  T2 <- as.numeric(T2)
  p.value <- pf(q=T2 * (v-p+1) / (v*p), df1=p, df2=v-p+1, lower.tail=FALSE)

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
              sp = sp))
}
#' @importFrom stats pf qf
two_mean_vector_test_james <- function(xbar1, s1, n1,
                                       xbar2, s2, n2,
                                       delta0=NULL, alpha=0.05) {
  p <- ncol(s1)
  xbar1 <- matrix(xbar1, ncol=1)
  xbar2 <- matrix(xbar2, ncol=1)
  S1 <- s1/n1    # Represents S1 tilde, uppercase
  S2 <- s2/n2    # Represents S2 tilde, uppercase
  S  <- S1 + S2  # Represents S tilde

  T2 <- t(xbar1-xbar2) %*% solve(S) %*% (xbar1-xbar2)
  T2 <- as.numeric(T2)

  Sinv <- solve(S)
  b1 <- Sinv %*% S1
  b2 <- Sinv %*% S2
  trb1 <- sum(diag(b1))
  trb2 <- sum(diag(b2))
  A <- 1 + ( trb1^2/(n1 - 1) + trb2^2/(n2 - 1) ) / (2 * p)
  B <- (sum(diag(b1%*%b1))/(n1 - 1) + sum(diag(b2%*%b2))/(n2 - 1) +
          0.5 * trb1^2/(n1 - 1) + 0.5 * trb2^2/(n2 - 1)) / (p * (p + 2))
  x2 <- qchisq(p=1 - alpha, df=p)
  delta <- A + B * x2
  critic_value <- x2 * delta
  p.value <- NULL
  method <- 'James test for two mean vectors'
  statistic <- c(T2, critic_value)
  names(statistic) <- c('T2', 'critic_value')
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
              data.name = data.name,
              critic_value = critic_value
              ))
}
#' @importFrom stats pf
two_mean_vector_test_yao <- function(xbar1, s1, n1,
                                     xbar2, s2, n2,
                                     delta0=NULL, alpha=0.05) {

  p <- ncol(s1)
  xbar1 <- matrix(xbar1, ncol=1)
  xbar2 <- matrix(xbar2, ncol=1)

  S1 <- s1/n1    # Represents S1 tilde
  S2 <- s2/n2    # Represents S2 tilde
  S  <- S1 + S2  # Represents S tilde

  T2 <- t(xbar1-xbar2) %*% solve(S) %*% (xbar1-xbar2)
  T2 <- as.numeric(T2)

  tr <- function(x) sum(diag(x)) # To obtain the trace easily

  # To obtain v
  b1 <- solve(S) %*% S1 %*% solve(S)  # An auxiliar element
  b2 <- solve(S) %*% S2 %*% solve(S)  # An auxiliar element
  v1 <- (t(xbar1-xbar2) %*% b1 %*% (xbar1-xbar2))/T2
  v1 <- v1^2/(n1-1)
  v2 <- (t(xbar1-xbar2) %*% b2 %*% (xbar1-xbar2))/T2
  v2 <- v2^2/(n2-1)
  v <- as.numeric(1/(v1+v2))

  p.value <- pf(q=T2*(v-p+1)/(v*p), df1=p, df2=v-p+1, lower.tail=FALSE)

  method <- 'Yao test for two mean vectors'
  statistic <- c(T2, T2*(v-p+1)/(v*p))
  names(statistic) <- c('T2', 'F')
  parameter <- c(p, v-p+1)
  names(parameter) <- c('df1', 'df2')
  alternative <- "mu1 is not equal to mu2 \n"
  estimate <- cbind(xbar1, xbar2)
  colnames(estimate) <- c('Sample 1', 'Sample 2')
  rownames(estimate) <- paste('xbar', 1:p, sep='_')
  data.name <- 'this test uses summarized data'

  return(list(statistic = statistic,
              parameter = parameter,
              p.value = p.value,
              estimate = estimate,
              alternative = alternative,
              method = method,
              data.name = data.name))
}
#' @importFrom stats pf
two_mean_vector_test_johansen <- function(xbar1, s1, n1,
                                          xbar2, s2, n2,
                                          delta0=NULL, alpha=0.05) {

  p <- ncol(s1)
  xbar1 <- matrix(xbar1, ncol=1)
  xbar2 <- matrix(xbar2, ncol=1)

  S1 <- s1/n1    # Represents S1 tilde
  S2 <- s2/n2    # Represents S2 tilde
  S  <- S1 + S2  # Represents S tilde

  T2 <- t(xbar1-xbar2) %*% solve(S) %*% (xbar1-xbar2)
  T2 <- as.numeric(T2)

  tr <- function(x) sum(diag(x)) # To obtain the trace easily

  # To obtain q, v and D
  b1 <- diag(p) - solve(solve(S1) + solve(S2)) %*% solve(S1)  # An auxiliar element
  b2 <- diag(p) - solve(solve(S1) + solve(S2)) %*% solve(S2)  # An auxiliar element
  d1 <- tr(b1 %*% b1) + (tr(b1))^2
  d1 <- d1/(n1-1)
  d2 <- tr(b2 %*% b2) + (tr(b2))^2
  d2 <- d2/(n2-1)
  D <- (d1 + d2)/2
  v <- (p*(p+2))/(3*D)
  q <- p + 2*D - (6*D)/(p*(p-1) + 2)
  p.value <- pf(q = T2/q, df1 = p, df2 = v, lower.tail = FALSE)

  method <- 'Johansen test for two mean vectors'
  statistic <- c(T2, T2/q)
  names(statistic) <- c('T2', 'F')
  parameter <- c(p, v-p+1)
  names(parameter) <- c('df1', 'df2')
  alternative <- "mu1 is not equal to mu2 \n"
  estimate <- cbind(xbar1, xbar2)
  colnames(estimate) <- c('Sample 1', 'Sample 2')
  rownames(estimate) <- paste('xbar', 1:p, sep='_')
  data.name <- 'this test uses summarized data'

  return(list(statistic = statistic,
              parameter = parameter,
              p.value = p.value,
              estimate = estimate,
              alternative = alternative,
              method = method,
              data.name = data.name))
}
#' @importFrom stats pf
two_mean_vector_test_nvm <- function(xbar1, s1, n1,
                                     xbar2, s2, n2,
                                     delta0=NULL, alpha=0.05) {

  p <- ncol(s1)
  xbar1 <- matrix(xbar1, ncol=1)
  xbar2 <- matrix(xbar2, ncol=1)

  S1 <- s1/n1    # Represents S1 tilde
  S2 <- s2/n2    # Represents S2 tilde
  S  <- S1 + S2  # Represents S tilde

  T2 <- t(xbar1-xbar2) %*% solve(S) %*% (xbar1-xbar2)
  T2 <- as.numeric(T2)

  tr <- function(x) sum(diag(x)) # To obtain the trace easily

  # To obtain v
  b1 <- S1 %*% solve(S) # An auxiliar element
  b2 <- S2 %*% solve(S) # An auxiliar element
  v1 <- tr(b1 %*% b1) + (tr(b1))^2
  v1 <- v1 / n1
  v2 <- tr(b2 %*% b2) + (tr(b2))^2
  v2 <- v2 / n2
  v <- (p + p^2) / (v1 + v2)

  p.value <- pf(q=T2*(v-p+1)/(v*p), df1=p, df2=v-p+1, lower.tail=FALSE)

  method <- 'Modified Nel and Van der Merwe test for two mean vectors'
  statistic <- c(T2, T2*(v-p+1)/(v*p))
  names(statistic) <- c('T2', 'F')
  parameter <- c(p, v-p+1)
  names(parameter) <- c('df1', 'df2')
  alternative <- "mu1 is not equal to mu2 \n"
  estimate <- cbind(xbar1, xbar2)
  colnames(estimate) <- c('Sample 1', 'Sample 2')
  rownames(estimate) <- paste('xbar', 1:p, sep='_')
  data.name <- 'this test uses summarized data'

  return(list(statistic = statistic,
              parameter = parameter,
              p.value = p.value,
              estimate = estimate,
              alternative = alternative,
              method = method,
              data.name = data.name))
}
