#' Tests for Equality of Two Normal Mean Vectors
#'
#' The function implements the test for \eqn{H_0: \mu_1 = \mu_2} versus \eqn{H_1: \mu_1} not = \eqn{\mu_2} when two random samples are obtained from two p-variate normal populations \eqn{Np(\mu_1, \Sigma_1)} and \eqn{Np(\mu_2, \Sigma_2)} respectively. By default, the function performs the Hotelling test for two normal mean vectors assumming equality in the covariance matrices. Other tests for the multivariate Behrens-Fisher problem are also implemented, see the argument \code{method} below.
#'
#' @param xbar1 a vector with the sample mean from population 1.
#' @param s1 a matrix with sample variances and covariances from population 1.
#' @param n1 sample size 1.
#' @param xbar2 a vector with the sample mean from population 2.
#' @param s2 a matrix with sample variances and covariances from population 2.
#' @param n2 sample size 2.
#' @param delta0 a number indicating the true value of the difference in means.
#' @param method a character string specifying the method, \code{"T2"} (default), \code{"james"} (James' first order test), please see the details section for other methods.
#' @param alpha the significance level only for method \code{"james"}, by default its value is 0.05.
#'
#' @details the \code{"method"} must be one of \code{"T2"} (default), \code{"james"} (James' first order test), \code{"yao"} (Yao's test), \code{"johansen"} (Johansen's test), \code{"nvm"} (Nel and Van der Merwe test), \code{"mnvm"} (modified Nel and Van der Merwe test), \code{"gamage"} (Gamage's test), \code{"yy"} (Yanagihara and Yuan test), \code{"byy"} (Bartlett Correction test), \code{"ks1"} (Second Order Procedure), \code{"ks2"} (Bias Correction Procedure). For James test the critic value is reported, we reject H0 if T2 > critic_value.
#'
#' @return A list with class \code{"htest"} containing the following components:
#' \item{statistic}{the value of the statistic.}
#' \item{parameter}{the degrees of freedom for the test.}
#' \item{p.value}{the p-value for the test.}
#' \item{estimate}{the estimated mean vectors.}
#' \item{method}{a character string indicating the type of test performed.}
#'
#' @author Freddy Hernandez, Jean Paul Piedrahita, Valentina Garcia.
#'
#' @example examples/examples_two_mean_vector_test.R
#'
#' @importFrom stats pf
#' @export
two_mean_vector_test <- function(xbar1, s1, n1, xbar2, s2, n2,
                                 delta0=NULL, method="T2", alpha=0.05) {

  if (! identical(dim(s1), dim(s2)))
    stop("The matrices s1 and s2 do not have same dimension")

  method <- match.arg(arg=method,
                      choices=c("T2", "james", "yao", "johansen",
                                "nvm", "mnvm", "gamage", "yy",
                                "byy", "mbyy", "ks1", "ks2"))

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
#' @importFrom stats qchisq
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
  v1 <- tr(S1 %*% S1) + (tr(S1))^2
  v1 <- v1 / (n1-1)
  v2 <- tr(S2 %*% S2) + (tr(S2))^2
  v2 <- v2 / (n2-1)
  v <- (tr(S %*% S) + tr(S)^2) / (v1 + v2)

  p.value <- pf(q=T2*(v-p+1)/(v*p), df1=p, df2=v-p+1, lower.tail=FALSE)

  method <- 'Nel and Van der Merwe test for two mean vectors'
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
two_mean_vector_test_mnvm <- function(xbar1, s1, n1,
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
  v1 <- v1 / (n1-1)
  v2 <- tr(b2 %*% b2) + (tr(b2))^2
  v2 <- v2 / (n2-1)
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
#' @importFrom stats rchisq
two_mean_vector_test_gamage <- function(xbar1, s1, n1,
                                        xbar2, s2, n2,
                                        nrep=2000,
                                        delta0=NULL, alpha=0.05) {

  p <- ncol(s1)
  xbar1 <- matrix(xbar1, ncol=1)
  xbar2 <- matrix(xbar2, ncol=1)
  S1 <- s1/n1    # Represents S1 tilde, uppercase
  S2 <- s2/n2    # Represents S2 tilde, uppercase
  S  <- S1 + S2  # Represents S tilde

  T2 <- t(xbar1-xbar2) %*% solve(S) %*% (xbar1-xbar2)
  T2 <- as.numeric(T2)

  # The next function computes the square root of the positive definite matrix
  my_sqrt <- function(A) {
    L <- diag(eigen(A)$values)
    P <- eigen(A)$vectors
    P %*% L %*% t(P)  # Este producto da como resultado la matriz A
    A_raiz_cuadrada <- P %*% L^0.5 %*% t(P)
    A_raiz_cuadrada
  }

  # The next function implements pvalue given in expression 3.9 of Gamage (2004)
  one_gen_pvalue <- function(v1, n1, n2, p) {
    d <- eigen(v1)$values
    z0i2 <- rchisq(n=p, df=1)
    Q1 <- rchisq(n=1, df=n1-p)
    Q2 <- rchisq(n=1, df=n2-p)
    T1 <- sum(d * z0i2) / Q1 + sum((1-d/(n1-1)) * z0i2) * (n2-1) / Q2
    T1
  }

  v1 <- ((n1-1)/n1) * solve(my_sqrt(S)) %*% S1 %*% solve(my_sqrt(S))
  T1 <- replicate(n=2000, one_gen_pvalue(v1, n1, n2, p))
  p.value <- mean(T1 > T2)

  method <- 'Gamage test for two mean vectors'
  statistic <- c(T2)
  names(statistic) <- c('T2')
  parameter <- NULL
  names(parameter) <- NULL
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
              data.name = data.name,
              T1 = T1))
}
#' @importFrom stats pf
two_mean_vector_test_yy <- function(xbar1, s1, n1,
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

  n <- n1+n2
  N <- n-2

  Sl1 <- (n2/n) * s1      # Represents S1 line
  Sl2 <- (n1/n) * s2      # Represents S2 line
  Sl <- Sl1 + Sl2         # Represents S line

  # Auxiliar elements
  a1 <- (n2^2 * (n-2))/(n^2 * (n1-1))
  a2 <- (n1^2 * (n-2))/(n^2 * (n2-1))
  b1 <- s1 %*% solve(Sl)
  b2 <- s2 %*% solve(Sl)

  # To obtain psi1 and psi2
  psi1 <- a1 * tr(b1)^2 + a2 * tr(b2)^2
  psi2 <- a1 * tr(b1%*%b1) + a2 * tr(b2%*%b2)

  # To obtain theta1 and theta2
  theta1 <- p * psi1 + (p-2)*psi2
  theta1 <- theta1/(p*(p+2))
  theta2 <- psi1 + 2*psi2
  theta2 <- theta2/(p*(p+2))

  # To obtain the v
  v1 <- (n-2-theta1)^2
  v2 <- (n-2)*theta2 - theta1
  v <- v1/v2

  p.value <- pf(q=T2*(n-2-theta1)/((n-2)*p), df1=p, df2=v, lower.tail=FALSE)

  method <- 'Yanagihara and Yuan test for two mean vectors'
  statistic <- c(T2, T2*(n-2-theta1)/((n-2)*p))
  names(statistic) <- c('T2', 'F')
  parameter <- c(p, v)
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
#' @importFrom stats pchisq
two_mean_vector_test_byy <- function(xbar1, s1, n1,
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

  n <- n1+n2
  N <- n-2

  Sl1 <- (n2/n) * s1      # Represents S1 line
  Sl2 <- (n1/n) * s2      # Represents S2 line
  Sl <- Sl1 + Sl2         # Represents S line

  # Auxiliar elements
  a1 <- (n2^2 * (n-2))/(n^2 * (n1-1))
  a2 <- (n1^2 * (n-2))/(n^2 * (n2-1))
  b1 <- s1 %*% solve(Sl)
  b2 <- s2 %*% solve(Sl)

  # To obtain psi1 and psi2
  psi1 <- a1 * tr(b1)^2 + a2 * tr(b2)^2
  psi2 <- a1 * tr(b1%*%b1) + a2 * tr(b2%*%b2)

  # To obtain c1 and c2
  c1 <- (psi1+psi2)/p
  c2 <- 2*(p+3)*psi1 + 2*(p+4)*psi2
  c2 <- c2/(p*(p+2))

  p.value <- pchisq(q=T2*(N-c1)/N, df=p, lower.tail=FALSE)

  method <- 'Bartlett Correction test for two mean vectors'
  statistic <- c(T2, T2*(N-c1)/N)
  names(statistic) <- c('T2', 'X-squared')
  parameter <- p
  names(parameter) <- 'df'
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
#' @importFrom stats pchisq
two_mean_vector_test_mbyy <- function(xbar1, s1, n1,
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

  n <- n1+n2
  N <- n-2

  Sl1 <- (n2/n) * s1      # Represents S1 line
  Sl2 <- (n1/n) * s2      # Represents S2 line
  Sl <- Sl1 + Sl2         # Represents S line

  # Auxiliar elements
  a1 <- (n2^2 * (n-2))/(n^2 * (n1-1))
  a2 <- (n1^2 * (n-2))/(n^2 * (n2-1))
  b1 <- s1 %*% solve(Sl)
  b2 <- s2 %*% solve(Sl)

  # To obtain psi1 and psi2
  psi1 <- a1 * tr(b1)^2 + a2 * tr(b2)^2
  psi2 <- a1 * tr(b1%*%b1) + a2 * tr(b2%*%b2)

  # To obtain c1 and c2
  c1 <- (psi1+psi2)/p
  c2 <- 2*(p+3)*psi1 + 2*(p+4)*psi2
  c2 <- c2/(p*(p+2))

  # To obtain beta1 y beta2
  beta1 <- 2/(c2 - 2*c1)
  beta2 <- (p+2)*c2 - 2*(p+4)*c1
  beta2 <- beta2/(2*(c2 - 2*c1))

  p.value <- pchisq(q=(N*beta1 + beta2)*log(1 + T2/(N*beta1)), df=p, lower.tail=FALSE)

  method <- 'Modified Bartlett Correction test for two mean vectors'
  statistic <- c(T2, (N*beta1 + beta2)*log(1 + T2/(N*beta1)))
  names(statistic) <- c('T2', 'X-squared')
  parameter <- p
  names(parameter) <- 'df'
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
#' @importFrom expm %^%
two_mean_vector_test_ks1 <- function(xbar1, s1, n1,
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

  tr <- function(x) sum(diag(x))   # To obtain the trace easily

  n <- n1+n2
  N <- n-2

  Sl1 <- (n2/n) * s1      # Represents S1 line
  Sl2 <- (n1/n) * s2      # Represents S2 line
  Sl <- Sl1 + Sl2         # Represents S line

  # Auxiliar elements
  x1 <- s1 %*% solve(Sl)
  x2 <- s2 %*% solve(Sl)

  # Function a (return a number)
  a <- function(i, l){
    if (i == 1){
      res <- tr(x1)^l
    }

    else if (i == 2){
      res <- tr(x2)^l
    }

    return(res)
  }

  # Function b (return a number)
  b <- function(q, r, s){
    #library(expm)
    res <- tr((x1 %^% q) %*% (x2 %^% r))
    res <- res^s
    return(res)
  }

  # c values (numbers)
  c1 <- (n - n1)^2 * (n-2)
  c1 <- c1/(n^2 * (n1-1))

  c2 <- (n - n2)^2 * (n-2)
  c2 <- c2/(n^2 * (n2-1))

  # d values (numbers)
  d1 <- (n - n1)^3 * (n-2)^2
  d1 <- d1/(n^3 * (n1-1)^2)

  d2 <- (n - n2)^3 * (n-2)^2
  d2 <- d2/(n^3 * (n2-1)^2)

  # psi values
  psi1 <- c1 * c2 * (a(1,1)*b(1,2,1) + a(2,1)*b(2,1,1))
  psi2 <- c1 * c2 * b(2,2,1)
  psi3 <- c1 * c2 * a(1,1) * a(2,1) * b(1,1,1)
  psi4 <- c1 * c2 * a(1,2) * a(2,2)
  psi5 <- c1 * c2 * (a(1,2)*a(2,1)^2 + a(1,1)^2*a(2,2))
  psi6 <- c1 * c2 * b(1,1,1)^2
  psi7 <- c1 * c2 * a(1,1)^2 * a(2,1)^2
  psi8 <- c1 * c2 * b(1,1,2)

  # --------------------------------------------------------------------------
  # theta values (numbers)

  # theta1
  theta1 <- c1 * (p*a(1,1)^2 + (p-2)*a(1,2)) + c2 * (p*a(2,1)^2 + (p-2)*a(2,2))
  theta1 <- theta1/(p*(p+2))

  # theta2
  th2_1 <- d1 * (4*p^2 * a(1,3) + (p-2)*(3*p + 4)*a(1,1)*a(1,2) + p*(p+2)*a(1,1)^3)
  th2_2 <- d2 * (4*p^2 * a(2,3) + (p-2)*(3*p + 4)*a(2,1)*a(2,2) + p*(p+2)*a(2,1)^3)
  theta2 <- (th2_1 + th2_2)/(p*(p+2)*(p+4))

  # theta3
  th3_1 <- p^2*(5*p + 14)*a(1,4) + 4*(p+3)*(p+2)*(p-2)*a(1,1)*a(1,3)
  th3_1 <- th3_1 + p*(p+3)*(p-2)*a(1,2)^2 + 2*(p^3 + 5*p^2 + 7*p + 6)*a(1,2)*a(1,1)^2
  th3_1 <- th3_1 - p*(p+4)*a(1,1)^4
  th3_1 <- c1^2 * (th3_1)

  th3_2 <- p^2*(5*p + 14)*a(2,4) + 4*(p+3)*(p+2)*(p-2)*a(2,1)*a(2,3)
  th3_2 <- th3_2 + p*(p+3)*(p-2)*a(2,2)^2 + 2*(p^3 + 5*p^2 + 7*p + 6)*a(2,2)*a(2,1)^2
  th3_2 <- th3_2 - p*(p+4)*a(2,1)^4
  th3_2 <- c2^2 * (th3_2)

  aux3 <- 4*(p+3)*(p+2)*(p-2)*psi1 + 4*p*(p+2)*(p-2)*psi2
  aux3 <- aux3 + 4*p*(p+4)*(p+2)*psi3 - 2*p*(p-2)*psi4 - 2*(p+3)*(p-2)*psi5
  aux3 <- aux3 + 2*p*(p+4)*(p-2)*psi6 - 2*p*(p+4)*psi7 + 2*p*(p+4)*(3*p + 2)*psi8

  theta3 <- th3_1 + th3_2 + aux3
  theta3 <- theta3/(p*(p+2)*(p+4)*(p+6))

  # theta4
  theta4 <- c1 * (a(1,1)^2 + 2*a(1,2)) + c2 * (a(2,1)^2 + 2*a(2,2))
  theta4 <- theta4/(p*(p+2))

  # theta5
  th5_1 <- 4*(p^2 - 3*p + 4)*a(1,3) + 3*p*(p-4)*a(1,1)*a(1,2) + p^2*a(1,1)^3
  th5_1 <- d1*th5_1

  th5_2 <- 4*(p^2 - 3*p + 4)*a(2,3) + 3*p*(p-4)*a(2,1)*a(2,2) + p^2*a(2,1)^3
  th5_2 <- d2*th5_2

  theta5 <- (th5_1 + th5_2)/(p*(p+2)*(p+4))

  # theta6
  th6_1 <- 2*(p+1)*(5*p^2 - 14*p + 24)*a(1,4)
  th6_1 <- th6_1 + 4*(p-4)*(2*p^2 + 5*p + 6)*a(1,1)*a(1,3) + (p-2)*(p-4)*(2*p + 3)*a(1,2)^2
  th6_1 <- th6_1 + 2*(p+2)*(2*p^2 - p + 12)*a(1,2)*a(1,1)^2 - 3*(p^2 + 2*p - 4)*a(1,1)^4
  th6_1 <- c1^2 * th6_1

  th6_2 <- 2*(p+1)*(5*p^2 - 14*p + 24)*a(2,4)
  th6_2 <- th6_2 + 4*(p-4)*(2*p^2 + 5*p + 6)*a(2,1)*a(2,3) + (p-2)*(p-4)*(2*p + 3)*a(2,2)^2
  th6_2 <- th6_2 + 2*(p+2)*(2*p^2 - p + 12)*a(2,2)*a(2,1)^2 - 3*(p^2 + 2*p - 4)*a(2,1)^4
  th6_2 <- c2^2 * th6_2

  aux6 <- 4*(p-4)*(2*p^2 + 5*p + 6)*psi1 + 8*p*(p-2)*(p-4)*psi2 + 8*p*(p^2 + 4*p + 2)*psi3
  aux6 <- aux6 - 6*(p-2)*(p-4)*psi4 - 6*(p-4)*(p+2)*psi5 + 4*(p+3)*(p-2)*(p-4)*psi6
  aux6 <- aux6 - 6*(p^2 + 2*p - 4)*psi7 + 12*(p^3 + p^2 - 2*p + 8)*psi8

  theta6 <- th6_1 + th6_2 + aux6
  theta6 <- theta6/(p*(p+2)*(p+4)*(p+6))

  # To obtain Vs
  Vs1 <- 2*(N^2 - N*theta1 + theta2 - theta3)^2
  Vs2 <- N^2 * (N^2 - 2*N*theta1 + 2*N*theta4 + 2*theta5 - theta6)
  Vs3 <- (N^2 - N*theta1 + theta2 - theta3)^2
  Vs <- Vs1/(Vs2 - Vs3)

  # To obtain psi_s
  psi_s <- N^2 * Vs
  psi_s <- psi_s/(N^2 - N*theta1 + theta2 - theta3)

  p.value <- pf(q=T2*Vs/(p*psi_s), df1=p, df2=Vs, lower.tail=FALSE)

  method <- 'Kawasaki and Seo (Second order) test for two mean vectors'
  statistic <- c(T2, T2*Vs/(p*psi_s))
  names(statistic) <- c('T2', 'F')
  parameter <- c(p, Vs)
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
#' @importFrom expm %^%
# Kawasaki and Seo test 2 (Bias Correction procedure)
two_mean_vector_test_ks2 <- function(xbar1, s1, n1,
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

  tr <- function(x) sum(diag(x))   # To obtain the trace easily

  n <- n1+n2
  N <- n-2

  Sl1 <- (n2/n) * s1      # Represents S1 line
  Sl2 <- (n1/n) * s2      # Represents S2 line
  Sl <- Sl1 + Sl2         # Represents S line

  # Auxiliar elements
  x1 <- s1 %*% solve(Sl)
  x2 <- s2 %*% solve(Sl)

  # Function a (return a number)
  a <- function(i, l){
    if (i == 1){
      res <- tr(x1)^l
    }

    else if (i == 2){
      res <- tr(x2)^l
    }

    return(res)
  }

  # Function b (return a number)
  b <- function(q, r, s){
    #library(expm)
    res <- tr((x1 %^% q) %*% (x2 %^% r))
    res <- res^s
    return(res)
  }

  # c values (numbers)
  c1 <- (n - n1)^2 * (n-2)
  c1 <- c1/(n^2 * (n1-1))

  c2 <- (n - n2)^2 * (n-2)
  c2 <- c2/(n^2 * (n2-1))

  # d values (numbers)
  d1 <- (n - n1)^3 * (n-2)^2
  d1 <- d1/(n^3 * (n1-1)^2)

  d2 <- (n - n2)^3 * (n-2)^2
  d2 <- d2/(n^3 * (n2-1)^2)

  # rho values (numbers)
  rho1 <- sqrt((n1-1)/(n-2))
  rho2 <- sqrt((n2-1)/(n-2))

  # psi values
  psi1 <- c1 * c2 * (a(1,1)*b(1,2,1) + a(2,1)*b(2,1,1))
  psi2 <- c1 * c2 * b(2,2,1)
  psi3 <- c1 * c2 * a(1,1) * a(2,1) * b(1,1,1)
  psi4 <- c1 * c2 * a(1,2) * a(2,2)
  psi5 <- c1 * c2 * (a(1,2)*a(2,1)^2 + a(1,1)^2*a(2,2))
  psi6 <- c1 * c2 * b(1,1,1)^2
  psi7 <- c1 * c2 * a(1,1)^2 * a(2,1)^2
  psi8 <- c1 * c2 * b(1,1,2)

  # --------------------------------------------------------------------------
  # theta values (numbers)

  # theta1
  theta1 <- c1 * (p*a(1,1)^2 + (p-2)*a(1,2)) + c2 * (p*a(2,1)^2 + (p-2)*a(2,2))
  theta1 <- theta1/(p*(p+2))

  # theta2
  th2_1 <- d1 * (4*p^2 * a(1,3) + (p-2)*(3*p + 4)*a(1,1)*a(1,2) + p*(p+2)*a(1,1)^3)
  th2_2 <- d2 * (4*p^2 * a(2,3) + (p-2)*(3*p + 4)*a(2,1)*a(2,2) + p*(p+2)*a(2,1)^3)
  theta2 <- (th2_1 + th2_2)/(p*(p+2)*(p+4))

  # theta3
  th3_1 <- p^2*(5*p + 14)*a(1,4) + 4*(p+3)*(p+2)*(p-2)*a(1,1)*a(1,3)
  th3_1 <- th3_1 + p*(p+3)*(p-2)*a(1,2)^2 + 2*(p^3 + 5*p^2 + 7*p + 6)*a(1,2)*a(1,1)^2
  th3_1 <- th3_1 - p*(p+4)*a(1,1)^4
  th3_1 <- c1^2 * (th3_1)

  th3_2 <- p^2*(5*p + 14)*a(2,4) + 4*(p+3)*(p+2)*(p-2)*a(2,1)*a(2,3)
  th3_2 <- th3_2 + p*(p+3)*(p-2)*a(2,2)^2 + 2*(p^3 + 5*p^2 + 7*p + 6)*a(2,2)*a(2,1)^2
  th3_2 <- th3_2 - p*(p+4)*a(2,1)^4
  th3_2 <- c2^2 * (th3_2)

  aux3 <- 4*(p+3)*(p+2)*(p-2)*psi1 + 4*p*(p+2)*(p-2)*psi2
  aux3 <- aux3 + 4*p*(p+4)*(p+2)*psi3 - 2*p*(p-2)*psi4 - 2*(p+3)*(p-2)*psi5
  aux3 <- aux3 + 2*p*(p+4)*(p-2)*psi6 - 2*p*(p+4)*psi7 + 2*p*(p+4)*(3*p + 2)*psi8

  theta3 <- th3_1 + th3_2 + aux3
  theta3 <- theta3/(p*(p+2)*(p+4)*(p+6))

  # theta4
  theta4 <- c1 * (a(1,1)^2 + 2*a(1,2)) + c2 * (a(2,1)^2 + 2*a(2,2))
  theta4 <- theta4/(p*(p+2))

  # theta5
  th5_1 <- 4*(p^2 - 3*p + 4)*a(1,3) + 3*p*(p-4)*a(1,1)*a(1,2) + p^2*a(1,1)^3
  th5_1 <- d1*th5_1

  th5_2 <- 4*(p^2 - 3*p + 4)*a(2,3) + 3*p*(p-4)*a(2,1)*a(2,2) + p^2*a(2,1)^3
  th5_2 <- d2*th5_2

  theta5 <- (th5_1 + th5_2)/(p*(p+2)*(p+4))

  # theta6
  th6_1 <- 2*(p+1)*(5*p^2 - 14*p + 24)*a(1,4)
  th6_1 <- th6_1 + 4*(p-4)*(2*p^2 + 5*p + 6)*a(1,1)*a(1,3) + (p-2)*(p-4)*(2*p + 3)*a(1,2)^2
  th6_1 <- th6_1 + 2*(p+2)*(2*p^2 - p + 12)*a(1,2)*a(1,1)^2 - 3*(p^2 + 2*p - 4)*a(1,1)^4
  th6_1 <- c1^2 * th6_1

  th6_2 <- 2*(p+1)*(5*p^2 - 14*p + 24)*a(2,4)
  th6_2 <- th6_2 + 4*(p-4)*(2*p^2 + 5*p + 6)*a(2,1)*a(2,3) + (p-2)*(p-4)*(2*p + 3)*a(2,2)^2
  th6_2 <- th6_2 + 2*(p+2)*(2*p^2 - p + 12)*a(2,2)*a(2,1)^2 - 3*(p^2 + 2*p - 4)*a(2,1)^4
  th6_2 <- c2^2 * th6_2

  aux6 <- 4*(p-4)*(2*p^2 + 5*p + 6)*psi1 + 8*p*(p-2)*(p-4)*psi2 + 8*p*(p^2 + 4*p + 2)*psi3
  aux6 <- aux6 - 6*(p-2)*(p-4)*psi4 - 6*(p-4)*(p+2)*psi5 + 4*(p+3)*(p-2)*(p-4)*psi6
  aux6 <- aux6 - 6*(p^2 + 2*p - 4)*psi7 + 12*(p^3 + p^2 - 2*p + 8)*psi8

  theta6 <- th6_1 + th6_2 + aux6
  theta6 <- theta6/(p*(p+2)*(p+4)*(p+6))

  # --------------------------------------------------------------------------------
  # nu values (numbers)

  # nu1
  nu1_1 <- c1 * rho1^-2 * ((p-2)*a(1,1)^2 + (3*p - 2)*a(1,2))
  nu1_2 <- c2 * rho2^-2 * ((p-2)*a(2,1)^2 + (3*p - 2)*a(2,2))
  nu1 <- (nu1_1 + nu1_2)/(p*(p+2))

  # nu2
  nu2_1 <- d1 * ((7*p - 6)*a(1,3) + (5*p - 6)*a(1,1)*a(1,2) + 2*p*a(1,1)^3)
  nu2_2 <- d2 * ((7*p - 6)*a(2,3) + (5*p - 6)*a(2,1)*a(2,2) + 2*p*a(2,1)^3)
  nu2 <- (nu2_1 + nu2_2)/(p*(p+2))

  # nu3
  nu3_1 <- 4*(p-1)*a(1,4) + (3*p - 2)*a(1,1)*a(1,3) + (p-2)*a(1,2)^2
  nu3_1 <- c1^2 * (nu3_1 + 2*p * a(1,1)^2 * a(1,2))

  nu3_2 <- 4*(p-1)*a(2,4) + (3*p - 2)*a(2,1)*a(2,3) + (p-2)*a(2,2)^2
  nu3_2 <- c2^2 * (nu3_1 + 2*p * a(2,1)^2 * a(2,2))

  nu3_aux <- (3*p - 2)*psi1 + 2*(p-2)*psi2 + 4*p*psi3 + 2*(p-2)*psi6 + 2*(3*p - 2)*psi8

  nu3 <- (nu3_1 + nu3_2 + nu3_aux)/(p*(p+2))

  # nu4
  nu4_1 <- c1 * rho1^-2 *(2*a(1,1)^2 + 4*a(1,2))
  nu4_2 <- c2 * rho2^-2 *(2*a(2,1)^2 + 4*a(2,2))
  nu4 <- (nu4_1 + nu4_2)/(p*(p+2))

  # nu5
  nu5_1 <- d1 * (10*a(1,3) + 8*a(1,1)*a(1,2) + 2*a(1,1)^3)
  nu5_2 <- d2 * (10*a(2,3) + 8*a(2,1)*a(2,2) + 2*a(2,1)^3)
  nu5 <- (nu5_1 + nu5_2)/(p*(p+2))

  # nu6
  nu6_1 <- 6*a(1,4) + 4*a(1,1)*a(1,3) + 2*a(1,2)^2 + 2*a(1,1)^2 * a(1,2)
  nu6_1 <- c1^2 * nu6_1

  nu6_2 <- 6*a(2,4) + 4*a(2,1)*a(2,3) + 2*a(2,2)^2 + 2*a(2,1)^2 * a(2,2)
  nu6_2 <- c2^2 * nu6_2

  nu6_aux <- 4*(psi1 + psi2 + psi3 + psi6 + 2*psi8)

  nu6 <- (nu6_1 + nu6_2 + nu6_aux)/(p*(p+2))

  # To obtain theta 1 asterisk
  theta1_asterisk <- nu1 - nu2 + nu3

  # To obtain theta 4 asterisk
  theta4_asterisk <- nu4 - nu5 + nu6

  # To obtain Vbc
  Vbc1 <- 2*(N^2 - N*theta1 + theta2 - theta3 - theta1_asterisk)^2
  Vbc2 <- N^2 * (N^2 - 2*N*theta1 + 2*N*theta4 + 2*theta5 - theta6 - 2*theta1_asterisk + 2*theta4_asterisk)
  Vbc3 <- (N^2 - N*theta1 + theta2 - theta3 - theta1_asterisk)^2
  Vbc <- Vbc1/(Vbc2 - Vbc3)

  # To obtain psi_bc
  psi_bc <- N^2 * Vbc
  psi_bc <- psi_bc/(N^2 - N*theta1 + theta2 - theta3 - theta1_asterisk)

  p.value <- pf(q=T2*Vbc/(p*psi_bc), df1=p, df2=Vbc, lower.tail=FALSE)

  method <- 'Kawasaki and Seo (Bias Correction) test for two mean vectors'
  statistic <- c(T2, T2*Vbc/(p*psi_bc))
  names(statistic) <- c('T2', 'F')
  parameter <- c(p, Vbc)
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
