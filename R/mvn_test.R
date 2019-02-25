#' Modified Nel and Van der Merwe test
#'
#' This function implements the Modified Nel and Van der Merwe test the Behrens-Fisher problem \eqn{H0: \mu 1 = \mu 2} versus \eqn{H1: \mu1} not = \eqn{\mu2} without assuming anything about the covariance matrices.
#'
#' @param xbar1 a vector with the sample mean from population 1.
#' @param Sigma1 a matrix with sample variances and covariances from population 1.
#' @param n1 sample size 1.
#' @param xbar2 a vector with the sample mean from population 2.
#' @param Sigma2 a matrix with sample variances and covariances from population 2.
#' @param n2 sample size 2.
#'
#' @return A list with class \code{"htest"}.
#' @author Freddy Hernandez
#' @examples
#' # Example 4.1 from Nel and Van de Merwe (1986) page 3729
#' # Test H0: mu1 = mu2 versus H1: mu1 != mu2
#' n1 <- 45
#' xb1 <- c(204.4, 556.6)
#' s1 <- matrix(c(13825.3, 23823.4, 23823.4, 73107.4), ncol=2)
#'
#' n2 <- 55
#' xb2 <- c(130.0, 355.0)
#' s2 <- matrix(c(8632.0, 19616.7, 19616.7, 55964.5), ncol=2)
#'
#' mvn_test(xbar1=xb1, Sigma1=s1, n1=n1, xbar2=xb2, Sigma2=s2, n2=n2)
#'
#' \dontrun{
#' # Example using simulated data -----
#' # Parameters for the simulation
#' n1 <- n2 <- 500
#' mu1 <- c(0, 0)
#' sigma1 <- matrix(c(1, 0.5, 0.5, 1), ncol=2)
#' mu2 <- c(0, 0)
#' sigma2 <- matrix(c(1, -0.5, -0.5, 1), ncol=2)
#'
#' # Simulating the data
#' library(MASS)
#' dt1 <- mvrnorm(n1, mu=mu1, Sigma=sigma1)
#' dt2 <- mvrnorm(n2, mu=mu2, Sigma=sigma2)
#'
#' Compositional::james(dt1, dt2, R=2)
#' mvn_test(xbar1=colMeans(dt1), Sigma1=var(dt1), n1=nrow(dt1),
#'          xbar2=colMeans(dt2), Sigma2=var(dt2), n2=nrow(dt2))
#'  }
#'
#' @importFrom stats pf
#' @export
mvn_test <- function(xbar1, Sigma1, n1, xbar2, Sigma2, n2) {

  p <- ncol(Sigma1)
  xbar1 <- matrix(xbar1, ncol=1)
  xbar2 <- matrix(xbar2, ncol=1)

  S1 <- Sigma1/n1 # Represents S1 tilde
  S2 <- Sigma2/n2 # Represents S2 tilde
  S <- S1 + S2    # Represents S tilde

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

  parameter <- c(p, v-p+1)
  names(parameter) <- c('df1', 'df2')
  method <- 'Modified Nel and Van der Merwe test for two mean vectors'
  statistic <- c(T2, T2*(v-p+1)/(v*p))
  names(statistic) <- c('T2', 'F')

  alternative <- "mu1 is not equal to mu2 \n"
  estimate <- cbind(xbar1, xbar2)
  colnames(estimate) <- c('Sample 1', 'Sample 2')
  rownames(estimate) <- paste('xbar', 1:p, sep='_')
  data.name <- 'this test uses summarized data'

  res <- list(statistic=statistic,
              parameter=parameter,
              p.value=p.value,
              estimate=estimate,
              alternative=alternative,
              method=method,
              data.name=data.name
  )

  class(res) <- "htest"
  res

}
