#' Tests for homogeneity of covariances matrices
#'
#' The function implements the test for
#'  \eqn{H_0: \Sigma_1 = \Sigma_2 = ... = \Sigma_g} versus
#'  \eqn{H_1} at least one \eqn{\Sigma_i} is different.
#'
#' @param S a list with the sample covariance matrices.
#' @param N a list with the sample sizes.
#' @param method a character string specifying the method, "box" (default),
#' please see the details section for other methods.

#' @details the \code{"method"} must be one of \code{"box"} (default),
#' \code{"modified_LRT"}, \code{"wald_schott"}.
#'
#' To know in detail all tests implemented here the reader can visit
#' the vignette \url{https://fhernanb.github.io/stests/articles/Tests_Sigmas.html}
#'
#' @return A list with class \code{"htest"} containing the following components:
#' \item{statistic}{the value of the statistic.}
#' \item{parameter}{the degrees of freedom for the test.}
#' \item{p.value}{the p-value for the test.}
#' \item{estimate}{the estimated mean vectors.}
#' \item{method}{a character string indicating the type of test performed.}
#'
#' @references
#' Schott, J. R. (2001). Some tests for the equality of covariance matrices.
#' Journal of statistical planning and inference, 94(1), 25-36.
#'
#' Schott, J. R. (2007). A test for the equality of covariance matrices
#' when the dimension is large relative to the sample sizes.
#' Computational Statistics & Data Analysis, 51(12), 6535-6542.
#'
#' Mardia, K. V., Kent, J. T., & Bibby, J. M. (1979). Multivariate analysis.
#'
#' @example examples/examples_mult_var_matrices.R
#'
#' @author Freddy Hernandez.
#'
#' @importFrom stats pf
#' @export
mult_var_matrices_test <- function(S, N, method="box") {

  if (! var(unlist(lapply(S, dim))) == 0)
    stop("All matrices in list s must have identical dimension")

  if (! identical(length(S), length(N)))
    stop("The length of the lists S and N do not match")

  method <- match.arg(arg=method,
                      choices=c("box",
                                "modified_LRT",
                                "wald_schott"))

  # To generate the code for evaluating, without using cases
  my_code <- paste0("mult_var_matrices_test_", method,
                    "(S, N)")

  # To obtain the result
  result <- eval(parse(text=my_code))
  class(result) <- "htest"
  return(result)
}
#' @importFrom stats pchisq
mult_var_matrices_test_box <- function(S, N) {

  # Box test
  # Taken from Mardia et al., 1979, page 140

  # S: list with the var() output
  # N: list with the number of observations by group

  # Important convention
  # N_i: sample sizes for i-th group
  # n_i: N_i - 1

  p <- nrow(S[[1]])   # number of variables
  g <- length(S)      # number of groups
  v <- sum(unlist(N)) - g
  vg <- unlist(N) - 1
  aux1 <- sum(1/vg) - 1/v
  rho <- 1 - (2*p^2+3*p-1) * aux1 / (6*(p+1)*(g-1))
  # To obtain log(lambda)
  aux2 <- Map("*", S, vg)
  Sp <- Reduce("+", aux2) / v
  aux3 <- sum(log(unlist(lapply(S, det))) * vg)
  log_lambda <- (v * log(det(Sp)) - aux3) / (-2)
  # The statistic
  fi <- -2 * rho * log_lambda
  p.value <- pchisq(q=fi, df=p*(p+1)*(g-1)/2, lower.tail=FALSE)

  method <- "Box test for homogeneity of covariances"
  statistic <- fi
  names(statistic) <- c("phi")
  parameter <- c(p*(p+1)*(g-1)/2)
  names(parameter) <- c("df")
  alternative <- "at least one covariance matrix is different \n"
  estimate <- sprintf("Due to the high value of m, matrices S1, ..., S%d are not displayed.", g)
  data.name <- "this test uses summarized data"

  return(list(statistic = statistic,
              parameter = parameter,
              p.value = p.value,
              estimate = estimate,
              alternative = alternative,
              method = method,
              data.name = data.name))
}
#' @importFrom stats pchisq
mult_var_matrices_test_modified_LRT <- function(S, N) {

  # Bartlett's test or modified LRT
  # Taken from Schott (2001) page 26

  # S: list with the var() output
  # N: list with the number of observations by group

  # Important convention
  # N_i: sample sizes for i-th group
  # n_i: N_i - 1

  n <- unlist(N) - 1   # to obtain the n_i's
  m  <- nrow(S[[1]])   # number of variables
  g <- length(S)       # number of groups

  # Some check to avoid singular S_i
  if (any(unlist(N) < m))
    stop("With Ni < m, the sample covariance matrix Si is singular.
         See Schott (2007) A test for the equality of covariance matrices")

  # Obtaining S_pooled
  S_pooled <- 0
  for (i in 1:g) {
    S_pooled <- S_pooled + n[i] * S[[i]] / sum(n)
  }

  # The statistic and p-value
  M <- sum(n) * log(det(S_pooled))
  for (i in 1:g) {
    M <- M - n[i] * log(det(S[[i]]))
  }
  p.value <- pchisq(q=M, df=(g-1)*m*(m+1)/2, lower.tail=FALSE)

  # The usual way to report a test in R is:
  method <- "Modified likelihood ratio (or Bartlett's) test for homogeneity of covariances"
  statistic <- M
  names(statistic) <- c("M")
  parameter <- c((g-1)*m*(m+1)/2)
  names(parameter) <- c("df")
  alternative <- "at least one covariance matrix is different \n"
  estimate <- sprintf("Due to the high value of m, matrices S1, ..., S%d are not displayed.", g)
  data.name <- "this test uses summarized data"

  return(list(statistic = statistic,
              parameter = parameter,
              p.value = p.value,
              estimate = estimate,
              alternative = alternative,
              method = method,
              data.name = data.name))
}
#' @importFrom stats pchisq
mult_var_matrices_test_wald_schott <- function(S, N) {

  # Taken from Schott (2001) page 27

  # S: list with the var() output
  # N: list with the number of observations by group

  # Important convention
  # N_i: sample sizes for i-th group
  # n_i: N_i - 1

  n <- unlist(N) - 1
  m  <- nrow(S[[1]])   # number of variables
  g <- length(S)       # number of groups

  # Obtaining S_pooled
  S_pooled <- 0
  for (i in 1:g) {
    S_pooled <- S_pooled + n[i] * S[[i]] / sum(n)
  }

  # Auxiliary function
  my_trace <- function(x) sum(diag(x))

  # To obtain the two parts to create the statistic W
  part1 <- 0
  part2 <- 0

  for (i in 1:g)
    part1 <- part1 + my_trace(S[[i]]%*%solve(S_pooled)%*%S[[i]]%*%solve(S_pooled)) * n[i] / sum(n)

  for (i in 1:g)
    for (j in 1: g)
      part2 <- my_trace(S[[i]]%*%solve(S_pooled)%*%S[[i]]%*%solve(S_pooled)) * n[i]*n[j]/(sum(n)*sum(n))

  # The statistic and p-value
  W <- sum(n) * (part1 - part2)
  p.value <- pchisq(q=W, df=(g-1)*m*(m+1)/2, lower.tail=FALSE)

  # The usual way to report a test in R is:
  method <- "Wald-Schott test for homogeneity of covariances"
  statistic <- W
  names(statistic) <- c("W")
  parameter <- c((g-1)*m*(m+1)/2)
  names(parameter) <- c("df")
  alternative <- "the covariance matrices are different \n"
  estimate <- sprintf("Due to the high value of m, matrices S1, ..., S%d are not displayed.", g)
  data.name <- "this test uses summarized data"

  return(list(statistic = statistic,
              parameter = parameter,
              p.value = p.value,
              estimate = estimate,
              alternative = alternative,
              method = method,
              data.name = data.name))
}
