#' Tests for homogeneity of covariances matrices
#'
#' The function implements the test for
#'  \eqn{H_0: \Sigma_1 = \Sigma_2 = ... = \Sigma_q} versus
#'  \eqn{H_1} at least one \eqn{\Sigma_i} is different.
#'
#' @param S a list with the sample covariance matrices.
#' @param N a list with the sample sizes.
#' @param method a character string specifying the method, "box" (default),
#' please see the details section for other methods.

#' @details the \code{"method"} must be one of \code{"box"} (default), \code{"xxx"}.
#'
#' @return A list with class \code{"htest"} containing the following components:
#' \item{statistic}{the value of the statistic.}
#' \item{parameter}{the degrees of freedom for the test.}
#' \item{p.value}{the p-value for the test.}
#' \item{estimate}{the estimated mean vectors.}
#' \item{method}{a character string indicating the type of test performed.}
#'
#' @author Freddy Hernandez.
#' @examples
#' # Example 5.2.3 from Diaz and Morales (2015) page 200
#' S1 <- matrix(c(12.65, -16.45,
#'                -16.45, 73.04), ncol=2, nrow=2)
#' S2 <- matrix(c(11.44, -27.77,
#'                -27.77, 100.64), ncol=2, nrow=2)
#' S3 <- matrix(c(14.46, -31.26,
#'                -31.26, 101.03), ncol=2, nrow=2)
#' N1 <- 26
#' N2 <- 23
#' N3 <- 25
#' S <- list(S1, S2, S3)
#' N <- list(N1, N2, N3)
#'
#' res <- mult_var_matrices_test(S, N, method="box")
#' res
#' plot(res, shade.col="tomato")
#'
#' # Example 5.3.4 from Mardia (1979) page 141
#' S1 <- matrix(c(132.99, 75.85, 35.82,
#'                75.85, 47.96, 20.75,
#'                35.82, 20.75, 10.79), ncol=3, nrow=3)
#' S2 <- matrix(c(432.58, 259.87, 161.67,
#'                259.87, 164.57, 98.99,
#'                161.67, 98.99, 63.87), ncol=3, nrow=3)
#' N1 <- 24
#' N2 <- 24
#' S <- list(S1, S2)
#' N <- list(N1, N2)
#'
#' res <- mult_var_matrices_test(S, N, method="box")
#' res
#' plot(res, from=20, to=30, shade.col="pink")
#'
#' @importFrom stats pf
#' @export
mult_var_matrices_test <- function(s, n, method="box") {

  if (! var(unlist(lapply(s, dim))) == 0)
    stop("All matrices in list s must have identical dimension")

  if (! identical(length(s), length(n)))
    stop("The length of the lists S and N do not match")

  method <- match.arg(arg=method,
                      choices=c("box",
                                "modified_LRT"))

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

  # Bartlett’s test or modified LRT
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
  method <- "Modified likelihood ratio (or Bartlett’s) test for homogeneity of covariances"
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
