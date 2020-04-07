#' Tests for homogeneity of covariances matrices
#'
#' The function implements the test for \eqn{H_0: \Sigma_1 = \Sigma_2 = ... = \Sigma_q} versus \eqn{H_1} at least one \eqn{\Sigma_i} is different.
#'
#' @param s a list with the sample covariance matrices.
#' @param n a list with the sample sizes.
#' @param method a character string specifying the method, "box" (default), please see the details section for other methods.

#' @details the \code{"method"} must be one of \code{"box"} (default), \code{"xxx"} (pronto otros metodos).
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
#' s1 <- matrix(c(12.65, -16.45,
#'                -16.45, 73.04), ncol=2, nrow=2)
#' s2 <- matrix(c(11.44, -27.77,
#'                -27.77, 100.64), ncol=2, nrow=2)
#' s3 <- matrix(c(14.46, -31.26,
#'                -31.26, 101.03), ncol=2, nrow=2)
#' n1 <- 26
#' n2 <- 23
#' n3 <- 25
#' s <- list(s1, s2, s3)
#' n <- list(n1, n2, n3)
#'
#' res <- mult_var_matrices_test(s, n, method="box")
#' res
#' plot(res, shade.col='tomato')
#'
#' # Example 5.3.4 from Mardia (1979) page 141
#' s1 <- matrix(c(132.99, 75.85, 35.82,
#'                75.85, 47.96, 20.75,
#'                35.82, 20.75, 10.79), ncol=3, nrow=3)
#' s2 <- matrix(c(432.58, 259.87, 161.67,
#'                259.87, 164.57, 98.99,
#'                161.67, 98.99, 63.87), ncol=3, nrow=3)
#' n1 <- 24
#' n2 <- 24
#' s <- list(s1, s2)
#' n <- list(n1, n2)
#'
#' res <- mult_var_matrices_test(s, n, method="box")
#' res
#' plot(res, from=20, to=30, shade.col='pink')
#'
#' @importFrom stats pf
#' @export
mult_var_matrices_test <- function(s, n, method='box') {

  if (! var(unlist(lapply(s, dim))) == 0)
    stop("All matrices in list s must have identical dimension")

  if (! identical(length(s), length(n)))
    stop("The length of the lists s and n do not match")

  method <- match.arg(arg=method,
                      choices=c("box"))

  # To generate the code for evaluating, without using cases
  my_code <- paste0("mult_var_matrices_test_", method,
                    "(s, n)")

  # To obtain the result
  result <- eval(parse(text=my_code))
  class(result) <- "htest"
  return(result)
}
#' @importFrom stats pchisq
mult_var_matrices_test_box <- function(s, n) {
  p <- nrow(s[[1]])   # number of variables
  g <- length(s)      # number of groups
  N <- sum(unlist(n))
  v <- N - g
  vg <- unlist(n) - 1
  aux1 <- sum(1/vg) - 1/v
  rho <- 1 - (2*p^2+3*p-1) * aux1 / (6*(p+1)*(g-1))
  # To obtian log(lambda)
  aux2 <- Map("*", s, vg)
  Sp <- Reduce("+", aux2) / v
  aux3 <- sum(log(unlist(lapply(s, det))) * vg)
  log_lambda <- (v * log(det(Sp)) - aux3) / (-2)
  # The statistic
  fi <- -2 * rho * log_lambda
  p.value <- pchisq(q=fi, df=p*(p+1)*(g-1)/2, lower.tail=FALSE)

  method <- 'Box test for homogeneity of covariances'
  statistic <- fi
  names(statistic) <- c('phi')
  parameter <- c(p*(p+1)*(g-1)/2)
  names(parameter) <- c('df')
  alternative <- "at least one covariance matrix is different \n"
  estimate <- NULL
  data.name <- 'this test uses summarized data'

  return(list(statistic = statistic,
              parameter = parameter,
              p.value = p.value,
              estimate = estimate,
              alternative = alternative,
              method = method,
              data.name = data.name))
}
