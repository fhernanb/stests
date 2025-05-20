#' Permutation test using SD
#'
#' @author Freddy Hernandez, \email{fhernanb@unal.edu.co}
#'
#' This function performs test for dependent data in samples
#' with a small number of subjects.
#'
#' @param x numeric vector of data values.
#' @param x_sd numeric vector with standard deviations.
#' @param y numeric vector of data values.
#' @param y_sd numeric vector with standard deviations.
#' @param alternative a character string specifying the alternative
#'  hypothesis, must be one of \code{"two.sided"} (default),
#'  \code{"greater"} or \code{"less"}. You can specify just the initial letter.
#'
#' @references
#' Ermolinskiy, P. (2024). An Extension of the Permutation Test for
#' Dependent Data in Samples with a Small Number of Subjects.
#' https://hal.science/hal-04558513/
#'
#' @details
#' This function is based in the Python function proposed
#' by Ermolinskiy (2024).
#'
#' @return A list with class \code{"htest"} containing the following
#'  components:
#' \item{statistic}{the value of the statistic.}
#' \item{p.value}{the p-value for the test.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{method}{a character string indicating the type of test performed.}
#'
#' @examples
#' # Example 0
#' x <- c(5, 6)
#' x_sd <- c(1, 2)
#' y <- c(7, 9)
#' y_sd <- c(1, 2)
#'
#' res <- permu_test_SD(x, x_sd, y, y_sd, alternative="less")
#' res
#'
#' # Example 1 of Ermolinskiy (2024).
#' x <- c(20, 33, 20, 24)
#' x_sd <- c(2, 2, 2, 2)
#' y <- c(21, 34, 22, 25)
#' y_sd <- c(2, 2, 2, 2)
#'
#' res <- permu_test_SD(x, x_sd, y, y_sd, alternative="less")
#' res
#'
#' @export
permu_test_SD <- function(x, x_sd, y, y_sd, alternative = "two.sided") {

  # To ensure that the inputs are vector with the same length
  error_handling(x, x_sd, y, y_sd)

  alternative <- match.arg(arg=alternative,
                           choices=c("two.sided", "greater", "less"))

  # This is the average difference in percent between
  # the second and the first group
  statistic <- mean(y / x * 100 - 100)
  names(statistic) <- "Percentage_difference"

  # All possible combinations of values using standard deviation
  x_combination <- combinations_mean_and_sd(x, x_sd)
  y_combination <- combinations_mean_and_sd(y, y_sd)

  # generate all possible combinations
  n <- nrow(x_combination)
  rep_indices_x <- rep(1:n, each=n)
  rep_indices_y <- rep(1:n, times=n)
  comb1 <- x_combination[rep_indices_x, ]
  comb2 <- y_combination[rep_indices_y, ]

  # creating the combinations of -1 and 1 to multiply percentage differences
  m <- ncol(comb1)
  multi <- expand.grid(replicate(m, c(-1, 1), simplify = FALSE))

  # Auxiliar function to obtain all combinatios of the percentage differences
  # varying the signs
  aux_fun <- function(x) rowMeans(t(x * t(multi)))

  # To obtain the percentage differences
  res <- (comb2 - comb1) / comb1
  res2 <- apply(res, 1, aux_fun)
  percent_diff <- 100 * res2

  # calculate the p-value
  if (alternative == "two.sided") {
    p.value <- mean(abs(percent_diff) < abs(statistic))
  } else if (alternative == "greater") {
    p.value <- mean(percent_diff < statistic)
  } else if (alternative == "less") {
    p.value <- mean(percent_diff > statistic)
  }

  method <- "Paired permutation with SD"

  res <- list(statistic = statistic,
              p.value=p.value,
              alternative=alternative,
              method=method)
  class(res) <- "htest"
  res
}
combinations_mean_and_sd <- function(x, x_sd) {
  # Generate all possible combinations
  sign_options <- c(-1, 1, 0)

  sign_combinations <- expand.grid(replicate(length(x),
                                             sign_options,
                                             simplify=FALSE))

  combinations <- t(t(sign_combinations) * x_sd + x)
  return(combinations)
}
error_handling <- function(x, x_sd, y, y_sd) {
  if (!is.vector(x) || !is.vector(x_sd) || !is.vector(y) || !is.vector(y_sd)) {
    stop("The input must be lists")
  }
  if (length(x) != length(y) || length(x_sd) != length(y_sd) || length(x) != length(x_sd)) {
    stop("All lists must have the same length")
  }
}
