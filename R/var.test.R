#' Variance test using vectors
#'
#' This function performs the test for a single variance or two variances given the vectors. This function is a generalization of \code{var.test} function from \code{stats} package.
#'
#' @param x a (non-empty) numeric vector of data values.
#' @param y an optional (non-empty) numeric vector of data values.
#' @param alternative a character string specifying the alternative
#'  hypothesis, must be one of \code{two.sided} (default),
#'  \code{greater} or \code{less}. You can specify just the initial letter.
#' @param null.value the hypothesized number (variance or ratio of the variances) in the null hypothesis.
#' @param conf.level confidence level of the interval, by default its value is 0.95.
#'
#' @return A list with class \code{htest} containing the following
#'  components:
#' \item{statistic}{the value of the statistic.}
#' \item{p.value}{the p-value for the test.}
#' \item{conf.int}{a confidence interval for the variance.}
#' \item{estimate}{the sample variance (or ratio of the sample variances)}
#' \item{null.value}{the specified hypothesized value for alternative hypothesis.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{method}{a character string indicating the type of test performed.}
#' \item{data.name}{a character string giving the name of the data.}
#'
#' @examples
#' # One sample -----
#'
#' # Interval confidence
#' duration <- c(1470, 1510, 1690, 1740, 1900, 2000, 2030,
#'               2010, 2190, 2200, 2290, 2380, 2390, 2480,
#'               2500, 2580, 2700)
#' var.test(x=duration, conf.level=0.95)
#'
#' # Hypothesis testing
#' #   H0: sigma2 = 100
#' #   H1: sigma2 > 100
#' weight <- c(775, 780, 781, 795, 803, 810, 823)
#' res1 <- var.test(x=weight, alternative='greater', null.value=100)
#' res1
#' # Using the plot function
#' plot(res1)
#'
#' # Two samples -----
#'
#' # Hypothesis testing
#' #   H0: sigma1/sigma2 = 1
#' #   H1: sigma1/sigma2 != 1
#' x1 <- rnorm(50, mean = 0, sd = 2)
#' x2 <- rnorm(30, mean = 1, sd = 1)
#' res2 <- var.test(x1, x2)
#' res2
#' plot(res2, from=0, to=10)
#'
#' @importFrom stats var
#' @export
var.test <- function(x, y=NULL,
                     alternative='two.sided',
                     null.value=1, conf.level=0.95) {

  # Checking if the information is correct
  if (! is.numeric(x))
    stop(paste("The x vector must be numeric", "\n", ""))
  if (length(x) <= 1)
    stop(paste("not enough 'x' observations", "\n", ""))

  if (! is.null(y)) {
    if (! is.numeric(y))
      stop(paste("The y vector must be numeric", "\n", ""))
    if (length(y) <= 1)
      stop(paste("not enough 'y' observations", "\n", ""))
  }

  # To check if the null.value is positive
  if(null.value <= 0)
    stop(paste("The null value must be positive", "\n", ""))

  # Argument Verification Using Partial Matching
  alternative <- match.arg(arg=alternative,
                           choices=c("two.sided","greater","less"))

  # The next variable is used to indicate if we have raw o summarized data
  raw.data <- TRUE
  name_x <- deparse(substitute(x))
  name_y <- deparse(substitute(y))

  if (is.null(y))
    res <- var_test_one(var(x), length(x), alternative,
                        conf.level, null.value,
                        raw.data, name_x)
  else
    res <- var_test_two(var(x), length(x), var(y), length(y),
                        alternative, conf.level, null.value,
                        raw.data, name_x, name_y)

  class(res) <- "htest"
  res
}
