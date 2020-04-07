#' Mean test using vector
#'
#' This function performs the mean test using raw data (vector).
#'
#' @param x a (non-empty) numeric vector of data values.
#' @param sigma2 population variance which is known.
#' @param alternative a character string specifying the alternative
#'  hypothesis, must be one of \code{two.sided} (default),
#'  \code{greater} or \code{less}. You can specify just the initial letter.
#' @param mu the hypothesized number in the null hypothesis.
#' @param conf.level confidence level of the interval, by default its value is 0.95.
#' @return A list with class \code{htest} containing the following
#'  components:
#' \item{statistic}{the value of the statistic.}
#' \item{p.value}{the p-value for the test.}
#' \item{conf.int}{a confidence interval for the variance.}
#' \item{estimate}{the estimated mean.}
#' \item{null.value}{the specified hypothesized value for alternative hypothesis.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{method}{a character string indicating the type of test performed.}
#' \item{data.name}{a character string giving the name of the data.}
#'
#' @examples
#' # Example 1
#' # H0: mu = 350
#' # Ha: mu < 350 with sigma2=25
#' content <- c(355, 353, 352, 346, 345, 345, 353, 353, 344, 350)
#' res1 <- z.test(x=content, mu=350, sigma2=25, alternative='less')
#' res1
#' plot(res1, shade.col='deepskyblue', col='deepskyblue')
#'
#' # Example 2
#' # H0: mu = 170
#' # Ha: mu != 170 with sigma2=25
#' x <- rnorm(n=80, mean=171, sd=5)
#' res2 <- z.test(x=x, mu=170, sigma2=25, alternative='two.sided')
#' res2
#' plot(res2, shade.col='indianred1', col='indianred1')
#'
#' @importFrom stats pnorm qnorm
#' @export
z.test <- function(x, sigma2,
                   alternative='two.sided',
                   mu=0, conf.level=0.95) {

  # Checking if the information is correct
  if (! is.numeric(x))
    stop(paste("The x vector must be numeric", "\n", ""))
  if (length(x) <= 1)
    stop(paste("not enough 'x' observations", "\n", ""))

  # Checking the sigma2 value
  if (is.null(sigma2) | !is.numeric(sigma2))
    stop(paste("Check the sigma2 value", "\n", ""))
  if (sigma2 <= 0)
    stop("The variance sigma2 must be positive")

  # Argument Verification Using Partial Matching
  alternative <- match.arg(arg=alternative,
                           choices=c("two.sided","greater","less"))

  meanx <- mean(x)
  nx <- length(x)
  alpha <- 1 - conf.level

  # Alternative two.sided
  if (alternative == 'two.sided') {
    statistic <- (meanx - mu) / sqrt(sigma2 / nx)
    p.value <- 2 * pnorm(q=abs(statistic), lower.tail=FALSE)
    quantiles <- c(-qnorm(p=alpha/2, lower.tail=FALSE),
                   qnorm(p=alpha/2, lower.tail=FALSE))
    conf.int <- meanx + quantiles * sqrt(sigma2 / nx)
  }
  # Alternative less
  if (alternative == 'less') {
    statistic <- (meanx - mu) / sqrt(sigma2 / nx)
    p.value <- pnorm(q=statistic, lower.tail=TRUE)
    conf.int <- c(-Inf,
                  meanx + qnorm(p=1-alpha) * sqrt(sigma2 / nx))
  }
  # Alternative greater
  if (alternative == 'greater') {
    statistic <- (meanx - mu) / sqrt(sigma2 / nx)
    p.value <- pnorm(q=statistic, lower.tail=FALSE)
    conf.int <- c(meanx + qnorm(p=alpha) * sqrt(sigma2 / nx),
                  Inf)
  }

  # To ensure that the output values are in the correct form
  names(statistic) <- 'Z'
  attr(conf.int, 'conf.level') <- conf.level
  estimate <- meanx
  names(estimate) <- paste("mean of", deparse(substitute(x)))
  null.value <- mu
  names(null.value) <- 'mean'
  method <- 'One Sample z-test'
  data.name <- deparse(substitute(x))

  res <- list(statistic=statistic,
              p.value=p.value,
              conf.int=conf.int,
              estimate=estimate,
              null.value=null.value,
              alternative=alternative,
              method=method,
              data.name=data.name)

  class(res) <- "htest"
  res
}
