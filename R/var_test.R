#' Variance test using values
#'
#' This function performs the test for a single variance or two variances using values, not the vectors.
#'
#' @param varx sample variance for sample x.
#' @param nx sample size for sample x.
#' @param vary sample variance for sample y.
#' @param ny sample size for sample y.
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
#'
#' @examples
#' # Examples with ONE sample
#'
#' # Example 7.7.1 from Wayne (2013), http://tinyurl.com/y6z49hrw
#' var_test(varx=670.81, nx=16, null.value=600, alternative='two.sided')
#'
#' # Exercise 7.7.5 from Wayne (2013), http://tinyurl.com/y6z49hrw
#' var_test(varx=30, nx=25, null.value=25, alternative='greater')
#'
#' # Using the plot to illustrate Hypothesis Test
#' mytest1 <- var_test(varx=30, nx=25, null.value=25, alternative='greater')
#' mytest1
#' plot(mytest1)
#'
#' # Examples with TWO samples
#'
#' # Example 7.8 from Montgomery (1996)
#' var_test(varx=5.1^2, nx=12, vary=4.7^2, ny=15, conf.level=0.90)
#'
#' # Example 8.17 from Montgomery (1996)
#' mytest2 <- var_test(varx=3.84, nx=20, vary=4.54, ny=20)
#' mytest2
#' plot(mytest2)
#'
#' @export
var_test <- function(varx, nx, vary=NULL, ny=NULL,
                     alternative='two.sided',
                     null.value=1, conf.level=0.95) {

  # Checking if the information is correct

  # To check if the information about x sample is correct
  if (varx <= 0)
    stop(paste("The variance x must be positive", "\n", ""))
  if (nx <= 0)
    stop(paste("The sample size nx must be positive", "\n", ""))
  if (nx %% 1 != 0)
    stop(paste("The sample size nx must be integer", "\n", ""))

  # To check if the user provided information about y sample
  if (xor(is.null(vary), is.null(ny)))
    stop("Some information about sample y is missing", "\n", "")

  # To check if the information about x sample is correct
  if (! is.null(vary) & ! is.null(ny)) {
    if (vary <= 0)
      stop(paste("The variance y must be positive", "\n", ""))
    if (ny <= 0)
      stop(paste("The sample size ny must be positive", "\n", ""))
    if (ny %% 1 != 0)
      stop(paste("The sample size ny must be integer", "\n", ""))
  }

  # To check if the conf.level is a number between 1 and 0
  if(conf.level <= 0 || conf.level >= 1)
    stop("The conf.level argument must be > 0 and < 1", "\n", "" )
  if(conf.level < 0.5)
    warning("Confidence levels are often close to 1")

  # To check if the null.value is positive
  if(null.value <= 0)
    stop(paste("The null value must be positive", "\n", ""))

  # Argument Verification Using Partial Matching
  alternative <- match.arg(arg=alternative,
                           choices=c("two.sided","greater","less"))
  # The next variable is used to indicate if we have raw o summarized data
  raw.data <- FALSE

  if (is.null(vary))
    res <- var_test_one(varx, nx, alternative,
                        conf.level, null.value, raw.data)
  else
    res <- var_test_two(varx, nx, vary, ny,
                        alternative, conf.level, null.value, raw.data)

  class(res) <- "htest"
  res
}
#' @importFrom stats pchisq qchisq
var_test_one <- function(varx, nx, alternative, conf.level,
                         null.value, raw.data, name_x=NULL) {
  alpha <- 1 - conf.level
  # Alternative two.sided
  if (alternative == 'two.sided') {
    quantiles <- c(qchisq(p=alpha/2, df=nx-1, lower.tail=F),
                   qchisq(p=1-alpha/2, df=nx-1, lower.tail=F))
    conf.int <- (nx-1) * varx / quantiles
    statistic <- (nx-1) * varx / null.value
    p.value <- 2 * min(c(pchisq(statistic, nx-1, lower.tail=F),
                         pchisq(statistic, nx-1, lower.tail=T)))
  }
  # Alternative less
  if (alternative == 'less') {
    quantiles <- c(qchisq(p=conf.level, df=nx-1, lower.tail=T),
                   0)
    conf.int <- (nx-1) * varx / quantiles
    statistic <- (nx-1) * varx / null.value
    p.value <- pchisq(statistic, nx-1)
  }
  # Alternative greater
  if (alternative == 'greater') {
    quantiles <- c(Inf,
                   qchisq(p=conf.level, df=nx-1, lower.tail=F))
    conf.int <- (nx-1) * varx / quantiles
    statistic <- (nx-1) * varx / null.value
    p.value <- pchisq(statistic, nx-1, lower.tail=F)
  }

  # To ensure that the output values are in the correct form
  names(statistic) <- 'X-squared'
  parameter <- nx - 1
  names(parameter) <- 'df'
  attr(conf.int, 'conf.level') <- conf.level
  estimate <- varx
  names(estimate) <- 'variance of x'
  null.value <- null.value
  names(null.value) <- 'variance'
  method <- 'X-squared test for variance'
  if (raw.data) data.name <- name_x
  else data.name <- paste('varx =', varx, 'and nx =', nx)

  res <- list(statistic=statistic,
              parameter=parameter,
              p.value=p.value,
              conf.int=conf.int,
              estimate=estimate,
              null.value=null.value,
              alternative=alternative,
              method=method,
              data.name=data.name)
  return(res)
}
#' @importFrom stats pf qf
var_test_two <- function(varx, nx, vary, ny, alternative, conf.level,
                         null.value, raw.data, name_x=NULL, name_y=NULL) {
  alpha <- 1 - conf.level
  # Alternative two.sided
  if (alternative == 'two.sided') {
    quantiles <- c(qf(p=alpha/2,   df1=nx-1, df2=ny-1, lower.tail=F),
                   qf(p=1-alpha/2, df1=nx-1, df2=ny-1, lower.tail=F))
    conf.int <- (varx / vary) / quantiles
    statistic <- (varx / vary) / null.value
    p.value <- 2 * min(c(pf(statistic, nx-1, ny-1, lower.tail=F),
                         pf(statistic, nx-1, ny-1, lower.tail=T)))
  }
  # Alternative less
  if (alternative == 'less') {
    quantiles <- c(Inf,
                   qf(p=conf.level, df1=nx-1, df2=ny-1, lower.tail=F))
    conf.int <- (varx / vary) / quantiles
    statistic <- (varx / vary) / null.value
    p.value <- pf(q=statistic, df1=nx-1, df2=ny-1, lower.tail=T)
  }
  # Alternative greater
  if (alternative == 'greater') {
    quantiles <- c(qf(p=conf.level, df1=nx-1, df2=ny-1, lower.tail=T),
                   0)
    conf.int <- (varx / vary) / quantiles
    statistic <- (varx / vary) / null.value
    p.value <- pf(q=statistic, df1=nx-1, df2=ny-1, lower.tail=F)
  }

  # To ensure that the output values are in the correct form
  names(statistic) <- 'F'
  parameter <- c(nx-1, ny-1)
  names(parameter) <- c('num df', 'denom df')
  attr(conf.int, 'conf.level') <- conf.level
  estimate <- varx / vary
  names(estimate) <- 'ratio of variances'
  null.value <- null.value
  names(null.value) <- 'ratio of variances'
  method <- 'F test to compare two variances'
  if (raw.data) data.name <- paste(name_x, ' and ', name_y)
  else data.name <- data.name <- paste('varx =', varx, ', nx =', nx,
                                       ', vary =', vary, 'and ny =', ny)

  res <- list(statistic=statistic,
              parameter=parameter,
              p.value=p.value,
              conf.int=conf.int,
              estimate=estimate,
              null.value=null.value,
              alternative=alternative,
              method=method,
              data.name=data.name)
  return(res)
}
