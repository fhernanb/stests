#' Mean test using values
#' 
#' This function performs the mean test using summarized data (values).
#' 
#' @param meanx sample mean for sample x.
#' @param nx sample size for sample x.
#' @param sigma2 population variance which is known.
#' @param alternative a character string specifying the alternative 
#'  hypothesis, must be one of \code{two.sided} (default), 
#'  \code{greater} or \code{less}. You can specify just the initial letter.
#' @param mu the hypothesized number in the null hypothesis. 
#' @param conf.level confidence level of the interval, by default its value is 0.95.
#'  
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
#' # Example 13.1 from Freund et al. (2000)
#' res1 <- z_test(meanx=8.091, nx=25, mu=8, sigma2=0.16^2, alternative='two.sided')
#' res1
#' plot(res1)
#' 
#' # Example 13.2 from Freund et al. (2000)
#' res2 <- z_test(meanx=21819, nx=100, mu=22000, sigma2=1295^2, alternative='less')
#' res2
#' plot(res2)
#' 
#' @importFrom stats pnorm qnorm
#' @export
z_test <- function(meanx, nx, sigma2=NULL,
                   alternative='two.sided',
                   mu=0, conf.level=0.95) {
  
  # Checking if the information is correct
  if (nx <= 0) 
    stop(paste("The sample size nx must be positive", "\n", ""))
  if (nx %% 1 != 0) 
    stop(paste("The sample size nx must be integer", "\n", ""))

  # Checking the sigma2 value
  if (is.null(sigma2) | !is.numeric(sigma2))
    stop(paste("Check the sigma2 value", "\n", ""))
  if (sigma2 <= 0)
    stop("The variance sigma2 must be positive")
  
  # Argument Verification Using Partial Matching
  alternative <- match.arg(arg=alternative,
                           choices=c("two.sided","greater","less"))
  
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
  names(estimate) <- 'mean of x'
  method <- 'Z test for mean'
  data.name <- deparse(substitute(x))
  
  res <- list(statistic=statistic,
              p.value=p.value,
              conf.int=conf.int,
              estimate=estimate,
              null.value=mu,
              alternative=alternative,
              method=method,
              data.name=data.name)

  class(res) <- "htest"
  res
}


