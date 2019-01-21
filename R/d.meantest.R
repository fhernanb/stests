#' D-value for hypothesis test using vectors
#' 
#' This function shows the D-value of a hypothesis test for two samples using non-summarized data, that is, vectors.
#' 
#' @param x vector sample data x.
#' @param y vector sample data y.
#' @param alternative a character string specifying the alternative hypothesis, must be one of \code{"less"} (default), \code{"greater"} or \code{"two.sided"}. You can specify just the initial letter.
#' @param nrep number of repetions to obtain the D-value, by default is 1000000
#'
#' @return A list with class \code{"htest"} containing the following 
#'  components:
#' \item{d.value}{the d-value for the test.}
#' \item{estimate}{the estimated mean or difference in means depending
#'       on whether it was a one-sample test or a two-sample test.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{method}{a character string indicating the type of test performed.}
#' 
#' @details In this case, it is important to explain that the sample size is assumed as n = 1 
#' because the D-value provides an individual analysis.
#' 
#' @examples
#' d.meantest(x=c(2,3,4,5,6), y=c(3,4,5,6,7), alternative="two.sided")
#' 
#' @importFrom stats rnorm
#' @export
d.meantest <- function(x, y, alternative='less', nrep=1000000) {
  meanx <- mean(x)
  meany <- mean(y)
  varx <- var(x)
  vary <- var(y)
  x <- rnorm(n=nrep, mean=meanx, sd=sqrt(varx))
  y <- rnorm(n=nrep, mean=meany, sd=sqrt(vary))
  if (alternative == 'less') 
    dvalue <-  mean(y > x)
  if (alternative == 'greater')
    dvalue <-  mean(y < x)
  if (alternative == 'two.sided') {
    dvalue1 <-  mean(y > x)
    dvalue2 <-  mean(y < x)
    dvalue <- 2 * min(c(dvalue1, dvalue2))
  }
  estimate <- c(meanx, meany)
  names(estimate) <- c('mean of x', 'mean of y')
  method <- 'Two Sample d.test'
  data.name <- paste('The D-value was calculated using','nrep=', nrep)
  statistic <- dvalue 
  names(statistic) <- 'd.value'
  
  res <- list(d.value=dvalue,
              statistic= statistic,
              estimate = estimate,
              alternative = alternative,
              method = method,
              data.name = data.name)
  
  class(res) <- 'htest'
  
  return(res)
}