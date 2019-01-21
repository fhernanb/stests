#' D-value for hypothesis test using values
#' 
#' This function shows the D-value of a hypothesis test for two samples using summarized values, not the vectors.
#' 
#' @param meanx sample mean for sample x.
#' @param varx sample variance for sample x.
#' @param meany sample mean for sample y.
#' @param vary sample variance for sample y. 
#' @param alternative a character string specifying the alternative hypothesis, must be one of \code{"less"} (default), \code{"greater"} or \code{"two.sided"}. You can specify just the initial letter.
##' @param nrep number of repetions to obtain the D-value, by default is 1000000
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
#' # --- Examples with TWO-SAMPLES and equal variances ---
#' d_meantest(meanx=250, varx=20^2, meany=249, vary=20^2, alternative='greater')
#'  
#' @importFrom stats rnorm
#' @export
d_meantest <- function(meanx, varx, meany, vary,
                    alternative='less', nrep=1000000) {
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
  data.name <- paste('The D-value was calculated using', 'nrep=', nrep)
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
