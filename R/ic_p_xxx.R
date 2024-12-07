#' Wald confidence interval for Binomial proportion
#'
#' @author Olga Bustos, \email{oabustos@unal.edu.co}
#'
#' @description
#' This function obtains the confidence interval for a proportion. It is vectorized, so the user can evaluate it using single values or a vector.
#'
#' @param x a number or a vector with the number of successes.
#' @param n a number or a vector with the number of trials.
#' @param conf.level confidence level for the returned confidence interval. By default is 0.95.
#'
#' @seealso \link{ic_p}.
#'
#' @return A matrix with the lower and upper limits.
#'
#' @example examples/examples_ic_p_xxx.R
#' @export
#'
ic_p_wald <- function(x, n, conf.level=0.95) {
  q <- qnorm(p = (1 + conf.level)/2)
  p <- x/n
  lower <- max(p - q*sqrt(p*(1 - p)/n), 0)
  upper <- min(p + q*sqrt(p*(1 - p)/n), 1)
  return(c(lower, upper))
}
ic_p_wald <- Vectorize(ic_p_wald)
