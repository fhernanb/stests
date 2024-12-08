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
#' @references
#' \insertRef{Wald1943}{stests}
#'
#' @seealso \link{ci_p}.
#'
#' @details
#' The expression to obtain the confidence interval is given below:
#'
#' \eqn{\hat{p} - z_{\alpha/2} \sqrt{\frac{\hat{p}(1 - \hat{p})}{n}} \leq p \leq \hat{p} + z_{\alpha/2} \sqrt{\frac{\hat{p}(1 - \hat{p})}{n}}},
#'
#' where \eqn{\hat{p} = \frac{x}{n}} is the sample proportion, \eqn{x} the number of observed successes in the sample with size \eqn{n}. The value \eqn{z_{\alpha/2}} is the \eqn{1-\alpha/2} percentile of the standard normal distribution (e.g., \eqn{z_{0.025} = 1.96} for a 95\% confidence interval).
#'
#' @return A matrix with the lower and upper limits.
#'
#' @example examples/examples_ci_p_xxx.R
#' @export
#'
ci_p_wald <- function(x, n, conf.level=0.95) {
  q <- qnorm(p = (1 + conf.level)/2)
  p <- x/n
  lower <- max(p - q*sqrt(p*(1 - p)/n), 0)
  upper <- min(p + q*sqrt(p*(1 - p)/n), 1)
  return(c(lower, upper))
}
ci_p_wald <- Vectorize(ci_p_wald)
