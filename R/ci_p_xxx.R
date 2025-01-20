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
#' @seealso \link{ci_p}.
#'
#' @details
#' The expression to obtain the confidence interval is given below:
#'
#' \eqn{\hat{p} - z_{\alpha/2} \sqrt{\frac{\hat{p}(1 - \hat{p})}{n}} \leq p \leq \hat{p} + z_{\alpha/2} \sqrt{\frac{\hat{p}(1 - \hat{p})}{n}}},
#'
#' where \eqn{\hat{p} = \frac{x}{n}} is the sample proportion, \eqn{x} the
#' number of observed successes in the sample with size \eqn{n}. The
#' value \eqn{z_{\alpha/2}} is the \eqn{1-\alpha/2} percentile of the
#' standard normal distribution (e.g., \eqn{z_{0.025} = 1.96} for a 95\%
#' confidence interval).
#'
#' @return A matrix with the lower and upper limits.
#'
#' @examples
#' ci_p_wald(x=15, n=50, conf.level=0.95)
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
#'
#'
#'
#' Agresti-Coull confidence interval for Binomial proportion
#'
#' @author Omar David Mercado Turizo, \email{omercado@unal.edu.co}
#'
#' @description
#' This function calculates the Agresti-Coull confidence interval for a Binomial proportion. It is vectorized, allowing the evaluation of single values or vectors.
#'
#' @param x A number or a vector with the number of successes.
#' @param n A number or a vector with the number of trials.
#' @param conf.level Confidence level for the returned confidence interval. By default, it is 0.95.
#'
#' @seealso \link{ci_p}.
#'
#' @details
#' The Agresti-Coull interval is an approximate confidence interval for the Binomial proportion \eqn{p}.
#' The limits are calculated based on an adjusted proportion \eqn{\tilde{p}} and its standard error. The mathematical definitions are as follows:
#'  Adjusted proportion: \eqn{\tilde{p} = \frac{x + 2}{n + 4}};
#'  Adjusted standard error: \eqn{se = \sqrt{\frac{\tilde{p}(1 - \tilde{p})}{n + 4}}};
#'  Confidence limits: \eqn{\tilde{p} \pm z_{\alpha/2} \cdot se},
#'
#' where \eqn{z_{\alpha/2}} is the critical value of the standard normal distribution. The limits are truncated to the range \eqn{[0, 1]}.
#'
#' @return A vector with the lower and upper limits.
#'
#' @examples
#' ci_p_agresti_coull(x=15, n=50, conf.level=0.95)
#'
#' @export
#'
ci_p_agresti_coull <- function(x, n, conf.level = 0.95) {
  alpha <- 1 - conf.level
  z_alpha <- qnorm(1 - alpha / 2) # Critic value
  # Estimated proportion
  pi_tilde <- (x + 2) / (n + 4)
  # Standard error
  se <- sqrt((pi_tilde * (1 - pi_tilde)) / (n + 4))
  lower <- max(0, pi_tilde - z_alpha * se)
  upper <- min(1, pi_tilde + z_alpha * se)
  return(c(lower, upper))
}

ci_p_agresti_coull <- Vectorize(ci_p_agresti_coull)
