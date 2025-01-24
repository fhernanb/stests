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
#' @references
#' Wald, A. (1949). Statistical decision functions. The Annals of Mathematical Statistics, 165-205.
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
#'
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
#' @references
#' Agresti, A., & Coull, B. A. (1998). Approximate is better than “exact” for interval estimation of binomial proportions. The American Statistician, 52(2), 119-126.
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
ci_p_agresti_coull <- function(x, n, conf.level=0.95) {
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
#'
#'
#'
#' Rindskopf  confidence interval for Binomial proportion
#'
#' @author David Esteban Cartagena Mejía, \email{dcartagena@unal.edu.co}
#'
#' @description
#' This function calculates the Rindskopf confidence interval for a Binomial proportion. It is vectorized, allowing the evaluation of single values or vectors.
#'
#' @param x A number or a vector with the number of successes.
#' @param n A number or a vector with the number of trials.
#' @param conf.level Confidence level for the returned confidence interval. By default, it is 0.95.
#'
#' @seealso \link{ci_p}.
#'
#' @references
#' Rindskopf, D. (2000). Commentary: Approximate is better than “exact” for interval
#' estimation of binomial proportions. The American Statistician, 54, 88.
#'
#' @details
#' The expression to calculate the confidence interval according to the Rindskopf approach is given by:
#'
#' \eqn{\phi = \text{logit}(\pi) = \log\left(\frac{\pi}{1 - \pi}\right)},
#'
#' where the maximum likelihood estimator for \eqn{\phi} is:
#'
#' \eqn{\hat{\phi}_{ML} = \log\left(\frac{x + 0.5}{n - x + 0.5}\right)},
#'
#' and its standard error is:
#'
#' \eqn{\text{se}(\hat{\phi}_{ML}) = \sqrt{\frac{1}{x + 0.5} + \frac{1}{n - x + 0.5}}}.
#'
#' The adjustment of adding 0.5 successes and non-successes ensures that intervals can also be computed for the cases where \eqn{x = 0} or \eqn{x = n} (where otherwise the maximum likelihood estimator and standard error would be infinite).
#'
#' Since the scale of \eqn{\phi} is \eqn{(- \infty, \infty)}, this interval respects the boundary constraints. Back-transformation to the scale of \eqn{\pi} is performed using the inverse logit function:
#'
#' \eqn{\pi = \text{expit}(\phi) = \frac{\exp(\phi)}{1 + \exp(\phi)}}.
#'
#' Thus, the confidence interval for \eqn{\pi} in the original scale is the Rindskopf confidence interval, as proposed by Rindskopf.
#'
#' @return A vector with the lower and upper limits.
#'
#' @examples
#' ci_p_rindskopf(x=15, n=50, conf.level=0.95)
#'
#' @export
#'
ci_p_rindskopf <- function(x, n, conf.level=0.95) {
  num <- x + 0.5
  denom <- n - x + 0.5
  phi_hat <- log(num / denom)
  se_phi <- sqrt(1 / num + 1 / denom)
  z <- qnorm(1 - (1 - conf.level) / 2)
  lp <- phi_hat - z * se_phi
  up <- phi_hat + z * se_phi
  alpha <- 1-conf.level

  # Límites del intervalo
  if (x == 0) {
    upper <- (alpha / 2)^(1 / n)
    lower <- 0
  } else if (x == n) {
    lower <- 1 - (alpha / 2)^(1 / n)
    upper <- 1
  } else {
    # Límites del intervalo
    if (x == 0) {
      upper <- (alpha / 2)^(1 / n)
      lower <- 0
    } else if (x == n) {
      lower <- 1 - (alpha / 2)^(1 / n)
      upper <- 1
    } else {
      lower<- exp(lp) / (1 + exp(lp))
      upper <- exp(up) / (1 + exp(up))
    }
  }
  return(c(lower,upper))
}
ci_p_rindskopf <- Vectorize(ci_p_rindskopf)
#'
#'
#'
#' Clopper-Pearson confidence interval for Binomial proportion
#'
#' @author Omar David Mercado Turizo, \email{omercado@unal.edu.co}
#'
#' @description
#' This function calculates the exact Clopper-Pearson confidence interval for a Binomial proportion. It is vectorized, allowing the evaluation of single values or vectors.
#'
#' @param x A number or a vector with the number of successes.
#' @param n A number or a vector with the number of trials.
#' @param conf.level Confidence level for the returned confidence interval. By default, it is 0.95.
#'
#' @seealso \link{ci_p}.
#'
#' @references
#' Clopper, C. J., & Pearson, E. S. (1934). The use of confidence or fiducial limits illustrated in the case of the binomial. Biometrika, 26(4), 404-413.
#'
#' @details
#' The Clopper-Pearson interval is an exact confidence interval for the Binomial proportion \eqn{p}.
#' The limits of the interval are derived based on the Beta distribution. For the special cases where \eqn{x = 0} or \eqn{x = n}, the limits are calculated directly.
#'
#' The mathematical definitions are as follows. If \eqn{x = 0}, the lower limit is \eqn{0}, and the upper limit is \eqn{1 - (\alpha / 2)^{1/n}}.
#' If \eqn{x = n}, the lower limit is \eqn{(\alpha / 2)^{1/n}}, and the upper limit is \eqn{1}.
#' Otherwise, the limits are given by the \eqn{\alpha / 2} and \eqn{1 - \alpha / 2} quantiles of the Beta distribution with parameters adjusted by the data,
#' that is, \eqn{B_{x+1/2,n-x+1/2,1-\alpha/2} \leq p \leq B_{x+1/2,n-x+1/2,\alpha/2}}.
#'
#' @return A vector with the lower and upper limits.
#'
#' @examples
#' ci_p_clopper_pearson(x=15, n=50, conf.level=0.95)
#'
#' @export
#'
ci_p_clopper_pearson <- function(x, n, conf.level=0.95) {
  alpha <- 1 - conf.level

  if (x == 0) {
    lower <- 0
    upper <- 1 - (alpha / 2)^(1 / n)
  } else if (x == n) {
    lower <- (alpha / 2)^(1 / n)
    upper <- 1
  } else {
    lower <- qbeta(alpha / 2, x, n - x + 1)
    upper <- qbeta(1 - alpha / 2, x + 1, n - x)
  }

  return(c(lower, upper))
}
ci_p_clopper_pearson <- Vectorize(ci_p_clopper_pearson)
#'
#'
#'
#' Add-4 Wald-t Confidence Interval for Binomial Proportion
#'
#' @author David Esteban Cartagena Mejía, \email{dcartagena@unal.edu.co}
#'
#' @description
#' This function calculates the Add-4 Wald-t confidence interval for a Binomial proportion using the method XX proposed by Pan (2002).
#' It modifies the classical Wald interval by adding four pseudo-observations and uses a t-distribution with adjusted degrees of freedom \eqn{\nu}, as specified in equation (2.9).
#' The function is vectorized, allowing for evaluation of single values or vectors.
#'
#' @param x A number or a vector with the number of successes.
#' @param n A number or a vector with the number of trials.
#' @param conf.level Confidence level for the returned confidence interval. By default, it is 0.95.
#'
#' @seealso \link{ci_p}.
#'
#' @references
#' Falta la referencia.
#'
#' @details
#' The Add-4 Wald-t confidence interval improves the performance of the Wald interval by adding 2 successes and 2 failures to the observed data, effectively modifying the estimated proportion:
#'
#' \deqn{\tilde{\pi} = \frac{x + 2}{n + 4}.}
#'
#' The variance \eqn{V(\tilde{\pi}, n+4)} is given by:
#'
#' \deqn{V(\tilde{\pi}, n+4) = \frac{\tilde{\pi}(1 - \tilde{\pi})}{n + 4}.}
#'
#' The degrees of freedom \eqn{\nu} are calculated using equation (2.9):
#'
#' \deqn{\nu = \frac{2 V(\tilde{\pi(2)}, n+4)^2}{\Omega(\tilde{\pi(2)}, n+4)},}
#'
#' where \eqn{\Omega(\tilde{\pi(2)}, n+4)} is defined as:
#'
#' \deqn{\Omega(p, n) = \frac{p - p^2}{n^3} + \frac{p + (6n - 7)p^2 + 4(n - 1)(n - 3)p^3 - 2(n - 1)(2n - 3)p^4}{n^5} - \frac{2(p + (2n - 3)p^2 - 2(n - 1)p^3)}{n^4}.}
#'
#' The confidence interval is then calculated as:
#'
#' \deqn{\text{Lower} = \tilde{\pi} - t \sqrt{\frac{\tilde{\pi}(1 - \tilde{\pi})}{n+4}},}
#' \deqn{\text{Upper} = \tilde{\pi} + t \sqrt{\frac{\tilde{\pi}(1 - \tilde{\pi})}{n+4}},}
#'
#' where \eqn{t} is the critical value from the t-distribution with \eqn{\nu} degrees of freedom.
#'
#' @return A vector with the lower and upper limits of the confidence interval.
#'
#' @examples
#' ci_p_add_4(x=15, n=50, conf.level=0.95)
#' ci_p_add_4(x=0,  n=50, conf.level=0.95)
#' ci_p_add_4(x=50, n=50, conf.level=0.95)
#'
#' @export
#'
ci_p_add_4 <- function(x, n, conf.level=0.95) {
  # Proporcion ajustada
  pi_tilde <- (x + 2) / (n + 4)

  # Funcion para calcular Omega
  Omega <- function(p, n) {
    term1 <- (p - p^2) / n^3
    term2 <- (p + (6 * n - 7) * p^2 + 4 * (n - 1) * (n - 3) * p^3 - 2 * (n - 1) * (2 * n - 3) * p^4) / n^5
    term3 <- (2 * (p + (2 * n - 3) * p^2 - 2 * (n - 1) * p^3)) / n^4
    return(term1 + term2 - term3)
  }

  # Varianza de p ajustada
  V <- function(p, n) {
    return(p * (1 - p) / n)
  }

  # Grados de libertad (nu)
  nu <- (2 * V(pi_tilde, n + 4)^2) / Omega(pi_tilde, n + 4)

  # Valor critico de la t-distribucion
  alpha <- 1 - conf.level
  t <- qt(1 - alpha / 2, df = nu)

  # Error estandar
  se <- sqrt(V(pi_tilde, n + 4))

  # Calculo de los límites
  lower <- max(0, pi_tilde - t * se)
  upper <- min(1, pi_tilde + t * se)

  return(c(lower = lower, upper = upper))
}
ci_p_add_4 <- Vectorize(ci_p_add_4)
#'
#'
#'
#' Arcsine Wald Confidence Interval with Continuity Correction for Binomial Proportion
#'
#' @author David Esteban Cartagena Mejía, \email{dcartagena@unal.edu.co}
#'
#' @description
#' This function calculates the Arcsine Wald confidence interval with continuity correction for a Binomial proportion.
#' It transforms the proportion to the arcsine scale, applies the Wald interval with continuity correction, and then back-transforms the result to the original scale.
#' The method is vectorized, allowing for evaluation of single values or vectors.
#'
#' @param x A number or a vector with the number of successes.
#' @param n A number or a vector with the number of trials.
#' @param conf.level Confidence level for the returned confidence interval. By default, it is 0.95.
#'
#' @seealso \link{ci_p}.
#'
#' @references
#' Me falta la referencia.
#'
#' @details
#' The Arcsine Wald confidence interval with continuity correction adjusts the classical Wald interval by transforming the proportion \eqn{\pi} to the arcsine scale.
#' The parameter \eqn{\pi} represents the true proportion of successes in a Binomial experiment, defined as:
#'
#' \deqn{\pi = \frac{x \pm 0.5}{n},}
#'
#' where \eqn{x} is the number of successes and \eqn{n} is the number of trials.
#'
#' On the arcsine scale, the transformed parameter is given by:
#'
#' \deqn{\phi = \arcsin(\sqrt{\pi}).}
#'
#' The standard error on the arcsine scale is constant:
#'
#' \deqn{\text{se}(\phi) = \frac{1}{\sqrt{4n}},}
#'
#' where \eqn{n} is the number of trials. The confidence interval on the arcsine scale is:
#'
#' \deqn{\text{Lower}(\phi) = \max\left(0, \phi - z \cdot (\phi)\right),}
#' \deqn{\text{Upper}(\phi) = \min\left(\frac{\pi}{2}, \phi + z \cdot (\phi)\right),}
#'
#' where \eqn{z} is the critical value from the standard normal distribution at the specified confidence level.
#'
#' Back-transforming the limits to the original scale gives:
#'
#' \deqn{\text{Lower}(\pi) = \sin^2(\text{Lower}(\phi)),}
#' \deqn{\text{Upper}(\pi) = \sin^2(\text{Upper}(\phi)).}
#'
#' Special cases are handled explicitly:
#'
#' - If \eqn{x = 0}, the lower limit is 0, and the upper limit is calculated as \eqn{(\alpha / 2)^{1/n}}.
#'
#' - If \eqn{x = n}, the upper limit is 1, and the lower limit is calculated as \eqn{1 - (\alpha / 2)^{1/n}}.
#'
#' These adjustments ensure that the confidence interval is valid and well-behaved, even at the boundaries of the parameter space.
#'
#' @return A vector with the lower and upper limits of the confidence interval.
#'
#' @examples
#' arcsine_wald_cc_ci(x=15, n=50, conf.level=0.95)
#' arcsine_wald_cc_ci(x=0,  n=50, conf.level=0.95)
#' arcsine_wald_cc_ci(x=50, n=50, conf.level=0.95)
#'
#' @export
#'
ci_p_arcsine_cc <- function(x, n, conf.level=0.95) {
  alpha <- 1 - conf.level
  pi_hat1 <- (x - 0.5) / n
  pi_hat2 <- (x + 0.5) / n
  pi_hat1 <- max(0, pi_hat1)
  pi_hat2 <- min(1, pi_hat2)

  phi_hat1 <- asin(sqrt(pi_hat1))
  phi_hat2 <- asin(sqrt(pi_hat2))

  se_phi <- 1 / sqrt(4 * n)
  z <- qnorm(1 - (1 - conf.level) / 2)

  phi_lower <- max(0, phi_hat1 - z * se_phi)
  phi_upper <- min(pi / 2, phi_hat2 + z * se_phi)

  if (x == 0) {
    pi_lower <- 0
    pi_upper <- (alpha / 2)^(1 / n)
  } else if (x == n) {
    pi_lower <- 1 - (alpha / 2)^(1 / n)
    pi_upper <- 1
  } else {
    pi_lower <- sin(phi_lower)^2
    pi_upper <- sin(phi_upper)^2
  }

  return(c(lower=pi_lower, upper=pi_upper))
}
ci_p_arcsine_cc <- Vectorize(ci_p_arcsine_cc)
