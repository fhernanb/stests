#' Wald confidence interval for binomial proportion
#'
#' @author Olga Bustos, \email{oabustos@unal.edu.co}
#'
#' @description
#' This function calculates the confidence interval for a proportion.
#' It is vectorized, allowing users to evaluate it using either
#' single values or vectors.
#'
#' @param x a number or a vector with the number of successes.
#' @param n a number or a vector with the number of trials.
#' @param conf.level confidence level for the returned confidence interval.
#' By default is 0.95.
#'
#' @seealso \link{ci_p}.
#'
#' @references
#' Wald, A. (1949). Statistical decision functions. The Annals of
#' Mathematical Statistics, 165-205.
#'
#' @details
#' The expression to obtain the confidence interval is given below:
#'
#' \eqn{\hat{p} - z_{\alpha/2} \sqrt{\frac{\hat{p}(1 - \hat{p})}{n}} \leq p \leq \hat{p} + z_{\alpha/2} \sqrt{\frac{\hat{p}(1 - \hat{p})}{n}}},
#'
#' where \eqn{\hat{p}=\frac{x}{n}} is the sample proportion, \eqn{x} the
#' number of observed successes in the sample with size \eqn{n}. The
#' value \eqn{z_{\alpha/2}} is the \eqn{1-\alpha/2} percentile of the
#' standard normal distribution (e.g., \eqn{z_{0.025}=1.96} for a 95\%
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
  q <- qnorm(p=(1 + conf.level)/2)
  p <- x/n
  lower <- max(p - q*sqrt(p*(1 - p)/n), 0)
  upper <- min(p + q*sqrt(p*(1 - p)/n), 1)
  return(c(lower, upper))
}
ci_p_wald <- Vectorize(ci_p_wald)
#'
#'
#'
#' Agresti-Coull confidence interval for binomial proportion
#'
#' @author Omar David Mercado Turizo, \email{omercado@unal.edu.co}
#'
#' @description
#' This function calculates the confidence interval for a proportion.
#' It is vectorized, allowing users to evaluate it using either
#' single values or vectors.
#'
#' @param x A number or a vector with the number of successes.
#' @param n A number or a vector with the number of trials.
#' @param conf.level Confidence level for the returned confidence interval.
#' By default, it is 0.95.
#'
#' @seealso \link{ci_p}.
#'
#' @references
#' Agresti, A., & Coull, B. A. (1998). Approximate is better than “exact” for
#' interval estimation of binomial proportions. The American Statistician,
#' 52(2), 119-126.
#'
#' @details
#' The Agresti-Coull interval is an approximate confidence interval for the
#' binomial proportion \eqn{p}.
#' The limits are calculated based on an adjusted proportion \eqn{\tilde{p}}
#' and its standard error. The mathematical definitions are as follows:
#'  Adjusted proportion: \eqn{\tilde{p}=\frac{x + 2}{n + 4}};
#'  Adjusted standard error: \eqn{se=\sqrt{\frac{\tilde{p}(1 - \tilde{p})}{n + 4}}};
#'  Confidence limits: \eqn{\tilde{p} \pm z_{\alpha/2} \cdot se},
#'
#' where \eqn{z_{\alpha/2}} is the critical value of the standard normal
#' distribution. The limits are truncated to the range \eqn{[0, 1]}.
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
#' Rindskopf  confidence interval for binomial proportion
#'
#' @author David Esteban Cartagena Mejía, \email{dcartagena@unal.edu.co}
#'
#' @description
#' This function calculates the confidence interval for a proportion.
#' It is vectorized, allowing users to evaluate it using either
#' single values or vectors.
#'
#' @param x A number or a vector with the number of successes.
#' @param n A number or a vector with the number of trials.
#' @param conf.level Confidence level for the returned confidence interval.
#' By default, it is 0.95.
#'
#' @seealso \link{ci_p}.
#'
#' @references
#' Rindskopf, D. (2000). Commentary: Approximate is better than “exact”
#' for interval
#' estimation of binomial proportions. The American Statistician, 54, 88.
#'
#' @details
#' The expression to calculate the confidence interval according to the
#' Rindskopf approach is given by:
#'
#' \eqn{\phi=\text{logit}(\pi)=\log\left(\frac{\pi}{1 - \pi}\right)},
#'
#' where the maximum likelihood estimator for \eqn{\phi} is:
#'
#' \eqn{\hat{\phi}_{ML}=\log\left(\frac{x + 0.5}{n - x + 0.5}\right)},
#'
#' and its standard error is:
#'
#' \eqn{\text{se}(\hat{\phi}_{ML})=\sqrt{\frac{1}{x + 0.5} + \frac{1}{n - x + 0.5}}}.
#'
#' The adjustment of adding 0.5 successes and non-successes ensures that
#' intervals can also be computed for the cases where \eqn{x=0} or \eqn{x=n}
#' (where otherwise the maximum likelihood estimator and standard error would
#' be infinite).
#'
#' Since the scale of \eqn{\phi} is \eqn{(- \infty, \infty)}, this interval
#' respects the boundary constraints. Back-transformation to the scale
#' of \eqn{\pi} is performed using the inverse logit function:
#'
#' \eqn{\pi=\text{expit}(\phi)=\frac{\exp(\phi)}{1 + \exp(\phi)}}.
#'
#' Thus, the confidence interval for \eqn{\pi} in the original scale is the
#' Rindskopf confidence interval, as proposed by Rindskopf.
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
#' Clopper-Pearson confidence interval for binomial proportion
#'
#' @author Omar David Mercado Turizo, \email{omercado@unal.edu.co}
#'
#' @description
#' This function calculates the confidence interval for a proportion.
#' It is vectorized, allowing users to evaluate it using either
#' single values or vectors.
#'
#' @param x A number or a vector with the number of successes.
#' @param n A number or a vector with the number of trials.
#' @param conf.level Confidence level for the returned confidence interval.
#' By default, it is 0.95.
#'
#' @seealso \link{ci_p}.
#'
#' @references
#' Clopper, C. J., & Pearson, E. S. (1934). The use of confidence or fiducial
#' limits illustrated in the case of the binomial. Biometrika, 26(4), 404-413.
#'
#' @details
#' The Clopper-Pearson interval is an exact confidence interval for the
#' binomial proportion \eqn{p}.
#' The limits of the interval are derived based on the Beta distribution.
#' For the special cases where \eqn{x=0} or \eqn{x=n}, the limits are
#' calculated directly.
#'
#' The mathematical definitions are as follows:
#'
#' - If \eqn{x=0}, the lower limit is \eqn{0}, and the upper limit is \eqn{1 - (\alpha / 2)^{1/n}}.
#'
#' - If \eqn{x=n}, the lower limit is \eqn{(\alpha / 2)^{1/n}}, and the
#' upper limit is \eqn{1}.
#'
#' Otherwise, the limits are given by:
#'
#' \deqn{\text{Lower Limit}=B_{1-\alpha/2, x, n-x+1}}
#'
#' \deqn{\text{Upper Limit}=B_{\alpha/2, x+1, n-x}}
#'
#' where \eqn{B_{\omega, a, b}} is the \eqn{100\%(1-\omega)}
#' percentile of the Beta distribution with parameters
#' \eqn{a} and \eqn{b}.
#'
#' Due to the relationship between Beta and F distributions,
#' the limits can be written as:
#'
#'  \deqn{\text{Lower Limit}=\frac{1}{1+\frac{n-x+1}{x}F_{\alpha/2, \, 2(n-x+1), \, 2x}}}
#'
#'  \deqn{\text{Upper Limit}=\frac{\frac{x+1}{n-x} F_{\alpha/2, \, 2(x+1), \, 2(n-x)}}{1+\frac{x+1}{n-x} F_{\alpha/2, \, 2(x+1), \, 2(n-x)}}}
#'
#' where \eqn{F_{\omega, a, b}} is the \eqn{100\%(1-\omega)}
#' percentile of the F distribution with parameters
#' \eqn{a} and \eqn{b}.
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
    lower <- qbeta(p=1-alpha/2, shape1=x, shape2=n-x+1, lower.tail=F)
    upper <- qbeta(p=alpha/2, shape1=x+1, shape2=n-x, lower.tail=F)
  }

  return(c(lower, upper))
}
ci_p_clopper_pearson <- Vectorize(ci_p_clopper_pearson)
#'
#'
#'
#' Add-4 Wald-t Confidence Interval for binomial proportion
#'
#' @author David Esteban Cartagena Mejía, \email{dcartagena@unal.edu.co}
#'
#' @description
#' This function calculates the confidence interval for a proportion.
#' It is vectorized, allowing users to evaluate it using either
#' single values or vectors.
#'
#' @param x A number or a vector with the number of successes.
#' @param n A number or a vector with the number of trials.
#' @param conf.level Confidence level for the returned confidence interval. By default, it is 0.95.
#'
#' @seealso \link{ci_p}.
#'
#' @references
#' Pan, W. (2002). Approximate confidence intervals for one proportion and
#' difference of two proportions. Computational Statistics and Data Analysis,
#' 40(1), 143–157
#'
#' @details
#' The Add-4 Wald-t confidence interval improves the performance of the Wald
#' interval by adding 2 successes and 2 failures to the observed data,
#' effectively modifying the estimated proportion:
#'
#' \deqn{\hat{p}=\frac{x + 2}{n + 4}.}
#'
#' The variance \eqn{V(\hat{p}, n+4)} is given by:
#'
#' \deqn{V(\hat{p}, n+4)=\frac{\hat{p}(1 - \hat{p})}{n + 4}.}
#'
#' The degrees of freedom \eqn{\nu} are calculated using equation (2.9):
#'
#' \deqn{\nu=\frac{2 V(\hat{p}, n+4)^2}{\Omega(\hat{p}, n+4)},}
#'
#' where \eqn{\Omega(\hat{p}, n+4)} is defined as:
#'
#' \deqn{\Omega(p, n)=\frac{p - p^2}{n^3} + \frac{p + (6n - 7)p^2 + 4(n - 1)(n - 3)p^3 - 2(n - 1)(2n - 3)p^4}{n^5} - \frac{2(p + (2n - 3)p^2 - 2(n - 1)p^3)}{n^4}.}
#'
#' The confidence interval is then calculated as:
#'
#' \deqn{\text{Lower}=\hat{p} - t \cdot \sqrt{\frac{\tilde{\pi}(1 - \hat{p})}{n+4}},}
#' \deqn{\text{Upper}=\hat{p} + t \cdot \sqrt{\frac{\tilde{\pi}(1 - \hat{p})}{n+4}},}
#'
#' where \eqn{t} is the critical value from the t-distribution with \eqn{\nu}
#' degrees of freedom.
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
  # Proporción ajustada
  pi_tilde <- (x + 2) / (n + 4)

  # Función para calcular Omega
  Omega <- function(p, n) {
    term1 <- (p - p^2) / n^3
    term2 <- (p + (6 * n - 7) * p^2 + 4 * (n - 1) * (n - 3) * p^3 - 2 * (n - 1) * (2 * n - 3) * p^4) / n^5
    term3 <- (2 * (p + (2 * n - 3) * p^2 - 2 * (n - 1) * p^3)) / n^4
    return(term1 + term2 - term3)
  }

  # Variancia de p ajustada
  V <- function(p, n) {
    return(p * (1 - p) / n)
  }

  # Grados de libertad (nu)
  nu <- (2 * V(pi_tilde, n + 4)^2) / Omega(pi_tilde, n + 4)

  # Valor crítico de la t-distribución
  alpha <- 1 - conf.level
  t <- qt(1 - alpha / 2, df=nu)

  # Error estándar
  se <- sqrt(V(pi_tilde, n + 4))

  # Cálculo de los límites
  lower <- max(0, pi_tilde - t * se)
  upper <- min(1, pi_tilde + t * se)

  return(c(lower, upper))
}
ci_p_add_4 <- Vectorize(ci_p_add_4)
#'
#'
#'
#' ArcSine confidence interval for binomial proportion
#'
#' @author Victor David Usuga Duque, \email{vusuga@unal.edu.co}
#'
#' @description
#' This function calculates the confidence interval for a proportion.
#' It is vectorized, allowing users to evaluate it using either
#' single values or vectors.
#'
#' @param x a number or a vector with the number of successes.
#' @param n a number or a vector with the number of trials.
#' @param conf.level confidence level for the returned confidence interval.
#' By default is 0.95.
#'
#' @seealso \link{ci_p}.
#'
#' @references
#' No reference.
#'
#' @details
#' The expression to obtain the confidence interval is given below:
#'
#' \eqn{\sin^2(\arcsin(\sqrt{\tilde{p}}) \mp \frac{z_{\alpha/2}}{2\sqrt{n}})},
#'
#' where \eqn{\tilde{p} = x/n}.
#' The value \eqn{z_{\alpha/2}} is the \eqn{1-\alpha/2} percentile of the
#' standard normal distribution (e.g., \eqn{z_{0.025}=1.96} for a 95\%
#' confidence interval).
#'
#' @return A matrix with the lower and upper limits.
#'
#' @examples
#' ci_p_arcsine(x= 0, n=50, conf.level=0.95)
#' ci_p_arcsine(x=15, n=50, conf.level=0.95)
#' ci_p_arcsine(x=50, n=50, conf.level=0.95)
#'
#' @export
#'
ci_p_arcsine <- function(x, n, conf.level=0.95) {
  pi_hat <- x / n # Estimador de maxima verosimilitud para π
  phi_hat <- asin(sqrt(pi_hat)) # Transformacion arcoseno
  # Desviacion estandar aproximada para phi
  se_phi <- 1 / sqrt(4*n)
  # Valor critico z para el nivel de confianza deseado
  z <- qnorm(1 - (1 - conf.level) / 2)
  # Intervalo en la escala de phi
  phi_lower <- max(0, phi_hat - z * se_phi) # Truncar a 0 si es necesario
  phi_upper <- min(pi / 2, phi_hat + z * se_phi) # Truncar a π/2 si es necesario
  lower <- sin(phi_lower)^2
  upper <- sin(phi_upper)^2
  return(c(lower, upper))
}
ci_p_arcsine <- Vectorize(ci_p_arcsine)
#'
#'
#'
#' Arcsine Wald Confidence Interval with Continuity Correction Anscombe for
#' binomial proportion
#'
#' @author David Esteban Cartagena Mejía, \email{dcartagena@unal.edu.co}
#'
#' @description
#' This function calculates the confidence interval for a proportion.
#' It is vectorized, allowing users to evaluate it using either
#' single values or vectors.
#'
#' @param x A number or a vector with the number of successes.
#' @param n A number or a vector with the number of trials.
#' @param conf.level Confidence level for the returned confidence interval.
#' By default, it is 0.95.
#'
#' @seealso \link{ci_p}.
#'
#' @references
#' Anscombe, F. J. (1948). Transformations of Poisson, binomial and
#' negative-binomial data. Biometrika, 35, 246–254.
#'
#' @details
#' The expression to obtain the confidence interval is given below:
#'
#' \eqn{\sin^2(\arcsin(\sqrt{\tilde{p}}) \mp \frac{z_{\alpha/2}}{2\sqrt{n}})},
#'
#' where \eqn{\tilde{p} = \frac{x+3/8}{n+3/4}}.
#' The value \eqn{z_{\alpha/2}} is the \eqn{1-\alpha/2} percentile of the
#' standard normal distribution (e.g., \eqn{z_{0.025}=1.96} for a 95\%
#' confidence interval).
#'
#' @return A vector with the lower and upper limits of the confidence interval.
#'
#' @examples
#' ci_p_arcsine_anscombe(x= 0, n=50, conf.level=0.95)
#' ci_p_arcsine_anscombe(x=15, n=50, conf.level=0.95)
#' ci_p_arcsine_anscombe(x=50, n=50, conf.level=0.95)
#'
#' @export
#'
ci_p_arcsine_anscombe <- function(x, n, conf.level=0.95) {
  pi_hat <- (x+3/8) / (n + 3/4) # Estimador de maxima verosimilitud para π
  phi_hat <- asin(sqrt(pi_hat)) # Transformacion arcoseno
  # Desviacion estandar aproximada para phi
  se_phi <- 1 / sqrt(4*n)
  # Valor critico z para el nivel de confianza deseado
  z <- qnorm(1 - (1 - conf.level) / 2)
  # Intervalo en la escala de phi
  phi_lower <- max(0, phi_hat - z * se_phi) # Truncar a 0 si es necesario
  phi_upper <- min(pi / 2, phi_hat + z * se_phi) # Truncar a π/2 si es necesario
  lower <- sin(phi_lower)^2
  upper <- sin(phi_upper)^2
  return(c(lower, upper))
}
ci_p_arcsine_anscombe <- Vectorize(ci_p_arcsine_anscombe)
#'
#'
#'
#' Wilson confidence interval for binomial proportion
#'
#' @author Victor David Usuga, \email{vusuga@unal.edu.co}
#'
#' @description
#' This function calculates the confidence interval for a proportion.
#' It is vectorized, allowing users to evaluate it using either
#' single values or vectors.
#'
#' @param x a number or a vector with the number of successes.
#' @param n a number or a vector with the number of trials.
#' @param conf.level confidence level for the returned confidence interval.
#' By default is 0.95.
#'
#' @seealso \link{ci_p}.
#'
#' @references
#' Wilson, E. B. (1927). Probable inference, the law of succession, and
#' statistical inference. Journal of the American Statistical Association,
#' 22(158), 209-212.
#'
#' @details
#' The expression to obtain the confidence interval is given below:
#'
#' \eqn{\frac{\hat{p}+ \frac{z_{\alpha/2}^2}{2n}}{\widetilde{n}} \pm \frac{z_{\alpha/2}^2}{\widetilde{n}}  \sqrt{ (\hat{p}(1 - \hat{p}) + \frac{z_{\alpha/2}^2}{4n} )/n}},
#'
#' where \eqn{\hat{p}=\frac{x}{n}} is the sample proportion, \eqn{\widetilde{n}=1 + \frac{ z_{\alpha/2}^2}{n}}, \eqn{x} the
#' number of observed successes in the sample with size \eqn{n}.
#' The value \eqn{z_{\alpha/2}} is the \eqn{1-\alpha/2} percentile of the
#' standard normal distribution (e.g., \eqn{z_{0.025}=1.96} for a 95\%
#' confidence interval).
#'
#' @return A matrix with the lower and upper limits.
#'
#' @examples
#' ci_p_wilson(x= 0, n=50, conf.level=0.95)
#' ci_p_wilson(x=15, n=50, conf.level=0.95)
#' ci_p_wilson(x=50, n=50, conf.level=0.95)
#'
#' @export
#'
ci_p_wilson <- function(x, n, conf.level=0.95) {
  pi_ml <- x / n
  q_alpha <- qnorm(1 - (1 - conf.level) / 2)
  numerator <- x + (q_alpha^2) / 2
  denominator <- n + q_alpha^2
  adjustment <- (q_alpha * sqrt(n) / denominator) * sqrt(pi_ml * (1 - pi_ml) + (q_alpha^2) / (4 * n))
  lower <- numerator / denominator - adjustment
  upper <- numerator / denominator + adjustment
  return(c(lower, upper))
}
ci_p_wilson <- Vectorize(ci_p_wilson)
#'
#'
#'
#' Bayesian confidence interval for binomial proportion using Jeffreys
#' prior (non-informative prior).
#'
#' @author Rusvelt Jose Meza San Martin, \email{rmezas@unal.edu.co}
#'
#' @description
#' This function calculates the confidence interval for a proportion.
#' It is vectorized, allowing users to evaluate it using either
#' single values or vectors.
#'
#' @param x A number or a vector with the number of successes.
#' @param n A number or a vector with the number of trials.
#' @param conf.level Confidence level for the returned confidence interval.
#' By default, it is 0.95.
#'
#' @seealso \link{ci_p}.
#'
#' @details
#' The Jeffreys prior is a non-informative prior (\eqn{Beta(0.5, 0.5)})
#' based on the Fisher information. The posterior distribution is Beta:
#'
#' \deqn{p | x \sim \text{Beta}(x+0.5, n-x+0.5)}
#'
#' The limits of the Bayesian confidence interval are derived from the
#' quantiles of the posterior distribution:
#'
#' \deqn{\text{Lower Limit}=B_{1-\alpha/2, x+0.5, n-x+0.5}}
#'
#' \deqn{\text{Upper Limit}=B_{\alpha/2, x+0.5, n-x+0.5}}
#'
#' where \eqn{B_{\omega, a, b}} is the \eqn{100\%(1-\omega)}
#' percentile of the Beta distribution with parameters
#' \eqn{a} and \eqn{b}.
#'
#' @return A vector with the lower and upper limits.
#'
#' @examples
#' # Example with a single value
#' ci_p_jeffreys(x=5, n=20, conf.level=0.95)
#'
#' # Example with vectors
#' x_values <- c(5, 10, 15)
#' n_values <- c(20, 30, 40)
#' ci_p_jeffreys(x=x_values, n=n_values, conf.level=0.95)
#'
#' @export
ci_p_jeffreys <- function(x, n, conf.level=0.95) {
  # x: número de éxitos
  # n: tamaño de la muestra
  # conf.level: nivel de confianza

  alpha <- 1 - conf.level

  # Limite inferior (cuantil inferior de la distribución beta posterior)
  lower <- qbeta(alpha / 2, 0.5 + x, 0.5 + n - x)

  # Limite superior (cuantil superior de la distribución beta posterior)
  upper <- qbeta(1 - alpha / 2, 0.5 + x, 0.5 + n - x)

  # Retornar los límites
  return(c(lower, upper))
}
# Vectorizar la función
ci_p_jeffreys <- Vectorize(ci_p_jeffreys)
#'
#'
#'
#' Highest Posterior Density (HPD) interval for binomial proportion using
#' Jeffreys prior
#'
#' @author Rusvelt Jose Meza San Martín, \email{rmezas@unal.edu.co}
#'
#' @description
#' This function calculates the Highest Posterior Density (HPD) interval
#' for a Binomial proportion using Jeffreys prior. It is vectorized,
#' allowing the evaluation of single values or vectors.
#'
#' @param x A number or a vector with the number of successes.
#' @param n A number or a vector with the number of trials.
#' @param conf.level Confidence level for the returned confidence interval.
#' By default, it is 0.95.
#'
#' @seealso \link{ci_p}.
#'
#' @details
#' The Jeffreys prior is a non-informative prior (\eqn{Beta(0.5, 0.5)})
#' based on the Fisher information. The posterior distribution is Beta:
#'
#' \deqn{p | x \sim \text{Beta}(0.5 + x, 0.5 + n - x)}
#'
#' The HPD interval is calculated using the posterior samples from the
#' Beta distribution:
#'
#' \deqn{	\text{HPD Interval}=[L, U]}
#'
#' where \eqn{L} and \eqn{U} are the bounds of the smallest interval
#' containing \eqn{1 - \alpha} posterior probability.
#'
#' @return A vector with the lower and upper limits of the HPD interval.
#'
#' @examples
#' # Example with a single value
#' ci_p_hpd_jeffreys(x=5, n=20, conf.level=0.95)
#'
#' # Example with vectors
#' x_values <- c(5, 10, 15)
#' n_values <- c(20, 30, 40)
#' ci_p_hpd_jeffreys(x=x_values, n=n_values, conf.level=0.95)
#'
#' @export
#' @importFrom HDInterval hdi
ci_p_hpd_jeffreys <- function(x, n, conf.level=0.95) {

  res <- HDInterval::hdi(qbeta, credMass=conf.level,
                         shape1=0.5+x, shape2=0.5+n-x)
  return(as.matrix(res))
}
# Vectorizar la función
ci_p_hpd_jeffreys <- Vectorize(ci_p_hpd_jeffreys)
#'
#'
#'
#'
#'
#'
#' Likelihood Ratio confidence interval for binomial proportion
#'
#' @author David Esteban Cartagena Mejía, \email{dcartagena@unal.edu.co}
#'
#' @description
#' This function calculates the Likelihood Ratio (LRT) confidence interval for
#' a binomial proportion.
#' It is vectorized, allowing the evaluation of single values or vectors.
#'
#' @param x A number or a vector with the number of successes.
#' @param n A number or a vector with the number of trials.
#' @param conf.level Confidence level for the returned confidence interval.
#' By default, it is 0.95.
#'
#' @seealso \link{ci_p}.
#'
#' @references
#' Somerville, M. C., & Brown, R. S. (2013). Exact likelihood ratio
#' and score confidence intervals for the binomial proportion.
#' Pharmaceutical statistics, 12(3), 120-128.
#'
#' @details
#' This function computes the confidence interval for a binomial proportion
#' based on the Likelihood Ratio Test (LRT). The confidence interval is
#' defined as the set of values for \eqn{p} that satisfy the following
#' condition:
#'
#' \deqn{-2 \log \left(\frac{L(p)}{L(\hat{p}_{ML})}\right) \leq \chi^2_{\gamma}(1),}
#'
#' where \eqn{L(p)} is the likelihood function for the binomial model,
#' \eqn{\hat{p}_{ML} = x / n} is the maximum likelihood estimator
#' for \eqn{p}, and \eqn{\chi^2_{\gamma}(1)} is the \eqn{1-\gamma}
#' quantile of the chi-square distribution with 1
#' degree of freedom.
#'
#' The confidence limits are calculated numerically using the \code{uniroot}
#' function in R. Special care is taken
#' to handle edge cases where \eqn{x = 0} or \eqn{x = n}, setting the
#' limits to 0 or 1, respectively.
#'
#' @return A vector with the lower and upper limits.
#'
#' @examples
#' ci_p_LRT(x = 15, n = 50, conf.level = 0.95)
#'
#' @export
#'
ci_p_LRT <- function(x, n, conf.level = 0.95) {
  chi_crit <- qchisq(conf.level, df = 1)
  alpha <- 1 - conf.level

  BinLLR <- function(x, n, p) {
    phat <- x / n
    if (p == 0 || p == 1) return(Inf) # Evita problemas de logaritmos
    -2 * (x * log(p / phat) + (n - x) * log((1 - p) / (1 - phat)))
  }

  LRT_function <- function(p) {
    BinLLR(x, n, p) - chi_crit
  }

  # Limites del intervalo
  if (x == 0) {
    upper <- (alpha / 2)^(1 / n)
    lower <- 0
  } else if (x == n) {
    lower <- 1 - (alpha / 2)^(1 / n)
    upper <- 1
  } else {
    lower <- tryCatch(
      uniroot(LRT_function,
              lower = 1e-10,
              upper = x / n - 1e-5, tol = 1e-8)$root,
      error = function(e) 0
    )
    upper <- tryCatch(
      uniroot(LRT_function,
              lower = x / n + 1e-5,
              upper = 1 - 1e-10, tol = 1e-8)$root,
      error = function(e) 1
    )
  }
  return(c(lower, upper))
}
ci_p_LRT <- Vectorize(ci_p_LRT)
#'
#'
#'
#' Mid-p confidence interval for binomial proportion
#'
#' @author Rusvelt Jose Meza San Martín, \email{rmezas@unal.edu.co}
#'
#' @description
#' This function calculates the mid-p confidence interval for a
#' binomial proportion. It is vectorized, allowing the evaluation
#' of single values or vectors.
#'
#' @param x A number or a vector with the number of successes.
#' @param n A number or a vector with the number of trials.
#' @param conf.level confidence level for the returned confidence interval.
#' By default is 0.95.
#'
#' @seealso \link{ci_p}.
#'
#' @references
#' Berry, G., & Armitage, P. (1995). Mid-P confidence intervals:
#' a brief review. Journal of the Royal Statistical Society Series
#' D: The Statistician, 44(4), 417-423.
#'
#' @details
#' The mid-p confidence interval adjusts the exact binomial
#' interval by averaging the tail probabilities of the observed value.
#' The limits are found by solving the following equations:
#'
#' \deqn{\sum_{j=0}^{x-1} P(X = j) + 0.5 \cdot P(X = x) = \alpha / 2}
#'
#' \deqn{\sum_{j=x+1}^{n} P(X = j) + 0.5 \cdot P(X = x) = \alpha / 2}
#'
#' where \eqn{P(X = j)} is the binomial probability.
#'
#' @return A vector with the lower and upper limits of the confidence
#' interval.
#'
#' @examples
#' # Example with a single value
#' ci_p_mid_p(x = 5, n = 20, conf.level = 0.95)
#'
#' # Example with vectors
#' x_values <- c(5, 10, 15)
#' n_values <- c(20, 30, 40)
#' ci_p_mid_p(x = x_values, n = n_values, conf.level = 0.95)
#'
#' @export
ci_p_mid_p <- function(x, n, conf.level = 0.95) {

  alpha <- 1 - conf.level

  # Upper limit
  upper_limit <- uniroot(
    function(p) {
      sum(sapply(0:(x - 1), function(j) dbinom(x=j, prob=p, size=n))) +
        0.5 * dbinom(x=x, prob=p, size=n) - alpha / 2
    },
    interval = c(0, 1),
    tol = 1e-8
  )$root

  # Lower limit
  lower_limit <- uniroot(
    function(p) {
      sum(sapply((x + 1):n, function(j) dbinom(x=j, prob=p, size=n))) +
        0.5 * dbinom(x=x, prob=p, size=n) - alpha / 2
    },
    interval = c(0, 1),
    tol = 1e-8
  )$root

  # Retornar los límites
  return(c(lower_limit, upper_limit))
}
ci_p_mid_p <- Vectorize(ci_p_mid_p)
#'
#'
#'
#' Agresti-Caffo confidence interval for binomial proportion
#'
#' @author David Esteban Cartagena Mejía, \email{dcartagena@unal.edu.co}
#'
#' @description
#' This function calculates the Agresti-Caffo confidence interval for
#' a Binomial
#' proportion, which includes a Bayesian adjustment by adding 2 successes
#' and 2 failures to stabilize the estimate. It is vectorized,
#' allowing the evaluation of single values or vectors.
#'
#' @param x A number or a vector with the number of successes.
#' @param n A number or a vector with the number of trials.
#' @param conf.level Confidence level for the returned confidence interval.
#' By default, it is 0.95.
#'
#' @seealso \link{ci_p}.
#'
#' @references
#' Agresti, Alan, and Brian Caffo. "Simple and effective confidence
#' intervals for proportions and differences of proportions result
#' from adding two successes and two failures." The American
#' Statistician 54.4 (2000): 280-288.
#'
#' @details
#' The Agresti-Caffo confidence interval incorporates a simple
#' Bayesian adjustment by adding 2 successes and 2 failures
#' to the data. The adjusted proportion is calculated as:
#'
#' \deqn{\hat{p} = \frac{x + 2}{n + 4}}
#'
#' The standard error for the adjusted proportion is given by:
#'
#' \deqn{\text{se} = \sqrt{\frac{\hat{p} (1 - \hat{p})}{n + 4}}.}
#'
#' The confidence interval is then constructed using the \eqn{t}-distribution
#' with \eqn{n - 1} degrees of freedom:
#'
#' \deqn{\text{CI} = \hat{p} \pm t_{n-1, \alpha/2} \cdot \text{se},}
#'
#' where \eqn{t_{n-1, \alpha/2}} is the critical value of the
#' \eqn{t}-distribution at a two-tailed significance level of
#' \eqn{\alpha = 1 - \text{conf.level}}.
#'
#' @return A vector with the lower and upper limits of the confidence interval.
#'
#' @examples
#' # Example with a single value
#' ci_p_agresti_caffo(x = 15, n = 50, conf.level = 0.95)
#'
#'
#' @export
#'
ci_p_agresti_caffo <- function(x, n, conf.level = 0.95) {
  p_hat <- (x + 2) / (n + 4)
  se <- sqrt((p_hat * (1 - p_hat)) / (n + 4))
  t <- qt(1 - (1 - conf.level) / 2, n - 1)
  alpha <- 1-conf.level
  # Limits
  if (x == 0) {
    upper <- (alpha / 2)^(1 / n)
    lower <- 0
  } else if (x == n) {
    lower <- 1 - (alpha / 2)^(1 / n)
    upper <- 1
  } else {
    # Limits
    if (x == 0) {
      upper <- (alpha / 2)^(1 / n)
      lower <- 0
    } else if (x == n) {
      lower <- 1 - (alpha / 2)^(1 / n)
      upper <- 1
    } else {
      lower <- p_hat - t * se
      upper <- p_hat + t * se
    }
  }

  return(c(lower, upper))
}
ci_p_agresti_caffo <- Vectorize(ci_p_agresti_caffo)
#'
#'
#'
#' Score confidence interval with continuity correction for binomial proportion
#'
#' @author Omar David Mercado Turizo, \email{omercado@unal.edu.co}
#'
#' @description
#' This function calculates the score confidence interval with continuity
#' correction for a binomial proportion. It is vectorized, allowing the
#' evaluation of single values or vectors.
#'
#' @param x A number or a vector with the number of successes.
#' @param n A number or a vector with the number of trials.
#' @param conf.level Confidence level for the returned confidence interval.
#' By default, it is 0.95.
#'
#' @references
#' Missing reference.
#'
#' @seealso \link{ci_p}.
#'
#' @details
#' The score confidence interval with continuity correction is an
#' adjusted interval for the binomial proportion \eqn{p}.
#'
#' The mathematical definitions are as follows:
#'
#' Lower limit: \eqn{\frac{2np + z^2 - 1 - z \sqrt{z^2 - 2 - \frac{1}{n} + 4p(nq + 1)}}{2n + 2z^2}}.
#'
#' Upper limit: \eqn{\frac{2np + z^2 + 1 + z \sqrt{z^2 + 2 - \frac{1}{n} + 4p(nq - 1)}}{2n + 2z^2}}.
#'
#' Where \eqn{p = x / n} is the sample proportion, \eqn{q = 1 - p} its complement, and \eqn{z} is the critical value of the standard normal distribution.
#'
#' The limits are truncated to the range \eqn{[0, 1]}.
#'
#' @return A vector with the lower and upper limits.
#'
#' @examples
#' # Example with a single value
#' ci_p_score_cc(x = 15, n = 50, conf.level = 0.95)
#' @export
#'
ci_p_score_cc <- function(x, n, conf.level = 0.95) {
  z <- qnorm((1 + conf.level) / 2) # Valor crítico
  p <- x / n                       # Proporción muestral
  q <- 1 - p                       # Complemento de la proporción

  # Cálculo del límite inferior
  lower <- (2 * n * p + z^2 - 1 - z * sqrt(z^2 - 2 - (1 / n) + 4 * p * (n * q + 1))) /
    (2 * n + 2 * z^2)
  lower <- max(lower, 0) # Garantizar que no sea menor a 0

  # Cálculo del límite superior
  upper <- (2 * n * p + z^2 + 1 + z * sqrt(z^2 + 2 - (1 / n) + 4 * p * (n * q - 1))) /
    (2 * n + 2 * z^2)
  upper <- min(upper, 1) # Garantizar que no sea mayor a 1

  return(c(lower, upper))
}
ci_p_score_cc <- Vectorize(ci_p_score_cc)
#'
#'
#'
#' Highest Posterior Density (HPD) interval for binomial proportion
#'
#' @author Omar David Mercado Turizo, \email{omercado@unal.edu.co}
#'
#' @description
#' This function calculates the Highest Posterior Density (HPD) interval
#' for a binomial proportion using a Bayesian approach. It is vectorized,
#' allowing the evaluation of single values or vectors.
#'
#' @param x A number or a vector with the number of successes.
#' @param n A number or a vector with the number of trials.
#' @param conf.level Confidence level for the returned confidence interval.
#' By default, it is 0.95.
#' @param prior The prior distribution to use. Options are
#' "uniform" (default) or "jeffreys".
#'
#' @references
#' Missing reference.
#'
#' @seealso \link{ci_p}.
#'
#' @details
#' The HPD interval is a Bayesian credible interval for the Binomial
#' proportion \eqn{p}. The posterior distribution is calculated based
#' on the Beta prior:
#'
#' - "uniform": \eqn{\text{Beta}(1, 1)}.
#'
#' - "jeffreys": \eqn{\text{Beta}(0.5, 0.5)}.
#'
#' The limits of the interval are computed using the quantiles
#' of the Beta posterior distribution:
#'
#' - Lower limit: \eqn{\text{qbeta}((1 - \text{conf.level}) / 2, \alpha + x, \beta + n - x)}.
#'
#' - Upper limit: \eqn{\text{qbeta}(1 - (1 - \text{conf.level}) / 2, \alpha + x, \beta + n - x)}.
#'
#' @return A vector with the lower and upper limits.
#'
#' @examples
#' # Example with a single value
#' ci_p_hpd(x = 15, n = 50, conf.level = 0.95)
#' @export
#'
ci_p_hpd <- function(x, n, conf.level = 0.95, prior = "uniform") {
  # Prior's parameters
  if (prior == "uniform") {
    alpha_prior <- 1
    beta_prior <- 1
  } else if (prior == "jeffreys") {
    alpha_prior <- 0.5
    beta_prior <- 0.5
  } else {
    stop("Unkwon prior, please use 'uniform' or 'jeffreys'.")
  }

  # Posterior Beta
  alpha_post <- alpha_prior + x
  beta_post  <- beta_prior + n - x

  # Limits HPD
  lower <- qbeta((1 - conf.level) / 2, alpha_post, beta_post)
  upper <- qbeta(1 - (1 - conf.level) / 2, alpha_post, beta_post)

  return(c(lower, upper))
}
ci_p_hpd <- Vectorize(ci_p_hpd)
#'
#'
#'
#' Wald Continuity-Corrected confidence interval for binomial proportion
#'
#' @author David Esteban Cartagena Mejía, \email{dcartagena@unal.edu.co}
#'
#' @description
#' This function calculates the Wald continuity-corrected confidence interval
#' for a binomial proportion.
#' It is vectorized, allowing the evaluation of single values or vectors.
#'
#' @param x A number or a vector with the number of successes.
#' @param n A number or a vector with the number of trials.
#' @param conf.level Confidence level for the returned confidence interval.
#' By default, it is 0.95.
#'
#' @references
#' Pires, Ana M., and Conceiçao Amado. "Interval estimators for a binomial
#' proportion: Comparison of twenty methods".
#' REVSTAT-Statistical Journal 6.2 (2008): 165-197.
#'
#' @seealso \link{ci_p}.
#'
#' @details
#' The Wald continuity-corrected confidence interval adjusts the standard Wald interval for small sample sizes or when the proportion
#' \eqn{\hat{p}} is near 0 or 1. It incorporates a continuity correction to improve accuracy.
#'
#' The estimated proportion is given by:
#'
#' \deqn{\hat{p} = \frac{x}{n},}
#'
#' and its complement is:
#'
#' \deqn{\hat{q} = 1 - \hat{p}.}
#'
#' The continuity-corrected interval incorporates the critical value \eqn{z} from the standard normal distribution:
#'
#' \deqn{\text{Lower} = \hat{p} - z \sqrt{\frac{\hat{p}\hat{q}}{n}} - \frac{1}{2n},}
#' \deqn{\text{Upper} = \hat{p} + z \sqrt{\frac{\hat{p}\hat{q}}{n}} + \frac{1}{2n}.}
#'
#'
#' These adjustments ensure that the confidence interval is valid even
#' at the boundaries of the parameter space.
#'
#' @return A vector with the lower and upper limits of the confidence interval.
#'
#' @examples
#' ci_p_wald_cc(x = 15, n = 50, conf.level = 0.95)
#' ci_p_wald_cc(x = 0, n = 50, conf.level = 0.95)
#' ci_p_wald_cc(x = 50, n = 50, conf.level = 0.95)
#'
#' @export
#'
ci_p_wald_cc <- function(x, n, conf.level = 0.95) {
  alpha <- 1 - conf.level
  z <- qnorm(1 - alpha / 2)
  p_hat <- x / n
  q_hat <- 1 - p_hat

  # Limites del intervalo
  if (x == 0) {
    upper <- (alpha / 2)^(1 / n)
    lower <- 0
  } else if (x == n) {
    lower <- 1 - (alpha / 2)^(1 / n)
    upper <- 1
  } else {
    lower <-  p_hat - z * sqrt(p_hat * q_hat / n) - 1 / (2 * n)
    upper <-  p_hat + z * sqrt(p_hat * q_hat / n) + 1 / (2 * n)
  }

  return(c(lower,upper))
}
ci_p_wald_cc <- Vectorize(ci_p_wald_cc)
#'
#'
#'
#' Recentered Wald Interval for binomial proportion
#'
#' @author David Esteban Cartagena Mejía, \email{dcartagena@unal.edu.co}
#'
#' @description
#' This function calculates the recentered Wald confidence interval for a
#' binomial proportion.
#' It adjusts the classical Wald interval to improve accuracy near the
#' boundaries of the parameter space.
#' The method is vectorized, allowing for evaluation of single values or vectors.
#'
#' @param x A number or a vector with the number of successes.
#' @param n A number or a vector with the number of trials.
#' @param conf.level Confidence level for the returned confidence interval.
#' By default, it is 0.95.
#'
#' @references
#' Pires, Ana M., and Conceiçao Amado. "Interval estimators for a binomial
#' proportion: Comparison of twenty methods".
#' REVSTAT-Statistical Journal 6.2 (2008): 165-197.
#'
#' @seealso \link{ci_p}.
#'
#' @details
#' The recentered Wald interval modifies the classical Wald interval by
#' incorporating a recentering term and applying bounds to ensure that the
#' interval remains within the parameter space \eqn{[0, 1]}.
#'
#' The critical value \eqn{z} is obtained from the standard normal
#' distribution for the specified confidence level:
#'
#' \deqn{z = \Phi^{-1}(1 - \alpha / 2),}
#'
#' where \eqn{\alpha = 1 - \text{conf.level}}.
#'
#' The confidence limits are calculated as:
#'
#' \deqn{\text{Lower} = \max\left(\frac{x + z^2 / 2}{n + z^2} - z \sqrt{\frac{x}{n^2} \left(1 - \frac{x}{n}\right)}, 0\right),}
#' \deqn{\text{Upper} = \min\left(\frac{x + z^2 / 2}{n + z^2} + z \sqrt{\frac{x}{n^2} \left(1 - \frac{x}{n}\right)}, 1\right).}
#'
#' Special cases are handled explicitly:
#'
#' - If \eqn{x = 0}, the lower limit is 0, and the upper limit is
#' calculated as \eqn{(\alpha / 2)^{1/n}}.
#'
#' - If \eqn{x = n}, the upper limit is 1, and the lower limit is
#' calculated as \eqn{1 - (\alpha / 2)^{1/n}}.
#'
#' @return A vector with the lower and upper limits of the confidence interval.
#'
#' @examples
#' ci_p_wald_recentered(x =  0, n = 50, conf.level = 0.95)
#' ci_p_wald_recentered(x = 22, n = 50, conf.level = 0.95)
#' ci_p_wald_recentered(x = 50, n = 50, conf.level = 0.95)
#'
#' @export
#'
ci_p_wald_recentered <- function(x, n, conf.level = 0.95) {
  alpha <- 1 - conf.level
  c2 <- qnorm(1 - alpha / 2)

  if (x == 0) {
    lower <- 0
    upper <- (alpha / 2)^(1 / n)
  } else if (x == n) {
    lower <- 1 - (alpha / 2)^(1 / n)
    upper <- 1
  } else {
    lower <- max((x + c2^2 / 2) / (n + c2^2) - c2 * sqrt((x / n^2) * (1 - x / n)), 0)
    upper <- min((x + c2^2 / 2) / (n + c2^2) + c2 * sqrt((x / n^2) * (1 - x / n)), 1)
  }

  return(c(lower, upper))
}
ci_p_wald_recentered <- Vectorize(ci_p_wald_recentered)
#'
#'
#'
#' Recentered Wald Interval with Continuity Correction for binomial proportion
#'
#' @author David Esteban Cartagena Mejía, \email{dcartagena@unal.edu.co}
#'
#' @description
#' This function calculates the recentered Wald confidence interval with
#' continuity correction for a binomial proportion.
#' It adjusts the classical Wald interval by introducing a recentering term
#' and a continuity correction, improving accuracy for small sample sizes
#' and boundary cases.
#' The method is vectorized, allowing for evaluation of single
#' values or vectors.
#'
#' @param x A number or a vector with the number of successes.
#' @param n A number or a vector with the number of trials.
#' @param conf.level Confidence level for the returned confidence interval.
#' By default, it is 0.95.
#'
#' @references
#' Pires, Ana M., and Conceiçao Amado. "Interval estimators for a binomial
#' proportion: Comparison of twenty methods".
#' REVSTAT-Statistical Journal 6.2 (2008): 165-197.
#'
#' @seealso \link{ci_p}.
#'
#' @details
#' The recentered Wald interval with continuity correction adjusts the
#' classical Wald interval by incorporating a recentering term and a
#' continuity correction to account for the discreteness of the
#' binomial distribution.
#'
#' The critical value \eqn{z} is obtained from the standard normal
#' distribution for the specified confidence level:
#'
#' \deqn{z = \Phi^{-1}(1 - \alpha / 2),}
#'
#' where \eqn{\alpha = 1 - \text{conf.level}}.
#'
#' The confidence limits are calculated as:
#'
#' \deqn{\text{Lower} = \max\left(\frac{x + z^2 / 2}{n + z^2} - \left[z \sqrt{\frac{x}{n^2} \left(1 - \frac{x}{n}\right)} + \frac{1}{2n}\right], 0\right),}
#' \deqn{\text{Upper} = \min\left(\frac{x + z^2 / 2}{n + z^2} + \left[z \sqrt{\frac{x}{n^2} \left(1 - \frac{x}{n}\right)} + \frac{1}{2n}\right], 1\right).}
#'
#' Special cases are handled explicitly:
#'
#' - If \eqn{x = 0}, the lower limit is 0, and the upper limit is
#' calculated as \eqn{(\alpha / 2)^{1/n}}.
#'
#' - If \eqn{x = n}, the upper limit is 1, and the lower limit is
#' calculated as \eqn{1 - (\alpha / 2)^{1/n}}.
#'
#' These adjustments ensure that the confidence interval is valid and
#' well-behaved, even at the boundaries of the parameter space.
#'
#' @return A vector with the lower and upper limits of the confidence interval.
#'
#' @examples
#' ci_p_wald_recentered_cc(x =  0, n = 50, conf.level = 0.95)
#' ci_p_wald_recentered_cc(x = 25, n = 50, conf.level = 0.95)
#' ci_p_wald_recentered_cc(x = 50, n = 50, conf.level = 0.95)
#'
#' @export
#'
ci_p_wald_recentered_cc <- function(x, n, conf.level = 0.95) {
  alpha <- 1 - conf.level
  c1 <- qnorm(1 - alpha / 2)  # Valor crítico Z

  if (x == 0) {
    lower <- 0
    upper <- (alpha / 2)^(1 / n)
  } else if (x == n) {
    lower <- 1 - (alpha / 2)^(1 / n)
    upper <- 1
  } else {
    error_term <- c1 * sqrt((x / n^2) * (1 - (x / n))) + 1 / (2 * n)

    lower <- max((x + c1^2 / 2) / (n + c1^2) - error_term, 0)
    upper <- min((x + c1^2 / 2) / (n + c1^2) + error_term, 1)
  }

  return(c(lower,upper))
}
ci_p_wald_recentered_cc <- Vectorize(ci_p_wald_recentered_cc)
#'
#'
#'
#' Wald-t Confidence Interval for binomial proportion
#'
#' @author David Esteban Cartagena Mejía, \email{dcartagena@unal.edu.co}
#'
#' @description
#' This function calculates the Wald-t confidence interval for a Binomial
#' proportion using the method XIX proposed by Pan (2002).
#' It incorporates the use of a t-distribution with adjusted degrees of
#' freedom \eqn{\nu}, as specified in equation (2.8).
#' The function is vectorized, allowing for evaluation of single
#' values or vectors.
#'
#' @param x A number or a vector with the number of successes.
#' @param n A number or a vector with the number of trials.
#' @param conf.level Confidence level for the returned confidence interval.
#' By default, it is 0.95.
#'
#' @references
#' Pires, Ana M., and Conceiçao Amado. "Interval estimators for a binomial
#' proportion: Comparison of twenty methods".
#' REVSTAT-Statistical Journal 6.2 (2008): 165-197.
#'
#' @seealso \link{ci_p}.
#'
#' @details
#' The Wald-t confidence interval is a modification of the classical
#' Wald interval that uses the t-distribution instead of the normal
#' distribution for greater accuracy in small samples. The estimated
#' proportion \eqn{\hat{p}} is adjusted as:
#'
#' \deqn{\hat{p} = \frac{x + 2}{n + 4},}
#'
#' which reduces bias in the interval estimation. The variance
#'  \eqn{V(\hat{p}, n)} is given by:
#'
#' \deqn{V(\hat{p}, n) = \frac{\hat{p}(1 - \hat{p})}{n}.}
#'
#' The degrees of freedom \eqn{\nu} are calculated using equation (2.8):
#'
#' \deqn{\nu = \frac{2 V(\hat{p}, n)^2}{\Omega(\hat{p}, n)},}
#'
#' where \eqn{\Omega(\hat{p}, n)} is defined as:
#'
#' \deqn{\Omega(\hat{p}, n) = \frac{\hat{p} - \hat{p}^2}{n^3} + \frac{\hat{p} + (6n - 7)\hat{p}^2 + 4(n - 1)(n - 3)\hat{p}^3 - 2(n - 1)(2n - 3)\hat{p}^4}{n^5} - \frac{2(\hat{p} + (2n - 3)\hat{p}^2 - 2(n - 1)\hat{p}^3)}{n^4}.}
#'
#' The confidence interval is then calculated as:
#'
#' \deqn{\text{Lower} = \max\left(0, \hat{p} - t \cdot \sqrt{V(\hat{p}, n)}\right),}
#' \deqn{\text{Upper} = \min\left(1, \hat{p} + t \cdot \sqrt{V(\hat{p}, n)}\right),}
#'
#' where \eqn{t} is the critical value from the t-distribution
#' with \eqn{\nu} degrees of freedom.
#'
#' @return A vector with the lower and upper limits of the confidence interval.
#'
#' @examples
#' ci_p_wald_t(x =  0, n = 50, conf.level = 0.95)
#' ci_p_wald_t(x = 15, n = 50, conf.level = 0.95)
#' ci_p_wald_t(x = 50, n = 50, conf.level = 0.95)
#'
#' @export
#'
ci_p_wald_t <- function(x, n, conf.level = 0.95) {
  p_hat <- (x + 2) / (n + 4)
  Omega <- function(p, n) {
    term1 <- (p - p^2) / n^3
    term2 <- (p + (6 * n - 7) * p^2 + 4 * (n - 1) * (n - 3) * p^3 - 2 * (n - 1) * (2 * n - 3) * p^4) / n^5
    term3 <- (2 * (p + (2 * n - 3) * p^2 - 2 * (n - 1) * p^3)) / n^4
    return(term1 + term2 - term3)
  }
  V <- function(p, n) {
    return(p * (1 - p) / n)
  }
  nu <- (2 * V(p_hat, n)^2) / Omega(p_hat, n)

  # Cálculo del error estándar
  se <- sqrt(p_hat * (1 - p_hat) / n)
  t <- qt(conf.level + (1 - conf.level) / 2, df = nu)
  lower <- max(p_hat - t * se,0)
  upper <- min(p_hat + t * se,1)

  return(c(lower, upper))
}
ci_p_wald_t <- Vectorize(ci_p_wald_t)
#'
#'
#'
#' Wald with Blyth and Still approximation Binomial Score confidence
#' interval for binomial proportion
#'
#' @author David Esteban Cartagena Mejía, \email{dcartagena@unal.edu.co}
#'
#' @description
#' This function calculates the Wald Binomial Score confidence interval for
#' a binomial proportion.
#' It is vectorized, allowing the evaluation of single values or vectors.
#'
#' @param x A number or a vector with the number of successes.
#' @param n A number or a vector with the number of trials.
#' @param conf.level Confidence level for the returned confidence interval.
#' By default, it is 0.95.
#'
#' @references
#' Blyth, C.R. and Still, H.A. (1983). Binomial confidence intervals,
#' Journal of the American Statistical Association, 78, 108–116.
#'
#' @seealso \link{ci_p}.
#'
#' @details
#' The Wald Binomial Score confidence interval is an adjusted version of
#' the Wald interval,
#' designed to improve accuracy in small samples and near the boundaries
#' of the parameter space.
#'
#' The bounds for the interval confidence are:
#'
#' \eqn{\hat{p} \mp \left[ \frac{Z_{\alpha/2} \sqrt{\hat{p} (1-\hat{p})}}{\sqrt{n - Z_{\alpha/2}^2 - \frac{2Z_{\alpha/2}}{\sqrt{n}} - \frac{1}{n}}} + \frac{1}{2n} \right]}
#'
#' This interval is particularly useful when \eqn{n} is small or when
#' the proportion \eqn{\hat{p}} is close to 0 or 1.
#'
#' @return A vector with the lower and upper limits of the confidence interval.
#'
#' @examples
#' ci_p_wald_bs(x = 15, n = 50, conf.level = 0.95)
#'
#' @export
#'
ci_p_wald_bs <- function(x, n, conf.level = 0.95) {
  p_hat <- x / n
  q_hat <- 1 - p_hat
  z <- qnorm(1 - (1 - conf.level) / 2)
  # Adjustment factor
  wn_z <- sqrt(n - z^2 - (2 * z / sqrt(n)) - (1 / n))
  term1 <- 1 / wn_z
  term2 <- z * sqrt((p_hat * q_hat) + 1 / (2 * n))
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
      lower <-  p_hat - z * sqrt(p_hat * q_hat / n) - 1 / (2 * n)
      upper <-  p_hat + z * sqrt(p_hat * q_hat / n) + 1 / (2 * n)
    }
  }
  return(c(lower, upper))
}
ci_p_wald_bs <- Vectorize(ci_p_wald_bs)
