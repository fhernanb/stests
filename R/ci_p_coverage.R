#' Coverage for interval confidence p
#'
#' @author David Esteban Cartagena Mejía, \email{dcartagena@unal.edu.co}
#'
#' @description
#' This function calculates the true coverage for any interval confidence for p.
#'
#' @param n number of trials.
#' @param p true value for p.
#' @param conf.level nominal confidence level for the returned confidence interval. By default is 0.95.
#' @param intervalType type of confidence interval, possible choices are listed in \link{ci_p}.
#'
#' @seealso \link{ci_p}.
#'
#' @references
#' Park, H., & Leemis, L. M. (2019). Ensemble confidence intervals for binomial proportions. Statistics in Medicine, 38(18), 3460-3475.
#'
#' @details
#' This function was inspired by the binomTestCoveragePlot() function from conf package and Park & Leemis (2019).
#'
#' @return A dataframe with Method, n, p and true coverage.
#'
#' @examples
#' ci_p_coverage(n=10, p=0.45, intervalType="wald")
#' ci_p_coverage(n=10, p=0.45, intervalType="wilson")
#' ci_p_coverage(n=10, p=c(0.10, 0.25, 0.40), intervalType="wilson")
#'
#' @export
#' @importFrom stats aggregate
ci_p_coverage <- function(n, p, conf.level=0.95, intervalType="wald") {
  res <- ci_p_coverage_single(n, p, conf.level, intervalType)
  res <- as.data.frame(t(res))
  return(res)
}
#'
#'
#'
ci_p_coverage_single <- function(n, p, conf.level=0.95, intervalType="wald") {
  # To select the interval type
  intervalType <- paste0("ci_p_", intervalType)

  x <- rev(0:n) # Generar valores de x
  n_rep <- rep(n, length(x)) # Repetir n para cada valor de x

  # Evaluar el metodo para calcular los intervalos de confianza
  ci <- eval(parse(text=paste0(intervalType, "(x, n_rep, conf.level=", conf.level, ")")))

  # Crear un marco de datos con los resultados del intervalo
  data <- data.frame(Method=rep(intervalType, length(x)),
                     x=x, n=n_rep,lower=ci[1, ], upper=ci[2, ])

  # Combinar con p mediante merge
  results <- merge(data, data.frame(p=p))

  # Calcular coverage
  results$coverage <- with(results, (p >= lower & p <= upper) * dbinom(x, n[1], p))

  coverage <- with(results, (p >= lower & p <= upper) * dbinom(x, n[1], p))
  coverage <- sum(coverage)

  # Agregar resultados
  results <- aggregate(results["coverage"], results[c("Method", "n", "p")], sum)

  res <- data.frame(intervalType, n, p, coverage)
  return(res)
}
ci_p_coverage_single <- Vectorize(ci_p_coverage_single)
#'
#'
#'
#'
#' Plot for confidence interval p - coverage vs p for fixed n
#'
#' @author David Esteban Cartagena Mejía, \email{dcartagena@unal.edu.co}
#'
#' @description
#' This function plots the coverage for any confidence interval for p.
#'
#' @param n number of trials.
#' @param conf.level nominal confidence level for the returned confidence interval. By default is 0.95.
#' @param intervalType type of confidence interval, possible choices are listed in \link{ci_p}.
#' @param seq_p sequence for the values of \eqn{n}. By default is \code{seq(from=0.01, to=0.99, length.out=50)}.
#' @param plot logical value to obtain the plot, TRUE by default.
#' @param col color for the coverage curve.
#' @param linecolor color for the line representing the conf.level.
#' @param ... further arguments and graphical parameters passed to plot function.
#'
#' @seealso \link{ci_p}.
#'
#' @references
#' Park, H., & Leemis, L. M. (2019). Ensemble confidence intervals for binomial proportions. Statistics in Medicine, 38(18), 3460-3475.
#'
#' @details
#' This function was inspired by the binomTestCoveragePlot() function from conf package and Park & Leemis (2019).
#'
#' @return A dataframe with Method, n, p and true coverage and the plot.
#'
#' @example examples/examples_ci_p_coverage_plot.R
#'
#' @export
#' @importFrom graphics abline grid
ci_p_coverage_plot <- function(n, conf.level=0.95,
                               intervalType="wald", plot=TRUE,
                               seq_p=seq(from=0.01, to=0.99, length.out=50),
                               col="deepskyblue2", linecolor="tomato", ...) {

  # Calcular la cobertura utilizando la funcion ci_p_coverage
  cover <- ci_p_coverage(n=n,
                         p=seq_p,
                         conf.level=conf.level,
                         intervalType=intervalType)

  if (plot) {
    plot(x=seq_p,
         y=cover$coverage,
         type="l", col=col, lwd=2,
         xlab="p", ylab="Coverage",
         main=intervalType, ...)
    abline(h=conf.level, col=linecolor)
    grid()
  }

  return(cover)
}
#'
#'
#'
#'
#' Plot for confidence interval p - coverage vs n for fixed p
#'
#' @author Laura Juliana Insignares Montes, \email{linsignares@unal.edu.co}
#'
#' @description
#' This function plots the coverage for any confidence interval for p.
#'
#' @param p true proportion \eqn{p}.
#' @param conf.level nominal confidence level for the returned confidence interval. By default is 0.95.
#' @param intervalType type of confidence interval, possible choices are listed in \link{ci_p}.
#' @param seq_n sequence with the values of sample size \eqn{n}. By default is the sequence 10, 20, 30, ..., 480, 490, 500.
#' @param plot logical value to obtain the plot, TRUE by default.
#' @param col color for the coverage curve.
#' @param linecolor color for the line representing the conf.level.
#' @param ... further arguments and graphical parameters passed to plot function.
#'
#' @seealso \link{ci_p}.
#'
#' @references
#' Park, H., & Leemis, L. M. (2019). Ensemble confidence intervals for binomial proportions. Statistics in Medicine, 38(18), 3460-3475.
#'
#' @details
#' This function was inspired by the binomTestCoveragePlot() function from conf package and Park & Leemis (2019).
#'
#' @return A dataframe with Method, n, p and true coverage and the plot.
#'
#' @example examples/examples_ci_p_coverage_plot2.R
#'
#' @export
#' @importFrom graphics abline grid
ci_p_coverage_plot2 <- function(p=0.5,
                                seq_n=seq(from=10, to=500, by=10),
                                conf.level=0.95,
                                intervalType="wald",
                                plot=TRUE,
                                col="deepskyblue2", linecolor="tomato", ...) {

  # Calcular la cobertura utilizando la funcion ci_p_coverage
  cover <- ci_p_coverage(n=seq_n,
                       p=p,
                       conf.level=conf.level,
                       intervalType=intervalType)

  # Graficar la cobertura en funcion del tamano de la muestra
  if (plot) {
    plot(x=seq_n,
         y=cover$coverage,
         type="l", col=col, lwd=2, las=1,
         xlab="n", ylab="Coverage",
         main=intervalType, ...)
    abline(h=conf.level, col=linecolor)
    grid()
  }

  return(cover)
}
