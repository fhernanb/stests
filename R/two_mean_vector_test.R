#' prueba de hipotesis para
#'
#' This function implements the test for \eqn{H0: \mu 1 = \mu 2} versus \eqn{H1: \mu1} not = \eqn{\mu2}. Muchachos, por favor revisar si los argumentos estan bien escritos.
#'
#' @param xbar1 a vector with the sample mean from population 1.
#' @param Sigma1 a matrix with sample variances and covariances from population 1.
#' @param n1 sample size 1.
#' @param xbar2 a vector with the sample mean from population 2.
#' @param Sigma2 a matrix with sample variances and covariances from population 2.
#' @param n2 sample size 2.
#' @param delta0 que es esto
#' @param alpha yo creo que este parametro lo debemos quita, que opinan?
#' @param method tipo de prueba de hipotesis a usar. Las opciones disponibles son T2, bla bla. Por defecto es T2.
#'
#' @return por ahora retorna una lista de cosas, debemos hacer que la salida sea htest y que entregue lo que se necesita.
#'
#' @author Freddy Hernandez
#' @examples
#' # Ejemplo 3.5.2 de Diaz y Morales (2005).
#' # Valentina y Jean Paul, debemos buscar mejor un ejemplo
#' # de un libro en ingles.
#' x1_barra <- c(15.97, 15.91, 27.19, 22.75)
#' x2_barra <- c(12.34, 13.91, 16.59, 21.94)
#' s1 <- matrix(c(5.192, 4.545, 6.522, 5.25, 4.545, 13.184, 6.76, 6.266, 6.522,
#'                6.76, 28.67, 14.468, 5.25, 6.266, 14.468, 16.645), ncol = 4)
#' s2 <- matrix(c(9.136, 7.549, 5.531, 4.151, 7.549, 18.6, 5.446, 5.446, 5.531,
#'                5.446, 13.55, 13.55, 4.151, 5.446, 13.55, 28), ncol = 4)
#'
#' two_mean_vector_test(xbar1=x1_barra, Sigma1=s1, n1=32,
#'                      xbar2=x2_barra, Sigma2=s2, n2=32)
#'
#' @importFrom stats pf
#' @export
two_mean_vector_test <- function(xbar1, Sigma1, n1,
                                 xbar2, Sigma2, n2,
                                 delta0=NULL, alpha=0.05,
                                 method="T2") {

  if (! identical(dim(Sigma1), dim(Sigma2)))
    stop("The matrices Sigma1 and Sigma2 do not have same dimension")

  method <- match.arg(arg=method,
                      choices=c("T2", "james"))

  # To generate the code for evaluating, without using cases
  my_code <- paste0("two_mean_vector_test_", method,
                    "(xbar1, Sigma1, n1, xbar2, Sigma2, n2, delta0=delta0, alpha=alpha)")

  # To obtain the result
  result <- eval(parse(text=my_code))
  return(result)
}
#' @importFrom stats pf qf
two_mean_vector_test_T2 <- function(xbar1, Sigma1, n1,
                                    xbar2, Sigma2, n2,
                                    delta0=NULL, alpha=0.05) {
  p <- ncol(Sigma1)
  v <- n1 + n2 - 2
  if (is.null(delta0)) delta0 <- matrix(0, ncol=1, nrow=p)
  xbar1 <- matrix(xbar1, ncol=1)
  xbar2 <- matrix(xbar2, ncol=1)

  sp <- ((n1 - 1) * Sigma1 + (n2 - 1) * Sigma2) / (n1 + n2 - 2)
  T2 <- t(xbar1-xbar2-delta0) %*% solve(sp * (1/n1+1/n2)) %*% (xbar1-xbar2-delta0)
  T2 <- as.numeric(T2)
  p.value <- pf(q=T2 * (v-p+1) / (v*p), df1=p, df2=v-p+1, lower.tail=FALSE)
  critic.value <- qf(p=alpha, df1=p, df2=v-p+1, lower.tail=F) * (v*p) / (v-p+1)
  return(list(sp=sp, T2=T2, p.value=p.value,
              critic.value=critic.value))
}
#' @importFrom stats pf qf
two_mean_vector_test_james <- function(xbar1, Sigma1, n1,
                                       xbar2, Sigma2, n2,
                                       delta0=NULL, alpha=0.05) {
  p <- ncol(Sigma1)
  xbar1 <- matrix(xbar1, ncol=1)
  xbar2 <- matrix(xbar2, ncol=1)

  S1 <- Sigma1/n1 # Represents S1 tilde
  S2 <- Sigma2/n2 # Represents S2 tilde
  S <- S1 + S2    # Represents S tilde

  T2 <- t(xbar1-xbar2) %*% solve(S) %*% (xbar1-xbar2)
  T2 <- as.numeric(T2)

  # To obtain delta
  Sinv <- solve(S)
  b1 <- Sinv %*% S1
  b2 <- Sinv %*% S2
  trb1 <- sum(diag(b1))
  trb2 <- sum(diag(b2))
  A <- 1 + (trb1^2/(n1 - 1) + trb2^2/(n2 - 1))/(2 * p)
  B <- (sum(b1^2)/(n1 - 1) + sum(b2^2)/(n2 - 1) + 0.5 *
          trb1^2/(n1 - 1) + 0.5 * trb2^2/(n2 - 1))/(p * (p + 2))
  x2 <- qchisq(1 - alpha, p)
  delta <- (A + B * x2)
  twoha <- x2 * delta
  p.value <- pchisq(q=T2/delta, df=p, lower.tail=FALSE)

  return(list(T2=T2, p.value=p.value))
}
