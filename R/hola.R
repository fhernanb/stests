#' Histograma rosado
#'
#' Sirve para crear un histograma de una muestra aleatoria de una normal.
#'
#' @param media es la media de la poblacion.
#' @param varianza la varianza de la poblacion.
#' @param n sample size.
#'
#' @author Freddy Hernandez
#' @examples
#' hola(media=170, varianza=10, n=59)
#' hola(media=-50, varianza=1, n=23)
#'
#' @importFrom stats rnorm
#' @importFrom graphics hist
#' @export
#'
hola <- function(media, varianza, n) {
  x <- rnorm(n=n, mean=media, sd=sqrt(varianza))
  hist(x, main="Mi histograma", col="pink")
}
