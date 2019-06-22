#' Crea un histograma rosado
#'
#' Sirve para construir un histograma de una muestra de una poblacion normal.
#'
#' @param n tamano de muestra.
#' @param media es la media de la poblacion.
#' @param varianza es la varianza de la poblacion.
#'
#' @examples
#' mi_hist(n=10, media=15, varianza=3)
#' mi_hist(n=100, media=45, varianza=12)
#'
#' @importFrom stats rnorm
#' @importFrom graphics hist
#' @export
mi_hist <- function(n, media, varianza) {
  x <- rnorm(n=n, mean=media, sd=sqrt(varianza))
  hist(x, col="pink", main="Mi histograma")
}
