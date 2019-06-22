#' Grafica una densidad
#'
#' sirve para graficar la densidad de una muestra de una normal.
#'
#' @param n tamano de muestra.
#' @param media media de la poblacion.
#' @param varianza varianza de la poblacion
#'
#' @examples
#' mi_densidad(n=10,media=0,varianza = 1)
#' mi_densidad(n=100,media = 10,varianza= 4)
#'
#' @importFrom stats rnorm density
#' @importFrom graphics plot
#' @export
mi_densidad <- function(n,media,varianza){
  x <- rnorm(n,media,varianza)
  plot(density(x), col='red', main = 'Mi densidad')
}
