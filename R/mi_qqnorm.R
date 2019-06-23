#' Crea un grafico qqnorm
#'
#' Sirve para explorar la normalidad de una muestra mediante un gráfico cuantil cuantil
#'
#' @param n tamano de muestra
#' @param media media de la poblacion
#' @param varianza varianza de la poblacion
#'
#' @examples
#' mi_qqnorm(n = 30, media = 10, varianza = 4)
#' mi_qqnorm(n = 100, media = 25, varianza = 12)
#'
#' @importFrom stats rnorm
#' @importFrom stats qqnorm
#' @importFrom stats qqline
#' @export
mi_qqnorm <- function(n, media, varianza){
  x <- rnorm(n = n, mean = media, sd = sqrt(varianza))
  qqnorm(x)
  qqline(y = x, col = 'blue', lwd = 2, lty = 2)
}
