#' suma
#'
#' Esta funcion sirve para obtener la suma de dos numeros reales.
#'
#' @param x A number.
#' @param y A number.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' suma(1, 1)
#' suma(10, 1)
#' @export
suma <- function(x, y) {
  res <- x + y
  class(res) <- 'sumita'
  res
}
