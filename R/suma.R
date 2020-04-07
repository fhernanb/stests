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
#'
#' Plot components from ETS model
#'
#' Produces a plot of the level, slope and seasonal components from an ETS
#' model.
#'
#' \code{plot} will produce an equivalent plot as a ggplot object.
#'
#' @param x Object of class \dQuote{sumita}.
#' @param col color for the plot.
#' @param las value for the plot.
#' @param ... Other plotting parameters to affect the plot.
#' @return None. Function produces a plot
#' @author Rob J Hyndman & Mitchell O'Hara-Wild
#' @seealso \code{\link{suma}}
#' @keywords hplot
#' @examples
#'
#' fit <- suma(10, 5)
#' plot(fit)
#'
#' @importFrom graphics plot
#' @export
plot.sumita <- function(x, col='red', las=1, ...) {
  plot(x=1:5, y=5:1, col=col, las=las)
}
