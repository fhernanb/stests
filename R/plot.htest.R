#' To plot p-value area
#'
#' This is a generic function to depict the p-value area for a object with class \code{htest}.
#'
#' @param x Object of class \dQuote{htest}.
#' @param col color for the observed statistic.
#' @param shade.col color for the shaded area.
#' @param cex A numerical value giving the amount by which plotting the p-value.
#' @param from inicio
#' @param to final
#' @param ... Other plotting parameters to affect the plot.
#' @return None. Function produces a plot.
#' @author Freddy Hernandez
#' @seealso \code{\link{z.test}}
#' @keywords plot.htest
#' @examples
#'
#' # Example
#' # H0: mu = 170
#' # Ha: mu != 170 with sigma2=25
#' x <- rnorm(n=80, mean=171, sd=5)
#' res <- z.test(x=x, mu=170, sigma2=25, alternative='two.sided')
#' res
#' plot(res, col='blue', shade.col='tomato')
#'
#' @importFrom graphics mtext title
#' @export
plot.htest <- function(x, col='red', shade.col='red', cex=0.8,
                       from, to, ...) {

  # Z test
  if (x$method %in% c('One Sample z-test',
                      'Z test for mean')) {
    if (x$alternative == 'less')
      shade.norm(x=x$statistic, tail='lower',
                 las=1, shade.col=shade.col, cex=cex)
    if (x$alternative == 'greater')
      shade.norm(x=x$statistic, tail='upper',
                 las=1, shade.col=shade.col, cex=cex)
    if (x$alternative == 'two.sided')
      shade.norm(x=x$statistic, tail="two",
                 las=1, shade.col=shade.col, cex=cex)
    title(main='Shaded area corresponds to p-value')
    mtext(text=round(x$statistic, digits=4),
          side=1, at=x$statistic,
          col=col, cex=cex, adj=0.5)
  }

  # for 1 prop test
  if (x$method %in% c('1-sample proportions test with continuity correction',
                      '1-sample proportions test without continuity correction')) {

    x$statistic <- sign(x$estimate - x$null.value) * sqrt(x$statistic)

    if (x$alternative == 'less')
      shade.norm(x=x$statistic, tail='lower',
                 las=1, shade.col=shade.col, cex=cex)
    if (x$alternative == 'greater')
      shade.norm(x=x$statistic, tail='upper',
                 las=1, shade.col=shade.col, cex=cex)
    if (x$alternative == 'two.sided')
      shade.norm(x=x$statistic, tail="two",
                 las=1, shade.col=shade.col, cex=cex)
  }

  # for 2 prop test
  if (x$method %in% c('2-sample test for equality of proportions with continuity correction',
                      '2-sample test for equality of proportions without continuity correction')) {

    x$statistic <- sign(x$estimate[1] - x$estimate[2]) * sqrt(x$statistic)

    if (x$alternative == 'less')
      shade.norm(x=x$statistic, tail='lower',
                 las=1, shade.col=shade.col, cex=cex)
    if (x$alternative == 'greater')
      shade.norm(x=x$statistic, tail='upper',
                 las=1, shade.col=shade.col, cex=cex)
    if (x$alternative == 'two.sided')
      shade.norm(x=x$statistic, tail="two",
                 las=1, shade.col=shade.col, cex=cex)
  }

  # t test
  if (x$method %in% c('One Sample t-test',
                      'Welch Two Sample t-test',
                      ' Two Sample t-test',
                      'Two Sample t-test',
                      'Paired t-test')) {
    if (x$alternative == 'less')
      shade.t(x=x$statistic, tail='lower',
              nu=x$parameter,
              las=1, shade.col=shade.col, cex=cex)
    if (x$alternative == 'greater')
      shade.t(x=x$statistic, tail='upper',
              nu=x$parameter,
              las=1, shade.col=shade.col, cex=cex)
    if (x$alternative == 'two.sided')
      shade.t(x=x$statistic, tail="two",
              nu=x$parameter,
              las=1, shade.col=shade.col, cex=cex)
    title(main='Shaded area corresponds to p-value')
    mtext(text=round(x$statistic, digits=4),
          side=1, at=x$statistic,
          col=col, cex=cex, adj=0.5)
  }

  # Chi squared test
  if (x$method %in% c('X-squared test for variance')) {
    if (x$alternative == 'less')
      shade.chi(x=x$statistic, nu=x$parameter,
                tail='lower', las=1,
                shade.col=shade.col, cex=cex)
    if (x$alternative == 'greater')
      shade.chi(x=x$statistic, nu=x$parameter,
                tail='upper', las=1,
                shade.col=shade.col, cex=cex)
    if (x$alternative == 'two.sided')
      shade.chi(nu=x$parameter, tail="two",
                las=1, shade.col=shade.col, cex=cex,
                prob.to.each.tail=x$p.value/2)
    title(main='Shaded area corresponds to p-value')
    mtext(text=round(x$statistic, digits=4),
          side=1, at=x$statistic,
          col=col, cex=cex, adj=0.5)
  }

  # F test
  if (x$method %in% c('F test to compare two variances')) {
    if (x$alternative == 'less')
      shade.F(x=x$statistic,
              nu1=x$parameter[1],
              nu2=x$parameter[2],
              tail='lower', las=1,
              shade.col=shade.col, cex=cex, ...)
    if (x$alternative == 'greater')
      shade.F(x=x$statistic,
              nu1=x$parameter[1],
              nu2=x$parameter[2],
              tail='upper', las=1,
              shade.col=shade.col, cex=cex, ...)
    if (x$alternative == 'two.sided')
      shade.F(nu1=x$parameter[1],
              nu2=x$parameter[2],
              tail="two",
              shade.col=shade.col, cex=cex,
              prob.to.each.tail=x$p.value/2, las=1, ...)
    title(main='Shaded area corresponds to p-value')
    mtext(text=round(x$statistic, digits=4),
          side=1, at=x$statistic,
          col=col, cex=cex, adj=0.5)
  }

  # X2 test for mean vector
  if (x$method %in% c('X2 test for mean vector')) {
    shade.chi(x=x$statistic, nu=x$parameter,
              tail='upper', las=1,
              shade.col=shade.col, cex=cex)
    title(main='Shaded area corresponds to p-value')
    mtext(text=round(x$statistic, digits=4),
          side=1, at=x$statistic,
          col=col, cex=cex, adj=0.5)
  }

  # T2 test for mean vector
  if (x$method %in% c('T2 test for mean vector')) {
    #shade.F(x=x$statistic[2],
    #        nu1=x$parameter[1],
    #        nu2=x$parameter[2],
    #        tail='upper', las=1,
    #        shade.col=shade.col, cex=cex)

    shadow.dist(dist='df',
                param=list(df1=x$parameter[1], df2=x$parameter[2]),
                a=x$statistic[2], type='upper',
                from=from, to=to, col.shadow=shade.col, ...)

    title(main='Shaded area corresponds to p-value')
    mtext(text=round(x$statistic[2], digits=4),
          side=1, at=x$statistic[2],
          col=col, cex=cex, adj=0.5)

    df1 <- x$parameter[1]
    df2 <- x$parameter[2]
    legend("top",
           legend=bquote(paste(italic(X)," ~ ",italic(F),"(",.(df1),",",.(df2),")", sep = "")),
           bty="n", adj=0.5)
  }

}

