#' Shading functions for interpretation of pdf probabilities
#'
#' Creates plots with lower, upper, two-tailed, and middle of the distribution shading for popular pdfs.
#'
#' @param x A quantile, i.e. \eqn{X = x}, or if \code{tail = "two.custom"} ins \code{shade.norm}, a two element vector specifying the upper bound of the lower tail and the lower bound of the upper tail.
#' @param from To be used with \code{tail = "middle"}; the value \emph{a} in \eqn{P(a < X < b)}.
#' @param to To be used with \code{tail = "middle"}; the value \emph{b} in \eqn{P(a < X < b)}.
#' @param sigma Standard deviation for the nomral distribution.
#' @param mu Mean of the normal distribution.
#' @param tail One of four possibilities: \code{"lower"} provides lower tail shading, \code{"upper"} provides upper tail shading, \code{"two"} provides two tail shading, and \code{"middle"} provide shading in the middle of the pdf, between \code{"from"} and \code{"to"}.  The additional option \code{"two.custom"} is allowed for \code{shade.norm}. This allows calculation of asymmetric two tailed probabilities.  It requires that the argument \code{x} is a two element vector with elements denoting the upper bound of the lower tail and the lower bound of the upper tail.  For discrete pdfs (binomial and Poisson) the possibility \code{"X=x"} is also allowed, and will be equivalent to the density. Two tailed probability is not implemented for \code{shade.poi}.
#' @param show.p Logical; indicating whether probabilities are to be shown.
#' @param show.d Logical; indicating whether densities are to be shown.
#' @param show.dist Logical; indicating whether parameters for the distribution are to be shown.
#' @param nu Degrees of freedom.
#' @param nu1 Numerator degrees of freedom for the \emph{F}-distribution.
#' @param nu2 Denominator degrees of freedom for the \emph{F}-distribution.
#' @param prob.to.each.tail Probability to be apportioned to each tail in the \emph{F} and Chi-square distributions if \code{tail = "two"}.
#' @param digits Number of digits to be reported in probsabilities and densities.
#' @param legend.cex Character expansion for legends in plots.
#' @param shade.col Color of probability shading.
#' @param ... Arguments to be passed to methods, such as graphical parameters (see par).
#'
#' @author Ken Aho
#-------------------------- normal distribution--------------------------------#
#' @importFrom graphics legend curve polygon
#' @importFrom stats pnorm dnorm df
#' @export
shade.norm <- function(x=NULL, from=NULL, to=NULL, sigma=1, mu=0,
                       tail="lower", show.p=TRUE, show.d=FALSE,
                       show.dist=TRUE, digits=5, legend.cex=.9,
                       shade.col="gray", ...) {
  xv<-seq(mu-4*sigma,mu+4*sigma,sigma/1000)
  yv<-dnorm(xv,mean=mu,sd=sigma)
  curve(dnorm(x,mu,sigma),from=mu-4*sigma,to=mu+4*sigma,ylab=expression(paste(italic(f),"(",italic(x),")", sep = "")),xlab = expression(italic(x)),...)

  if(tail=="lower"){
    polygon(c(xv[xv<=x],x),c(yv[xv<=x],yv[xv==mu-4*sigma]),col=shade.col)
    p<-round(pnorm(x,mu,sigma,lower.tail=TRUE),digits)
    d<-round(dnorm(x,mean=mu,sd=sigma),digits)
    if(show.p==TRUE&show.d==FALSE)legend("topright",legend=bquote(paste(italic(P),"(",italic(X)<=.(x),") = ",.(p), sep = "")),bty="n",cex=legend.cex)
    if(show.d==TRUE&show.p==FALSE)legend("topright",legend=bquote(paste("",italic(f),"(",.(x),") = ",.(d), sep = "")),bty="n",cex=legend.cex)
    if(show.d==TRUE&show.p==TRUE)legend("topright",legend=bquote(paste(italic(P),"(",italic(X)<=.(x),") = ",.(p),",  ",italic(f),"(",.(x),") = ",.(d), sep = "")),bty="n",cex=legend.cex)
    if(show.dist==TRUE)legend("topleft",legend=bquote(paste(italic(X)," ~ ",italic(N),"(",.(mu)," , ",.(sigma^2),")", sep = "")),bty="n",cex=legend.cex)
  }

  if(tail=="upper"){
    polygon(c(x,xv[xv>=x]),c(yv[xv==mu+4*sigma],yv[xv>=x]),col=shade.col)
    p<-round(pnorm(x,mu,sigma,lower.tail=FALSE),digits)
    d<-round(dnorm(x,mean=mu,sd=sigma),digits)
    if(show.p==TRUE&show.d==FALSE)legend("topright",legend=bquote(paste(italic(P),"(",italic(X)>=.(x),") = ",.(p), sep = "")),bty="n",cex=legend.cex)
    if(show.d==TRUE&show.p==FALSE)legend("topright",legend=bquote(paste("",italic(f),"(",.(x),") = ",.(d), sep = "")),bty="n",cex=legend.cex)
    if(show.d==TRUE&show.p==TRUE)legend("topright",legend=bquote(paste(italic(P),"(",italic(X)>=.(x),") = ",.(p),",  ",italic(f),"(",.(x),") = ",.(d), sep = "")),bty="n",cex=legend.cex)
    if(show.dist==TRUE)legend("topleft",legend=bquote(paste(italic(X)," ~ ",italic(N),"(",.(mu)," , ",.(sigma^2),")", sep = "")),bty="n",cex=legend.cex)
  }

  if(tail=="two"){
    polygon(c(xv[xv<=-abs(x)],-abs(x)),c(yv[xv<=-abs(x)],yv[xv==mu-4*sigma]),col=shade.col)
    polygon(c(abs(x),xv[xv>=abs(x)]),c(yv[xv==mu+4*sigma],yv[xv>=abs(x)]),col=shade.col)
    p<-round(2*pnorm(abs(x),mu,sigma,lower.tail=FALSE),digits)
    if(show.p==TRUE)legend("topright",bty="n",cex=legend.cex,legend=bquote(paste(2%*%italic(P),"(", italic(X)>="|",.(x),"|) = ",.(p), sep = "")))
    if(show.dist==TRUE)legend("topleft",legend=bquote(paste(italic(X)," ~ ",italic(N),"(",.(mu)," , ",.(sigma^2),")", sep = "")),bty="n",cex=legend.cex)
  }

  if(tail=="two.custom"){
    polygon(c(xv[xv<=x[1]],x[1]),c(yv[xv<=x[1]],yv[xv==mu-4*sigma]),col=shade.col)
    polygon(c(x[2],xv[xv>=x[2]]),c(yv[xv==mu+4*sigma],yv[xv>=x[2]]),col=shade.col)
    p<-round(pnorm(x[1],mu,sigma),digits) + round(pnorm(x[2],mu,sigma, lower.tail = FALSE),digits)
    if(show.p==TRUE)legend("topright",bty="n",cex=legend.cex,legend=bquote(paste(italic(P),"(",.(x[1]) >= "", italic(X)>=.(x[2]),") = ",.(p), sep = "")))
    if(show.dist==TRUE)legend("topleft",legend=bquote(paste(italic(X)," ~ ",italic(N),"(",.(mu)," , ",.(sigma^2),")", sep = "")),bty="n",cex=legend.cex)
  }

  if(tail=="middle"){
    polygon(c(xv[xv<=mu+4*sigma],mu+4*sigma),c(yv[xv<=mu+4*sigma],yv[xv==mu-4*sigma]),col=shade.col)
    polygon(c(xv[xv<=from],from),c(yv[xv<=from],yv[xv==mu-4*sigma]),col="white")
    polygon(c(to,xv[xv>=to]),c(yv[xv==mu+4*sigma],yv[xv>=to]),col="white")
    p<-round(pnorm(to,mu,sigma)-pnorm(from,mu,sigma),digits)
    if(show.p==TRUE)legend("topright",legend=bquote(paste(italic(P),"(",.(from)<="", italic(X)<=.(to),") = ",.(p), sep = "")),bty="n",cex=legend.cex)
    if(show.dist==TRUE)legend("topleft",legend=bquote(paste(italic(X)," ~ ",italic(N),"(",.(mu)," , ",.(sigma^2),")", sep = "")),bty="n",cex=legend.cex)
  }
}
#-------------------------- t-distribution --------------------------------#
#' @importFrom stats dt pt qt
#' @export
#' @rdname shade.norm
shade.t <- function(x=NULL, from=NULL, to=NULL, nu=3, tail="lower",
                    show.p=TRUE, show.d=FALSE, show.dist=TRUE,
                    digits=5, legend.cex=.9, shade.col="gray", ...){
  sigma<-qt(.975,nu)
  xv<-seq(-4*sigma,4*sigma,sigma/1000)
  yv<-dt(xv,df=nu)
  curve(dt(x,nu),from=-4*sigma,to=4*sigma,xlab=expression(italic(x)),ylab=expression(paste(italic(f),"(",italic(x),")", sep = "")),...)

  if(tail=="lower"){
    polygon(c(xv[xv<=x],x),c(yv[xv<=x],yv[xv==-4*sigma]),col=shade.col)
    p<-round(pt(x,nu,lower.tail=TRUE),digits)
    d<-round(dt(x,nu),digits)
    if(show.p==TRUE&show.d==FALSE)legend("topright",legend=bquote(paste(italic(P),"(",italic(X)<=.(x),") = ",.(p), sep = "")),bty="n",cex=legend.cex)
    if(show.d==TRUE&show.p==FALSE)legend("topright",legend=bquote(paste("",italic(f),"(",.(x),") = ",.(d), sep = "")),bty="n",cex=legend.cex)
    if(show.d==TRUE&show.p==TRUE)legend("topright",legend=bquote(paste(italic(P),"(",italic(X)<=.(x),") = ",.(p),",  ",italic(f),"(",.(x),") = ",.(d), sep = "")),bty="n",cex=legend.cex)
    if(show.dist==TRUE)legend("topleft",legend=bquote(paste(italic(X)," ~ ",italic(t),"(",.(nu),")", sep = "")),bty="n",cex=legend.cex)
  }

  if(tail=="upper"){
    polygon(c(x,xv[xv>=x]),c(yv[xv==4*sigma],yv[xv>=x]),col=shade.col)
    p<-round(pt(x,nu,lower.tail=FALSE),digits)
    d<-round(dt(x,nu),digits)
    if(show.p==TRUE&show.d==FALSE)legend("topright",legend=bquote(paste(italic(P),"(",italic(X)>=.(x),") = ",.(p), sep = "")),bty="n",cex=legend.cex)
    if(show.d==TRUE&show.p==FALSE)legend("topright",legend=bquote(paste("",italic(f),"(",.(x),") = ",.(d), sep = "")),bty="n",cex=legend.cex)
    if(show.d==TRUE&show.p==TRUE)legend("topright",legend=bquote(paste(italic(P),"(",italic(X)>=.(x),") = ",.(p),",  ",italic(f),"(",.(x),") = ",.(d), sep = "")),bty="n",cex=legend.cex)
    if(show.dist==TRUE)legend("topleft",legend=bquote(paste(italic(X)," ~ ",italic(t),"(",.(nu),")", sep = "")),bty="n",cex=legend.cex)
  }

  if(tail=="two"){
    polygon(c(xv[xv<=-abs(x)],-abs(x)),c(yv[xv<=-abs(x)],yv[xv==4*sigma]),col=shade.col)
    polygon(c(abs(x),xv[xv>=abs(x)]),c(yv[xv==4*sigma],yv[xv>=abs(x)]),col=shade.col)
    p<-round(2*pt(abs(x),nu,lower.tail=FALSE),digits)
    if(show.p==TRUE)legend("topright",bty="n",cex=legend.cex,legend=bquote(paste(2%*%italic(P),"(", italic(X)>="|",.(x),"|) = ",.(p), sep = "")))
    if(show.dist==TRUE)legend("topleft",legend=bquote(paste(italic(X)," ~ ",italic(t),"(",.(nu),")", sep = "")),bty="n",cex=legend.cex)
  }

  if(tail=="middle"){
    polygon(c(xv[xv<=4*sigma],4*sigma),c(yv[xv<=4*sigma],yv[xv==-4*sigma]),col=shade.col)
    polygon(c(xv[xv<=from],from),c(yv[xv<=from],yv[xv==-4*sigma]),col="white")
    polygon(c(to,xv[xv>=to]),c(yv[xv==4*sigma],yv[xv>=to]),col="white")
    p<-round(pt(to,nu)-pt(from,nu),digits)
    if(show.p==TRUE)legend("topright",legend=bquote(paste(italic(P),"(",.(from)<="", italic(X)<=.(to),") = ",.(p), sep = "")),bty="n",cex=legend.cex)
    if(show.dist==TRUE)legend("topleft",legend=bquote(paste(italic(X)," ~ ",italic(t),"(",.(nu),")", sep = "")),bty="n",cex=legend.cex)
  }
}
#------------------------F distribution---------------------------#
#' @importFrom stats pf
#' @export
#' @rdname shade.norm
shade.F<-function(x=NULL, from=NULL, to=NULL, nu1=1, nu2=5, tail="lower",
                  show.p=TRUE, show.d=FALSE, show.dist=TRUE,
                  prob.to.each.tail=0.025, digits=5, legend.cex=.9,
                  shade.col="gray", ...){

  sigma<-qf(.9999,nu1,nu2)
  xv<-seq(0,sigma,sigma/1000)
  yv<-df(xv,nu1,nu2)

  curve(df(x,nu1,nu2),from=0,to=sigma,xlab=expression(italic(x)),ylab=expression(paste(italic(f),"(",italic(x),")", sep = "")), ...)

  if(tail=="lower"){
    polygon(c(xv[xv<=x],x),c(yv[xv<=x],yv[xv==0]),col=shade.col)
    p<-round(pf(x,nu1,nu2,lower.tail=TRUE),digits)
    d<-round(df(x,nu1,nu2),digits)
    if(show.p==TRUE&show.d==FALSE)legend("topright",legend=bquote(paste(italic(P),"(",italic(X)<=.(x),") = ",.(p), sep = "")),bty="n",cex=legend.cex)
    if(show.d==TRUE&show.p==FALSE)legend("topright",legend=bquote(paste("",italic(f),"(",.(x),") = ",.(d), sep = "")),bty="n",cex=legend.cex)
    if(show.d==TRUE&show.p==TRUE)legend("topright",legend=bquote(paste(italic(P),"(",italic(X)<=.(x),") = ",.(p),",  ",italic(f),"(",.(x),") = ",.(d), sep = "")),bty="n",cex=legend.cex)
    if(show.dist==TRUE)legend("top",legend=bquote(paste(italic(X)," ~ ",italic(F),"(",.(nu1),",",.(nu2),")", sep = "")),bty="n",cex=legend.cex,adj=0.5)
  }
  if(tail=="upper"){
    polygon(c(x,xv[xv>=x]),c(yv[xv==sigma],yv[xv>=x]),col=shade.col)
    p<-round(pf(x,nu1,nu2,lower.tail=FALSE),digits)
    d<-round(df(x,nu1,nu2),digits)
    if(show.p==TRUE&show.d==FALSE)legend("topright",legend=bquote(paste(italic(P),"(",italic(X)>=.(x),") = ",.(p), sep = "")),bty="n",cex=legend.cex)
    if(show.d==TRUE&show.p==FALSE)legend("topright",legend=bquote(paste("",italic(f),"(",.(x),") = ",.(d), sep = "")),bty="n",cex=legend.cex)
    if(show.d==TRUE&show.p==TRUE)legend("topright",legend=bquote(paste(italic(P),"(",italic(X)>=.(x),") = ",.(p),",  ",italic(f),"(",.(x),") = ",.(d), sep = "")),bty="n",cex=legend.cex)
    if(show.dist==TRUE)legend("top",legend=bquote(paste(italic(X)," ~ ",italic(F),"(",.(nu1),",",.(nu2),")", sep = "")),bty="n",cex=legend.cex,adj=0.5)
  }
  if(tail=="two"){
    qs<-round(qf(c(prob.to.each.tail,1-prob.to.each.tail),nu1,nu2),1)
    polygon(c(xv[xv<=qs[1]],qs[1]),c(yv[xv<=qs[1]],yv[xv==sigma]),col=shade.col)
    polygon(c(qs[2],xv[xv>=qs[2]]),c(yv[xv==sigma],yv[xv>=qs[2]]),col=shade.col)
    p<-round(prob.to.each.tail*2,digits)
    if(show.p==TRUE)legend("topright",bty="n",cex=legend.cex,legend=bquote(paste(italic(P),"(",.(qs[1])<=italic(X), " and ", italic(X)>=.(qs[2]),") = ",.(p))))
    if(show.dist==TRUE)legend("top",legend=bquote(paste(italic(X)," ~ ",italic(F),"(",.(nu1),",",.(nu2),")", sep = "")),bty="n",cex=legend.cex,adj=0.5)
  }
  if(tail=="middle"){
    polygon(c(xv[xv<=sigma],sigma),c(yv[xv<=sigma],yv[xv==0]),col=shade.col)
    polygon(c(xv[xv<=from],from),c(yv[xv<=from],yv[xv==0]),col="white")
    polygon(c(to,xv[xv>=to]),c(yv[xv==sigma],yv[xv>=to]),col="white")
    p<-round(pf(to,nu1,nu2)-pf(from,nu1,nu2),digits)
    if(show.p==TRUE)legend("topright",legend=bquote(paste(italic(P),"(",.(from)<="", italic(X)<=.(to),") = ",.(p), sep = "")),bty="n",cex=legend.cex)
    if(show.dist==TRUE)legend("top",legend=bquote(paste(italic(X)," ~ ",italic(F),"(",.(nu1),",",.(nu2),")", sep = "")),bty="n",cex=legend.cex,adj=0.5)
  }
}
#-------------------------Chi-square distribution------------------------------#
#' @importFrom stats pchisq dchisq
#' @export
#' @rdname shade.norm
shade.chi <- function(x=NULL, from=NULL, to=NULL, nu=1, tail="lower",
                      show.p=TRUE, show.d=FALSE, show.dist=TRUE,
                      prob.to.each.tail=0.025, digits=5, legend.cex=.9,
                      shade.col="gray", ...){
  sigma<-qchisq(.9999,nu)
  xv<-seq(0,sigma,sigma/10000)
  yv<-dchisq(xv,nu)
  curve(dchisq(x,nu),from=0,to=sigma,xlab=expression(italic(x)),ylab=expression(paste(italic(f),"(",italic(x),")", sep ="")), ...)

  if(tail=="lower"){
    if(nu<3){
      ax <- abs(xv-x)
      xn <- ax==min(ax)
      polygon(c(0, x, x, sort(xv[xv <= x], decreasing = T)), c(0, 0, yv[xn], sort(yv[xv <= x])), col = shade.col)
    }
    if(nu>=3)polygon(c(xv[xv <= x], x), c(yv[xv <= x], yv[xv == 0]), col = shade.col)

    d<-round(dchisq(x,nu),digits)
    p<-round(pchisq(x,nu,lower.tail=TRUE), digits)
    if(show.p==TRUE&show.d==FALSE)legend("topright",legend=bquote(paste(italic(P),"(",italic(X)<=.(x),") = ",.(p), sep = "")),bty="n",cex=legend.cex)
    if(show.d==TRUE&show.p==FALSE)legend("topright",legend=bquote(paste("",italic(f),"(",.(x),") = ",.(d), sep = "")),bty="n",cex=legend.cex)
    if(show.d==TRUE&show.p==TRUE)legend("topright",legend=bquote(paste(italic(P),"(",italic(X)<=.(x),") = ",.(p),",  ",italic(f),"(",.(x),") = ",.(d), sep = "")),bty="n",cex=legend.cex)
    if(show.dist==TRUE)legend("top",legend=bquote(paste(italic(X)," ~ ",chi^2,"(",.(nu),")")),bty="n",cex=legend.cex,adj=0.5)
  }

  if(tail=="upper"){
    polygon(c(x,xv[xv>=x]),c(yv[xv==sigma],yv[xv>=x]),col=shade.col)
    p<-round(pchisq(x,nu,lower.tail=FALSE),digits)
    d<-round(dchisq(x,nu),digits)
    if(show.p==TRUE&show.d==FALSE)legend("topright",legend=bquote(paste(italic(P),"(",italic(X)>=.(x),") = ",.(p), sep = "")),bty="n",cex=legend.cex)
    if(show.d==TRUE&show.p==FALSE)legend("topright",legend=bquote(paste("",italic(f),"(",.(x),") = ",.(d), sep = "")),bty="n",cex=legend.cex)
    if(show.d==TRUE&show.p==TRUE)legend("topright",legend=bquote(paste(italic(P),"(",italic(X)>=.(x),") = ",.(p),",  ",italic(f),"(",.(x),") = ",.(d), sep = "")),bty="n",cex=legend.cex)
    if(show.dist==TRUE)legend("top",legend=bquote(paste(italic(X)," ~ ",chi^2,"(",.(nu),")")),bty="n",cex=legend.cex,adj=0.5)
  }
  if(tail=="two"){
    qs<-round(qchisq(c(prob.to.each.tail,1-prob.to.each.tail),nu),1)
    if(nu >=3){
      polygon(c(xv[xv<=qs[1]],qs[1]),c(yv[xv<=qs[1]],yv[xv==sigma]),col=shade.col)
      polygon(c(qs[2],xv[xv>=qs[2]]),c(yv[xv==sigma],yv[xv>=qs[2]]),col=shade.col)
    }
    p<-round(prob.to.each.tail*2,digits)
    if(show.p==TRUE)legend("topright",bty="n",cex=legend.cex,legend=bquote(paste(italic(P),"(",.(qs[1])<=italic(X), " and ", italic(X)>=.(qs[2]),") = ",.(p))))
    if(show.dist==TRUE)legend("top",legend=bquote(paste(italic(X)," ~ ",chi^2,"(",.(nu),")")),bty="n",cex=legend.cex,adj=0.5)
  }
  if(tail=="middle"){
    if(nu >= 3){
      polygon(c(xv[xv<=sigma],sigma),c(yv[xv<=sigma],yv[xv==0]),col=shade.col)
      polygon(c(xv[xv<=from],from),c(yv[xv<=from],yv[xv==0]),col="white")
      polygon(c(to,xv[xv>=to]),c(yv[xv==sigma],yv[xv>=to]),col="white")
    }
    p<-round(pchisq(to,nu)-pchisq(from,nu),digits)
    if(show.p==TRUE)legend("topright",legend=bquote(paste(italic(P),"(",.(from)<="",italic(X)<=.(to),") = ",.(p), sep = "")),bty="n",cex=legend.cex)
    if(show.dist==TRUE)legend("top",legend=bquote(paste(italic(X)," ~ ",chi^2,"(",.(nu),")")),bty="n",cex=legend.cex,adj=0.5)
  }
}
