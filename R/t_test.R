#' Student's t-Test using values
#' 
#' This function performs the t-test for one and two samples using summarized values, not the vectors.
#' 
#' @param meanx sample mean for sample x.
#' @param varx sample variance for sample x.
#' @param nx sample size for sample x.
#' @param meany sample mean for sample y.
#' @param vary sample variance for sample y. 
#' @param ny sample size for sample y.
#' @param alternative a character string specifying the alternative 
#'  hypothesis, must be one of \code{"two.sided"} (default), 
#'  \code{"greater"} or \code{"less"}. You can specify just the initial letter.
#' @param mu the hypothesized number (mean) in the null hypothesis. 
#' @param conf.level confidence level of the interval, by default its value is 0.95.
#' @param var.equal a logical variable indicating whether to treat the
#'  two variances as being equal. If \code{TRUE} then the pooled variance
#'  is used to estimate the variance otherwise the Welch (or Satterthwaite)
#'  approximation to the degrees of freedom is used. 
#'
#' @return A list with class \code{"htest"} containing the following 
#'  components:
#' \item{statistic}{the value of the statistic.}
#' \item{parameter}{the degrees of freedom for the t-statistic.}
#' \item{p.value}{the p-value for the test.}
#' \item{conf.int}{a confidence interval for the mean appropiate to the
#'       alternative hypothesis.}
#' \item{estimate}{the estimated mean or difference in means depending
#'       on whether it was a one-sample test or a two-sample test.}
#' \item{null.value}{the specified hypothesized value for alternative hypothesis
#'       value for the mean or mean difference depending on whether it was
#'       a one-sample test or a two-sample test.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{method}{a character string indicating the type of test performed.}
#' 
#' @examples
#' # Examples with ONE sample
#' 
#' # Example 9.2 from Walpole
#' t_test(meanx=42, varx=11.9^2, nx=12, 
#'        mu=50, alternative='less')
#'
#' # Example 11.5 from A. F. Siegel & C. Morgan 
#' t_test(meanx=100, varx=12^2, nx=9,
#'        mu=83, alternative='two.sided')
#'        
#' # Example 11.6 from Murray        
#' t_test(meanx=0.053, varx=0.003^2, nx=10,
#'        mu=0.050, alternative='two.sided') 
#'        
#' # --- Examples with TWO-SAMPLES and equal variances ---
#' 
#' # Example 9.3 From Walpole
#' t_test(meanx=85, varx=4^2, nx=12,
#'        meany=81, vary=5^2, ny=10, 
#'        alternative='two.sided', mu=0, var.equal=TRUE)
#'       
#' # Example 13.5 from J. Freund's
#' t_test(meanx=546, varx=31^2, nx=4,
#'        meany=492, vary=26^2, ny=4, 
#'        alternative='two.sided', mu=0, var.equal=TRUE)  
#'         
#' # --- Examples with TWO-SAMPLES and different variances ---
#' 
#' # Example from
#' t_test(meanx =6.8, varx=1.8^2, nx=13,
#'        meany=5.3, vary=1.6^2,  ny=15, 
#'        alternative='two.sided', mu=0, var.equal=FALSE)   
#'     
#' @export
t_test <- function(meanx, varx, nx, 
                   meany=NULL, vary=NULL, ny=NULL, 
                   alternative='two.sided', mu=0, 
                   conf.level=0.95, var.equal=FALSE){
  
  # Checking if the information is correct
  
  # To check if the information about x sample is correct
  if (varx <= 0) 
    stop(paste("The variance x must be positive", "\n", ""))
  if (nx <= 0) 
    stop(paste("The sample size nx must be positive", "\n", ""))
  if (nx %% 1 != 0) 
    stop(paste("The sample size nx must be integer", "\n", ""))
  
  # To check if the user provided information about y sample
  if (xor(is.null(vary), is.null(ny))) 
    stop("Some information about sample y is missing", "\n", "")
  
  # To check if the information about y sample is corrected
  if (! is.null(vary) & ! is.null(ny) & ! is.null(meany)) {
    if (vary <= 0) 
      stop(paste("The variance y must be positive", "\n", ""))
    if (ny <= 0) 
      stop(paste("The sample size ny must be positive", "\n", ""))
    if (ny %% 1 != 0) 
      stop(paste("The sample size ny must be integer", "\n", ""))
    if (is.null(meany) == TRUE)
      stop(paste("Some information about y sample is missing", "\n",""))
  }
  
  # To check if the conf.level it's a number between 1 and 0
  if(conf.level <= 0 || conf.level >= 1)
    stop("The conf.level argument must be > 0 and < 1", "\n", "" )
  if(conf.level < 0.5)
    warning("Confidence levels are often close to 1, eg. 0.95")
  
  # Argument Verification Using Partial Matching
  alternative <- match.arg(arg=alternative,
                           choices=c("two.sided","greater","less"))
  
  # To perform the test
  if (is.null(meany) & is.null(vary) & is.null(ny))
    res <- t_test_one(meanx, varx, nx, alternative, mu, conf.level)
  
  else {
    if (var.equal == TRUE)
      res <- t_test_two_equal(meanx, varx, nx, meany, vary, ny, 
                              alternative, mu, conf.level, var.equal)
    
    if (var.equal == FALSE)
      res <- t_test_two_difer(meanx, varx, nx, meany, vary, ny, 
                              alternative, mu, conf.level, var.equal)
  }
  
  class(res) <- "htest"
  res   
  
}
#' @importFrom stats pt qt
t_test_one <- function(meanx, varx, nx, alternative, mu,
                       conf.level) {
  
  alpha <- 1 - conf.level
  
  if (alternative == 'two.sided') {
    statistic <- (meanx - mu) / sqrt(varx / nx)
    p.value <- 2 * pt(q=abs(statistic), df=nx-1, lower.tail=FALSE)
    quantiles <- c(-qt(p=alpha/2, df=nx-1, lower.tail=FALSE),
                   qt(p=alpha/2, df=nx-1, lower.tail=FALSE))
    conf.int <- meanx + quantiles * sqrt(varx / nx)
  }
  
  if (alternative == 'less') {
    statistic <- (meanx - mu) / sqrt(varx / nx)
    p.value <- pt(q=statistic, df=nx-1, lower.tail=TRUE)
    conf.int <- c(-Inf,
                  meanx + qt(p=1-alpha, df=nx-1) * sqrt(varx / nx))
  }
  
  if (alternative == 'greater') {
    statistic <- (meanx - mu) / sqrt(varx / nx)
    p.value <- pt(q=statistic, df=nx-1, lower.tail=FALSE)
    conf.int <- c(meanx + qt(p=alpha, df=nx-1) * sqrt(varx / nx),
                  Inf)
  }
  
  # To ensure that the output values are in the correct form
  names(statistic) <- 't'
  parameter <- nx - 1
  names(parameter) <- 'df'
  attr(conf.int, 'conf.level') <- conf.level
  estimate <- meanx
  names(estimate) <- 'mean of x'
  null.value <- mu
  names(null.value) <- 'mean'
  method <- 'One Sample t-test'
  data.name <- paste('meanx = ', meanx, ', var = ', varx, ' and nx = ', nx, sep='')
  
  res <- list(statistic=statistic,
              parameter=parameter,
              p.value=p.value,
              conf.int=conf.int,
              estimate=estimate,
              null.value=null.value,
              alternative=alternative,
              method=method,
              data.name=data.name)
  return(res)
}
#' @importFrom stats pt qt
t_test_two_difer <- function(meanx, varx, nx,
                             meany, vary, ny,
                             alternative, mu,
                             conf.level, var.equal) {
  
  
  alpha <- 1 - conf.level
  
  df <- (varx/nx + vary/ny)^2 / ((varx/nx)^2 / (nx-1) + (vary/ny)^2 / (ny-1))
  se <- sqrt(varx/nx + vary/ny)
  statistic <- (meanx - meany - mu) / se
  
  if (alternative == 'two.sided') {
    p.value <- 2 * pt(q=abs(statistic), df=df, lower.tail=FALSE)
    quantiles <- c(-qt(p=alpha/2, df=df, lower.tail=FALSE),
                   qt(p=alpha/2, df=df, lower.tail=FALSE))
    conf.int <- (meanx-meany) + quantiles * se
  }
  
  if (alternative == 'less') {
    p.value <- pt(q=statistic, df=df, lower.tail=TRUE)
    conf.int <- c(-Inf,
                  (meanx-meany) + qt(p=alpha, df=df, lower.tail=F) * se)
  }
  
  if (alternative == 'greater') {
    p.value <- pt(q=statistic, df=df, lower.tail=FALSE)
    conf.int <- c((meanx-meany) - qt(p=alpha, df=df, lower.tail=F) * se,
                  Inf)
  }
  # To ensure that the output values are in the correct form
  names(statistic) <- 't'
  parameter <- df
  names(parameter) <- 'df'
  attr(conf.int, 'conf.level') <- conf.level
  estimate <- c(meanx, meany)
  names(estimate) <- c('mean of x', 'mean of y')
  null.value <- mu
  names(null.value) <- 'difference in means'
  method <- 'Welch Two Sample t-test'
  data.name <- paste('meanx = ', meanx, ', nx = ', nx,
                     ', meany = ', meany, ' and ny = ', ny, sep='')
  
  res <- list(statistic=statistic,
              parameter=parameter,
              p.value=p.value,
              conf.int=conf.int,
              estimate=estimate,
              null.value=null.value,
              alternative=alternative,
              method=method,
              data.name=data.name)
  return(res)
}
#' @importFrom stats pt qt
t_test_two_equal <- function(meanx, varx, nx,
                             meany, vary, ny,
                             alternative, mu,
                             conf.level, var.equal) {
  
  alpha <- 1 - conf.level
  
  df <- nx + ny - 2
  sp2 <- ((nx-1) * varx + (ny-1) * vary) / df
  se <- sqrt(sp2/nx + sp2/ny)
  statistic <- (meanx - meany - mu) / se
  
  if (alternative == 'two.sided') {
    mu <- as.numeric(mu)
    alt <- paste("true difference in means is not equal to", mu)
    p.value <- 2 * pt(q=abs(statistic), df=df, lower.tail=FALSE)
    quantiles <- c(-qt(p=alpha/2, df=df, lower.tail=FALSE),
                   qt(p=alpha/2, df=df, lower.tail=FALSE))
    conf.int <- (meanx-meany) + quantiles * se
  }
  
  if (alternative == 'less') {
    mu <- as.numeric(mu)
    alt <- paste("true difference in means is less than", mu)
    p.value <- pt(q=statistic, df=df, lower.tail=TRUE)
    conf.int <- c(-Inf,
                  (meanx - meany) + qt(p=alpha, df=df, lower.tail=F) * se)
  }
  
  if (alternative == 'greater') {
    mu <- as.numeric(mu)
    alt<- paste("true difference in means is greater than", mu)
    p.value <- pt(q=statistic, df=df, lower.tail=FALSE)
    conf.int <- c((meanx-meany) - qt(p=alpha, df=df, lower.tail=F) * se,
                  Inf)
  }
  
  # To ensure that the output values are in the correct form
  names(statistic) <- 't'
  parameter <- df
  names(parameter) <- 'df'
  attr(conf.int, 'conf.level') <- conf.level
  estimate <- c(meanx, meany)
  names(estimate) <- c('mean of x', 'mean of y')
  null.value <- mu
  names(null.value) <- 'difference in means'
  method <- 'Two Sample t-test'
  data.name <- paste('meanx =', meanx, ', nx =', nx,
                     ', meany =', meany, 'and ny =', ny)
  
  res <- list(statistic = statistic,
              parameter = parameter,
              p.value = p.value,
              conf.int = conf.int,
              estimate = estimate,
              null.value = null.value,
              alternative = alternative,
              method = method,
              data.name = data.name)
  return(res)
}
