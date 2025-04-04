% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/z.test.R
\name{z.test}
\alias{z.test}
\title{Mean test using vector}
\usage{
z.test(x, sigma2, alternative = "two.sided", mu = 0, conf.level = 0.95)
}
\arguments{
\item{x}{a (non-empty) numeric vector of data values.}

\item{sigma2}{population variance which is known.}

\item{alternative}{a character string specifying the alternative
hypothesis, must be one of \code{two.sided} (default),
\code{greater} or \code{less}. You can specify just the initial letter.}

\item{mu}{the hypothesized number in the null hypothesis.}

\item{conf.level}{confidence level of the interval, by default its value is 0.95.}
}
\value{
A list with class \code{htest} containing the following
 components:
\item{statistic}{the value of the statistic.}
\item{p.value}{the p-value for the test.}
\item{conf.int}{a confidence interval for the variance.}
\item{estimate}{the estimated mean.}
\item{null.value}{the specified hypothesized value for alternative hypothesis.}
\item{alternative}{a character string describing the alternative hypothesis.}
\item{method}{a character string indicating the type of test performed.}
\item{data.name}{a character string giving the name of the data.}
}
\description{
This function performs the mean test using raw data (vector).
}
\examples{
# Example 1
# H0: mu = 350
# Ha: mu < 350 with sigma2=25
content <- c(355, 353, 352, 346, 345, 345, 353, 353, 344, 350)
res1 <- z.test(x=content, mu=350, sigma2=25, alternative='less')
res1
plot(res1, shade.col='deepskyblue', col='deepskyblue')

# Example 2
# H0: mu = 170
# Ha: mu != 170 with sigma2=25
x <- rnorm(n=80, mean=171, sd=5)
res2 <- z.test(x=x, mu=170, sigma2=25, alternative='two.sided')
res2
plot(res2, shade.col='indianred1', col='indianred1')

}
