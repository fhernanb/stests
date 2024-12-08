% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ci_p.R
\name{ci_p}
\alias{ci_p}
\title{Confidence intervals for Binomial proportions}
\usage{
ci_p(x, n, conf.level = 0.95, intervalType = "wald")
}
\arguments{
\item{x}{a number or a vector with the number of successes.}

\item{n}{a number or a vector with the number of trials.}

\item{conf.level}{confidence level for the returned confidence interval. By default is 0.95.}

\item{intervalType}{type of confidence interval, possible choices are: "wald".}
}
\value{
A dataframe with the input information and the confidence interval.
}
\description{
This function obtains the confidence interval for a proportion.
}
\examples{
# Examples for single values
ci_p(x=5, n=15, intervalType="wald", conf.level=0.90)
ci_p(x=7, n=17, intervalType="wald", conf.level=0.95)
ci_p(x=9, n=19, intervalType="wald", conf.level=0.97)

# Examples for multiple values
ci_p(x=c(5, 7, 9),
     n=c(15, 17, 19),
     intervalType="wald",
     conf.level=c(0.90, 0.95, 0.97))
}
\author{
Freddy Hernandez, \email{fhernanb@unal.edu.co}
}
