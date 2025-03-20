#' Confidence intervals for Binomial proportions
#'
#' @author Freddy Hernandez, \email{fhernanb@unal.edu.co}
#'
#' @description
#' This function obtains the confidence interval for a proportion.
#'
#' @param x a number or a vector with the number of successes.
#' @param n a number or a vector with the number of trials.
#' @param conf.level confidence level for the returned confidence interval. By default is 0.95.
#' @param intervalType type of confidence interval, possible choices are:
#' "wald", "agresti_coull", "rindskopf",
#' "clopper_pearson", "add_4", "arcsine_cc",
#' "arcsine", "arcsine_ac", "wilson",
#' "ci_p_jeffreys", "hpd_jeffreys", "LRT",
#' "mid_p", "agresti_caffo",
#'
#' @return A dataframe with the input information and the confidence interval.
#'
#' @seealso \link{ci_p_wald},
#' \link{ci_p_agresti_coull},
#' \link{ci_p_rindskopf},
#' \link{ci_p_clopper_pearson},
#' \link{ci_p_add_4},
#' \link{ci_p_arcsine_cc},
#' \link{ci_p_arcsine},
#' \link{ci_p_arcsine_ac},
#' \link{ci_p_wilson},
#' \link{ci_p_jeffreys},
#' \link{ci_p_hpd_jeffreys}
#' \link{ci_p_mid_p},
#' \link{ci_p_agresti_caffo},
#'
#' @example examples/examples_ci_p.R
#' @export
#'
ci_p <- function(x, n, conf.level=0.95, intervalType="wald") {
  # Checks
  if (any(x < 0)) stop("x must not be negative.")
  if (any(n < 0)) stop("n must not be negative.")
  if (any(x > n)) stop("x must be lower than n.")
  if (missing(x)) stop("argument 'x' is missing, with no default")
  if (missing(n)) stop("argument 'n' is missing, with no default")
  if (any(conf.level < 0 | conf.level > 1))
    stop("conf.level must be between 0 and 1.")
  # To ensure integer values
  if (max(abs(x - round(x))) > 1e-07)
    stop("x must be nonnegative and integer")
  if (max(abs(n - round(n))) > 1e-07)
    stop("n must be nonnegative and integer")

  # To select the interval type
  method <- match.arg(arg=intervalType,
                      choices=c("wald",
                                "agresti_coull",
                                "rindskopf",
                                "clopper_pearson",
                                "add_4",
                                "arcsine_cc",
                                "arcsine",
                                "arcsine_ac",
                                "wilson",
                                "jeffreys",
                                "hpd_jeffreys",
                                "LRT",
                                "mid_p",
                                "agresti_caffo"))

  # To generate the code for evaluating, without using cases
  my_code <- paste0("ci_p_", method,
                    "(x=x, n=n, conf.level=conf.level)")

  # To obtain the result
  result <- eval(parse(text=my_code))

  # To extract the lower and upper limits
  lower <- result[1, ]
  upper <- result[2, ]

  # To format the result
  result <- data.frame(x=x, n=n, proportion=x/n, lower=lower,
                       upper=upper, conf.level=conf.level)

  return(result)
}

