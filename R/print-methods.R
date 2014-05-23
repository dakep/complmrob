#' Print the object
#' 
#' Print information about the models returned by \code{\link{complmrob}} or \code{\link{bootcoefs}}.
#' For a detailed description see the help on \code{\link[=summary.complmrob]{summary}}.
#' 
#' @param x the object to be printed.
#' @param conf.level the confidence level for the confidence interval.
#' @param conf.type the type of the printed confidence interval.
#' @param ... ignored.
#' 
#' @seealso \code{\link{summary.complmrob}}
#' @describeIn print for robust linear regression models with compositional data
#' @export
print.complmrob <- function(x, conf.level = 0.95, ...) {
    print(summary(x, conf.level = conf.level))
}

#' @describeIn print for bootstrapped robust linear regression models with or without compositional data
#' @export
print.bootcoefs <- function(x, conf.level = 0.95, conf.type = "perc", ...) {
    print(summary(x, conf.level = conf.level, conf.type = conf.type))
}