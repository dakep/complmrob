#' @export
print.complmrob <- function(x, conf.level = 0.95, conf.type = "perc", ...) {
    print(summary(x, conf.level = conf.level, conf.type = conf.type))
}

#' @export
print.bootcoefs <- function(x, conf.level = 0.95, conf.type = "perc", ...) {
    print(summary(x, conf.level = conf.level, conf.type = conf.type))
}