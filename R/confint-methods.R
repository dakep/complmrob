#' Extract the confidence intervals of the coefficients from the compositional lmrob object
#' 
#' @param object the returned object from a call to bootcoefs
#' @param param a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level the confidence level required.
#' @param type the type of interval required (see the type argument of \code{\link{boot.ci}}).
#' @param ... Currently ignored
#' 
#' @importFrom boot boot.ci
#' @import robustbase
#' @export
confint.bccomplmrob <- function(object, param, level = 0.95, type = c("bca", "perc", "norm", "basic", "stud"), ...) {
    type = match.arg(type);
    
    outtype <- switch(type,
        bca = "bca",
        perc = "percent",
        norm = "normal",
        basic = "basic",
        stud = "student",
        "percent"
    )
    ci <- lapply(object$bootres, function(boot.out) {
        boot::boot.ci(boot.out, conf = level, type = type)[[outtype]][1 , c(4, 5), drop = TRUE];
    });

    ci <- do.call(rbind, ci);
    colnames(ci) <- stats:::format.perc((1 + c(-1, 1) * level) / 2, 3);
    
    return(ci);
}


#' Extract the confidence intervals of the coefficients from the bootstrapped lmrob object
#' 
#' @param object the returned object from a call to bootcoefs
#' @param param a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level the confidence level required.
#' @param type the type of interval required (see the type argument of \code{\link{boot.ci}}).
#' @param ... Currently ignored
#' 
#' @importFrom boot boot.ci
#' @import robustbase
#' @export
confint.bclmrob <- function(object, param, level = 0.95, type = c("bca", "perc", "norm", "basic", "stud"), ...) {
    type = match.arg(type);
    
    outtype <- switch(type,
        bca = "bca",
        perc = "percent",
        norm = "normal",
        basic = "basic",
        stud = "student",
        "percent"
    )

    ci <- lapply(seq_len(ncol(object$bootres$t)), function(i) {
        boot.ci(object$bootres, conf = level, type = type, index = i)[[outtype]][1, 4:5, drop = TRUE]
    });
    
    ci <- do.call(rbind, ci);
    colnames(ci) <- stats:::format.perc((1 + c(-1, 1) * level) / 2, 3);
    
    return(ci);
}