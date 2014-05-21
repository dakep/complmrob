#' Extract the confidence intervals of the coefficients from the compositional lmrob object
#' 
#' @param object the returned object from a call to complmrob.
#' @param param a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level the confidence level required.
#' @param type the type of interval required (see the type argument of \code{\link{boot.ci}}). If \code{theoretical}, the
#'     	theoretical confidence interval from are returned
#' @param ... Currently ignored
#' 
#' @importFrom boot boot.ci
#' @import robustbase
#' @export
confint.complmrob <- function(object, param, level = 0.95, type = c("bca", "perc", "norm", "basic", "stud", "theoretical"), ...) {
    type = match.arg(type);
    
    outtype <- switch(type,
        bca = "bca",
        perc = "percent",
        norm = "normal",
        basic = "basic",
        stud = "student",
        "theoretical"
    )
    
    if(is.numeric(object$R) && type != "theoretical") {		
        ci <- lapply(object$boot.out, function(boot.out) {
            boot::boot.ci(boot.out, conf = level, type = type)[[outtype]][1 , c(4, 5), drop = TRUE];
        })
    } else {
        ci <- list();
        
        if(object$intercept == TRUE) {
            ci <- list("(Intercept)" = confint(object$models[[1]], level = level)[1, ]);	
        }
        
        ci <- c(ci, lapply(object$models, function(m) {
            confint(m, level = level)[object$coefind, ];
        }));
    }
    
    ci <- do.call(rbind, ci);
    colnames(ci) <- stats:::format.perc((1 + c(-1, 1) * level) / 2, 3);
    
    return(ci);
}
