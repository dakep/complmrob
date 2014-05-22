#' Function for robust linear regression using compositional predictors
#' 
#' 
#' 
#' @param formula The formula for the regression model
#' @param data The data.frame to use
#' @return A list of type \code{complmrob} with fields
#' @import robustbase
#' @importFrom boot boot
#' @import parallel
#' @export
complmrob <- function(formula, data) {
    #
    # Initialize auxiliary variables
    #
    mf <- match.call(expand.dots = FALSE);
    
    m <- match(c("formula", "data"), names(mf), 0);
    mf <- mf[c(1, m)];
    mf$drop.unused.levels <- TRUE;
    mf[[1]] <- as.name("model.frame");
    mf <- eval(mf, parent.frame());
    mt <- attr(mf, "terms");
    y <- model.response(mf, "numeric");
    
    compPred <- attr(mt, "term.labels");
    npred <- length(compPred);
    
    int <- (attr(mt, "intercept") == 1);
    coefind <- 1L + int;
    
    #
    # Initialize return object
    #
    ret <- list(
        coefficients = numeric(npred),
        models = vector("list", npred),
        npred = npred,
        predictors = compPred,
        coefind = coefind,
        call = match.call(),
        intercept = int
    );
    
    class(ret) <- "complmrob";
    
    names(ret$coefficients) <- compPred;
    names(ret$models) <- compPred;
 
    for(predInd in seq_along(compPred)) {
        tmpPred <- isomLR(as.matrix(mf[ , compPred]), predInd);
        tmp <- data.frame(y = y, tmpPred);
        
        tmpFormula <- as.formula(sprintf(" y ~ %s", paste(colnames(tmpPred), collapse = " + ")))
        
        if(int == FALSE) {
            tmpFormula <- update(tmpFormula, . ~ . - 1)
        }
        
        m <- robustbase::lmrob(tmpFormula, data = tmp);
        
        ret$models[[predInd]] <- m;
        ret$coefficients[predInd] <- coef(m)[coefind];
    }
    
    # If the intercept is included in the model, evaluate it as well, but only once!
    # The variable m holds the last fit model!
    if(int == TRUE) {
        ret$coefficients <- c(coef(m)["(Intercept)"], ret$coefficients);
    }

    return(ret);
}