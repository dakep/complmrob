#' Function for robust linear regression using compositional predictors
#' 
#' 
#' 
#' @param formula The formula for the regression model
#' @param data The data.frame to use
#' @param boottype One of \code{fast}, \code{residuals} or \code{cases} (see Details for more information)
#' @param R The number of bootstrap replicates to generate to estimate the standard error and confidence
#'     	intervals. If NULL or missing, the theoretical standard error and confidence intervals are returned.
#' @param boot.params a list of parameters additionally passed to the boot function for bootstrapping
#' 		(see \code{\link{boot}} for available arguments). Especially interesting if parallel execution
#' 		is required.
#' @return A list of type \code{complmrob} with fields
#' @import robustbase
#' @importFrom boot boot
#' @import parallel
#' @export
complmrob <- function(formula, data, boottype = c("fast", "residuals", "cases"), R = 100, ncpus = NULL, cl = NULL) {
    #
    # Initialize auxiliary variables
    #
    boottype <- match.arg(boottype);
    
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
        std.errors = numeric(npred),
        models = vector("list", npred),
        boot.out = vector("list", npred),
        R = R,
        npred = npred,
        predictors = compPred,
        coefind = coefind,
        formula = formula,
        intercept = int,
        boottype = boottype
    );
    
    class(ret) <- "complmrob";
    
    names(ret$coefficients) <- compPred;
    names(ret$std.errors) <- compPred;
    names(ret$models) <- compPred;
    names(ret$boot.out) <- compPred;
    
    ##
    ## Setup clusters (if any)
    ##
    needToShutdownCluster <- FALSE;
    parallel <- "no";
    if(!is.null(ncpus) && is.null(cl)) {
        cl <- parallel::makeForkCluster(nnodes = ncpus);
        needToShutdownCluster <- TRUE;
    }
    
    if(!is.null(cl)) {
        parallel <- "snow";
        tryCatch({
            parallel::clusterEvalQ(cl, {
                suppressPackageStartupMessages(library(robustbase));
            });
            parallel::clusterExport(cl, varlist = c("isomLR"));
        }, error = function(e) {
            if(needToShutdownCluster == TRUE) {
                parallel::stopCluster(cl);
                parallel <<- "no";
                needToShutdownCluster <<- FALSE;
            }
            cl <<- NULL;
        }, finally = function(...) {})
    }
    
    ##
    ## Do the bootstrapping
    ##
    bootParams <- list(
        R = quote(R),
        coefind = quote(coefind),
        parallel = quote(parallel),
        ncpus = quote(length(cl)),
        cl = quote(cl)
    )
    
    tryCatch({
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
            
            if(is.numeric(R)) {
                if(boottype == "fast") {
                    bootParams$data <- quote(tmp);
                    bootParams$control <- quote(bootStatFastControl(m));
                    bootParams$statistic <- quote(bootStatFast);
                } else if(boottype == "residuals") {
                    bootParams$data <- quote(data.frame(tmpPred, fit = fitted(m), resid = residuals(m)));
                    bootParams$weights = quote(m$rweights);
                    bootParams$statistic <- quote(bootStatResiduals);
                    bootParams$intercept <- quote(int);
                } else { #boottype == "cases"
                    bootParams$data <- quote(tmp);
                    bootParams$weights = quote(m$rweights);
                    bootParams$statistic <- quote(bootStatCases);
                    bootParams$intercept <- quote(int);
                }
                
                bc <- as.call(c(list(expression(boot::boot)[[1]]), bootParams))
                br <- eval(bc);
                ret$std.errors[predInd] <- sd(br$t, na.rm = TRUE);
                ret$boot.out[[predInd]] <- br;
            } else {
                ret$std.errors[predInd] <- sqrt(m$cov[coefind, coefind]);
            }
        }
        
        # If the intercept is included in the model, evaluate it as well, but only once!
        # The variable m holds the last fit model!
        if(int == TRUE) {
            ret$coefficients <- c(coef(m)["(Intercept)"], ret$coefficients);
            
            if(is.numeric(R)) {
                bootParams$coefind <- 1;
                if(boottype == "fast") {
                    bootParams$data <- quote(tmp);
                    bootParams$control <- quote(bootStatFastControl(m));
                    bootParams$statistic <- quote(bootStatFast);
                } else if(boottype == "residuals") {
                    bootParams$data <- quote(data.frame(tmpPred, fit = fitted(m), resid = residuals(m)));
                    bootParams$statistic <- quote(bootStatResiduals);
                    bootParams$intercept <- quote(int);
                } else { #boottype == "cases"
                    bootParams$data <- quote(tmp);
                    bootParams$statistic <- quote(bootStatCases);
                    bootParams$intercept <- quote(int);
                }
                
                bc <- as.call(c(list(expression(boot::boot)[[1]]), bootParams))
                br <- eval(bc);
                ret$std.errors <- c("(Intercept)" = sd(br$t, na.rm = TRUE), ret$std.errors);
                ret$boot.out[["(Intercept)"]] <- br;
                # Reorder list
                ret$boot.out <- ret$boot.out[c(length(ret$boot.out), seq_len(length(ret$boot.out) - 1))];
            } else {
                ret$std.errors <- c("(Intercept)" = sqrt(m$cov[1L, 1L]), ret$std.errors);
            }
        }
    },
        error = function(e) {
            print(e);
        }, finally = {
            if(needToShutdownCluster == TRUE) {
                parallel::stopCluster(cl);
            }
        });
    
    return(ret);
}