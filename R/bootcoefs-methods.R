#' Get the bootstrap distribution of the parameter estimates
#'
#' @param object the object to bootstrap the coefficients from
#' @param R the number of replicates
#' @param method one of \code{"frb"} for fast and robust bootstrap, \code{"residuals"} to resample
#'      the residuals or \code{"cases"} to resample the cases
#' @param ncpus the number of CPUs to utilize for bootstrapping
#' @param cl a snow or parallel cluster to use for bootstrapping
#' @param ... further arguments needed by the specializations
#' @export
bootcoefs <- function(object, R = 999, method = c("frb", "residuals", "cases"), ncpus = NULL, cl = NULL, ...) {
    UseMethod("bootcoefs", object);
}

#' Bootstrap the coefficients from the compositional robust regression model
#' 
#' @param object the complmrob model for which the estimates are bootstrapped
#' @param R the number of replicates
#' @param method one of \code{"frb"} for fast and robust bootstrap, \code{"residuals"} to resample
#'      the residuals or \code{"cases"} to resample the cases
#' @param ncpus the number of CPUs to utilize for bootstrapping
#' @param cl a snow or parallel cluster to use for bootstrapping
#' @param ... currently ignored
#' @export
#' @importFrom boot boot
bootcoefs.complmrob <- function(object, R = 999, method = c("frb", "residuals", "cases"), ncpus = NULL, cl = NULL, ...) {
    #
    # Initialize auxiliary variables
    #
    method <- match.arg(method);
    
    clSetup <- setupCluster(ncpus, cl);
    
    bootParams <- list(
        R = quote(R),
        coefind = quote(object$coefind),
        parallel = quote(clSetup$parallel),
        ncpus = quote(length(clSetup$cl)),
        cl = quote(clSetup$cl)
    )

    bootres <- list();
    
    tryCatch({
        if(method == "frb") {
            bootres <- lapply(object$models, function(m, bootParams) {
                bootParams$data <- quote(model.frame(m));
                bootParams$control <- quote(bootStatFastControl(m));
                bootParams$statistic <- quote(bootStatFast);
                
                bc <- as.call(c(list(expression(boot::boot)[[1]]), bootParams))
                return(eval(bc));
            }, bootParams);
            
            if(object$intercept == TRUE) {
                m <- object$models[[1]];
                bootParams$coefind <- 1;
                bootParams$data <- quote(model.frame(m));
                bootParams$control <- quote(bootStatFastControl(m));
                bootParams$statistic <- quote(bootStatFast);
                
                bc <- as.call(c(list(expression(boot::boot)[[1]]), bootParams))
                bootres[["(Intercept)"]] <- eval(bc);
            }
        } else if(method == "residuals") {
            bootres <- lapply(object$models, function(m, bootParams) {
                respInd <- attr(m$terms, "response");
                
                bootParams$data <- quote(data.frame(model.frame(m)[ , -respInd, drop = FALSE],
                    fit = fitted(m), resid = residuals(m)));
                bootParams$weights = quote(m$rweights);
                bootParams$statistic <- quote(bootStatResiduals);
                bootParams$intercept <- quote(object$intercept);
                
                bc <- as.call(c(list(expression(boot::boot)[[1]]), bootParams))
                return(eval(bc));
            }, bootParams);
            
            if(object$intercept == TRUE) {
                m <- object$models[[1]];
                respInd <- attr(m$terms, "response");
                bootParams$coefind <- 1;
                bootParams$data <- quote(data.frame(model.frame(m)[ , -respInd, drop = FALSE],
                    fit = fitted(m), resid = residuals(m)));
                bootParams$weights = quote(m$rweights);
                bootParams$statistic <- quote(bootStatResiduals);
                bootParams$intercept <- quote(object$intercept);
                
                bc <- as.call(c(list(expression(boot::boot)[[1]]), bootParams))
                bootres[["(Intercept)"]] <- eval(bc);
            }
        } else {
            bootres <- lapply(object$models, function(m, bootParams) {
                bootParams$data <- quote(model.frame(m));
                bootParams$formula <- quote(formula(m$terms));
                bootParams$weights = quote(m$rweights);
                bootParams$statistic <- quote(bootStatCases);
                
                bc <- as.call(c(list(expression(boot::boot)[[1]]), bootParams))
                return(eval(bc));
            }, bootParams);
            
            if(object$intercept == TRUE) {
                m <- object$models[[1]];
                bootParams$coefind <- 1;
                bootParams$data <- quote(model.frame(m));
                bootParams$formula <- quote(formula(m$terms));
                bootParams$weights = quote(m$rweights);
                bootParams$statistic <- quote(bootStatCases);
                
                bc <- as.call(c(list(expression(boot::boot)[[1]]), bootParams))
                bootres[["(Intercept)"]] <- eval(bc);
            }
        }
        
        if(object$intercept == TRUE) {
            #reorder bootres so that the intercept is first
            bootres <- bootres[c(length(bootres), seq_len(length(bootres) - 1))];
        }
    }, error = function(e) {
        print(e);
    }, finally = {
        if(clSetup$needToShutdownCluster == TRUE) {
            parallel::stopCluster(clSetup$cl);
        }
    });

    ret <- list(
        bootres = bootres,
        model = object,
        R = R
    );
    
    class(ret) <- c("bootcoefs", "bccomplmrob");
    return(ret);
}

#' Bootstrap the coefficients from a single robust regression model
#' 
#' @param object the lmrob model object for which the estimates are bootstrapped
#' @param R the number of replicates
#' @param method one of \code{"frb"} for fast and robust bootstrap, \code{"residuals"} to resample
#'      the residuals or \code{"cases"} to resample the cases
#' @param ncpus the number of CPUs to utilize for bootstrapping
#' @param cl a snow or parallel cluster to use for bootstrapping
#' @param ... currently ignored
#' @export
#' @importFrom boot boot
bootcoefs.lmrob <- function(object, R = 999, method = c("frb", "residuals", "cases"), ncpus = NULL, cl = NULL, ...) {
    #
    # Initialize auxiliary variables
    #
    method <- match.arg(method);

    clSetup <- setupCluster(ncpus, cl);

    bootParams <- list(
        R = quote(R),
        coefind = quote(object$coefind),
        parallel = quote(clSetup$parallel),
        ncpus = quote(length(clSetup$cl)),
        cl = quote(clSetup$cl)
    )

    bootres = NULL;
    
    tryCatch({
        if(method == "frb") {            
            bootres <- boot::boot(data = model.frame(object), statistic = bootStatFast,
                R = R, parallel = clSetup$parallel, ncpus = length(clSetup$cl), cl = clSetup$cl,
                coefind = seq_along(coef(object)), control = bootStatFastControl(object));
        } else if(method == "residuals") {
            respInd <- attr(object$terms, "response");
            tmpData <- data.frame(model.frame(object)[ , -respInd, drop = FALSE], fit = fitted(object), resid = residuals(object));

            bootres <- boot::boot(data = tmpData, statistic = bootStatResiduals, weigths = object$rweights,
                R = R, parallel = clSetup$parallel, ncpus = length(clSetup$cl), cl = clSetup$cl,
                intercept = (attr(object$terms, "intercept") == 1), coefind = seq_along(coef(object)));
        } else {
            bootres <- boot::boot(data = model.frame(object), statistic = bootStatCases, weights = object$rweights,
                R = R, parallel = clSetup$parallel, ncpus = length(clSetup$cl), cl = clSetup$cl,
                coefind = seq_along(coef(object)), formula = formula(object$terms));
        }
    }, error = function(e) {
        print(e);
    }, finally = {
        if(clSetup$needToShutdownCluster == TRUE) {
            parallel::stopCluster(clSetup$cl);
        }
    });

    ret <- list(
        bootres = bootres,
        model = object,
        R = R
    );
    
    class(ret) <- c("bootcoefs", "bclmrob");
    return(ret);
}

#' @import parallel
setupCluster <- function(ncpus, cl) {
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
                loadNamespace("robustbase");
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

    return(list(
        needToShutdownCluster = needToShutdownCluster,
        parallel = parallel,
        cl = cl
    ));
}