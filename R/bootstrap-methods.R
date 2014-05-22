#' Function to calculate the coefficient(s) of the robust linear regression model
#' from a bootstrapped sample
#'
#' @param dat The original data set (with the columns fit, resid and the predictor variables)
#' @param inds The resampled indices
#' @param coefind The index of the coefficient to extract
#' @param intercept If the model includes an intercept term
#' @param maxTries The maximum number of tries to increase the maxit control arguments for the S estimator
#' 
#' @import robustbase
bootStatResiduals <- function(dat, inds, coefind, intercept = TRUE, maxTries = 4L) {
    formula <- fit + resid[inds] ~ .;
    
    if(intercept == FALSE) {
        formula <- fit + resid[inds] ~ . - 1;
    }
    
    suppressWarnings(m <- robustbase::lmrob(formula, data = dat));
    
    # Ensure convergence!
    itcount <- 0L;
    while(!m$converged && itcount < maxTries) {
        itcount <- itcount + 1L;
        suppressWarnings(m <- update(m, maxit.scale = m$control$maxit.scale + 200, max.it = m$control$max.it + 50))
    }
    if(itcount == maxTries) {
        return(NA_real_);
    }
    return(coef(m)[coefind])
}

#' Function to calculate the coefficient(s) of the robust linear regression model
#' from a bootstrapped sample
#'
#' @param dat The original data set
#' @param inds The resampled indices
#' @param coefind The index of the coefficient to extract
#' @param formula the formula to fit the model
#' @param maxTries The maximum number of tries to increase the maxit control arguments for the S estimator
#' 
#' @import robustbase
bootStatCases <- function(dat, inds, coefind, formula, maxTries = 4L) {
    suppressWarnings(m <- robustbase::lmrob(formula, data = dat[inds, ]));
    
    # Ensure convergence!
    itcount <- 0L;
    while(!m$converged && itcount < maxTries) {
        itcount <- itcount + 1L;
        suppressWarnings(m <- update(m, maxit.scale = m$control$maxit.scale + 200, max.it = m$control$max.it + 50))
    }
    if(itcount == maxTries) {
        return(NA_real_);
    }
    return(coef(m)[coefind])
}

#' Function to extract the necessary data from the model and to calculate the correction
#' factors for the fast and robust bootstrap
#'
#' @param model The lmrob model
#' @import robustbase
#' @references Salibian-Barrera M., Zamar R.: Bootstrapping robust estimates of regression, 2002
bootStatFastControl <- function(model) {
    X <- model.matrix(model);
    n <- attr(X, "dim")[1];
    p <- attr(X, "dim")[2];
    vfactorinv <- (n - p) * model$init.S$control$bb;
#     vfactor <- 1 / (model$init.S$control$bb * n);
    
    scaledRes <- model$residuals / model$scale;
    scaledResS <- model$init.S$residuals / model$scale;
    
    w <- model$rweights / model$scale;
    v <- (model$scale / vfactorinv) *
        robustbase::Mchi(scaledResS, deriv = 0, cc = model$init.S$control$tuning.chi, psi = model$init.S$control$psi);
    
    wp2 <- robustbase::Mpsi(scaledRes, deriv = 1, cc = model$control$tuning.psi, psi = model$control$psi);
    
    sf <- crossprod(X, diag(wp2)) %*% X;
    sfinv <- NULL;

    tryCatch({
        sfinv <- solve(sf);
    }, error = function(e) {
        warning("The data is (almost) singular, results will not be reliable.");
        svddecomp <- svd(sf);
        sfinv <<- tcrossprod(svddecomp$v %*% diag(1/svddecomp$d), svddecomp$u);
    });
    M <- model$scale * sfinv %*% crossprod(X, diag(w)) %*% X
    
    chideriv <- robustbase::Mchi(scaledResS, deriv = 1, cc = model$init.S$control$tuning.chi, psi = model$init.S$control$psi);
    
    a <- (chideriv %*% scaledResS)[1, 1, drop = TRUE];
    d <- (sfinv %*% crossprod(X, wp2 * model$residuals)) * vfactorinv / (a);# * model$scale);
    
    ret <- list(
        M = M,
        d = d,
        coefficients = coef(model),
        v = v,
        wsqrt = sqrt(w),
        scale = model$scale,
        residualsS = model$init.S$residuals,
        terms = terms(model)
    )
    
    class(ret) <- "bootStatFastControl";
    return(ret);
}

#' Function to calculate the bootstrapped coefficient of a robust linear regression model
#'
#' @param data The original data frame
#' @param inds The resampled indices
#' @param control The control object as returned by \code{\link{bootStatFastControl}}.
#' @param coefind The index of the coefficient to extract
#' 
#' @references Salibian-Barrera M., Zamar R.: Bootstrapping robust estimates of regression, 2002
#' @import robustbase
bootStatFast <- function(data, inds, control, coefind) {
    bsResidS <- control$residualsS[inds];
    
    mf <- model.frame(control$terms, data[inds, , drop = FALSE]);
    X <- model.matrix(control$terms, mf);
    y <- model.response(mf, "numeric");
    
    wts <- control$wsqrt[inds];
    
    bsCoef <- .lm.fit(X * wts, y * wts)$coefficients;
    bsScale <- sum(control$v[inds]);

    bootCoefs <- control$coefficients + control$M %*% (bsCoef - control$coefficients) + control$d * (bsScale - control$scale);
    
    return(bootCoefs[coefind])
}