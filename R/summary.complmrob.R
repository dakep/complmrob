#' @import robustbase
#' @export
summary.complmrob <- function(x, type = c("bootstrapped", "theoretical"), conf.level = 0.95, conf.type = "perc", ...) {    
    type <- match.arg(type);
    
    if(!is.numeric(x$R)) {
        type <- "theoretical";
    }
    
    ret <- list(
        obj = x,
        ci = NULL,
        stats = NULL,
        r.squared = NULL,
        adj.r.squared = NULL,
        scale = NULL,
        type = type
    );
    
    if(type == "theoretical") {
        intercSe <- sqrt(x$models[[1]]$cov[x$coefind, x$coefind]);
        intercTval <- x$models[[1]]$coefficients[1] / intercSe;
        
        thParams <- list();
        
        if(x$intercept == TRUE) {
            thParams <- list("(Intercept)" = c(
                se = intercSe,
                tval = intercTval,
                pval = 2 * pt(abs(intercTval), x$models[[1]]$df.residual, lower.tail = FALSE)
            ));
        }
        
        thParams <- c(thParams, lapply(x$models, function(m) {
            se <- sqrt(m$cov[x$coefind, x$coefind]);
            tval <- m$coefficients[x$coefind] / se;
            return(c(
                se = se,
                tval <- tval,
                pval = 2 * pt(abs(tval), m$df.residual, lower.tail = FALSE)
            ));
        }));
        
        ret$stats <- as.data.frame(do.call(rbind, thParams))
        colnames(ret$stats) <- c("Std. Error", "t value", "Pr(>|t|)");
        
        ret$ci <- list();
        
        if(x$intercept == TRUE) {
            ret$ci <- list("(Intercept)" = confint(x$models[[1]], level = conf.level)[1L, ]);
        }
        
        ret$ci <- c(ret$ci, lapply(x$models, function(m) {
            confint(m, level = conf.level)[x$coefind, ]
        }));
        
        ret$ci <- do.call(rbind, ret$ci);
        
    } else {
        statBs <- do.call(rbind, lapply(x$boot.out, function(bo) {
            bias <- mean(bo$t, na.rm = TRUE) - bo$t0;
            pval <- (sum((sign(bo$t0) * bo$t) < 0, na.rm = TRUE) + 1) / (sum(!is.na(bo$t)) + 1);
            c(bias = bias, pval = pval)
        }));
        ret$stats <- cbind("bias" = statBs[ , 1L], "Std. Error" = x$std.errors, "Pr(b<>0)" = statBs[ , 2L]);
        
        ret$ci <- confint(x, level = conf.level, type = conf.type);
    }
    
    sm <- summary(x$models[[1]]);
    sm$residuals
    ret$r.squared <- sm$r.squared;
    ret$adj.r.squared <- sm$adj.r.squared;
    ret$scale = x$models[[1]]$scale;
    
    class(ret) <- "summary.complmrob";
    return(ret);
}

#' @export
print.summary.complmrob <- function(x, digits = max(3, getOption("digits") - 3), signif.stars = getOption("show.signif.stars"), ...) {
    cat("Robust linear regression with compositional covariates\n");
    print(x$obj$formula);
    if(x$type == "bootstrapped") {
        cat("\nStandard errors and derived statistics are base on ", x$obj$R, " replications\n", sep = "");
    }
    cat("\n");
    
    prdf <- cbind("Estimate" = coef(x$obj), x$stats);
    printCoefmat(prdf, digits = digits, signif.stars = signif.stars, cs.ind = seq_len(ncol(prdf) - 1), tst.ind = NULL)
    
    cat("\nConfidence intervals:\n");
    printCoefmat(cbind("Estimate" = coef(x$obj), x$ci), P.values = FALSE, cs.ind = seq_len(3), tst.ind = NULL);
    
    cat("\nRobust residual standard error:", format(signif(x$scale, digits)), "\n");
    cat("Multiple R-squared: ", formatC(x$r.squared, digits = digits), sep = "");
    cat("\tAdjusted R-squared: ", formatC(x$adj.r.squared, digits = digits), "\n", sep = "");
}
