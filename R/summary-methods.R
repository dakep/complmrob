#' @import robustbase
#' @export
summary.bccomplmrob <- function(x, conf.level = 0.95, conf.type = "perc", ...) {    
    ret <- list(
        obj = x$model,
        R = x$R,
        ci = NULL,
        stats = NULL,
        r.squared = NULL,
        adj.r.squared = NULL,
        scale = NULL,
        type = "bootstrapped"
    );
    
    statBs <- do.call(rbind, lapply(x$bootres, function(bo) {
        bias <- mean(bo$t, na.rm = TRUE) - bo$t0;
        pval <- (sum((sign(bo$t0) * bo$t) < 0, na.rm = TRUE) + 1) / (sum(!is.na(bo$t)) + 1);
        c(bias = bias, se = sd(bo$t), pval = pval)
    }));
    
    ret$stats <- cbind("bias" = statBs[ , 1L], "Std. Error" = statBs[ , 2L], "Pr(b<>0)" = statBs[ , 3L]);
    
    ret$ci <- confint(x, level = conf.level, type = conf.type);
    
    sm <- summary(x$model$models[[1]]);
    ret$r.squared <- sm$r.squared;
    ret$adj.r.squared <- sm$adj.r.squared;
    ret$scale = x$model$models[[1]]$scale;

    class(ret) <- "summary.complmrob";
    return(ret);
}

#' @import robustbase
#' @importFrom boot boot.ci
#' @export
summary.bclmrob <- function(x, conf.level = 0.95, conf.type = "perc", ...) {    
    ret <- list(
        obj = x$model,
        R = x$R,
        ci = NULL,
        stats = NULL,
        r.squared = NULL,
        adj.r.squared = NULL,
        scale = NULL,
        type = "bootlmrob"
    );

    nc <- ncol(x$bootres$t);
    bsl <- split(x$bootres$t, rep.int(seq_len(nc), times = rep.int(nrow(x$bootres$t), nc)));

    statBs <- do.call(rbind, mapply(function(t, t0) {
        bias <- mean(t, na.rm = TRUE) - t0;
        pval <- (sum((sign(t0) * t) < 0, na.rm = TRUE) + 1) / (sum(!is.na(t)) + 1);
        c(bias = bias, se = sd(t), pval = pval)
    }, bsl, x$bootres$t0, SIMPLIFY = FALSE));
    
    ret$stats <- cbind("bias" = statBs[ , 1L], "Std. Error" = statBs[ , 2L], "Pr(b<>0)" = statBs[ , 3L]);

    ret$ci <- confint(x, level = conf.level, type = conf.type);
    
    sm <- summary(x$model);
    ret$r.squared <- sm$r.squared;
    ret$adj.r.squared <- sm$adj.r.squared;
    ret$scale = x$model$scale;
    
    class(ret) <- "summary.complmrob";
    return(ret);
}

#' @import robustbase
#' @export
summary.complmrob <- function(x, conf.level = 0.95, ...) {
    ret <- list(
        obj = x,
        ci = NULL,
        stats = NULL,
        r.squared = NULL,
        adj.r.squared = NULL,
        scale = NULL,
        type = "theoretical"
    );
    
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

    sm <- summary(x$models[[1]]);
    ret$r.squared <- sm$r.squared;
    ret$adj.r.squared <- sm$adj.r.squared;
    ret$scale = x$models[[1]]$scale;
    
    class(ret) <- "summary.complmrob";
    return(ret);
}

#' @export
print.summary.complmrob <- function(x, digits = max(3, getOption("digits") - 3), signif.stars = getOption("show.signif.stars"), ...) {
    cat("Robust linear regression with compositional covariates\n");
    print(x$obj$call);
    if(x$type == "bootstrapped" || x$type == "bootlmrob") {
        cat("\nStandard errors and derived statistics are base on ", x$R, " bootstrap replications\n", sep = "");
    }
    cat("\nCoefficients:\n");
    
    prdf <- cbind("Estimate" = coef(x$obj), x$stats);
    printCoefmat(prdf, digits = digits, signif.stars = signif.stars, cs.ind = seq_len(ncol(prdf) - 1), tst.ind = NULL)
    
    cat("\nConfidence intervals:\n");
    printCoefmat(cbind("Estimate" = coef(x$obj), x$ci), P.values = FALSE, cs.ind = seq_len(3), tst.ind = NULL);
    
    cat("\nRobust residual standard error:", format(signif(x$scale, digits)), "\n");
    cat("Multiple R-squared: ", formatC(x$r.squared, digits = digits), sep = "");
    cat("\tAdjusted R-squared: ", formatC(x$adj.r.squared, digits = digits), "\n", sep = "");
}
