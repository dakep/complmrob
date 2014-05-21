#' Diagnostic plots for the robust regression model with compositional covariats
#'
#' @param x The complmrob object to plot
#' @param type One of \code{response}, \code{model} or \code{bootstrap}
#' @param se Should the confidence interval be shown in the response plot
#' @param conf.level If the confidence interval is shown in the response plot, this parameter sets
#'      the level of the confidence interval
#' 
#' @importFrom reshape2 melt
#' @import ggplot2
#' @export
plot.complmrob <- function(x, type = c("response", "model", "bootstrap"), se = TRUE, conf.level = 0.95, ...) {    
    type <- match.arg(type);
    
    if(type == "model") {
        plot(x$models[[1]]);
    } else if(type == "bootstrap") {
        resamples <- reshape2::melt(lapply(x$boot.out, function(bo) {
            bo$t[ , 1, drop = TRUE]
        }));
        
        colnames(resamples) <- c("value", "coef")
        trueCoefs <- data.frame(coef = names(x$coefficients), value = unname(x$coefficients));
        
        p <- ggplot2::ggplot(resamples, aes(x = value)) +
            ggplot2::stat_density(aes(ymax = ..density..), geom = "line") +
            ggplot2::geom_vline(data = trueCoefs, aes(xintercept = value), linetype = "dashed", size = 1) +
            ggplot2::facet_wrap(~ coef, scales = "free") +
            ggplot2::theme_bw();
        
        return(p);
    } else {
        respVar <- as.character(m$formula[[2]]);
        
        y <- unname(model.response(model.frame(x$models[[1]]), "numeric"));
        compParts <- reshape2::melt(lapply(x$models, function(m) {
            return(model.frame(m)[ , 2L])
        }))
        
        colnames(compParts) <- c("value", "part");
        compParts$part <- factor(compParts$part, levels = names(x$models))
        
        X <- cbind(y, compParts);
        
        p <- ggplot2::ggplot(X, ggplot2::aes(x = value, y = y)) +
            ggplot2::geom_point() +
            ggplot2::stat_smooth(method = complmrob.wrapper, complmrob.model = x, se = se, level = conf.level) +
            ggplot2::facet_grid(. ~ part, scales = "fixed") +
            ggplot2::ylab(respVar) +
            ggplot2::theme_bw();
        
        return(p);
    }
}


complmrob.wrapper <- function(formula, data, weights = NULL, complmrob.model) {    
    part <- names(complmrob.model$models)[as.integer(data$PANEL[1])]
    
    x <- complmrob.model$models[[part]]
    x$part <- part;
    class(x) <- "complmrob.part";
    return(x)
}

#' @import robustbase
predictdf.complmrob.part <- function(model, xseq, se, level) {
    getSeq <- function(x, length = 100) {
        # 		r <- range(x)
        # 		seq(from = r[1], to = r[2], length = length)
        rep.int(median(x), times = length)
    };
    
    origdata <- model$x[ , -which(colnames(model$x) == model$part), drop = FALSE];
    if(attr(model$terms, "intercept") == 1) {
        origdata <- origdata[ , -1, drop = FALSE];
    }
    preddata <- data.frame(part = xseq, apply(origdata, 2, getSeq, length = length(xseq)));
    
    colnames(preddata)[1] <- model$part;
    class(model) <- "lmrob";
    
    pred <- predict(model, newdata = preddata, interval = ifelse(se, "confidence", "none"), level = level);
    if(se == TRUE) {
        colnames(pred) <- c("y", "ymin", "ymax");
    } else {
        pred <- data.frame(y = pred);
    }
    return(data.frame(x = xseq, pred))
}