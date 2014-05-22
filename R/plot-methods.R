#' Diagnostic plots for the robust regression model with compositional covariats
#'
#' @param x The complmrob object to plot
#' @param y ignored
#' @param type one of \code{"response"} to plot the response or \code{"model"} to get the standard
#'      lmrob model diagnostic plots.
#' @param se should the confidence interval be shown in the response plot
#' @param conf.level if the confidence interval is shown in the response plot, this parameter sets
#'      the level of the confidence interval
#' @param ... ignored
#' 
#' @import ggplot2
#' @method plot complmrob
#' @export
plot.complmrob <- function(x, y = NULL, type = c("response", "model"), se = TRUE, conf.level = 0.95, ...) {
    type <- match.arg(type);

    if(type == "model") {
        plot(x$models[[1]]);
    } else {
        respVar <- as.character(x$formula[[2]]);

        y <- unname(model.response(model.frame(x$models[[1]]), "numeric"));
        compParts <- lapply(x$models, function(m) {
            return(model.frame(m)[ , 2L, drop = TRUE])
        });

        X <- data.frame(y = y, value = do.call(c, compParts),
            part = factor(rep.int(names(x$models), rep.int(length(y), length(x$models))), names(x$models)));

        p <- ggplot2::ggplot(X, ggplot2::aes(x = value, y = y)) +
            ggplot2::geom_point() +
            ggplot2::stat_smooth(method = complmrob.wrapper, complmrob.model = x, se = se, level = conf.level) +
            ggplot2::facet_grid(. ~ part, scales = "fixed") +
            ggplot2::ylab(respVar) +
            ggplot2::theme_bw();

        return(p);
    }
}

#' Plot the distribution of the bootstrap estimates
#'
#' @param x the object to plot, as returned by the \code{\link{bootcoefs}} function
#' @param y ignored
#' @param conf.level the level of the confidence interval that is plotted as shaded region under the
#'      density estimate.
#' @param conf.type the confidence interval type, see \code{\link{boot.ci}} for details
#' @param kernel the kernel used for density estimation, see \code{\link{density}} for details.
#' @param adjust see \code{\link{density}} for details.
#' @param ... ignored
#' 
#' @import ggplot2
#' @method plot bootcoefs
#' @export
plot.bootcoefs <- function(x, y = NULL, conf.level = 0.95, conf.type = "perc", kernel = "gaussian", adjust = 1, ...) {
    replicates <- list();
    if(class(x$bootres) == "boot") {
        nc <- ncol(x$bootres$t);
        replicates <- split(x$bootres$t, rep.int(seq_len(nc), rep.int(nrow(x$bootres$t), nc)));
        names(replicates) <- names(x$model$coefficients);
    } else {
        replicates <- lapply(x$bootres, function(bo) {
            bo$t[ , 1, drop = TRUE]
        });
    }

    replicatesLong <- na.omit(data.frame(x = do.call(c, replicates),
        coef = rep.int(names(replicates), sapply(replicates, length))));
    
    ci <- confint(x, level = conf.level, type = conf.type);
    ci <- split(ci, rep.int(seq_len(nrow(ci)), ncol(ci)));
    
    replicatesDE <- mapply(function(t, ci) {
        de <- density(t, from = ci[1], to = ci[2], kernel = kernel, adjust = adjust, na.rm = TRUE);
        return(data.frame(x = de$x, y = de$y));
    }, replicates, ci, SIMPLIFY = FALSE);

    replicatesDens <- do.call(rbind, replicatesDE);
    replicatesDens$coef <- factor(rep.int(names(replicatesDE), times = sapply(replicatesDE, nrow)), levels = names(replicatesDE));

    trueCoefs <- data.frame(coef = names(x$model$coefficients), x = unname(x$model$coefficients));
    
    p <- ggplot2::ggplot(replicatesDens, aes(x = x)) +
        ggplot2::geom_area(mapping = aes(y = y), fill = "#56B4E9", alpha = 0.4) +
        ggplot2::stat_density(data = replicatesLong, mapping = aes(ymax = ..density..), geom = "line",
            kernel = kernel, adjust = adjust) +
        ggplot2::geom_vline(data = trueCoefs, aes(xintercept = x), linetype = "dashed", size = 1) +
        ggplot2::facet_wrap(~ coef, scales = "free") +
        ggplot2::scale_y_continuous(expand = c(0, 0)) +
        ggplot2::ggtitle(sprintf("Distribution of %d bootstrap estimates with %s confidence interval", x$R,
            format.perc(conf.level, 2))) +
        ggplot2::xlab(NULL) +
        ggplot2::ylab(NULL) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            axis.text.y = ggplot2::element_blank(), 
            axis.ticks.y = ggplot2::element_blank());
    
    return(p);
}

complmrob.wrapper <- function(formula, data, weights = NULL, complmrob.model) {
    part <- names(complmrob.model$models)[as.integer(data$PANEL[1])]

    x <- complmrob.model$models[[part]]
    x$part <- part;
    class(x) <- "complmrob.part";
    return(x)
}

#' This function is used by ggplot2 to predict the values for a \code{complmrob} model
#' 
#' The function is exported solely because the ggplot2 method \code{predictdf} is not exported
#' and thus the generics system would not be working otherwise.
#' 
#' @param model the complmrob.part model the prediction should be done for
#' @param xseq the sequence of x values to predict for
#' @param se should the confidence interval be returned as well
#' @param level the level of the confidence interval (if any)
#' @import robustbase
#' @export
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
