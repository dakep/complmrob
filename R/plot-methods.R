globalVariables(c("..density..", "value", "ymin", "ymax"));

#' Diagnostic plots for the robust regression model with compositional covariates
#'
#' Plot the response or the model diagnostic plots for robust linear regression model with compositional
#' data
#'
#' The response plot shows the value on the first component of the orthonormal basis versus the response
#' and the fitted values. For the fitted values, the other components are set to the median of the values
#' in that direction. This usually causes aberrant predictions when plotting on the *percent* scale.
#'
#' For the model diagnostic plots see the details in the help file for \code{\link[robustbase]{plot.lmrob}}.
#' The model diagnostic plots are the same for all sub-models fit to the data transformed with the different
#' orthonormal basis.
#'
#' @param x the object returned by \code{\link{complmrob}}.
#' @param y ignored.
#' @param type one of \code{"response"} to plot the response or \code{"model"} to get the standard
#'      \code{\link[robustbase]{lmrob}} model diagnostic plots. Partial matching is performed, so any unique
#'      abbreviation of the two possible values is excepted (e.g., \code{"r"} for the response plot).
#' @param se should the confidence interval be shown in the response plot.
#' @param conf.level if the confidence interval is shown in the response plot, this parameter sets
#'      the level of the confidence interval.
#' @param scale should the x-axis in the response plot be in percentage or in the ILR-transformed scale?
#' @param theme the ggplot2 theme to use for the response plot.
#' @param pointStyle a list with style parameters for the points in the response plot (possible entries
#'      are \code{color}, \code{size}, \code{alpha}, and \code{shape}). If \code{color} and/or \code{shape} is a vector
#'      of length equal to the number of observations in the model, the points will be colored/shaped according
#'      to this vector.
#' @param lineStyle  list with style parameters for the smoothing lines in the response plot (possible entries
#'      are \code{color}, \code{width}, and \code{linetype})
#' @param seBandStyle a list with style parameters (\code{color} and \code{alpha}) for the confidence band (if \code{se} is \code{TRUE})
#' @param stack how the facets are laid out in the response plot. \code{"horizontal"} for side by side and \code{"vertical"}
#'      for on top of each other.
#' @param ... further arguments to the model diagnostic plot method (see \code{\link[robustbase]{plot.lmrob}} for details).
#'
#' @import ggplot2
#' @importFrom ggplot2 ggplot aes geom_point ylab xlab facet_grid rel theme_bw stat_smooth scale_x_continuous
#' @importFrom scales percent
#' @method plot complmrob
#' @export
#' @examples
#' data <- data.frame(lifeExp = state.x77[, "Life Exp"], USArrests[ , -3])
#' mUSArr <- complmrob(lifeExp ~ ., data = data)
#' plot(mUSArr)
#' plot(mUSArr, type = "model") # for the model diagnostic plots
plot.complmrob <- function(x, y = NULL, type = c("response", "model"), se = TRUE, conf.level = 0.95,
    scale = c("ilr", "percent"), theme = theme_bw(),
    pointStyle = list(color = "black", size = rel(1), alpha = 1, shape = 19),
    lineStyle = list(color = "grey20", width = rel(1), linetype = "solid"),
    seBandStyle = list(color = "gray80", alpha = 0.5),
    stack = c("horizontal", "vertical"), ...) {

    type <- match.arg(type);
    stack <- match.arg(stack);
    scale <- match.arg(scale);

    if(type == "model") {
        return(plot(x$models[[1]], ...));
    }
    respVar <- as.character(x$formula[[2]]);

    y <- unname(model.response(model.frame(x$models[[1]]), "numeric"));
    compParts <- lapply(x$models, function(m) {
        trX <- as.matrix(model.frame(model.frame(m)[ , -1L, drop = FALSE]));
        if(scale == "percent") {
            return(isomLRinv(trX, perc = TRUE)[ , 1L, drop = TRUE]);
        }
        return(trX[ , 1L]);
    });

    pointStyle <- c(pointStyle, list(color = "black", size = rel(1), alpha = 1, shape = 19));
    lineStyle <- c(lineStyle, list(color = "grey20", width = rel(1), linetype = "solid"))
    seBandStyle <- c(seBandStyle, list(color = "gray80", alpha = 0.5));

    pointStyle <- pointStyle[!duplicated(names(pointStyle))]
    lineStyle <- lineStyle[!duplicated(names(lineStyle))]
    seBandStyle <- seBandStyle[!duplicated(names(seBandStyle))]

    pointStyle$data <- data.frame(y = y, value = do.call(c, compParts),
                                  part = factor(rep.int(names(x$models),
                                                        rep.int(length(y), length(x$models))),
                                                names(x$models)));

    pointMapping <- list();
    shapeAndColorIdent <- identical(pointStyle$shape, pointStyle$color);

    if(length(pointStyle$color) == length(y)) {
        pointStyle$data$color <- rep.int(pointStyle$color, length(x$models));
        pointMapping$color <- pointMapping$fill <- quote(color);

        # Overriding has to be done twice as the list actually has two elements with the same name
        # (as with the above assignment the old item is just obfuscated by the new one)
        pointStyle$color <- NULL;
        pointStyle$color <- NULL;
    }

    if(length(pointStyle$shape) == length(y)) {
        if(shapeAndColorIdent) {
            pointMapping$shape <- quote(color);
        } else {
            pointStyle$data$shape <- rep.int(pointStyle$shape, length(x$models));
            pointMapping$shape <- quote(shape);
        }

        # Overriding has to be done twice as the list actually has two elements with the same name
        # (as with the above assignment the old item is just obfuscated by the new one)
        pointStyle$shape <- NULL;
        pointStyle$shape <- NULL;
    }

    pointStyle$mapping = switch(1L + length(pointMapping),
                                NULL,
                                do.call(aes, pointMapping));

    # Build prediction data
    nPred <- 1000L

    partPred <- lapply(x$models, function (model) {
        part <- attr(terms(model), 'term.labels')[[1L]]
        partColumn <- which(colnames(model$x) == part)

        origdata <- model$x[ , -partColumn, drop = FALSE]
        if(attr(model$terms, "intercept") == 1) {
            origdata <- origdata[ , -1L, drop = FALSE];
        }

        xRange <- range(model$x[ , part])
        xRange <- diff(xRange) * c(-0.05, 0.05) + xRange; # Expand limits by 5 percent
        xseq <- seq.int(xRange[[1]], xRange[[2]], length.out = nPred)

        predDataMatrix <- matrix(rep(colMeans(origdata), each = nPred),
                                 nrow = nPred)

        suppressWarnings(preddata <- data.frame(part = xseq, predDataMatrix));

        colnames(preddata)[[1]] <- part;
        colnames(preddata)[-1] <- colnames(origdata);
        class(model) <- "lmrob";

        pred <- predict(model, newdata = preddata, level = conf.level,
                        interval = ifelse(se, "confidence", "none"));
        if(se == TRUE) {
            colnames(pred) <- c("y", "ymin", "ymax");
        } else {
            pred <- data.frame(y = pred);
        }
        if(scale == 'percent') {
            Zinv <- isomLRinv(as.matrix(preddata));
            return(data.frame(part = part, value = Zinv[ , 1], pred))
        } else {
            return(data.frame(part = part, value = xseq, pred))
        }
    })

    partPred <- do.call(rbind, partPred)

    p <- ggplot(data = partPred, aes(x = value, y = y)) +
        do.call(geom_point, pointStyle, quote = TRUE) +
        geom_line(size = lineStyle$width, color = lineStyle$color) +
        geom_ribbon(mapping = aes(ymin = ymin, ymax = ymax),
                    fill = seBandStyle$color, alpha = seBandStyle$alpha) +
        ylab(respVar) +
        xlab("Share") +
        theme +
        switch(stack,
               vertical = facet_grid(part ~ ., scales = "fixed"),
               facet_grid(. ~ part, scales = "fixed"))

    if(scale == "percent") {
        p <- p + scale_x_continuous(labels = scales::percent)
    }

    return(p);
}

#' Plot the distribution of the bootstrap estimates
#'
#' Plot the distribution of the bootstrap estimates and the confidence intervals for the estimates
#'
#' @param x the object returned by a call to the \code{\link{bootcoefs}} function.
#' @param y ignored.
#' @param conf.level the level of the confidence interval that is plotted as shaded region under the
#'      density estimate.
#' @param conf.type the confidence interval type, see \code{\link[boot]{boot.ci}} for details.
#' @param kernel the kernel used for density estimation, see \code{\link[stats]{density}} for details.
#' @param adjust see \code{\link[stats]{density}} for details.
#' @param which which parameters to plot
#' @param theme the ggplot2 theme to use for the plot.
#' @param confStyle a list with style parameters for the confidence region below the density estimate (possible entries
#'      are \code{color}, and \code{alpha})
#' @param estLineStyle a list with style parameters for the line at the original parameter estimate (possible entries
#'      are \code{color}, \code{width}, \code{alpha}, and \code{linetype})
#' @param densityStyle a list with style parameters for the line of the density estimate (possible entries
#'      are \code{color}, \code{width}, \code{alpha}, and \code{linetype})
#' @param ... ignored
#'
#' @import ggplot2
#' @importFrom ggplot2 geom_area stat_density geom_segment scale_y_continuous ggtitle element_blank theme facet_wrap
#' @method plot bootcoefs
#' @seealso \code{\link[=confint.bccomplmrob]{confint}} to get the numerical values for the confidence intervals
#' @export
#' @examples
#' data <- data.frame(lifeExp = state.x77[, "Life Exp"], USArrests[ , -3])
#' mUSArr <- complmrob(lifeExp ~ ., data = data)
#' bc <- bootcoefs(mUSArr, R = 200) # this can take some time
#' plot(bc) # for the model diagnostic plots
plot.bootcoefs <- function(x, y = NULL, conf.level = 0.95, conf.type = "perc", kernel = "gaussian", adjust = 1,
    which = "all", theme = theme_bw(), confStyle = list(color = "#56B4E9", alpha = 0.4),
    estLineStyle = list(color = "black", width = rel(1), alpha = 1, linetype = "dashed"),
    densityStyle = list(color = "black", width = rel(0.5), alpha = 1, linetype = "solid"), ...) {

    allCoefNames <- names(coef(x$model));

    if(length(which) == 1 && which == "all") {
        which <- seq_along(allCoefNames);
    } else if(is.character(which)) {
        which <- na.omit(match(which, allCoefNames));
    }

    replicates <- list();
    if(class(x$bootres) == "boot") {
        tsub <- x$bootres$t[ , which, drop = FALSE];
        nc <- ncol(tsub);
        replicates <- split(tsub, rep.int(seq_len(nc), rep.int(nrow(tsub), nc)));
        names(replicates) <- allCoefNames[which];
    } else {
        replicates <- lapply(x$bootres[which], function(bo) {
            bo$t[ , 1, drop = TRUE]
        });
    }

    confStyle <- c(confStyle, list(color = "#56B4E9", alpha = 0.4));
    estLineStyle <- c(estLineStyle, list(color = "black", width = rel(1), alpha = 1, linetype = "dashed"));
    densityStyle <- c(densityStyle, list(color = "black", width = rel(0.5), alpha = 1, linetype = "solid"));

    replicatesLong <- na.omit(data.frame(x = do.call(c, replicates),
        coef = rep.int(names(replicates), sapply(replicates, length))));

    ci <- confint(x, level = conf.level, type = conf.type)[which, , drop = FALSE];
    ci <- split(ci, rep.int(seq_len(nrow(ci)), ncol(ci)));

    replicatesDE <- mapply(function(t, ci) {
        de <- density(t, from = ci[1], to = ci[2], kernel = kernel, adjust = adjust, na.rm = TRUE);
        return(data.frame(x = de$x, y = de$y));
    }, replicates, ci, SIMPLIFY = FALSE);

    replicatesDens <- do.call(rbind, replicatesDE);
    replicatesDens$coef <- factor(rep.int(names(replicatesDE), times = sapply(replicatesDE, nrow)), levels = names(replicatesDE));

    trueCoefs <- data.frame(coef = names(x$model$coefficients[which]), x = unname(x$model$coefficients[which]),
        y = c(by(replicatesDens, replicatesDens$coef, function(x) { return(max(x$y)) })));

    p <- ggplot(replicatesDens, aes(x = x)) +
        geom_area(mapping = aes(y = y), fill = confStyle$color, alpha = confStyle$alpha) +
        suppressWarnings(stat_density(data = replicatesLong,
                                      mapping = aes(ymax = ..density..), geom = "line",
                                      kernel = kernel, adjust = adjust,
                                      color = densityStyle$color,
                                      size = densityStyle$width,
                                      alpha = densityStyle$alpha,
                                      linetype = densityStyle$linetype)) +
        geom_segment(data = trueCoefs, aes(x = x, xend = x, y = 0, yend = 1.02 * y), linetype = estLineStyle$linetype,
            size = estLineStyle$width, color = estLineStyle$color, alpha = estLineStyle$alpha) +
        scale_y_continuous(expand = c(0, 0)) +
        ggtitle(sprintf("Distribution of %d bootstrap estimates with %s confidence interval", x$R,
            format.perc(conf.level, 2))) +
        xlab(NULL) +
        ylab(NULL) +
        theme +
        theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()
        );

    if(nlevels(replicatesDens$coef) > 1) {
        p <- p +  facet_wrap(~ coef, scales = "free");
    }

    return(p);
}
