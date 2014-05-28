#' Isometric log-ratio transformation for compositional data
#' 
#' Projects the p-dimensional compositional data on the (p-1)-dimensional simplex isometrically
#' 
#' @param x a numeric vector of length \code{p} or a numeric matrix with \code{p} columns
#' @param comp the component to use as the first compositional part
#' @return a numeric matrix with \code{(p-1)} columns with the transformed values. The name of the first
#'     	column is the name of the first part (the other names are according to the order of the
#' 		columns in the given matrix \code{x})
#'
#' @export
isomLR <- function(x, comp = 1) {
    if(is.character(comp)) {
        comp <- which(colnames(x) == comp);
    }
    
    ## assert x is numeric
    if(!is.numeric(x)) {
        stop("x must be a numeric");
    }
    
    if(!is.matrix(x)) {
        x <- matrix(x, nrow = 1)
    }
    
    x <- cbind(x[ , comp, drop = FALSE], x[ , -comp, drop = FALSE])
    
    D <- ncol(x);
    
    sc <- D - seq_len(D - 1);
    sc <- sqrt(sc / (sc + 1))
    
    ilr <- function(x) {
        denum <- rev(exp(cumsum(rev(log(x))) / seq_len(D)))[-1];
        sc * log(x[-D] / denum)
    };
    
    ret <- apply(x, 1, ilr);
    if(!is.matrix(ret)) {
        ret <- matrix(ret, ncol = 1);
    } else {
        ret <- t(ret);
    }
    colnames(ret) <- colnames(x)[-D]
    return(ret);
}

#' Inverse Isometric log-ratio transformation for compositional data
#' 
#' Projects the isometric log-ratio transformed (p-1) dimensional data back to the p-dimensional space
#' 
#' @param x a numeric vector of length \code{p-1} or a numeric matrix with \code{p-1} columns
#' @param perc should the result be a matrix with percentual shares (default \code{TRUE}).
#' @return a numeric matrix with \code{p} columns with the transformed values. The values in the matrix
#'      are not on the original scale, but the percentual shares are equal
#'
#' @export
isomLRinv <- function(x, perc = TRUE) {
    ## assert x is numeric
    if(!is.numeric(x)) {
        stop("x must be a numeric");
    }
    
    if(!is.matrix(x)) {
        x <- matrix(x, nrow = 1)
    }

    D <- as.integer(ncol(x) + 1);

    Zinv <- matrix(NA_real_, ncol = D, nrow = nrow(x));

    Zinv[ , 1] <- exp((sqrt(D - 1) / sqrt(D)) * x[ , 1, drop = TRUE]);

    normalizer <- D - seq_len(D - 1)
    norm <- -1 / sqrt((normalizer + 1) * normalizer);

    Y <- apply(x, 1, function(z) {
        norm * z
    });

    Zinv[ , -c(1, D)] <- do.call(cbind, lapply(1 + seq_len(D - 2), function(i) {
        exp(colSums(Y[seq_len(i - 1), , drop = FALSE]) + (sqrt(D - i) / sqrt(D - i + 1)) * x[, i])
    }));

    Zinv[ , D] <- exp(colSums(Y));
    
    if(perc == TRUE) {
        Zinv <- Zinv / rowSums(Zinv);
    }

    return(Zinv);
}
