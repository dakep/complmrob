#' Isometric log-ratio transformation for compositional data
#' 
#' Projects the p-dimensional compositional data on the (p-1)-dimensional simplex isometrically
#' 
#' @param x a numeric vector of length p or a numeric matrix with p columns
#' @param comp the component to use as the first compositional part
#' @return a numeric matrix with (p-1) columns with the transformed values. The name of the first
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