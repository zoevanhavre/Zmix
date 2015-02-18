 
#' Draw from Wishart 
#'
#' This function draws samples from a Wishart dist
#' @param v and s
#' @keywords Wishart
#' @export
#' @examples
#' rwish(2,1)


  rwish<-function (v, S) {
					    if (!is.matrix(S)) 
					        S <- matrix(S)
					    if (nrow(S) != ncol(S)) {
					        stop(message = "S not square in rwish().\n")
					    }
					    if (v < nrow(S)) {
					        stop(message = "v is less than the dimension of S in rwish().\n")
					    }
					    p <- nrow(S)
					    CC <- chol(S)
					    Z <- matrix(0, p, p)
					    diag(Z) <- sqrt(rchisq(p, v:(v - p + 1)))
					    if (p > 1) {
					        pseq <- 1:(p - 1)
					        Z[rep(p * pseq, pseq) + unlist(lapply(pseq, seq))] <- rnorm(p * 
					            (p - 1)/2)
					    }
					    return(crossprod(Z %*% CC))}
