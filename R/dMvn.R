 
#' Multivariate normal density, a bit faster
#'
#' This function draws samples from a Wishart dist
#' @param x mu Sigma
#' @keywords multivariate
#' @export
#' @examples
#' # dMvn

	#dMvn <- function(X,mu,Sigma) {
	#	    k <- length(X)
	#	    rooti <- backsolve(chol(Sigma),diag(k))
	#	    #    quads <- colSums((crossprod(rooti,(t(X)-mu)))^2)
	#	    #    quads <- colSums((crossprod(rooti,(X-mu)))^2)
	#	    quads <-tryCatch({ colSums((crossprod(rooti,(  X -mu)))^2)  }, error = function(e) {  colSums((tcrossprod(rooti,(  X -mu)))^2) })
	#	    return(exp(-(k/2)*log(2*pi) + sum(log(diag(rooti))) - .5*quads))
	#	    }
	dMvn <- function(X,mu,Sigma) {
		    k <- ncol(X)
			rooti <- backsolve(chol(Sigma),diag(k))
			quads <- colSums((crossprod(rooti,(t(X)-mu)))^2)
			return(exp(-(k/2)*log(2*pi) + sum(log(diag(rooti))) - .5*quads))}