 #' A dDirichlet Function
#'
#' density of dirichlet
#' @param x, Input parameter
#' @param alpha, Hyper-prior
#' @param log=False,  use log?
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))



 dDirichlet<-function (x, alpha, log = FALSE) {
				    dlog = lgamma(sum(alpha)) + sum((alpha - 1) * log(x)) - sum(lgamma(alpha))
				    result = ifelse(!log, exp(dlog), dlog)
				    return(result)
					}