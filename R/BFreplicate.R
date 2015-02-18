#' Bayes Factors Using Bayesm 
#'
#' This function draws samples from a Wishart dist
#' @param mydata, alpha, Krange
#' @keywords bayes factor, bayesm
#' @export
#' @examples
#' # see bayesm documentation
#' sim1_n100<-sim1func(n=300) 
#' BF<-getBF(sim1_n100$Y, 1, c(1:6))
#' BF

BFreplicate<-function(sim, n=100, numRep=20, iters=2000){ 
	require('parallel')
	yrep<-lapply(rep(n, numRep),  function(x) simMe( sim, x))
	dtc<-min(detectCores(), numRep)
	BFrun<-mclapply(yrep, getBFk0,runlength=iters, mc.cores=dtc) 
	return(unlist(BFrun) )}