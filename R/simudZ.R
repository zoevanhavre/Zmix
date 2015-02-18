	 
#' Simulate Univariate normal Mixture
#'
#' ...
#' @param  Means, covariance, Sample size, weights
#' @keywords Wishart
#' @export
#' @examples
#' #not run

simudZ <-function(n,mu,sig,p,k){
			Z = sample((1:k),n,prob = p, replace = T)    # randomly select a labelling (Z) from the number of mixtures k with probability p
			Y = rnorm(n,mean =mu[Z],sd=sqrt(sig[Z]))       # simulate random normal variable given mean of mixture chosen 
			list(Y = Y, Z = Z)   }