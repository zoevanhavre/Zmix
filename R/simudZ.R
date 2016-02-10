	 
#' Simulate data from a Univariate normal Mixture
#'
#' ...
#' @param  n: sample size to simulate
#' @param  mu: vector of mixture means
#' @param  sig: vector of mixture component variances
#' @param  p: vector of mixture weights
#' @param  k: Number of components.
#' @keywords Simulate, Gaussian, Univariate
#' @export
#' @examples simudZ(n=10, mu=c(-1,4,10), sig=c(1,1,2), p=c(.2, .7, .1), k=3)

simudZ <-function(n,mu,sig,p,k){
			Z = sample((1:k),n,prob = p, replace = T)    # randomly select a labelling (Z) from the number of mixtures k with probability p
			Y = rnorm(n,mean =mu[Z],sd=sqrt(sig[Z]))       # simulate random normal variable given mean of mixture chosen 
			list(Y = Y, Z = Z)   }