	 
#' Simulate multivariate normal Mixture
#'
#' ...
#' @param  Means, covariance, Sample size, weights
#' @keywords Wishart
#' @export
#' @examples
#' #not run
SimMVN<-function(Means, Covs, N, P){
			TrueK<-length(P)
			r<-length(Means[[1]])
			Ys<-matrix(0,nrow = N, ncol = r)  
			Zs<-matrix(0,nrow = N, ncol = 1)  
			for (.n in 1:N){
			#draw Z according to probability
			.z<-sample(c(1:TrueK), size=1, prob=P)
			Zs[.n,]<-.z
			#generate random sample
			Ys[.n,]<-rmvnorm(1, Means[[.z]], Covs[[.z]])
			}
			list(Y=Ys, Z=Zs)}