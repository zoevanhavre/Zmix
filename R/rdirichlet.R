#' A Dirichlet Function
#'
#' This function draws samples from q dirichlet distriution with hyperpars m.
#' @param m Vector of Alpha values (hyperparameter of dist)
#' @keywords dirichlet
#' @export
#' @examples
#' rdirichlet()
 

rdirichlet<-function(m=1,par){
			k=length(par)         # for two mixture model k=2
			mat=matrix(0,m,k)       # makes an empty 1x2 matrix to store values
			for (i in 1:m)          # at this stage m is 1, when is it not 1?
			{
			sim=rgamma(k,shape=par,scale=1)
			mat[i,]=sim/sum(sim)      # simulated gamma scaled to sum to 1 and stored in matrix 
			}
			mat               # matrix outputed
			}