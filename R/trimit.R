#' trim univariate gibbs output
#'
#' This function draws samples from a Wishart dist
#' @param stuff
#' @keywords stuff
#' @export
#' @examples
#' #nope


trimit<-function(Out=Out, nEnd=EndSize){
						yo<-length(Out$Bigmu)                                       #number of chains
						nmax<-length(Out$Loglike)
						mu<-Out$Bigmu[[yo]][c(nmax-nEnd+1):nmax,]    
						sig<-Out$Bigsigma[[yo]][c(nmax-nEnd+1):nmax,]
						ps<-Out$Bigp[[yo]][c(nmax-nEnd+1):nmax,]
						Loglike<-Out$Loglike[c(nmax-nEnd+1):nmax]
						zs<-Out$Zs[[yo]][,c(nmax-nEnd+1):nmax]
						SteadyScore<-Out$SteadyScore$K0[c(nmax-nEnd+1):nmax]
						list(Mu = mu,Sig=sig, Ps = ps, Loglike=Loglike, SteadyScore=SteadyScore, Zs=zs, YZ=Out$YZ)	
						}
		