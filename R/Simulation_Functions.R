#' Sim 1 function
#'
#' This function draws samples from a simulation
#' @param n
#' @keywords simulation
#' @export
#' @examples
#' #nope
sim1func<-function(n){
			S1<-list(mu = c(-1,10,4),sig=c(0.5,0.5,3),p=c(0.5,0.3,0.2),k=3) 
			simudZ(mu=S1$mu, sig=S1$sig, n=n, p=S1$p, k=S1$k) }
#' Sim 2 function
#'
#' This function draws samples from a simulation
#' @param n
#' @keywords simulation
#' @export
#' @examples
#' #nope			
		sim2func<-function(n){
			S2<-list(mu = c(1,1),sig=c(10,1),p=c(0.5,0.5),k=2) 
			simudZ(mu=S2$mu, sig=S2$sig, n=n, p=S2$p, k=S2$k)  }
#' Sim 3 function
#'
#' This function draws samples from a simulation
#' @param n
#' @keywords simulation
#' @export
#' @examples
#' #nope
		sim3func<-function(n){ 
			S3<-list(     "Means"=list( c(3,5), c(1,5)),
			"Covs"= list( matrix(c(10,0.1,0.1,0.5), nrow=2), 
			matrix(c(0.5,0.1,0.1,10), nrow=2)),
			"P"=c(0.7,0.3))
			SimMVN(Means=S3$Means, Covs=S3$Covs, N=n, P=S3$P)  }
#' Sim 4 function
#'
#' This function draws samples from a simulation
#' @param n
#' @keywords simulation
#' @export
#' @examples
#' #nope
		sim4func<-function(n){
			S4<-list( "Means"=list(   c( 5,2, 15,1), 
			c( 10,6, 5,10),               
			c( 10,-1,5,10)),
			"Covs"= list(matrix(c(10,1,1,1, 1,10,1,1, 1,1,10,1, 1,1,1,10), nrow=4), 
			matrix(c(1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1), nrow=4),    
			matrix(c(1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1), nrow=4)),
			"P"=  c(0.85, 0.1, 0.05)  )
			SimMVN(Means=S4$Means, Covs=S4$Covs, N=n, P=S4$P)  }
#' Sim 5 function
#'
#' This function draws samples from a simulation
#' @param n
#' @keywords simulation
#' @export
#' @examples
#' #nope
sim5func<-function(n){
			S5<-list(mu = c(6,10, 20),sig=c(1,1,0.5),p=c(0.6, 0.39,0.01),k=3) 
			simudZ(mu=S5$mu, sig=S5$sig, n=n, p=S5$p, k=S5$k) }


#' Sim 6 function
#'
#' This function draws samples from a simulation
#' @param n
#' @keywords simulation
#' @export
#' @examples
#' #nope


sim6func<-function(n){
			S6<-list(mu = c(15, 7, 1),sig=c(1, 1, 1),p=c(0.4,  0.35, 0.25),k=3) 
			simudZ(mu=S6$mu, sig=S6$sig, n=n, p=S6$p, k=S6$k) }



#' Sim Big func
#'
#' This function ...
#' @param n
#' @keywords simulation
#' @export
#' @examples
#' #nope


simMe<-function(whichOne=1, n){
	if(whichOne==1) return(sim6func(n))
	if(whichOne==2) return(sim1func(n))
	if(whichOne==3) return(sim2func(n))
	if(whichOne==4) return(sim5func(n))
	if(whichOne==5) return(sim2EASYfunc(n))
}



#' Sim 2 k EASY
#'
#' This function draws samples from a simulation
#' @param n
#' @keywords simulation
#' @export
#' @examples
#' #nope			
		sim2EASYfunc<-function(n){
			S2<-list(mu = c(-5,5),sig=c(1,2),p=c(0.6,0.4),k=2) 
			simudZ(mu=S2$mu, sig=S2$sig, n=n, p=S2$p, k=S2$k)  }