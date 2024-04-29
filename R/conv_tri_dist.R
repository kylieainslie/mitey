#' Convolution of the triangular distribution with the mixture component density (continuous case)
#'
#' We split the folded normal distribution for Primary-Secondary, Primary-Tertiary and Primary-Quaternary routes into two parts
#' + component 1: Co-Primary route
#' + component 2+3: Primary-Secondary route
#' + component 4+5: Primary-Tertiary route
#' + component 6+7: Primary-Quaternary route

conv_tri_dist <- function(x, sigma, r, mu, route){

  if(route == 1){
    f10 <- (2-2*x)*dhalfnorm(x,theta=sqrt(pi/2)/(sqrt(2)*sigma))
    f1lower<-	(x-r+1)*dhalfnorm(x,theta=sqrt(pi/2)/(sqrt(2)*sigma))
    f1upper<-	(r+1-x)*dhalfnorm(x,theta=sqrt(pi/2)/(sqrt(2)*sigma))
  }

  f20<-function(x,mu,sigma)      	(2-2*x)*dnorm(x,mean=mu,sd=sigma)
  f2lower<-function(x,r,mu,sigma) 	(x-r+1)*dnorm(x,mean=mu,sd=sigma)
  f2upper<-function(x,r,mu,sigma) 	(r+1-x)*dnorm(x,mean=mu,sd=sigma)

  f30<-function(x,mu,sigma)       	(2-2*x)*dnorm(x,mean=-mu,sd=sigma)
  f3lower<-function(x,r,mu,sigma) 	(x-r+1)*dnorm(x,mean=-mu,sd=sigma)
  f3upper<-function(x,r,mu,sigma) 	(r+1-x)*dnorm(x,mean=-mu,sd=sigma)

  f40<-function(x,mu,sigma)       	(2-2*x)*dnorm(x,mean=2*mu,sd=sqrt(2)*sigma)
  f4lower<-function(x,r,mu,sigma) 	(x-r+1)*dnorm(x,mean=2*mu,sd=sqrt(2)*sigma)
  f4upper<-function(x,r,mu,sigma) 	(r+1-x)*dnorm(x,mean=2*mu,sd=sqrt(2)*sigma)

  f50<-function(x,mu,sigma)       	(2-2*x)*dnorm(x,mean=-2*mu,sd=sqrt(2)*sigma)
  f5lower<-function(x,r,mu,sigma) 	(x-r+1)*dnorm(x,mean=-2*mu,sd=sqrt(2)*sigma)
  f5upper<-function(x,r,mu,sigma) 	(r+1-x)*dnorm(x,mean=-2*mu,sd=sqrt(2)*sigma)

  f60<-function(x,mu,sigma)       	(2-2*x)*dnorm(x,mean=3*mu,sd=sqrt(3)*sigma)
  f6lower<-function(x,r,mu,sigma) 	(x-r+1)*dnorm(x,mean=3*mu,sd=sqrt(3)*sigma)
  f6upper<-function(x,r,mu,sigma) 	(r+1-x)*dnorm(x,mean=3*mu,sd=sqrt(3)*sigma)

  f70<-function(x,mu,sigma)       	(2-2*x)*dnorm(x,mean=-3*mu,sd=sqrt(3)*sigma)
  f7lower<-function(x,r,mu,sigma) 	(x-r+1)*dnorm(x,mean=-3*mu,sd=sqrt(3)*sigma)
  f7upper<-function(x,r,mu,sigma) 	(r+1-x)*dnorm(x,mean=-3*mu,sd=sqrt(3)*sigma)
}
