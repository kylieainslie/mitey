# R-code for estimating the mean serial interval 
# method from Vink et al. 2014 AJE https://doi.org/10.1093/aje/kwu209

library(fdrtool)

# Note that you have to define your data in vector format
# data<-c()
# see Web Appendix 3 for simulated data and the data sets from the manuscript

# function for weighted variance, to be used in E-step
# weighted variance
weighted.var <- function(x, w, na.rm = FALSE) {
  if (na.rm) {
    w <- w[i <- !is.na(x)]
    x <- x[i]
  }
  sum.w <- sum(w)
  (sum.w*(sum(w*(x-weighted.mean(x,w))^2))) / (sum.w^2 - sum(w^2))
}

# mixture with 7 components
# we split the folded normal distribution for the PS, PT and PQ route into two parts
# component 1: CP route
# component 2+3: PS route
# component 4+5: PT route
# component 6+7: PQ route

# convolution of the triangular distribution with the mixture component density
# contiunous case
f10<-function(x,sigma)          	(2-2*x)*dhalfnorm(x,theta=sqrt(pi/2)/(sqrt(2)*sigma))
f1lower<-function(x,r,sigma)    	(x-r+1)*dhalfnorm(x,theta=sqrt(pi/2)/(sqrt(2)*sigma))
f1upper<-function(x,r,sigma)    	(r+1-x)*dhalfnorm(x,theta=sqrt(pi/2)/(sqrt(2)*sigma))

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


# discretization
p10<-function(d,sigma)        	integrate(f = f10, lower = d, upper = (d+1), sigma = sigma)
p1lower<-function(d,sigma)    	integrate(f=f1lower, lower=(d-1), upper=d, r = d, sigma = sigma)
p1upper<-function(d,sigma)   	 integrate(f=f1upper, lower=d, upper=(d+1), r = d, sigma = sigma)

p20<-function(d,mu,sigma)     	integrate(f = f20, lower = d, upper = (d+1), mu = mu, sigma = sigma)
p2lower<-function(d,mu,sigma) 	integrate(f=f2lower, lower=(d-1), upper=d, r = d, mu = mu, sigma = sigma)
p2upper<-function(d,mu,sigma) 	integrate(f=f2upper, lower=d, upper=(d+1), r = d, mu = mu, sigma = sigma)

p30<-function(d,mu,sigma)     	integrate(f = f30, lower = d, upper = (d+1), mu = mu, sigma = sigma)
p3lower<-function(d,mu,sigma) 	integrate(f=f3lower, lower=(d-1), upper=d, r = d, mu = mu, sigma = sigma)
p3upper<-function(d,mu,sigma) 	integrate(f=f3upper, lower=d, upper=(d+1), r = d, mu = mu, sigma = sigma)

p40<-function(d,mu,sigma)    	integrate(f = f40, lower = d, upper = (d+1), mu = mu, sigma = sigma)
p4lower<-function(d,mu,sigma) 	integrate(f=f4lower, lower=(d-1), upper=d, r = d, mu = mu, sigma = sigma)
p4upper<-function(d,mu,sigma) 	integrate(f=f4upper, lower=d, upper=(d+1), r = d, mu = mu, sigma = sigma)

p50<-function(d,mu,sigma)    	integrate(f = f50, lower = d, upper = (d+1), mu = mu, sigma = sigma)
p5lower<-function(d,mu,sigma) 	integrate(f=f5lower, lower=(d-1), upper=d, r = d, mu = mu, sigma = sigma)
p5upper<-function(d,mu,sigma) 	integrate(f=f5upper, lower=d, upper=(d+1), r = d, mu = mu, sigma = sigma)

p60<-function(d,mu,sigma)    	integrate(f = f60, lower = d, upper = (d+1), mu = mu, sigma = sigma)
p6lower<-function(d,mu,sigma) 	integrate(f=f6lower, lower=(d-1), upper=d, r = d, mu = mu, sigma = sigma)
p6upper<-function(d,mu,sigma) 	integrate(f=f6upper, lower=d, upper=(d+1), r = d, mu = mu, sigma = sigma)

p70<-function(d,mu,sigma)    	integrate(f = f70, lower = d, upper = (d+1), mu = mu, sigma = sigma)
p7lower<-function(d,mu,sigma) 	integrate(f=f7lower, lower=(d-1), upper=d, r = d, mu = mu, sigma = sigma)
p7upper<-function(d,mu,sigma) 	integrate(f=f7upper, lower=d, upper=(d+1), r = d, mu = mu, sigma = sigma)

j<-length(data) 

# EM algorithm: 
tau1<-numeric(j)    #mixture weights for (coprimary, coprimary) pairs
tau2<-numeric(j)    #mixture weights for (primary, secondary) pairs
tau3<-numeric(j)    #mixture weights for (secondary, primary) pairs 
tau4<-numeric(j)    #mixture weights for (primary, tertiary) pairs
tau5<-numeric(j)    #mixture weights for (tertiary, primary) pairs
tau6<-numeric(j)    #mixture weights for (primary, quaternary) pairs
tau7<-numeric(j)    #mixture weights for (quaternary, primary) pairs

# plausible starting values/educated guess
mu<-mean(data)
sigma<-sd(data)

# number of iterations
N<-50 

# E-step

# calculate the absolute probability of interval belonging to a component
for(k in 1:N){
  for(l in 1:j){
    if(data[l]==0){
      d1<-p10(data[l],sigma)[[1]]
      d2<-p20(data[l],mu,sigma)[[1]]
      d3<-p30(data[l],mu,sigma)[[1]]
      d4<-p40(data[l],mu,sigma)[[1]]
      d5<-p50(data[l],mu,sigma)[[1]]
      d6<-p60(data[l],mu,sigma)[[1]]
      d7<-p70(data[l],mu,sigma)[[1]]
    }
    else{
      d1<-p1lower(data[l],sigma)[[1]]+p1upper(data[l],sigma)[[1]]
      d2<-p2lower(data[l],mu,sigma)[[1]]+p2upper(data[l],mu,sigma)[[1]]
      d3<-p3lower(data[l],mu,sigma)[[1]]+p3upper(data[l],mu,sigma)[[1]]
      d4<-p4lower(data[l],mu,sigma)[[1]]+p4upper(data[l],mu,sigma)[[1]]
      d5<-p5lower(data[l],mu,sigma)[[1]]+p5upper(data[l],mu,sigma)[[1]]
      d6<-p6lower(data[l],mu,sigma)[[1]]+p6upper(data[l],mu,sigma)[[1]]
      d7<-p7lower(data[l],mu,sigma)[[1]]+p7upper(data[l],mu,sigma)[[1]]
    }
    
    # then calculate the relative probability of a data point to belong to one of
    # the components
    dummy<-d1+d2+d3+d4+d5+d6+d7
    tau1[l]<-d1/dummy
    tau2[l]<-d2/dummy
    tau3[l]<-d3/dummy
    tau4[l]<-d4/dummy
    tau5[l]<-d5/dummy
    tau6[l]<-d6/dummy
    tau7[l]<-d7/dummy
  }
  
  # now calculate the weights for each of the components
  w1<-sum(tau1)/j
  w2<-sum(tau2)/j
  w3<-sum(tau3)/j
  w4<-sum(tau4)/j
  w5<-sum(tau5)/j
  w6<-sum(tau6)/j
  w7<-sum(tau7)/j
  
  # M-step
  # estimates for the mean and standard deviation of the primary-secondary 
  # infection component can be calculated directly
  mu<-weighted.mean(data,tau2)
  sigma<-sqrt(weighted.var(data, tau2))
  print(c(mu,sigma))
}

