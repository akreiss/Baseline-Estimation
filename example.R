## In this example file it is assumed that you have the following variables
## available. Examples for these variables are provided in the file
## example_data.RData. If you load the data-file first, the following code
## should run without errors.
## covariates1,       : Two covariate lists each corresponding to one data-set
##  covariates2 
## events1,           : The two event lists corresponding to the same two 
##  events2             data-sets
## T                  : End of observation horizons for the two data-sets (we
##                      suppose that they have the same length).
## alpha_covariates1, : Covariates for the baselines on the two data-sets.
##  alpha_covariates2

## In the example below, we suppose that the time t is given in hours and that
## the weight function w cuts of the first and last hour.


## Load functions
source('./functions.R')

## Load Libraries
library(Matrix)

################################################################################
## Kernel, weight function and their integrals #################################
################################################################################

## Kernel used for Nelson estimator (triangular kernel)
myK <- function(x) {
  return(pmax(1-abs(x),0))
}

## Integral of myK
myKintegrate <- function(a) {
  out <- rep(0,length(a))
  
  ind <- which(a>-1 & a<=0)
  if(length(ind)>0) {
    out[ind] <- a[ind]+1+0.5*(a[ind]^2-1)
  }
  
  ind <- which(a>0 & a<=1)
  if(length(ind)>0) {
    out[ind] <- 0.5+a[ind]-0.5*a[ind]^2
  }
  
  ind <- which(a>1)
  if(length(ind)>0) {
    out[ind] <- 1
  }
  
  return(out)
}

## Weight function w (constant equal to one on the interval [1,T-1])
myw <- function(t) {
  ind <- which(t>=1 & t<=T-1)
  out <- rep(0,length(t))
  out[ind] <- 1
  return(out)
}

## Integral of weight function
mywIntegrate <- function(a) {
  out <- rep(0,length(a))
  
  ind <- which(a>1 & a<=T-1)
  if(length(ind)>0) {
    out[ind] <- a[ind]-1
  }
  
  ind <- which(a>T-1)
  if(length(ind)>0) {
    out[ind] <- T-2
  }
  
  return(out)
}

## Compute K^{(2)}
integrand <- function(u,v) {
  return(myK(u+v)*myK(u))
}
integral1 <- function(v) {
  out <- rep(0,length(v))
  for(i in 1:length(v)) {
    out[i] <- integrate(integrand,lower=-1,upper=1,v=v[i])$value^2
  }
  return(out)
}
Ksq <- integrate(integral1,lower=0,upper=2)$value
################################################################################

## Compute Partial Likelihood Estimator Based on data-set 1
epp1 <- events_per_pair(covariates1,events1)
beta_cpl_out <- nlm(cpl,runif(p),print.level=2,covariates=covariates1,
                    events=events1,T=T,epp=epp1)
beta_cpl <- beta_cpl_out$estimate

## Compute parametric Estimator based on data-set 1
param_out <- nlm(cfl,rnorm(p+d),print.level=2,iterlim=1000,betadim=p,
                 covariates=covariates1,alpha_covariates=alpha_covariates1,
                 events=events1,T=T,epp=epp1)
theta_hat <- param_out$estimate[(p+1):(p+d)]

## Compute Nelson-Aalen Estimator based on data-set 2
t <- seq(from=0,to=T,length.out=1000)
h <- 0.5
alpha_hat_NP <- NAE(t,h,beta_cpl,covariates2,events2,T,myK)

## Get parametric Baseline based on data-set 2
alpha_hat_P <- alpha_eval(t,alpha_covariates2,theta_hat,T)

## Plot the two estimators
dev.new()
plot(t/24,alpha_hat_NP,type="l",xlab="Time in Days",ylab="",main="Baseline Estimates")
lines(t/24,alpha_hat_P,lty=2)
legend("topleft",lty=c(1,2),legend=c("NP","P"))


## Compute the test statistic
test <- baseline_test(alpha_hat_NP,alpha_hat_P,myw,t,events2,covariates2,
              beta_cpl,myK,h,T,myKintegrate,alpha_covariates2,theta_hat,
              mywIntegrate,Ksq)
