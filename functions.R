## Computes the negative value of the Cox Partial Likelihood. This functions
## needs to be minimized in beta in order to obtain the partial maximum likeli-
## hood estimator
## Input:
##  beta       - Vector of beta values at which the likelihood is supposed to
##               be computed.
##  covariates - List of covariates, see documentation.pdf for details
##  events     - Matrix of events, see documentation.pdf for details
##  T          - Endpoint of observation interval [0,T]
##  epp        - List of events per pair, this can be computed once by using the
##               function events_per_pair below. If this is set to NA (the de-
##               fault), events_per_pair is automatically called. This takes
##               time and should not be done for repeated calls.
## Output:
##  The negative value of the partial likelihood together with its gradient and
##  the hessian matrix (as attributes "gradient" and "hessian", respectively).
cpl <- function(beta,covariates,events,T,epp=NA) {
  ## Output Variables
  LL <- 0
  grad <- rep(0,length(beta))
  hess <- matrix(0,ncol=length(beta),nrow=length(beta))
  
  ## Compute epp unless given
  if(is.na(epp)==TRUE) {
    epp <- events_per_pair(covariates,events)
  }
  
  ## Go through each covariate element
  for(k in 1:length(covariates$time)) {
    ## Find endpoints of segment
    begin_of_segment <- covariates$time[k]
    if(k==length(covariates$time)) {
      end_of_segment <- T
    } else {
      end_of_segment <- covariates$time[k+1]
    }
    
    ## Get number of events in the time frame of interest (OPTIMIEREN)
    noe <- epp[[k]]
    
    ## Compute first integral in partial likelihood over time segment
    I1 <- rowSums(mapply(function(a,b) {return(a*b)},covariates$covars[[k]],noe))
    first_integral <- beta%*%I1
    
    ## Compute second integral
    S <- as.numeric(Reduce('+',lapply(covariates$covars[[k]],function(a) {return(exp(beta%*%a))})))
    Sd <- Reduce('+',lapply(covariates$covars[[k]],function(a) {return(a*as.numeric(exp(beta%*%a)))}))
    Sh <- Reduce('+',lapply(covariates$covars[[k]],function(a) {return(a%*%t(a)*as.numeric(exp(beta%*%a)))}))
    second_integral <- log(S)*sum(noe)
    
    
    ## Update partial likelihood and derivatives
    LL <- LL+first_integral-second_integral
    grad <- grad+I1-Sd/S*sum(noe)
    hess <- hess-(Sh/S-Sd%*%t(Sd)/S^2)*sum(noe)
  }

  output <- -LL
  attr(output,"gradient") <- -grad
  attr(output,"hessian") <- -hess
  
  return(output)
}


## Computest the negative value of the full likelihood and its derivatives.
## Input:
##  par              - Parameters at which the likelihood shall be evaluated,
##                     the first betadim entries are the entries corresponding to
##                     beta, the remainig entries correspond to the parameter of
##                     the baseline.
##  betadim          - See par
##  covariates       - List of covariates, see documentation.pdf for details
##  alpha_covariates - List of baseline covariartes, see documentation.pdf for
##                     details
##  events           - Matrix of events, see documentation.pdf for details
##  T                - Endpoint of observation interval [0,T]
##  epp              - Same as in cpl above
## Output:
##  The negative value of the full likelihood together with its gradient and
##  the hessian matrix (as attributes "gradient" and "hessian", respectively)
cfl <- function(par,betadim,covariates,alpha_covariates,events,T,epp=NA) {
  print(par)
  p <- betadim
  q <- length(par)-p
  beta <-  par[   1 :p]
  theta <- par[(p+1):length(par)]

  out <- 0
  grad <- rep(0,p+q)
  hess <- matrix(0,ncol=p+q,nrow=p+q)
  
  ## Compute epp unless given
  if(is.na(epp)==TRUE) {
    epp <- events_per_pair(covariates,events)
  }
  
  ## Compute Contribution from \int log(alpha(theta,t))dN(t)
  for(k in 1:length(alpha_covariates$time)) {
    ## Find endpoints of segment
    begin_of_segment <- alpha_covariates$time[k]
    if(k==length(alpha_covariates$time)) {
      end_of_segment <- T
    } else {
      end_of_segment <- alpha_covariates$time[k+1]
    }
    
    noe <- length(which(events[,1]>=begin_of_segment & events[,1]<end_of_segment))
    ## Update likelihood and derivative
    out               <- out              +sum(alpha_covariates$covars[[k]]*theta)*noe
    grad[(p+1):(p+q)] <- grad[(p+1):(p+q)]+    alpha_covariates$covars[[k]]       *noe
  }
  
  ## Compute certain sums for later use
  S <- vector(mode="list",length=length(covariates$time))
  Sd <- vector(mode="list",length=length(covariates$time))
  Sh <- vector(mode="list",length=length(covariates$time))
  for(k in 1:length(covariates$time)) {
    S[[k]] <- as.numeric(Reduce('+',lapply(covariates$covars[[k]],function(a) {return(exp(beta%*%a))})))
    Sd[[k]] <- Reduce('+',lapply(covariates$covars[[k]],function(a) {return(a*as.numeric(exp(beta%*%a)))}))
    Sh[[k]] <- Reduce('+',lapply(covariates$covars[[k]],function(a) {return(a%*%t(a)*as.numeric(exp(beta%*%a)))}))
  }
  
  ## Compute contribution from \sum_{i,j} int log(Psi(X_{n,ij}(t);beta)) dN_{n,ij}(t)
  for(k in 1:length(covariates$time)) {
    ## Find endpoints of segment
    begin_of_segment <- covariates$time[k]
    if(k==length(covariates$time)) {
      end_of_segment <- T
    } else {
      end_of_segment <- covariates$time[k+1]
    }
    
    ## Get number of events in the time frame of interest (OPTIMIEREN)
    noe <- epp[[k]]
    
    ## Compute first integral in partial likelihood over time segment
    I1 <- rowSums(mapply(function(a,b) {return(a*b)},covariates$covars[[k]],noe))
    
    ## Update partial likelihood and derivatives
    out <- out+beta%*%I1
    grad[1:p] <- grad[1:p]+I1
  }
  

  
  ## Compute Contribution from \int alpha(theta,t)Xn(t,beta)dt
  for(k in 1:length(alpha_covariates$time)) {
    outer_start <- alpha_covariates$time[k]
    if(k==length(alpha_covariates$time)) {
      outer_end <- T
    } else {
      outer_end <- alpha_covariates$time[k+1]
    }
    for(l in 1:length(covariates$time)) {
      inner_start <- covariates$time[l]
      if(l==length(covariates$time)) {
        inner_end <- T
      } else {
        inner_end <- covariates$time[l+1]
      }
      
      ## Check if Time intervals overlap
      if(inner_end>outer_start & inner_start<outer_end) {
        ## Yes, there is overlap, so we compute something
        begin_of_segment <- max(c(outer_start,inner_start))
        end_of_segment   <- min(c(outer_end  ,inner_end))
        
        ## Update likelihood and derivative
        out <- out-(end_of_segment-begin_of_segment)*exp(sum(alpha_covariates$covars[[k]]*theta))*S[[l]]
        
        grad[1:p] <- grad[1:p]-(end_of_segment-begin_of_segment)*exp(sum(alpha_covariates$covars[[k]]*theta))*Sd[[l]]
        grad[(p+1):(p+q)] <- grad[(p+1):(p+q)]-(end_of_segment-begin_of_segment)*alpha_covariates$covars[[k]]*exp(sum(alpha_covariates$covars[[k]]*theta))*S[[l]]
        
        hess[1:p,1:p] <- hess[1:p,1:p]-(end_of_segment-begin_of_segment)*exp(sum(alpha_covariates$covars[[k]]*theta))*Sh[[l]]
        hess[(p+1):(p+q),(p+1):(p+q)] <- hess[(p+1):(p+q),(p+1):(p+q)]-(end_of_segment-begin_of_segment)*(alpha_covariates$covars[[k]]%*%t(alpha_covariates$covars[[k]]))*exp(sum(alpha_covariates$covars[[k]]*theta))*S[[l]]
        hess[1:p,(p+1):(p+q)] <- hess[1:p,(p+1):(p+q)]-(end_of_segment-begin_of_segment)*exp(sum(alpha_covariates$covars[[k]]*theta))*(Sd[[l]]%*%t(alpha_covariates$covars[[k]]))
        hess[(p+1):(p+q),1:p] <- hess[(p+1):(p+q),1:p]-(end_of_segment-begin_of_segment)*exp(sum(alpha_covariates$covars[[k]]*theta))*(alpha_covariates$covars[[k]]%*%t(Sd[[l]]))
      }
    }
  }
  
  output <- -out
  attr(output,"gradient") <- -grad
  attr(output,"hessian") <- -hess
  
  return(output)
}


## Computes the Nelson Aalen Estimator, the non-parametric baseline estimator.
## Input:
##  t          - Time points at which the estimator shall be evaluated
##  h          - Bandwidth to use
##  beta_tilde - Beta_tilde to be used
##  covariates - List of covariates, see documentation.pdf for details
##  events     - Matrix of events, see documentation.pdf for details
##  T          - Endpoint of observation interval [0,T]
##  kernelfunc - Kernel to be used, needs to be a function which takes vectors
##               of numbers as input and returns vectors of the same length.
## Output:
##  A vector of the same length as t which contains the Nelson-Aalen-Estimator
##  evaluated at the time points specified in t.
NAE <- function(t,h,beta_tilde,covariates,events,T,kernelfunc) {
  beta <- beta_tilde
  ## Set cumulated covariates evaluated at the jump times
  ccv <- rep(0,dim(events)[1])
  start <- 1
  for(k in 1:length(covariates$time)) {
    ## Find endpoints of segment
    begin_of_segment <- covariates$time[k]
    if(k==length(covariates$time)) {
      end_of_segment <- T
    } else {
      end_of_segment <- covariates$time[k+1]
    }
    
    ## Get relevant event indices
    rind <- which(events[,1]>=begin_of_segment & events[,1]<end_of_segment)
    
    if(length(rind)>0) {
      ccv[start:(start+length(rind)-1)] <- as.numeric(Reduce('+',lapply(covariates$covars[[k]],function(a) {return(exp(beta%*%a))})))
    }
  }

  ## compute the Nelson-Aalen Estimator
  eval_matrix <- matrix(events[,1],ncol=length(ccv),nrow=length(t),byrow = TRUE)-matrix(t,ncol=length(ccv),nrow=length(t))
  kernel_eval <- 1/h*kernelfunc(eval_matrix/h)
  NAE <- kernel_eval%*%(1/ccv)
  
  return(NAE)
}

## Compute number of events per pair
## Input:
##  covariates - List of covariates, see documentation.pdf for details
##  events     - Matrix of events, see documentation.pdf for details
## Output:
##  List of the same length as covariates$time. Each element of this list is a
##  vector. The i-th element of this vector contains the number of events which
##  happened in the corresponding covariate segment between the pair which cor-
##  responds to the i-th covariate in the covariate list.
events_per_pair <- function(covariates,events) {
  out <- vector(mode="list",length=length(covariates$time))
  ## Go through each covariate element
  for(k in 1:length(covariates$time)) {
    ## Find endpoints of segment
    begin_of_segment <- covariates$time[k]
    if(k==length(covariates$time)) {
      end_of_segment <- T
    } else {
      end_of_segment <- covariates$time[k+1]
    }
    
    ## Get relevant event indices
    rind <- which(events[,1]>=begin_of_segment & events[,1]<end_of_segment)
    
    ## Get number of events in the time frame of interest
    noe <- rep(0,length(covariates$covars[[k]]))
    for(i in 1:length(rind)) {
      noe[covariates$index[[k]][events[rind[i],2],events[rind[i],3]]] <- noe[covariates$index[[k]][events[rind[i],2],events[rind[i],3]]]+1
    }
    
    out[[k]] <- noe
  }
  return(out)
}

## Evaluates the Cox form of the baseline intensity for given parameters.
## Input:
##  t                - Time points at which the baseline shall be evaluated
##  alpha_covariates - Covariates used for the baseline, see documentation.pdf
##  theta            - Parameter to be used
##  T                - Endpoint of observation interval [0,T]
alpha_eval <- function(t,alpha_covariates,theta,T) {
  out <- rep(-Inf,length(t))
  
  for(k in 1:length(alpha_covariates$time)) {
    ## Find endpoints of segment
    begin_of_segment <- alpha_covariates$time[k]
    if(k==length(alpha_covariates$time)) {
      end_of_segment <- T
    } else {
      end_of_segment <- alpha_covariates$time[k+1]
    }
    
    ind <- which(t>=begin_of_segment & t<end_of_segment)
    
    out[ind] <- sum(alpha_covariates$covars[[k]]*theta)
  }
  
  return(exp(out))
}



## This functions computes the test statistic and the appropriate scaling.
## Input:
##  alpha_NP         - Vector containing the values of the non-parametric base-
##                     line estimator at the times t.
##  alpha_P          - Vector containing the values of the parametric baseline
##                     estimator at the times t.
##  wfunc            - Function w, must take vectors as input and returns a vec-
##                     tor of the same length
##  t                - See alpha_NP and alpha_P
##  events           - Matrix of events, see documentation.pdf for details
##  covariates       - List of covariates, see documentation.pdf for details
##  beta_hat         - Estimator for beta to be used
##  Kfunc            - Kernel to be used, needs to be a function which takes
##                     vectors of numbers as input and returns vectors of the
##                     same length.
##  h                - Bandwidth to be used
##  T                - Endpoint of observation interval [0,T]
##  Kintegrate       - Function f(x) which returns the integral of the kernel
##                     over (-infinity,x), x must be allowed to be a vector.
##  alpha_covariates - List of baseline covariartes, see documentation.pdf
##  theta_hat        - Value of theta to be used
##  wIntegrate       - Function f(x) which returns the integral of w over
##                     (-infinity,x), x must be allowed to be a vector.
##  Ksq              - Value of K^{(2)}
## Output: A list of the following elements
##  Tstat   - Estimated value of the test statistic
##  an      - Vector, estimate of a_n for each covariate segment
##  An      - Estimate for A_n
##  B       - Estimate for B
##  Tscaled - Scaled and centred test statistic (asymptotically N(0,1))
##  Xnbar   - Veector, Estimate for \bar{X}_n for each covariate segment
baseline_test <- function(alpha_NP,alpha_P,wfunc,t,events,covariates,beta_hat,Kfunc,h,T,Kintegrate,alpha_covariates,theta_hat,wIntegrate,Ksq) {
  ## Compute the midpoints of the intervals
  w <- wfunc(t)
  alpha_NP_midpoints <- (alpha_NP[1:(length(alpha_NP)-1)]+alpha_NP[2:length(alpha_NP)])/2
  alpha_P_midpoints  <- (alpha_P[ 1:(length(alpha_P) -1)]+alpha_P[ 2:length(alpha_P )])/2
  w_midpoints <- (w[1:(length(w)-1)]+w[2:length(w)])/2
  
  ## Compute the test statistic
  Tstat <- sum(diff(t)*(alpha_NP_midpoints-alpha_P_midpoints)^2*w_midpoints)
  
  ## Compute the scalings
  ## Add acolumn to the events stating which covariate index is relevant
  events <- add_covariate_index(events,covariates)
  
  ## Compute an estimate for \bar{X}_n and r_np_n
  barXn_hat <- rep(0,length=length(covariates$time))
  r_np_n <- rep(0,length(covariates$time))
  for(k in 1:length(covariates$time)) {
    r_np_n[k] <- length(covariates$covars[[k]])
    barXn_hat[k] <- as.numeric(Reduce('+',lapply(covariates$covars[[k]],function(a) {return(exp(beta_hat%*%a))})))/r_np_n[k]
  }
  
  ## an
  inti <- 0
  for(k in 1:length(covariates$time)) {
    ## Find endpoints of segment
    begin_of_segment <- covariates$time[k]
    if(k==length(covariates$time)) {
      end_of_segment <- T
    } else {
      end_of_segment <- covariates$time[k+1]
    }
    
    hfunc <- function(t) {
      return((Kintegrate((end_of_segment-t)/h)-Kintegrate((begin_of_segment-t)/h))*wfunc(t))
    }
    
    inti <- inti+integrate(hfunc,lower=0,upper=T)$value/r_np_n[k]
  }
  an <- inti^(-1/2)
  
  ## An
  An <- 0
  for(i in 1:length(events[,1])) {
    k <- events[i,5]
    An <- An+fn(events[i,1],events[i,1],Kfunc,wfunc,h,T)*(barXn_hat[k]^(-1)/r_np_n[k])^2
  }
  An <- an^2*An
  
  ## Gamma
  gamma <- rep(0,length(covariates$time))
  for(k in 1:length(covariates$time)) {
    gamma[k] <- an^4/r_np_n[k]^2
  }
  
  ## B
  B <- 0
  for(k in 1:length(alpha_covariates$time)) {
    outer_start <- alpha_covariates$time[k]
    if(k==length(alpha_covariates$time)) {
      outer_end <- T
    } else {
      outer_end <- alpha_covariates$time[k+1]
    }
    for(l in 1:length(covariates$time)) {
      inner_start <- covariates$time[l]
      if(l==length(covariates$time)) {
        inner_end <- T
      } else {
        inner_end <- covariates$time[l+1]
      }
      
      ## Check if Time intervals overlap
      if(inner_end>outer_start & inner_start<outer_end) {
        ## Yes, there is overlap, so we compute something
        begin_of_segment <- max(c(outer_start,inner_start))
        end_of_segment   <- min(c(outer_end  ,inner_end))
        
        B <- B+(end_of_segment-begin_of_segment)*barXn_hat[l]^(-2)*alpha_eval(0.5*(begin_of_segment+end_of_segment),alpha_covariates,theta_hat,T)^2*gamma[l]*(wIntegrate(end_of_segment)-wIntegrate(begin_of_segment))
      }
    }
  }
  B <- 4*Ksq*B                                                                                                                 

  
  return(list(Tstat=Tstat,an=an,An=An,B=B,Tscaled=(an^2*sqrt(h)*Tstat-An/sqrt(h))/sqrt(B),Xnbar=barXn_hat))
}


################################################################################
## Internal functions which, normally, need not be called by the user. #########
################################################################################


## This function adds a fifth column to the events matrix which contains the index
## of the covariate index which is active at time of the event.
#### Input ####
## "events" and "covariates" are both as in cpl.
#### Output ####
## A matrix which first four columns are identical to the input "events". Its fifth
## column contains the index of the covariate entry which is active at the time of
## the corresponding event.
add_covariate_index <- function(events,covariates) {
  events <- cbind(events,0)
  
  for(i in 1:(length(covariates$time)-1)) {
    indis <- which(events[,1]>=covariates$time[i] & events[,1]<covariates$time[i+1])
    events[indis,5] <- i
  }
  i <- length(covariates$time)
  indis <- which(events[,1]>=covariates$time[i])
  events[indis,5] <- i
  
  return(events)
}

## Compute fn
fn <- function(r,s,kernelfunc,w,h,T) {
  helper_function <- function(t,h,r,s,kernelfunc,w) {
    return(h*kernelfunc((s-t)/h)*kernelfunc((r-t)/h)*w(t))
  }
  return(integrate(helper_function,lower=max(c(min(c(r,s))-2*h,0)),upper=min(c(max(c(r,s))+2*h,T)),h=h,r=r,s=s,kernelfunc=kernelfunc,w=myw)$value)
}