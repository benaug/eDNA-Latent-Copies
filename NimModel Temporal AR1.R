NimModel <- nimbleCode({
  ###priors - Ecological Process Model###
  # sdlog.C.time ~ dunif(0,5)
  sdlog.C.sample ~ dunif(0,5)
  sdlog.C.time ~ T(dt(0, df = 7, sigma = 1),0,) #half-t prior
  # sdlog.C.sample ~ T(dt(0, df = 7, sigma = 1),0,) #half-t prior
  
  ###priors - Observation Model###
  #parameters relating copies to p.y
  gamma0 ~ dlogis(0,1)
  gamma1 ~ T(dnorm(0,sd=5),0,Inf)
  #copy number measurement error parameter
  alpha0 ~ dunif(0,3)
  
  ##Ecological Process Model - Concentration/Copy Number
  beta0 ~ dnorm(0,sd=10)
  #AR1 in continuous time, equivalent to spatial exponential
  #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4020183/#R2
  logit(rho) ~ dlogis(0,1)
  for(i in 1:M) {
    log.mu.time[i] <- beta0
    for(j in 1:M) {
      Sigma[i,j] <- sdlog.C.time^2*(1*(equals(i,j) + (1-equals(i,j))*pow(rho,lag[i,j])))
    }
  }
  log(C.time[1:M]) ~ dmnorm(log.mu.time[1:M],cov=Sigma[1:M,1:M])
  for(t in 1:M){ #times when samples recorded
    for(j in 1:J.time[t]){ #samples at time t
      #Sample Concentration at time t follows Normal RV with mean equal to expected C at time t
      log(C.sample[t,j]) ~ dnorm(log(C.time[t]),sd=sdlog.C.sample)
      # unit.sample[t,j] ~ dnorm(0,sd=1)
      # log(C.sample[t,j]) <- log(C.time[t]) + unit.sample[t,j]*sdlog.C.sample
      #Convert sample concentration to mean number of copies in sample
      lambda.rep[t,j] <- C.sample[t,j]*volume.rep[t,j] #expected N in sample
      #realized N in replicates
      for(k in 1:K.time[t,j]){
        N.rep[t,j,k] ~ dpois(lambda.rep[t,j])
        N.present[t,j,k] <- step( N.rep[t,j,k] - 0.1) #is at least 1 copy present in replicate?
      }
    }
  }
  ###Observation Model - detection and quantification
  for(t in 1:M){ #times when samples recorded
    for(j in 1:J.time[t]){ #samples at time t
      for(k in 1:K.time[t,j]){
        #logit-linear model for p.y - function of number of copies in each replicate.
        #subtracting 1 so gamma0 corresponds to 1 copy instead of 0
        logit(p.y[t,j,k]) <- gamma0 + gamma1*(N.rep[t,j,k]-1)
        sdlog.y[t,j,k] <- alpha0 #measurement sd
        #likelihood for detections, only if >0 true copy present in rep
        y.det[t,j,k] ~ dbern(p.y[t,j,k]*N.present[t,j,k])
        #likelihood for observations conditional on detection and N.rep
        y.time[t,j,k] ~ dObs(y.det = y.det[t,j,k],
                                       N.rep = N.rep[t,j,k],
                                       sdlog.y = sdlog.y[t,j,k],eta.q=eta.q)
      }
    }
  }
  ####Predictions
  for(t in 1:N.N.pred){
    sdlog.y.pred[t] <- alpha0
    log(y.pred[t]) ~ T(dnorm(log(N.pred[t]), sd = sdlog.y.pred[t]),log(eta.q),Inf)
  }
  for(t in 1:N.N.pred2){
    logit(p.y.pred[t]) <- gamma0 + gamma1*(N.pred2[t]-1)
  }
  
  #Variance decomposition, not including observation model
  #Decomposing variance in copies allocated to reps
  #marginal expectation considering volume.rep
  E.Marginal <-  exp(beta0 + log(volume.rep[1,1]) + (sdlog.C.time^2 + sdlog.C.sample^2)/2)
  V1 <- E.Marginal #level 1 variance
  V2 <- (E.Marginal^2)*exp(sdlog.C.time^2)*(exp(sdlog.C.sample^2)-1) #level 2 variance
  V3 <- (E.Marginal^2)*(exp(sdlog.C.time^2)-1) #level 3 variance
  Total.var <- V1 + V2 + V3
  VPC[1] <- V1/Total.var #variance allocating copies to reps
  VPC[2] <- V2/Total.var #variance in copies across samples
  VPC[3] <- V3/Total.var #variance in copies across time
})# end model


#Observation model distribution
dObs <- nimbleFunction(
  run = function(x = double(0),y.det = integer(0), N.rep = integer(0), sdlog.y = double(0),eta.q = double(0),log = integer(0)) {
    returnType(double(0))
    if(y.det>0){ #if there was a detection
      if(x>0){#did we observe a copy number?
        logProb <- dnorm(log(x),log(N.rep),sd=sdlog.y,log=TRUE)
      }else{#otherwise, censored
        logProb <- pnorm(log(eta.q),log(N.rep),sdlog.y,log=TRUE)
      }
    }else{ #otherwise, no detection
      logProb <- 0
    }
    if(log){
      return(logProb)
    }else{
      return(exp(logProb))
    }
  }
)
#dummy RNG not used. makes nimble happy.
rObs <- nimbleFunction(
  run = function(n = integer(0),y.det = integer(0), N.rep = integer(0), sdlog.y = double(0),eta.q=double(0)) {
    returnType(double(0))
    if(y.det>0){
      sim <- exp(rnorm(1,mean=log(N.rep),sd=sdlog.y))
      if(sim > eta.q){
        return(sim)
      }else{
        return(0)
      }
    }else{
      return(0)
    }
  }
)