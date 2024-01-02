NimModel <- nimbleCode({
  #Priors
  ######################################################################################
  ###priors - Concentration and Copy Number Model###
  #Concentration function parameters for C decline vs distance from source
  #Power law model
  p0 ~ dunif(0,50000)
  lambda.site ~ dunif(0,2)
  theta.site ~ dunif(0,1) #percent underestimation at first site, distance 0.5m
  #parameter for sample sd at a site
  zeta0 ~ dunif(0,5)
  
  ###priors - Inhibition Model###
  #inhibition probability parameters
  beta0.w ~ dlogis(0,1)
  beta1.w ~ T(dnorm(0,sd=5),-2,0)
  #thinning parameters
  theta0 ~ T(dlogis(0,1),-4,1.5)
  sd.w ~ T(dt(0,sigma=1,df=4),0,Inf)
  
  ###priors - Observation Model###
  #parameters relating copies to p.y
  gamma0 ~ dlogis(0,1)
  gamma1 ~ T(dnorm(0,sd=5),0,Inf)
  #copy number measurement error parameter
  alpha0 ~ dunif(0,0.6)

  #Likelihoods
  #########################################################################################
  ###Ecological Process Model - Concentration/Copy Number
  #thin first site, "plume effect"
  # C.site[1] <- (1/Q[1])*(p0*exp(-(dist[1]/(v*tau))^rho))*theta.site
  C.site[1] <- (1/Q[1])*p0*(dist[1]/v[1])^(-lambda.site)*theta.site
  for(i in 2:M){
    # C.site[i] <- (1/Q[i])*(p0*exp(-(dist[i]/(v*tau))^rho))
    C.site[i] <- (1/Q[i])*p0*(dist[i]/v[i])^(-lambda.site)
  }
  for(i in 1:M){
    sdlog.C.sample[i] <- zeta0 
    for(j in 1:J){
      log(C.sample[i,j]) ~ dnorm(log(C.site[i]),sd=sdlog.C.sample[i])
      #Convert sample concentration to expected number of copies in replicates for this sample
      lambda.rep[i,j] <- C.sample[i,j]*volume.rep[i,j] #expected N in replicates
      for(k in 1:K){
        N.rep[i,j,k] ~ dpois(lambda=lambda.rep[i,j]) #realized N in replicates - Poisson
      }
    }
  }
  ###Ecological Process Model - Inhibition Process
  for(i in 1:M){
    for(j in 1:J){
      unit.inhibit[i,j] ~ dnorm(0,sd=1) #prob of inhibition random effect, non-centered
      for(k in 1:K){
        #probability rep is inhibited|occupied by inhibitor
        logit(p.w[i,j,k]) <- beta0.w + beta1.w*(N.rep[i,j,k]-1) + unit.inhibit[i,j]*sd.w
        w[i,j,k] ~ dbern(p.w[i,j,k]) #rep-level inhibition state
        #Thinning process|inhibited
        logit(theta.thin[i,j,k]) <- theta0
        N.rep.thin[i,j,k] ~ dbinom(theta.thin[i,j,k],N.rep[i,j,k]) #thinned copies after inhibition
        #Which N.rep to use for each i,j,k depending on inhibition state
        N.rep.available[i,j,k] <- (1-w[i,j,k])*N.rep[i,j,k] + #if rep not inhibited
          w[i,j,k]*N.rep.thin[i,j,k] #if rep inhibited
        N.present[i,j,k] <- step( N.rep.available[i,j,k] - 0.1) #is at least 1 uninhibited copy present in replicate?
      }
    }
  }
  ###Observation Model - detection and quantification
  for(i in 1:M){
    for(j in 1:J){
      for(k in 1:K){
        #logit-linear model for p.y - function of number of copies in each replicate.
        #subtracting 1 so gamma0 corresponds to 1 copy instead of 0
        logit(p.y[i,j,k]) <- gamma0 + gamma1*(N.rep.available[i,j,k]-1) 
        sdlog.y[i,j,k] <- alpha0 #measurement sd
        #likelihood for detections, only if >0 true copy present in rep
        y.det[i,j,k] ~ dbern(p.y[i,j,k]*N.present[i,j,k])
        #likelihood for observations conditional on detection and N.rep
        y[i,j,k] ~ dObs(y.det = y.det[i,j,k],N.rep = N.rep.available[i,j,k],sdlog.y = sdlog.y[i,j,k],eta.q=eta.q)
      }
    }
  }
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