sim.data <- function(M=NA,J=NA,K=NA,dist=NA,volume.rep=NA,eta.q=NA,
                     p0=NA,lambda.site=NA,tau=NA,theta.site=NA,
                     gamma0=NA,gamma1=NA,
                     zeta0=NA,alpha0=NA,
                     Q=NA,v=NA,seed=NA,eDITHmod="PL"){
  if(!is.na(seed)){
    set.seed(seed)
  }
  
  if(eDITHmod=="PL"){
    if(is.na(lambda.site)){
      stop("if eDITHmod=PL, must provide lambda.site")
    }
  }else if(eDITHmod=="Exp"){
    if(is.na(tau)){
      stop("if eDITHmod=Exp, must provide tau")
    }
  }else{
    stop("eDITHmod must be PL or Exp")
  }
  
  ###Ecological Process Model - Concentration/Copy Number
  #Site concentration
  C.site <- rep(NA,M)
  if(eDITHmod=="PL"){
    C.site[1] <- (1/Q[1])*p0*(dist[1]/v[1])^(-lambda.site)*theta.site
    for(i in 2:M){
      C.site[i] <- (1/Q[i])*p0*(dist[i]/v[i])^(-lambda.site)
    }
  }else{
    C.site[1] <- (1/Q[1])*(p0*exp(-dist[1]/(v[1]*tau)))*theta.site
    for(i in 2:M){
      C.site[i] <- (1/Q[i])*(p0*exp(-dist[i]/(v[i]*tau)))
    }
  }
  #Sample concentration
  C.sample <- lambda.rep <- matrix(NA,M,J)
  N.rep <- N.present <- array(NA,dim=c(M,J,K))
  sdlog.C.sample <- rep(zeta0,M)
  for(i in 1:M){
    for(j in 1:J){
      C.sample[i,j] <- exp(rnorm(1,log(C.site[i]),sd=sdlog.C.sample[i]))
      #Convert sample concentration to expected number of copies in replicates for this sample
      lambda.rep[i,j] <- C.sample[i,j]*volume.rep[i,j] #expected N in replicates
      # p.rep[i,j] <- phi.rep/(phi.rep+lambda.rep[i,j]) #for NegBin copy model
      for(k in 1:K){
        N.rep[i,j,k] <- rpois(1,lambda=lambda.rep[i,j]) #realized N in replicates - Poisson
        # N.rep[i,j,k] ~ dnbinom(p=p.rep[i,j],size=phi.rep) #realized N in replicates - NegBin
        N.present[i,j,k] <- 1*(N.rep[i,j,k]>0) #is at least 1 copy present in replicate?
      }
    }
  }
  
  ###Observation Model - detection and quantification
  p.y <- sdlog.y <- y.det <- y <- array(NA,dim=c(M,J,K))
  for(i in 1:M){
    for(j in 1:J){
      for(k in 1:K){
        #logit-linear model for p.y - function of number of copies in each replicate.
        #subtracting 1 so gamma0 corresponds to 1 copy instead of 0
        p.y[i,j,k] <- plogis(gamma0 + gamma1*(N.rep[i,j,k]-1))
        sdlog.y[i,j,k] <- alpha0 #measurement sd
        #likelihood for detections, only if >0 true copy present in rep
        y.det[i,j,k] <- rbinom(1,1,p.y[i,j,k]*N.present[i,j,k])
        #likelihood for observations conditional on detection and N.rep
        if(y.det[i,j,k]>0){ #if there was a detection
          y[i,j,k] <- exp(rnorm(1,log(N.rep[i,j,k]),sd=sdlog.y[i,j,k]))
          if(y[i,j,k] < eta.q){ #is obs censored?
            y[i,j,k] <- 0
          }
        }else{ #otherwise, no detection
          y[i,j,k] <- 0
        }
      }
    }
  }
  return(list(y=y,y.det=y.det,
              C.site=C.site,C.sample=C.sample,
              N.rep=N.rep))
}