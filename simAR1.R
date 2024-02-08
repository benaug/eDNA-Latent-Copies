simAR1 <- function(beta0=NA,sdlog.C.time=NA,rho=NA,gamma0=NA,gamma1=NA,sdlog.y=NA,
                   times=NA,J.time=NA,K.time=NA,volume.rep=NA,eta.q=NA){
  M <- length(times)
  J <- max(J.time)
  K <- max(K.time)
  #get lag
  lag <- matrix(NA,M,M)
  for(i in 1:M){
    for(j in 1:M){
      lag[i,j] <- times[i] - times[j]
    }
  }
  lag <- abs(lag)
  #Temporal process
  library(mvtnorm)
  Sigma <- rho^lag
  diag(Sigma) <- 1
  Sigma <- Sigma*sdlog.C.time^2
  log.C.time <- rmvnorm(1,mean=rep(beta0,M),sigma=Sigma)[1,]
  # plot(log.C.time~times,type="l")
  C.time <- exp(log.C.time)
  
  #Sampling process
  log.C.sample <- matrix(NA,M,J)
  for(t in 1:M){
    log.C.sample[t,1:J.time[t]] <- rnorm(J.time[t],log.C.time[t],sdlog.C.sample)
  }
  C.sample <- exp(log.C.sample)
  
  #Replication process
  lambda.rep <- C.sample*volume.rep
  N.rep <- array(NA,dim=c(M,J,K))
  for(k in 1:K){
    N.rep[,,k] <- rpois(M*J,lambda.rep)
  }
  N.present <- 1*(N.rep>0)
  
  #Detection process
  p.y <- plogis(gamma0 + gamma1*(N.rep-1))
  y.det <- array(rbinom(M*J*K,1,p.y),dim=c(M,J,K))
  y.det[N.present==0] <- 0 #can't detect copies when not present
  
  #Quantification process
  log.y.time <- array(rnorm(M*J*K,log(N.rep),sd=sdlog.y),dim=c(M,J,K))
  y.time <- exp(log.y.time)
  y.time[y.det==0] <- 0 #can't quantify copies when not detected
  y.time[y.time<eta.q] <- 0
  # plot(log(y.time)~log(N.rep))
  return(list(y.time=y.time,y.det=y.det,log.C.time=log.C.time,log.C.sample=log.C.sample,
              lag=lag))
}
