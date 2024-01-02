#inhibitors, power law survival
source("sim.data.inhibitor.R")

#Data dimensions
M <- 13
J <- 10
K <- 5
# dist <- seq(0.5,1000,length.out=M) #distance of each site from source
dist <- c(0.5,2.5,5,10,20,40,80,125,200,300,400,500,1000) #distance of each site from source
if(length(dist)!=M)stop("dist must be of length M")

#eDITH Process Parameters
p0 <- 6192.59
lambda.site <- 0.59
theta.site <- 0.03

#Sample variation in concentration at a site
zeta0 <- 0.89

#Replicate Inhibition parameters
sd.w <- 3.27
beta0.w <- 1.85
beta1.w <- -0.03

#Inhibitor Thinning parameter
theta0 <- -2.29

#Detection parameters
gamma0 <-  -2.39
gamma1 <- 2.52

#Measurement error parameter
alpha0 <- 0.29


volume.rep <- matrix(0.01,M,J)
eta.q <- exp(28.5371-0.6411*50) #corresponds to Cq of 50 in real data


#hydrological inputs
mean.Q <- 0.01992154 #site-level discharge (m^3/s)
mean.v <- 0.245908 #mean stream velocity (m/s)

#using mean of each for all sites for simplicity
Q <- rep(mean.Q,M)
v <- rep(mean.v,M)

data <- sim.data.inhibitor(M=M,J=J,K=K,dist=dist,volume.rep=volume.rep,eta.q=eta.q,
                 p0=p0,lambda.site=lambda.site,theta.site=theta.site,
                 beta0.w=beta0.w,beta1.w=beta1.w,sd.w=sd.w,
                 theta0=theta0, gamma0=gamma0,gamma1=gamma1,
                 zeta0=zeta0,alpha0=alpha0,Q=Q,v=v,eDITHmod="PL")

plot(data$C.site~dist)
rowMeans(data$N.rep)

#extract data object
y <- data$y
y.det <- data$y.det
volume.rep <- matrix(volume.rep,M,J)

library(nimble)
library(coda)
source("Release NimModel Inhibit PL.R")
source("State Samplers.R")

##Initialize model using data
samp.means.copies <- apply(y+0.01,c(1,2),mean,na.rm=TRUE) #sample mean observed copies plus a small delta (0.01)
#may want to exclude zeros above if many non-detects
samp.means.concentration <- samp.means.copies/volume.rep #convert to concentration
nan.idx <- which(is.nan(samp.means.concentration)) #fill in any NAs for samples with no observations
if(length(nan.idx)>0){
  samp.means.concentration[nan.idx] <- mean(samp.means.concentration[-nan.idx])
}

C.site.init <- rowMeans(samp.means.concentration) 
C.sample.init <- samp.means.concentration
C0.site.init <- max(C.site.init)
##get eDITH inits that will give us finite starting likelihood
p0.init <- Q[1]*C0.site.init #good choice for p0, maximize likelihood for lambda.site given p0.init
lik <- function(par,sd=100){
  lambda.site <- par[1]
  mu <- p0.init/Q*(dist/v)^(-lambda.site)
  ll <- sum(dnorm(C.sample.init,mu,sd,log=TRUE))
  return(-ll)
}
fit <- optim(par=c(log(1)),fn=lik,method="Brent",lower=0,upper=2)
lambda.site.init <- fit$par

N.rep.init <- ceiling(y) #ceiling to round measurements <1 up to 1 so starting likelihood is finite
N.rep.init <- N.rep.init + 1 #then add 1 for cases where detected, but below eta.q


# minimal inits, takes longer to converge
Niminits <- list(p0=p0.init,lambda.site=lambda.site.init,
                 unit.inhibit=matrix(0,M,J),
                 C.site=C.site.init,log_C.sample=log(C.sample.init),
                 N.rep=N.rep.init,
                 w=array(0,dim=c(M,J,K))) #to ensure finite starting likelihood, init these at all 0 (no inhibition))


constants <- list(M=M,J=J,K=K,eta.q=eta.q,Q=Q,v=v,dist=dist)
Nimdata <- list(y=data$y,y.det=data$y.det,volume.rep=volume.rep)

# set parameters to monitor
parameters <- c('alpha0','zeta0','gamma0','gamma1','p0','theta0','lambda.site','sd.w',
                'beta0.w','theta.site','beta1.w')
parameters2 <- c('C.site','C.sample','N.rep','N.rep.available')

thin <- 5
thin2 <- 20

start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)
conf <- configureMCMC(Rmodel,monitors=parameters, thin=thin,monitors2=parameters2,thin2=thin2,useConjugacy = FALSE)

##add more samplers
#Sampler 1) jointly update w and N.rep - required
conf$removeSampler('w') #sampler below makes this redundant, can remove
#keep slice N.rep and N.rep.thin samplers for better mixing
# conf$removeSampler('N.rep')
# conf$removeSampler('N.rep.thin')
for(i in 1:M){
  for(j in 1:J){
    for(k in 1:K){
      calcNodes.w <- Rmodel$expandNodeNames(paste("w[",i,",",j,",",k,"]"))
      calcNodes.N.rep <- Rmodel$getDependencies(c(paste("N.rep.thin[",i,",",j,",",k,"]"),
                                                  Rmodel$getDependencies(paste("N.rep[",i,",",j,",",k,"]"))))
      calcNodes.p.w <- Rmodel$expandNodeNames(paste("p.w[",i,",",j,",",k,"]"))
      calcNodes.theta.thin <- Rmodel$expandNodeNames(paste("theta.thin[",i,",",j,",",k,"]"))
      calcNodes.all <- c(calcNodes.N.rep,calcNodes.w)
      conf$addSampler(target = paste("w[",i,",",j,",",k,"]"),
                      type = 'wSampler',control = list(calcNodes=calcNodes.all,
                                                       calcNodes.p.w=calcNodes.p.w,
                                                       calcNodes.theta.thin=calcNodes.theta.thin,
                                                       i=i,j=j,k=k),
                      silent = TRUE)
    }
  }
}
#Sampler 2) 2nd w update "shift up/down". Not required but I am using it, likely improves mixing, but not sure
for(i in 1:M){
  for(j in 1:J){
    for(k in 1:K){
      calcNodes.w <- Rmodel$expandNodeNames(paste("w[",i,",",j,",",k,"]"))
      calcNodes.N.rep <- Rmodel$getDependencies(c(paste("N.rep.thin[",i,",",j,",",k,"]"),
                                                  Rmodel$getDependencies(paste("N.rep[",i,",",j,",",k,"]"))))
      calcNodes.p.w <- Rmodel$expandNodeNames(paste("p.w[",i,",",j,",",k,"]"))
      calcNodes.theta.thin <- Rmodel$expandNodeNames(paste("theta.thin[",i,",",j,",",k,"]"))
      calcNodes.all <- c(calcNodes.N.rep,calcNodes.w) #can remove y nodes from here
      conf$addSampler(target = paste("w[",i,",",j,",",k,"]"),
                      type = 'wSampler2',control = list(calcNodes=calcNodes.all,
                                                        calcNodes.p.w=calcNodes.p.w,
                                                        calcNodes.theta.thin=calcNodes.theta.thin,
                                                        i=i,j=j,k=k),
                      silent = TRUE)
    }
  }
}

#jointly update C.sample/lambda.rep, w, and N.rep at sample level.
for(i in 1:M){
  for(j in 1:J){
    calcNodes.C.sample <- Rmodel$getDependencies(paste("log_C.sample[",i,",",j,"]"))
    calcNodes.unit.inhibit <- Rmodel$expandNodeNames(paste("unit.inhibit[",i,",",j,"]"))
    calcNodes.N.rep <- Rmodel$getDependencies(c(paste("N.rep[",i,",",j,",1:",K,"]"),
                                                paste("N.rep.thin[",i,",",j,",1:",K,"]")))
    calcNodes.w <- Rmodel$expandNodeNames(paste("w[",i,",",j,",1:",K,"]"))
    calcNodes.p.w <- Rmodel$expandNodeNames(paste("p.w[",i,",",j,",1:",K,"]"))
    calcNodes.theta.thin <- Rmodel$expandNodeNames(paste("theta.thin[",i,",",j,",1:",K,"]"))
    calcNodes.C.sample <- setdiff(calcNodes.C.sample,calcNodes.N.rep)
    calcNodes.z <- Rmodel$expandNodeNames(paste("z[",i,",",j,"]"))
    calcNodes1 <- c(calcNodes.C.sample,calcNodes.z,calcNodes.unit.inhibit)
    calcNodes2 <- c(calcNodes.N.rep)
    calcNodes.all <- c(calcNodes1,calcNodes2)
    conf$addSampler(target = paste("log_C.sample[",i,",",j,"]"),
                    type = 'C.sample.Sampler',control = list(calcNodes1=calcNodes1,
                                                             calcNodes2=calcNodes2,
                                                             calcNodes3=calcNodes.all,
                                                             calcNodes.p.w=calcNodes.p.w,
                                                             calcNodes.theta.thin=calcNodes.theta.thin,
                                                             i=i,j=j,K=K),
                    silent = TRUE)
  }
}

#keep independent samplers and use block sampler (independent ones help with convergence while block RW still tuning)
conf$addSampler(target = c("p0","lambda.site","theta.site"),
                type = 'RW_block',control = list(adaptive=TRUE),silent = TRUE)

conf$addSampler(target = c("beta0.w","beta1.w"),
                type = 'RW_block',control = list(adaptive=TRUE),silent = TRUE)
# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=10) #this will run in R, used for debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model
start.time2 <- Sys.time()
Cmcmc$run(5000,reset=FALSE) #short run for demonstration. Can run repeatedly to get more posterior samples.
end.time <- Sys.time()
end.time - start.time  # total time for compilation, replacing samplers, and fitting
end.time - start.time2 # post-compilation run time

mvSamples <- as.matrix(Cmcmc$mvSamples)
mvSamples2 <- as.matrix(Cmcmc$mvSamples2)

burnin <- 500
plot(mcmc(mvSamples[-c(1:burnin),]))

###Process posteriors in mvSamples2
burnin2 <- burnin/thin2
library(MCMCglmm)
#Site and sample concentrations
C.site.idx <- grep("C.site",colnames(mvSamples2))[1:M] #1:M if C.site.pred in model
C.sample.idx <- grep("C.sample",colnames(mvSamples2))
# C.site.log.ests <- colMeans(log(mvSamples2[-c(1:burnin2),C.site.idx]))
C.site.log.ests <- posterior.mode(mcmc(log(mvSamples2[-c(1:burnin2),C.site.idx])),adjust=1)
C.site.log.HPDs <- matrix(0,length(C.site.idx),2)
for(i in 1:M){
  C.site.log.HPDs[i,] <- HPDinterval(mcmc(log(mvSamples2[-c(1:burnin2),C.site.idx[i]])))
}
# C.sample.log.ests <- colMeans(log(mvSamples2[-c(1:burnin2),C.sample.idx]))
C.sample.log.ests <- posterior.mode(mcmc(log(mvSamples2[-c(1:burnin2),C.sample.idx])),adjust=1)
C.sample.log.ests2D <- matrix(C.sample.log.ests,M,J)
C.sample.log.HPDs <- array(0,dim=c(M,J,2))
idx <- 1
for(j in 1:J){
  for(i in 1:M){
    C.sample.log.HPDs[i,j,] <- HPDinterval(mcmc(log(mvSamples2[-c(1:burnin2),C.sample.idx[idx]])))
    idx <- idx+1
  }
}
##Note, all HPD intervals below don't account for multimodality. Use HDInterval R package for those.
N.rep.available.idx <- grep("N.rep.available",colnames(mvSamples2))
N.rep.idx <- grep("N.rep",colnames(mvSamples2))
N.rep.idx <- setdiff(N.rep.idx,N.rep.available.idx)
# plot(mcmc(mvSamples2[-c(1:burnin2),N.rep.idx[1:100]]))
# N.rep.est <- colMeans(mvSamples2[-c(1:burnin2),N.rep.idx])
N.rep.est <- round(posterior.mode(mcmc(mvSamples2[-c(1:burnin2),N.rep.idx]),adjust=1))
N.rep.HPDs <- HPDinterval(mcmc(mvSamples2[-c(1:burnin2),N.rep.idx]))
N.rep.est3D <- array(N.rep.est,dim=c(M,J,K))
N.rep.lower <- array(N.rep.HPDs[,1],dim=c(M,J,K))
N.rep.upper <- array(N.rep.HPDs[,2],dim=c(M,J,K))
N.rep.cover <- array(NA,dim=c(M,J,K))

N.rep.available.idx <- grep("N.rep.available",colnames(mvSamples2))
# N.rep.available.est <- colMeans(mcmc(mvSamples2[-c(1:burnin2),N.rep.available.idx])
N.rep.available.est <- round(posterior.mode(mcmc(mvSamples2[-c(1:burnin2),N.rep.available.idx]),adjust=1))
N.rep.available.HPDs <- HPDinterval(mcmc(mvSamples2[-c(1:burnin2),N.rep.available.idx]))
N.rep.available.est3D <- array(N.rep.available.est,dim=c(M,J,K))
N.rep.available.lower <- array(N.rep.available.HPDs[,1],dim=c(M,J,K))
N.rep.available.upper <- array(N.rep.available.HPDs[,2],dim=c(M,J,K))
N.rep.available.cover1 <- array(NA,dim=c(M,J,K))
N.rep.available.cover2 <- array(NA,dim=c(M,J,K))



###Color palette for plots
library(RColorBrewer)
cols <- brewer.pal(9,"Set3")

###Concentration Plot###
plot.log10 <- TRUE #can plot on log10 or natural log scale (as modeled)
offset.tmp <- seq(-0.04,0.04,length.out=J)
if(plot.log10){
  #log10 conversions for labels
  tmp <- range(C.sample.log.HPDs)/log(10,base=exp(1))
  tmp[1] <- floor(tmp[1])
  tmp[2] <- ceiling(tmp[2])
  plot.range1 <- seq(tmp[1],tmp[2],by=0.5)
  plot.range2 <- plot.range1/log(exp(1),base=10)
  tmp2 <- range(log(dist))/log(10,base=exp(1))
  tmp2[1] <- floor(tmp2[1])
  tmp2[2] <- ceiling(tmp2[2])
  plot.range3 <- seq(tmp2[1],tmp2[2],by=0.5)
  plot.range4 <- plot.range3/log(exp(1),base=8)
  plot(C.site.log.ests~log(dist),pch=16,ylim=range(C.sample.log.HPDs),xlab="Log Distance (Log10 Meters)",
       ylab="Log Site Concentration (Log10 Copies per Liter)",
       main="Estimated Site and Sample Concentration vs. Distance",xaxt='n',yaxt='n')
  axis(1,at=plot.range4,labels=round(plot.range3,2))
  axis(2,at=plot.range2,labels=plot.range1)
}else{
  plot(C.site.log.ests~log(dist),pch=16,ylim=range(C.sample.log.HPDs),xlab="Log Distance (Log Meters)",
       ylab="Log Site Concentration (Log Copies per Liter)",
       main="Estimated Site and Sample Concentration vs. Distance")
}
for(i in 1:M){
  for(j in 1:J){
    offset=sample(offset.tmp,replace=FALSE)
    points(C.sample.log.ests2D[i,j]~c(log(dist[i])+offset[j]),pch=16,col=cols[5])
    lines(x=rep(log(dist[i])+offset[j],2),y=C.sample.log.HPDs[i,j,1:2],col=cols[5])
  }
  lines(rep(log(dist[i]),2),C.site.log.HPDs[i,1:2],lwd=2)
}
points(C.site.log.ests~log(dist),pch=16)
legend("topright",legend=c("Site C","Sample C"),
       pch=c(16,16),lwd=c(NA,NA),col=c("black",cols[5]))


plot.log10 <- TRUE
if(plot.log10){
  plot(log(N.rep.est3D)~log(y),pch=16,main="Measured vs. True Log10 Copies per Replicate",
       xlab="Measured Log10 Copies",ylab="True Log10 Copies",
       xlim=c(0,max(log(N.rep.upper))),ylim=c(0,max(log(N.rep.upper))),
       xaxt='n',yaxt='n')
  tmp <- range(log(N.rep.HPDs))/log(10,base=exp(1))
  tmp[1] <- 0
  tmp[2] <- ceiling(tmp[2])
  plot.range1 <- seq(tmp[1],tmp[2],by=0.5)
  plot.range2 <- plot.range1/log(exp(1),base=10)
  axis(1,at=plot.range2,labels=round(plot.range1,2))
  axis(2,at=plot.range2,labels=plot.range1)
}else{
  plot(log(N.rep.est3D)~log(y),pch=16,main="Measured vs. True Log Copies per Replicate",
       xlab="Measured Log Copies",ylab="True Log Copies",
       xlim=c(0,max(log(N.rep.upper))),ylim=c(0,max(log(N.rep.upper))))
}
for(i in 1:M){
  for(j in 1:J){
    for(k in 1:K){
      N.rep.cover[i,j,k] <- 1*(N.rep.lower[i,j,k]<=data$N.rep[i,j,k]&N.rep.upper[i,j,k]>=data$N.rep[i,j,k])
      if(N.rep.cover[i,j,k]==1){
        lines(x=rep(log(y[i,j,k]),2),y=c(log(N.rep.lower[i,j,k]),log(N.rep.upper[i,j,k])))
      }else{
        lines(x=rep(log(y[i,j,k]),2),y=c(log(N.rep.lower[i,j,k]),log(N.rep.upper[i,j,k])),col="red")
      }
    }
  }
}
abline(0,1)
abline(h=log(eta.q),col="darkblue")



plot.log10 <- TRUE
if(plot.log10){
  plot(log(N.rep.available.est3D)~log(y),pch=16,main="Measured vs. Available Log10 Copies per Replicate",
       xlab="Measured Log10 Copies",ylab="Available Log10 Copies",ylim=c(0,max(log(N.rep.available.upper))),
       xaxt='n',yaxt='n')
  tmp <- range(log(N.rep.HPDs))/log(10,base=exp(1))
  tmp[1] <- 0
  tmp[2] <- ceiling(tmp[2])
  plot.range1 <- seq(tmp[1],tmp[2],by=0.5)
  plot.range2 <- plot.range1/log(exp(1),base=10)
  axis(1,at=plot.range2,labels=round(plot.range1,2))
  axis(2,at=plot.range2,labels=plot.range1)
}else{
  plot(log(N.rep.available.est3D)~log(y),pch=16,main="Measured vs. Available Log Copies per Replicate",
       xlab="Measured Log Copies",ylab="Available Log Copies",ylim=c(0,max(log(N.rep.available.upper))))
}
for(i in 1:M){
  for(j in 1:J){
    for(k in 1:K){
      #computing coverage here...
      N.rep.available.cover1[i,j,k] <- 1*(N.rep.available.lower[i,j,k]<=y[i,j,k]&N.rep.available.upper[i,j,k]>=y[i,j,k])
      N.rep.available.cover2[i,j,k] <- 1*(N.rep.available.lower[i,j,k]<=data$N.rep.available[i,j,k]&
                                            N.rep.available.upper[i,j,k]>=data$N.rep.available[i,j,k])
      if(N.rep.available.cover1[i,j,k]==1){
        lines(x=rep(log(y[i,j,k]),2),y=c(log(N.rep.available.lower[i,j,k]),log(N.rep.available.upper[i,j,k])))
      }else{
        lines(x=rep(log(y[i,j,k]),2),y=c(log(N.rep.available.lower[i,j,k]),log(N.rep.available.upper[i,j,k])),col="red")
      }
    }
  }
}
abline(0,1)
abline(h=log(eta.q),col="darkblue")


mean(N.rep.available.cover1[y>0]) #does N.rep.available cover y? May be low without enough data
#these should have close to 95% coverage
mean(N.rep.cover) #does N.rep cover true N.rep?
mean(N.rep.available.cover2) #does N.available cover true N.available?