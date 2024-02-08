#Load data simulator
source("simAR1.R")

#Temporal process parameters (continuous time AR1, equivalent to spatial exponential)
beta0 <- 6 #mean C (copies/liter) on natural log scale
sdlog.C.time <- 1 #temporal variation
rho <- 0.8 #temporal correlation (must be between 0 and 1)

#Sampling process parameters
sdlog.C.sample <- 0.75 #sampling variation

#Observation model parameters
gamma0 <- 1 #detection prob of 1 copy in replicate on logit link
gamma1 <- 4 #slope for relationship between detection prob and number of copies in replicate
sdlog.y <- 0.3 #copy measurement error given detection. (called alpha0 in nimble model)

times <- seq(0,250,2) #sampling times
M <- length(times) #total number of sampling times
J.time <- rep(3,M) #number of temporal replicates at each sampling time
J <- max(J.time) #max temporal reps
K.time <- matrix(5,M,J) #number of replicates per sample
K <- max(K.time) #maximum number of replicates per sample
volume.rep <- matrix(0.01,M,J) #sample volume associated with each rep (in liters)
eta.q <- 0 #minimum measurement allowed in quantification process on scale of absolute copies.
#can set this to limit false positives, or better meet model assumptions if they break down at very low copy number

data <- simAR1(beta0=beta0,sdlog.C.time=sdlog.C.time,rho=rho,
               gamma0=gamma0,gamma1=gamma1,sdlog.y=sdlog.y,
               times=times,J.time=J.time,K.time=K.time,
               volume.rep=volume.rep,eta.q=eta.q)

times2D <- matrix(times,M,J)
par(mfrow=c(1,1),ask=FALSE)
plot(NA,xlim=range(times),ylim=range(data$log.C.sample),xlab="Time",ylab="log Concentration",
     main="Concentration (natual log) through Time") #natural log
points(data$log.C.sample~times2D,pch=16)
lines(data$log.C.time~times,lwd=2,col="darkred")

mean(data$y.det)
hist(data$y.time)

#Fit model
y.time <- data$y.time
y.det <- data$y.det
lag <- data$lag

##Initialize model using data
samp.means.copies <- apply(y.time+0.01,c(1,2),mean,na.rm=TRUE) #sample mean observed copies plus a small delta (0.01)
#may want to exclude zeros above if many non-detects
samp.means.concentration <- samp.means.copies/volume.rep #convert to concentration
nan.idx <- which(is.nan(samp.means.concentration)) #fill in any NAs for samples with no observations
if(length(nan.idx)>0){
  samp.means.concentration[nan.idx] <- mean(samp.means.concentration[-nan.idx])
}

C.time.init <- rep(0,M)
C.time.init <- rowMeans(samp.means.concentration) 
C.sample.init <- matrix(0,M,J)
C.sample.init <- samp.means.concentration

library(nimble)
library(coda)
source("NimModel Temporal AR1.R")
N.rep.init <- ceiling(y.time) #ceiling to round measurements <1 up to 1 so starting likelihood is finite
N.rep.init <- N.rep.init + 1 #then add 1 for cases where detected, but below eta.q

##make objects used for model predictions
N.pred=seq(1:15)
N.N.pred=length(N.pred)

##make objects used for model predictions
N.pred <- round(c(1:5,exp(seq(log(6),6,length.out=25))))
N.N.pred <- length(N.pred)
N.pred2 <- 1:10
N.N.pred2 <- length(N.pred2)

Niminits <- list(log_C.time=log(C.time.init),log_C.sample=log(C.sample.init),
                 N.rep=N.rep.init,sdlog.C.sample=1,gamma0=5,gamma1=2,rho=0,beta0=mean(log(C.time.init)),
                 sdlog.C.time=1.5,alpha0=0.5) 
constants <- list(M=M,lag=lag,
                  J.time=J.time,K.time=K.time,eta.q=eta.q,
                  N.pred=N.pred,N.N.pred=N.N.pred,
                  N.pred2=N.pred2,N.N.pred2=N.N.pred2)

Nimdata <- list(y.time=y.time,y.det=y.det,volume.rep=volume.rep) 

# set parameters to monitor
parameters <- c('sdlog.C.time','alpha0','sdlog.C.sample','gamma0','gamma1','rho','beta0','VPC')
parameters2 <- c('C.time','C.sample','N.rep','y.pred','p.y.pred','log.mu.time',
                 'sdlog.C.time','sdlog.C.sample')

n.thin <- 5
n.thin2 <- 25

start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)
conf <- configureMCMC(Rmodel,monitors=parameters, thin=n.thin,
                      monitors2=parameters2,thin2=n.thin2,useConjugacy = FALSE)

#try blocking parameters if posteriors are strongly correlated
# conf$removeSampler(c("gamma0","gamma1"))
# conf$addSampler(target = c("gamma0","gamma1"),type = 'RW_block',control = list(adaptive=TRUE),silent = TRUE)

# conf$removeSampler(c("logit_rho","sdlog.C.time"))
# conf$addSampler(target = c("logit_rho","sdlog.C.time"),type = 'RW_block',control = list(adaptive=TRUE),silent = TRUE)

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=10) #this will run in R, used for debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model
start.time2 <- Sys.time()
Cmcmc$run(10000,reset=FALSE) #short run for demonstration. Can run repeatedly to get more posterior samples.1end.time <- Sys.time()
end.time <- Sys.time()
end.time - start.time  # total time for compilation, replacing samplers, and fitting
end.time - start.time2 # post-compilation run time

mvSamples = as.matrix(Cmcmc$mvSamples)
mvSamples2 <- as.matrix(Cmcmc$mvSamples2)

burnin <- 500
plot(coda::mcmc(mvSamples[-c(1:burnin),]))

#VPCs are proportion of variance in copies at rep, sample, and site-level (ignoring observation model)

par(mfrow=c(1,1),ask=FALSE)
burnin2 <- burnin*(n.thin/n.thin2)

###Process posteriors in mvSamples2
##Site and sample concentrations
C.time.idx <- grep("C.time",colnames(mvSamples2))
sdlog.C.time.idx <- grep("sdlog.C.time",colnames(mvSamples2))
C.time.idx <- C.time.idx[-which(C.time.idx==sdlog.C.time.idx)]
C.sample.idx <- grep("C.sample",colnames(mvSamples2))
sd.sample.idx <- grep("sdlog.C.sample",colnames(mvSamples2))
C.sample.idx <- setdiff(C.sample.idx,sd.sample.idx)
C.time.log.means <- (colMeans(log(mvSamples2[-c(1:burnin2),C.time.idx])))
C.time.log.HPDs <- matrix(0,length(C.time.idx),2)
for(t in 1:M){
  C.time.log.HPDs[t,] <- HPDinterval(mcmc(log(mvSamples2[-c(1:burnin2),C.time.idx[t]])))
}
C.sample.log.means <- colMeans(log(mvSamples2[-c(1:burnin2),C.sample.idx]))
C.sample.log.means2D <- matrix(C.sample.log.means,M,J)
C.sample.log.HPDs <- array(NA,dim=c(M,J,2))
idx <- 1
for(j in 1:max(J.time)){
  for(t in 1:M){
    if(J.time[t]>=j){
      C.sample.log.HPDs[t,j,] <- HPDinterval(mcmc(log(mvSamples2[-c(1:burnin2),C.sample.idx[idx]])))
    }
    idx <- idx+1
  }
}

log.mu.idx <- grep("log.mu.time",colnames(mvSamples2))
log.mu.time.est <- colMeans(mvSamples2[-c(1:burnin2),log.mu.idx])
log.mu.time.HPD <- HPDinterval(mcmc(mvSamples2[-c(1:burnin2),log.mu.idx]))

#Copies, true, thinned, available
N.rep.idx <- grep("N.rep",colnames(mvSamples2))
# plot(mcmc(mvSamples2[-c(1:burnin2),N.rep.idx]))
N.rep.est <- colMeans(coda::mcmc(mvSamples2[-c(1:burnin2),N.rep.idx]))
N.rep.est <- round(MCMCglmm::posterior.mode(coda::mcmc(mvSamples2[-c(1:burnin2),N.rep.idx]),adjust=1))
N.rep.HPDs <- HPDinterval(coda::mcmc(mvSamples2[-c(1:burnin2),N.rep.idx]))
N.rep.est3D <- array(N.rep.est,dim=c(M,J,K))
N.rep.lower <- array(N.rep.HPDs[,1],dim=c(M,J,K))
N.rep.upper <- array(N.rep.HPDs[,2],dim=c(M,J,K))


library(RColorBrewer)
cols <- brewer.pal(9,"Set3")

#raw data plot
# times3D <- array(times,dim=c(M,J,K))
# plot(log(y.time)~times3D,main="Raw Copy Number Measurements through Time",pch=16,col=cols[5],
#      xlab="Time (hours)",ylab="log10(Concentration)")


###Concentration Plot###
offset.tmp <- seq(-0.04,0.04,length.out=J)
#log10 conversions for labels
tmp <- range(C.sample.log.HPDs,na.rm=TRUE)/log(10,base=exp(1))
tmp[1] <- floor(tmp[1])
tmp[2] <- ceiling(tmp[2])
plot.range1 <- seq(tmp[1],tmp[2],by=0.5)
plot.range2 <- plot.range1/log(exp(1),base=10)
plot.range3 <- range(C.sample.log.HPDs,na.rm=TRUE)
plot.range3[2] <- plot.range3[2] + 0.1*diff(plot.range3) #make some room for legend
plot(log.mu.time.est~times,col=rgb(139,0,0,100,maxColorValue=255),type="l",lwd=2,
     ylab="Concentration (log 10 copies/liter)",ylim=plot.range3,
     xlab="Time (hours)",main="Site and Sample Concentration Through Time",yaxt='n')
axis(2,at=plot.range2,labels=plot.range1)
lines(log.mu.time.HPD[,1]~times,col=rgb(139,0,0,200,maxColorValue=255),lty=2)
lines(log.mu.time.HPD[,2]~times,col=rgb(139,0,0,200,maxColorValue=255),lty=2)
for(t in 1:M){
  lines(rep(times[t],2),C.time.log.HPDs[t,1:2],lwd=2)
  for(j in 1:J.time[t]){
    offset=sample(offset.tmp,replace=FALSE)
    points(C.sample.log.means2D[t,j]~c(times[t]+offset[j]),pch=16,col=cols[5])
    lines(x=rep(times[t]+offset[j],2),y=C.sample.log.HPDs[t,j,1:2],col=cols[5])
  }
}
points(C.time.log.means~times,pch=16)
legend("topleft",legend=c("Mean Site C","Temporal Site C","Sample C"),lwd=c(2,1,1),
       pch=c(NA,16,16),col=c(rgb(139,0,0,200,maxColorValue=255),"black",cols[5]))

#variance components plot
library(vioplot)
vioplot(mvSamples[-c(1:burnin),3],mvSamples[-c(1:burnin),2],mvSamples[-c(1:burnin),1],xaxt='n',
        main="Variance Decomposition of Copies Allocated to PCR Replicates",ylab="Proportion of Variance")
axis(1,at=1:3,labels=c("Temporal","Sampling","Pipetting"))
mtext("Source of Variance",side=1,line=3,at=2)


###N.rep plots###
N.rep.cover <- array(NA,dim=c(M,J,K))
plot(log(N.rep.est3D)~log(y.time),pch=16,main="Measured vs. True Log10 Copies per Replicate",
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

for(t in 1:M){
  for(j in 1:J.time[t]){
    for(k in 1:K.time[t,j]){
      N.rep.cover[t,j,k] <- 1*(N.rep.lower[t,j,k]<=y.time[t,j,k]&N.rep.upper[t,j,k]>=y.time[t,j,k])
      if(N.rep.cover[t,j,k]==1){
        lines(x=rep(log(y.time[t,j,k]),2),y=c(log(N.rep.lower[t,j,k]),log(N.rep.upper[t,j,k])))
      }else{
        lines(x=rep(log(y.time[t,j,k]),2),y=c(log(N.rep.lower[t,j,k]),log(N.rep.upper[t,j,k])),col="red")
      }
    }
  }
}
points(log(N.rep.est3D)~log(y.time),pch=16)
abline(0,1)
abline(h=log(eta.q),col="darkblue")

mean(N.rep.cover,na.rm=TRUE)


pred.idx <- grep("y.pred",colnames(mvSamples2))
pred.p.idx <- grep("p.y.pred",colnames(mvSamples2))
pred.sd.idx <- grep("sdlog.y.pred",colnames(mvSamples2))
pred.y.idx <- setdiff(pred.idx,c(pred.p.idx,pred.sd.idx))

#p.y.preds
p.y.pred.est <- colMeans(mvSamples2[-c(1:burnin2),pred.p.idx])
p.y.pred.HPD <- HPDinterval(mcmc(mvSamples2[-c(1:burnin2),pred.p.idx]))
plot(p.y.pred.est~N.pred2,pch=16,ylim=c(0,1),xlim=range(N.pred2),
     main='Replicate Detection Probability Vs. Available Copy Number',
     xlab="Available Copy Number",ylab="Replicate Detection Probability")
for(t in 1:N.N.pred2){
  lines(rep(N.pred2[t],2),p.y.pred.HPD[t,1:2])
}

#Plot measurement error distribution for each true copy number in rep
#just plotting results for first 8 N.preds below

library(vioplot)
vioplot(log(mvSamples2[-c(1:burnin2),pred.y.idx[1]]),
        log(mvSamples2[-c(1:burnin2),pred.y.idx[2]]),
        log(mvSamples2[-c(1:burnin2),pred.y.idx[3]]),
        log(mvSamples2[-c(1:burnin2),pred.y.idx[4]]),
        log(mvSamples2[-c(1:burnin2),pred.y.idx[5]]),
        log(mvSamples2[-c(1:burnin2),pred.y.idx[6]]),
        log(mvSamples2[-c(1:burnin2),pred.y.idx[7]]),
        log(mvSamples2[-c(1:burnin2),pred.y.idx[8]]))
abline(h=log(eta.q),lwd=2,col="darkblue")
vioplot((mvSamples2[-c(1:burnin2),pred.y.idx[1]]),
        (mvSamples2[-c(1:burnin2),pred.y.idx[2]]),
        (mvSamples2[-c(1:burnin2),pred.y.idx[3]]),
        (mvSamples2[-c(1:burnin2),pred.y.idx[4]]),
        (mvSamples2[-c(1:burnin2),pred.y.idx[5]]),
        (mvSamples2[-c(1:burnin2),pred.y.idx[6]]),
        (mvSamples2[-c(1:burnin2),pred.y.idx[7]]),
        (mvSamples2[-c(1:burnin2),pred.y.idx[8]]),
        main="Measurement Distribution for Replicate Copy Number",xlab="Replicate Copy Number",ylab="Measurement")
abline(h=eta.q,lwd=2,col="darkblue")
