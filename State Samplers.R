C.sample.Sampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    calcNodes1 <- control$calcNodes1
    calcNodes2 <- control$calcNodes2
    calcNodes3 <- control$calcNodes3
    calcNodes.p.w <- control$calcNodes.p.w
    # calcNodes.p.w.lifted <- control$calcNodes.p.w.lifted
    calcNodes.theta.thin <- control$calcNodes.theta.thin
    i <- control$i
    j <- control$j
    K <- control$K
  },
  run = function() {
    log.prop.for <- log.prop.back <- 0
    #get initial logProb for relevant nodes
    lp.initial <- model$calculate(calcNodes1) + model$calculate(calcNodes2)
    #get backprobs first
    log.prop.back <- log.prop.back + dnorm(model$log_C.sample[i,j],log(model$C.site[i]),model$sdlog.C.sample[i],log=TRUE) #prob of current C.sample
    log.prop.back <- log.prop.back + dnorm(model$unit.inhibit[i,j],0,sd=1,log=TRUE) #prob of current unit.inhibit
    for(k in 1:K){
      log.prop.back <- log.prop.back + dpois(model$N.rep[i,j,k],model$lambda.rep[i,j],log=TRUE) #prob of current N.rep|current lambda.rep
      log.prop.back <- log.prop.back + dbinom(model$w[i,j,k],1,model$p.w[i,j,k],log=TRUE) #prob of w|current p.w
      log.prop.back <- log.prop.back + dbinom(model$N.rep.thin[i,j,k],model$N.rep[i,j,k],model$theta.thin[i,j,k],log=TRUE) #prob of current N.rep.thin|N.rep
    }
    #propose new C.sample
    model$log_C.sample[i,j] <<- rnorm(1,log(model$C.site[i]),model$sdlog.C.sample[i])
    log.prop.for <- log.prop.for + dnorm(model$log_C.sample[i,j],log(model$C.site[i]),model$sdlog.C.sample[i],log=TRUE) #prob propose this C.sample
    #propose new unit.inhibit
    model$unit.inhibit[i,j] <<- rnorm(1,0,sd=1)
    log.prop.for <- log.prop.for + dnorm(model$unit.inhibit[i,j],0,sd=1,log=TRUE)
    #update lambda.rep (theta.thin updated below before N.thin)
    model$calculate(calcNodes1)
    #for each rep,
    for(k in 1:K){
      #propose new N.rep
      model$N.rep[i,j,k] <<- rpois(1,model$lambda.rep[i,j])
      #update p.w for new N.rep
      model$calculate(calcNodes.p.w[k])
      #update theta.thin for new N.rep (and unit thin)
      model$calculate(calcNodes.theta.thin[k])
      #propose new w state
      model$w[i,j,k] <<- rbinom(1,1,model$p.w[i,j,k])
      #propose new N.thin state from new N.rep state. 
      #Assuming thinning rate doesn't depend on N.rep. Or we'd need to propose that, too, before N.thin
      model$N.rep.thin[i,j,k] <<- rbinom(1,model$N.rep[i,j,k],model$theta.thin[i,j,k])
      log.prop.for <- log.prop.for + dpois(model$N.rep[i,j,k],model$lambda.rep[i,j],log=TRUE) # prob propose this N.rep|proposed lambda.rep
      log.prop.for <- log.prop.for + dbinom(model$w[i,j,k],1,model$p.w[i,j,k],log=TRUE) #prob propose this w|proposed w
      log.prop.for <- log.prop.for + dbinom(model$N.rep.thin[i,j,k],model$N.rep[i,j,k],model$theta.thin[i,j,k],log=TRUE) #prop propose this N.rep.thin|proposed N.rep
    }
    lp.proposed <- model$calculate(calcNodes2) + model$calculate(calcNodes1) #add on log_C.sample likelihood
    log_MH_ratio <- (lp.proposed + log.prop.back) - (lp.initial + log.prop.for)
    accept <- decide(log_MH_ratio)
    if(accept) {
      copy(from = model, to = mvSaved, row = 1, nodes = calcNodes3, logProb = TRUE)
    }else{
      copy(from = mvSaved, to = model, row = 1, nodes = calcNodes3, logProb = TRUE)    
    }   
  },
  methods = list( reset = function () {} )
)

wSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    calcNodes <- control$calcNodes
    calcNodes.p.w <- control$calcNodes.p.w
    calcNodes.theta.thin <- control$calcNodes.theta.thin
    i <- control$i
    j <- control$j
    k <- control$k
  },
  run = function() {
    #get initial logProb for relevant nodes
    lp.initial <- model$calculate(calcNodes)
    #get backwards proposal probs first
    log.prop.back <- dpois(model$N.rep[i,j,k],model$lambda.rep[i,j],log=TRUE) + #prob of current N.rep
      dbinom(model$w[i,j,k],1,model$p.w[i,j,k],log=TRUE) + #prob of proposing current w|current N.rep
      dbinom(model$N.rep.thin[i,j,k],model$N.rep[i,j,k],model$theta.thin[i,j,k],log=TRUE) #prob of current N.thin|current N.rep
    #propose new N.rep state first because p.w depends on N.rep
    model$N.rep[i,j,k] <<- rpois(1,model$lambda.rep[i,j])
    #update p.w for new N.rep
    model$calculate(calcNodes.p.w)
    #update theta.thin for new N.rep
    model$calculate(calcNodes.theta.thin)
    #propose new w state
    model$w[i,j,k] <<- rbinom(1,1,model$p.w[i,j,k])
    #propose new N.thin state from new N.rep state. 
    model$N.rep.thin[i,j,k] <<- rbinom(1,model$N.rep[i,j,k],model$theta.thin[i,j,k])
    lp.proposed <- model$calculate(calcNodes)
    #get forward proposal probs
    log.prop.for <- dpois(model$N.rep[i,j,k],model$lambda.rep[i,j],log=TRUE) +  #prob of proposing N.rep
      dbinom(model$w[i,j,k],1,model$p.w[i,j,k],log=TRUE) + #prob of proposing w|proposed N.rep that changes p.w
      dbinom(model$N.rep.thin[i,j,k],model$N.rep[i,j,k],model$theta.thin[i,j,k],log=TRUE) #prob of proposing N.thin|proposed N.rep
    
    log_MH_ratio <- (lp.proposed + log.prop.back) - (lp.initial + log.prop.for)
    accept <- decide(log_MH_ratio)
    if(accept) {
      copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    }else{
      copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
    }
  },
  methods = list( reset = function () {} )
)

#this one doesn't have to thread the needle from below (y measurement), but does from above (sample concentration)
wSampler2 <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    calcNodes <- control$calcNodes
    calcNodes.p.w <- control$calcNodes.p.w
    calcNodes.theta.thin <- control$calcNodes.theta.thin
    i <- control$i
    j <- control$j
    k <- control$k
  },
  run = function() {
    #get initial logProb for relevant nodes
    lp.initial <- model$calculate(calcNodes)
    #y measurement likelihood should not change!
    #1 if w=0, propose w=1 -> set N.rep.thin to N.rep, propose new N.rep(how exactly?), updates theta thin and p.w (shift up)
    #2. if w=1, propose w=0 -> set N.rep to N.rep.thin, updates theta thin and p.w, propose new N.rep.thin (shift down)
    if(model$w[i,j,k]==0){ #shift up
      #get backprobs first
      log.prop.back <- dbinom(model$N.rep.thin[i,j,k],model$N.rep[i,j,k],model$theta.thin[i,j,k],log=TRUE) #prob of proposing N.rep.thin
      #propose
      model$w[i,j,k] <<- 1
      model$N.rep.thin[i,j,k] <<- model$N.rep[i,j,k]
      #how to propose N.rep? Can propose from above
      model$N.rep[i,j,k] <<- rpois(1,model$lambda.rep[i,j]) #what if not larger than N.rep.thin prop? Always rejected. this is fine.
      log.prop.for <- dpois(model$N.rep[i,j,k],model$lambda.rep[i,j],log=TRUE) #prob of proposing N.rep
      log.prop.for <- dpois(model$N.rep[i,j,k],model$lambda.rep[i,j],log=TRUE) #prob of proposing N.rep
      model$calculate(calcNodes.p.w)
      model$calculate(calcNodes.theta.thin)
    }else{ #shift down
      #get backprobs first
      log.prop.back <- dpois(model$N.rep[i,j,k],model$lambda.rep[i,j],log=TRUE) #prob of proposing N.rep
      #propose
      model$w[i,j,k] <<- 0
      model$N.rep[i,j,k] <<- model$N.rep.thin[i,j,k]
      model$calculate(calcNodes.p.w)
      model$calculate(calcNodes.theta.thin)
      model$N.rep.thin[i,j,k] <<- rbinom(1,model$N.rep[i,j,k],model$theta.thin[i,j,k])
      log.prop.for <- dbinom(model$N.rep.thin[i,j,k],model$N.rep[i,j,k],model$theta.thin[i,j,k],log=TRUE) #prob of proposing N.rep.thin
    }
    
    lp.proposed <- model$calculate(calcNodes)
    log_MH_ratio <- (lp.proposed + log.prop.back) - (lp.initial + log.prop.for)
    accept <- decide(log_MH_ratio)
    if(accept) {
      copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    }else{
      copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
    }
  },
  methods = list( reset = function () {} )
)