library(nimble)
library(TeachingDemos)
library(MCMCglmm)
library(ggplot2)
library(qqplotr)


## Generalized Gama distribution
dGG <- function(t, scale, shape, tau, log = FALSE){
  val <- log(tau) - (shape*tau)*log(scale) - lgamma(shape) + ((shape*tau) - 1)*log(t) -
    (t/scale)^tau
  if(log) return(val) else return(exp(val))
}

pGG <- function(t, scale, shape, tau, log.p = FALSE){
  val <- pgamma( (t/scale)^tau, shape = shape, log.p = TRUE) 
  if(log.p) return(val) else return(exp(val))
}

rGG <- function(n, scale, shape, tau){
  p <- runif(n)
  out <- qgamma(p, shape = shape)^(1/tau)*scale
  return(as.vector(out))
}

## mixIGG distribution
dmixUGG <- function(x, kappa1=0.2, kappa2=0.2, gama=2, theta=0.5,tau=2, q=0.5, log = FALSE) {
  n_arg <- max(length(x), length(kappa1), length(kappa2), length(gama), length(theta), length(tau), length(q))
  x <- rep(x, length = n_arg)
  kappa1 <- rep(kappa1, length = n_arg)
  kappa2 <- rep(kappa2, length = n_arg)
  gama <- rep(gama, length = n_arg)
  theta <- rep(theta, length = n_arg)
  tau <- rep(tau, length = n_arg)
  q <- rep(q, length = n_arg)
  h1<- ((1-kappa1)/kappa1)*(1/qgamma(1-q, gama))^(1/tau)
  h2<- ((1-kappa2)/kappa2)*(1/qgamma(1-q, gama))^(1/tau)
  d <- theta *(1/x^2)*dGG((1-x)/x,scale=h1,shape=gama,tau=tau) +
    (1 - theta)*(1/x^2)*dGG((1-x)/x,scale=h2,shape=gama,tau=tau) 
  
  if (log== F) {
    return(d)
  }
  else {
    return(log(d))
  }
}

pmixUGG <- function(x, kappa1=0.2, kappa2=0.2, gama=2, theta=0.5,tau=2, q=0.5, lower.tail = TRUE, 
                      log.p = FALSE) {
  n_arg <- max(length(x), length(kappa1), length(kappa2), length(gama), length(theta), length(tau), length(q))
  x <- rep(x, length = n_arg)
  kappa1 <- rep(kappa1, length = n_arg)
  kappa2 <- rep(kappa2, length = n_arg)
  gama <- rep(gama, length = n_arg)
  theta <- rep(theta, length = n_arg)
  tau <- rep(tau, length = n_arg)
  if (!lower.tail) {
    x <- 1 - x
  }
  p <- as.vector(t(mapply(function(x, kappa1, kappa2, gama, theta, tau, q){
    integrate(function(x) dmixUGG(x,kappa1,kappa2,gama,theta,tau,q),0,x)$value}, x, kappa1, kappa2, gama, theta, tau, q)))
  if (!log.p) {
    return(p)
  }
  else {
    return(log(p))
  }
}

qmixUGG <- function(p, kappa1=0.2, kappa2=0.2, gama=2, theta=0.5,tau=2, q=0.5, lower.tail = TRUE, 
                      log.p = FALSE) {
  n_arg <- max(length(p), length(kappa1), length(kappa2), length(gama), length(theta), length(tau), length(q))
  p <- rep(p, length = n_arg)
  kappa1 <- rep(kappa1, length = n_arg)
  kappa2 <- rep(kappa2, length = n_arg)
  gama <- rep(gama, length = n_arg)
  theta <- rep(theta, length = n_arg)
  tau <- rep(tau, length = n_arg)
  p <- rep(p, length = n_arg)
  if(log.p) p <- exp(p)
  out <- as.vector(t(mapply(function(p, kappa1, kappa2, gama, theta, tau, q){
    stats::uniroot(f = function(x)
      (pmixUGG(x, kappa1, kappa2, gama, theta, tau, q, lower.tail) - p),
      interval = c(0.0001, 0.9999), maxiter = 100000, tol = 1e-30)$root},
    p, kappa1, kappa2, gama, theta, tau, q)))
  return(out)
}


rmixUGG <- function(n, kappa1=0.2, kappa2=0.2, gama=2, theta=0.5,tau=2, q=0.5) {
  kappa1 <- rep(kappa1, length = n)
  kappa2 <- rep(kappa2, length = n)
  theta <- rep(theta, length = n)
  gama <- rep(gama, length = n)
  tau <- rep(tau, length = n)
  q <- rep(q, length = n)
  h1<- ((1-kappa1)/kappa1)*(1/qgamma(1-q, gama))^(1/tau)
  h2<- ((1-kappa2)/kappa2)*(1/qgamma(1-q, gama))^(1/tau)
  u <- runif(n)
  out=0
  for(i in 1:n){
    if(u[i]<theta[i]){
      v1 <- rGG(1, h1[i],gama[i],tau[i])
      out[i]<-(1/(1+v1))
    }else{
      v2 <- rGG(1, h2[i],gama[i],tau[i])
      out[i]<-(1/(1+v2))
    }
  }
  return(out)
}


                          
dggamma <- nimbleFunction(
  run = function(x = double(), scale = double(), 
                 shape = double(), tau = double(), log = logical(0, default = 0)) {
    returnType(double())
    d <- log(tau) - (shape*tau)*log(scale) - lgamma(shape) + ((shape*tau) - 1)*log(x) -
      (x/scale)^tau
    if (log) return(d)
    else return(exp(d))
  })


registerDistributions(list(
  dggamma = list(
    BUGSdist = "dggamma(scale, shape, tau)",
    types = c('value = double()', 'scale = double()', 'shape = double()', 'tau = double()')
  )))

###################### Function to fit the mixUGG regression model #######################

## There are three available link functions: logit, probit and cloglog, use kappa.link
## to set this. With nchains you set the number of MCMC chains, with thin the user can
## specify how much the MCMC chains should be thinned out before storing them. niter
## is used to set the iteraction numbers and nburnin is used to set the burning-in.
## The argument q is to define which quantile will be modelled 0<q<1.                          
                          
mixUGG=function(kappa1.formula=formula,
                kappa2.formula=~1,data=NULL,
                kappa.link="logit",q=0.5, nchains = 1,thin=50, 
                niter = 90000, nburnin = 40000){
  
  mcall <- if(is.null(data)){ terms(kappa1.formula) 
  }else terms(kappa1.formula, data = data)
  X <- if(is.null(data)){ model.matrix(kappa1.formula) 
  }else model.matrix(kappa1.formula, data = data)  
  X2 <- if(is.null(data)){ model.matrix(kappa2.formula) 
  }else model.matrix(kappa2.formula, data = data)  
  y <- if(is.null(data)){ model.frame(mcall)[,1] 
  }else model.frame(mcall, data = data)[,1]
  
  if(ncol(X2)==1){X2=X}
  n=length(y)
  p=ncol(X)
  p2=ncol(X2)
  
  if(kappa.link=="logit"){
    pumpCode <- nimbleCode({
      for(k in 1:p){
        alpha[k] ~ dnorm(0, tau=1/1000)
      }
      
      for(k2 in 1:p2){
        beta[k2] ~ dnorm(0, tau=1/1000)
      }
      theta ~ dunif(0,1)
      tau ~ dinvgamma(0.1,0.1)
      gama ~ dinvgamma(0.1,0.1)
      for (i in 1:n){
        z[i] ~ dbern(theta)
        a[i] <- inprod(X[i,1:p],alpha[1:p])
        b[i] <- inprod(X2[i,1:p2],beta[1:p2])
        kappa1[i] <- ilogit(a[i])
        kappa2[i] <- ilogit(b[i])
        h1[i]<- ((1-kappa1[i])/kappa1[i])*(1/qgamma(1-q, gama))^(1/tau)
        h2[i]<- ((1-kappa2[i])/kappa2[i])*(1/qgamma(1-q, gama))^(1/tau)
        l1[i] <- log((1/y[i]^2))+dggamma((1-y[i])/y[i], shape = gama, scale=h1[i],tau=tau, log=TRUE)
        l2[i] <- log((1/y[i]^2))+dggamma((1-y[i])/y[i], shape = gama, scale=h2[i],tau=tau, log=TRUE) 
        zeros[i] ~ dpois(C - equals(z[i],1)*l1[i] - equals(z[i],0)*l2[i])
      }
    })
  }
  if(kappa.link=="probit"){
    pumpCode <- nimbleCode({
      for(k in 1:p){
        alpha[k] ~ dnorm(0, tau=1/1000)
      }
      
      for(k2 in 1:p2){
        beta[k2] ~ dnorm(0, tau=1/1000)
      }
      theta ~ dunif(0,1)
      tau ~ dinvgamma(0.1,0.1)
      gama ~ dinvgamma(0.1,0.1)
      for (i in 1:n){
        z[i] ~ dbern(theta)
        a[i] <- inprod(X[i,1:p],alpha[1:p])
        b[i] <- inprod(X2[i,1:p2],beta[1:p2])
        kappa1[i] <- iprobit(a[i])
        kappa2[i] <- iprobit(b[i])
        h1[i]<- ((1-kappa1[i])/kappa1[i])*(1/qgamma(1-q, gama))^(1/tau)
        h2[i]<- ((1-kappa2[i])/kappa2[i])*(1/qgamma(1-q, gama))^(1/tau)
        l1[i] <- log((1/y[i]^2))+dggamma((1-y[i])/y[i], shape = gama, scale=h1[i],tau=tau, log=TRUE)
        l2[i] <- log((1/y[i]^2))+dggamma((1-y[i])/y[i], shape = gama, scale=h2[i],tau=tau, log=TRUE) 
        zeros[i] ~ dpois(C - equals(z[i],1)*l1[i] - equals(z[i],0)*l2[i])
      }
    })
  }
  if(kappa.link=="cloglog"){
    
    pumpCode <- nimbleCode({
      for(k in 1:p){
        alpha[k] ~ dnorm(0, tau=1/1000)
      }
      
      for(k2 in 1:p2){
        beta[k2] ~ dnorm(0, tau=1/1000)
      }
      theta ~ dunif(0,1)
      tau ~ dinvgamma(0.1,0.1)
      gama ~ dinvgamma(0.1,0.1)
      for (i in 1:n){
        z[i] ~ dbern(theta)
        a[i] <- inprod(X[i,1:p],alpha[1:p])
        b[i] <- inprod(X2[i,1:p2],beta[1:p2])
        kappa1[i] <- icloglog(a[i])
        kappa2[i] <- icloglog(b[i])
        h1[i]<- ((1-kappa1[i])/kappa1[i])*(1/qgamma(1-q, gama))^(1/tau)
        h2[i]<- ((1-kappa2[i])/kappa2[i])*(1/qgamma(1-q, gama))^(1/tau)
        l1[i] <- log((1/y[i]^2))+dggamma((1-y[i])/y[i], shape = gama, scale=h1[i],tau=tau, log=TRUE)
        l2[i] <- log((1/y[i]^2))+dggamma((1-y[i])/y[i], shape = gama, scale=h2[i],tau=tau, log=TRUE) 
        zeros[i] ~ dpois(C - equals(z[i],1)*l1[i] - equals(z[i],0)*l2[i])
      }
    })
    
  }
  
  pumpConsts <- list(n = n, p = p, p2 = p2, q = q, 
                     C=10000)
  
  inicial<- list(beta = rep(0.1,p2), alpha = rep(0.1,p),
                 theta=0.1, z=rep(1,n), gama = 1, tau = 1)
  
  pumpData <- list(y = y, X=X, X2=X2, zeros=rep(0,n))
  
  mcmc.out <- nimbleMCMC(code = pumpCode, constants = pumpConsts,
                         data = pumpData, inits = inicial,
                         nchains = nchains,thin=thin, niter = niter, 
                         nburnin = nburnin,
                         monitors = c('alpha','beta','gama','tau','theta','z'),
                         summary = TRUE)
  
  return(list(mcmc.out=mcmc.out,X1=X1,X2=X2,
               y=y,q=q,kappa.link=kappa.link))
}


modelsummary=function(mod){
  p=ncol(mod$X)
  p2=ncol(mod$X2)
  prop=matrix(,nrow=p+p2+3,ncol=6)
  for(j in 1:(p+p2+3)){
    prop[j,1]=round(posterior.mode(as.mcmc(mod$mcmc.out$samples[,j])),2)
    prop[j,2]=round(median(mod$mcmc.out$samples[,j]),2)
    prop[j,3]=round(mean(mod$mcmc.out$samples[,j]),2)
    prop[j,4]=round(var(mod$mcmc.out$samples[,j]),4)
    prop[j,5]=length(which(mod$mcmc.out$samples[,j]>0))/length(mod$mcmc.out$samples[,1])
    prop[j,6]=paste("[",round(emp.hpd(mod$mcmc.out$samples[,j])[1],3),";",round(emp.hpd(mod$mcmc.out$samples[,j])[2],3),"]")
  }
  colnames(prop)<-c("Mode","Median","Mean","Var","P(par>0)","HPD")
  rownames(prop)<-c(colnames(mod$X),colnames(mod$X2),"delta","tau","theta")
  return(prop)
}

modelresiduals=function(mod,estimator="Mode",conf=0.95){
  X=mod$X
  X2=mod$X2
  q=mod$q
  kappa.link=mod$kappa.link
  if(estimator=="Mode"){
    if(kappa.link=="logit"){
      kappa1est=ilogit(X%*%posterior.mode(as.mcmc(mod$mcmc.out$samples[,1:p])))
      kappa2est=ilogit(X2%*%posterior.mode(as.mcmc(mod$mcmc.out$samples[,(p+1):(p+p2)])))
    }
    if(kappa.link=="probit"){
      kappa1est=iprobit(X%*%posterior.mode(as.mcmc(mod$mcmc.out$samples[,1:p])))
      kappa2est=iprobit(X2%*%posterior.mode(as.mcmc(mod$mcmc.out$samples[,(p+1):(p+p2)])))
    }
    if(kappa.link=="cloglog"){
      kappa1est=icloglog(X%*%posterior.mode(as.mcmc(mod$mcmc.out$samples[,1:p])))
      kappa2est=icloglog(X2%*%posterior.mode(as.mcmc(mod$mcmc.out$samples[,(p+1):(p+p2)])))
    }
    gamaest=posterior.mode(as.mcmc(mod$mcmc.out$samples[,(p+p2+1)]))
    thetaest=posterior.mode(as.mcmc(mod$mcmc.out$samples[,(p+p2+3)]))
    tauest=posterior.mode(as.mcmc(mod$mcmc.out$samples[,(p+p2+2)]))
    resiquan=qnorm(pmixUGG(y,kappa1est,kappa2est,
                           gamaest,thetaest,tauest,q=q))
  }
  if(estimator=="Median"){
    if(kappa.link=="logit"){
      kappa1est=ilogit(X%*%median(mod$mcmc.out$samples[,1:p]))
      kappa2est=ilogit(X2%*%median(mod$mcmc.out$samples[,(p+1):(p+p2)]))
    }
    if(kappa.link=="probit"){
      kappa1est=iprobit(X%*%median(mod$mcmc.out$samples[,1:p]))
      kappa2est=iprobit(X2%*%median(mod$mcmc.out$samples[,(p+1):(p+p2)]))
    }
    if(kappa.link=="cloglog"){
      kappa1est=icloglog(X%*%median(mod$mcmc.out$samples[,1:p]))
      kappa2est=icloglog(X2%*%median(mod$mcmc.out$samples[,(p+1):(p+p2)]))
    }
    gamaest=median(mod$mcmc.out$samples[,(p+p2+1)])
    thetaest=median(mod$mcmc.out$samples[,(p+p2+3)])
    tauest=median(mod$mcmc.out$samples[,(p+p2+2)])
    resiquan=qnorm(pmixUGG(y,kappa1est,kappa2est,
                           gamaest,thetaest,tauest,q=q))
  }
  if(estimator=="Mean"){
    if(kappa.link=="logit"){
      kappa1est=ilogit(X%*%mean(mod$mcmc.out$samples[,1:p]))
      kappa2est=ilogit(X2%*%mean(mod$mcmc.out$samples[,(p+1):(p+p2)]))
    }
    if(kappa.link=="probit"){
      kappa1est=iprobit(X%*%mean(mod$mcmc.out$samples[,1:p]))
      kappa2est=iprobit(X2%*%mean(mod$mcmc.out$samples[,(p+1):(p+p2)]))
    }
    if(kappa.link=="cloglog"){
      kappa1est=icloglog(X%*%mean(mod$mcmc.out$samples[,1:p]))
      kappa2est=icloglog(X2%*%mean(mod$mcmc.out$samples[,(p+1):(p+p2)]))
    }
    gamaest=mean(mod$mcmc.out$samples[,(p+p2+1)])
    thetaest=mean(mod$mcmc.out$samples[,(p+p2+3)])
    tauest=mean(mod$mcmc.out$samples[,(p+p2+2)])
    resiquan=qnorm(pmixUGG(y,kappa1est,kappa2est,
                           gamaest,thetaest,tauest,q=q))
  }
  smpbr <- data.frame(norm = resiquan)
  s1<-ggplot(data = smpbr, mapping = aes(sample = norm))  +
    stat_qq_band(conf = conf) +
    stat_qq_line() +
    stat_qq_point() +
    labs(x = "Theoretical Quantiles", y = "Quantile Residuals") +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "grey"),
      text=element_text(size=25,family="serif")
    )
  print(s1)
}
                          
