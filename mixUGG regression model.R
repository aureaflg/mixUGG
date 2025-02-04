library(nimble)
library(TeachingDemos)
library(MCMCglmm)
library(ggplot2)
library(qqplotr)
library(ggpubr)
library(ggrepel)


## Generalized Gamma distribution
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

## mixUGG distribution
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
## is used to set the iteraction numbers and nburnin is used to set the burn-in.
## The argument q defines which quantile will be modelled (0<q<1).                          
                          
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
  
  return(list(mcmc.out=mcmc.out,X=X,X2=X2,
               y=y,q=q,kappa.link=kappa.link))
}


## Example of how to use the mixUGG function to fit a mixUGG
## quantile regression model:  
## mod=mixUGG(y~x1+x2, kappa2.formula=~x1+x2+x3, q=0.25)    
                          
## In mod, y is the response variable, x1, x2, and x3 are covariates
## y~x1+x2 is the linear predictor related to the mixture component 1,
## while kappa2.formula=~x1+x2+x3 is related to the mixture component 2.
## If kappa2.formula was not set by the user, it is assumed the same linear
## predictor for both components. The default link function is the logit link,
## and the modelled quantile in this case is the first quartile (default is q=0.5).                          

  

  

# Model fit summary
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

# This function returns the proportion of observations calssified in each component
# as well as the classification of each observation according to the fitted model                          
classification=function(mod){
  p=ncol(mod$X)
  p2=ncol(mod$X2)
  n=length(mod$y)
  prop1=table(mod$mcmc.out$summary[(p+p2+4):(p+p2+3+n),2])[1]/n
  prop2=table(mod$mcmc.out$summary[(p+p2+4):(p+p2+3+n),2])[2]/n
  aux=data.frame('Proportion of obs coming from comp 1'=prop1,'Proportion of obs coming from comp 2'=prop2)
  yclass=mod$mcmc.out$summary[(p+p2+4):(p+p2+3+n),2]
  yclass[yclass==1]=2
  yclass[yclass==0]=1
  return(list(proportion=aux,"classification of observations"=yclass))
}
# Example: classification(mod)$proportion returns the proportion and
# classification(mod)$"classification of observations" returns the classification
# of each observation: 1 for component 1 and 2 for component 2.                         

# Quantile Residuals                          
modelresiduals=function(mod,estimator="Mode",conf=0.95){
  X=mod$X
  X2=mod$X2
  q=mod$q
  y=mod$y
  n=length(y)
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

# Information criteria
criteria=function(mod){
  X=mod$X
  X2=mod$X2
  q=mod$q
  p=ncol(X)
  p2=ncol(X2)
  y=mod$y
  n=length(y)
  m=length(mod$mcmc.out$samples[,1])
  aux = matrix(0,m,n)
  kappa1es=kappa2es=gamaes=thetaes=taues=lppd=pwaic=den=0
  
  if(mod$kappa.link=="logit"){
    for(j in 1:m){
      kappa1estj=ilogit(X%*%mod$mcmc.out$sample[j,1:p])
      kappa2estj=ilogit(X2%*%mod$mcmc.out$sample[j,(p+1):(p+p2)])
      gamaestj=mod$mcmc.out$sample[j,(p+p2+1)]
      thetaestj=mod$mcmc.out$sample[j,(p+p2+3)]
      tauestj=mod$mcmc.out$sample[j,(p+p2+2)]
      for(i in 1:n){
        aux[j,i] = dmixUGG(y[i],kappa1estj[i],kappa2estj[i],
                           gamaestj,thetaestj,tauestj,q)
      }
    }
    for(i in 1:n){
      for(j in 1:m){
        kappa1es[j]=ilogit(X[i,]%*%mod$mcmc.out$samples[j,1:p])
        kappa2es[j]=ilogit(X2[i,]%*%mod$mcmc.out$samples[j,(p+1):(p+p2)])
        gamaes[j]=mod$mcmc.out$samples[j,(p+p2+1)]
        thetaes[j]=mod$mcmc.out$samples[j,(p+p2+3)]
        taues[j]=mod$mcmc.out$samples[j,(p+p2+2)]
        den[j]=dmixUGG(y[i],kappa1es[j],kappa2es[j],gamaes[j],
                       thetaes[j],taues[j],q)
      }
      lppd[i]=log(mean(den))
      pwaic[i]=var(log(den))
    }
    kappa1est=ilogit(X%*%mod$mcmc.out$summary[1:p])
    kappa2est=ilogit(X2%*%mod$mcmc.out$summary[(p+1):(p+p2)])
  }
  if(mod$kappa.link=="probit"){
    for(j in 1:m){
      kappa1estj=iprobit(X%*%mod$mcmc.out$sample[j,1:p])
      kappa2estj=iprobit(X2%*%mod$mcmc.out$sample[j,(p+1):(p+p2)])
      gamaestj=mod$mcmc.out$sample[j,(p+p2+1)]
      thetaestj=mod$mcmc.out$sample[j,(p+p2+3)]
      tauestj=mod$mcmc.out$sample[j,(p+p2+2)]
      for(i in 1:n){
        aux[j,i] = dmixUGG(y[i],kappa1estj[i],kappa2estj[i],
                           gamaestj,thetaestj,tauestj,q)
      }
    }
    for(i in 1:n){
      for(j in 1:m){
        kappa1es[j]=iprobit(X[i,]%*%mod$mcmc.out$samples[j,1:p])
        kappa2es[j]=iprobit(X2[i,]%*%mod$mcmc.out$samples[j,(p+1):(p+p2)])
        gamaes[j]=mod$mcmc.out$samples[j,(p+p2+1)]
        thetaes[j]=mod$mcmc.out$samples[j,(p+p2+3)]
        taues[j]=mod$mcmc.out$samples[j,(p+p2+2)]
        den[j]=dmixUGG(y[i],kappa1es[j],kappa2es[j],gamaes[j],
                       thetaes[j],taues[j],q)
      }
      lppd[i]=log(mean(den))
      pwaic[i]=var(log(den))
    }
    kappa1est=iprobit(X%*%mod$mcmc.out$summary[1:p])
    kappa2est=iprobit(X2%*%mod$mcmc.out$summary[(p+1):(p+p2)])
  }
  if(mod$kappa.link=="cloglog"){
    for(j in 1:m){
      kappa1estj=icloglog(X%*%mod$mcmc.out$sample[j,1:p])
      kappa2estj=icloglog(X2%*%mod$mcmc.out$sample[j,(p+1):(p+p2)])
      gamaestj=mod$mcmc.out$sample[j,(p+p2+1)]
      thetaestj=mod$mcmc.out$sample[j,(p+p2+3)]
      tauestj=mod$mcmc.out$sample[j,(p+p2+2)]
      for(i in 1:n){
        aux[j,i] = dmixUGG(y[i],kappa1estj[i],kappa2estj[i],
                           gamaestj,thetaestj,tauestj,q)
      }
    }
    for(i in 1:n){
      for(j in 1:m){
        kappa1es[j]=icloglog(X[i,]%*%mod$mcmc.out$samples[j,1:p])
        kappa2es[j]=icloglog(X2[i,]%*%mod$mcmc.out$samples[j,(p+1):(p+p2)])
        gamaes[j]=mod$mcmc.out$samples[j,(p+p2+1)]
        thetaes[j]=mod$mcmc.out$samples[j,(p+p2+3)]
        taues[j]=mod$mcmc.out$samples[j,(p+p2+2)]
        den[j]=dmixUGG(y[i],kappa1es[j],kappa2es[j],gamaes[j],
                       thetaes[j],taues[j],q)
      }
      lppd[i]=log(mean(den))
      pwaic[i]=var(log(den))
    }
    kappa1est=icloglog(X%*%mod$mcmc.out$summary[1:p])
    kappa2est=icloglog(X2%*%mod$mcmc.out$summary[(p+1):(p+p2)])
  }
  
  cpo = 1/colMeans(1/aux)
  LPML=sum(log(cpo)) 
  gamaest=mod$mcmc.out$summary[(p+p2+1),1]
  thetaest=mod$mcmc.out$summary[(p+p2+3)]
  tauest=mod$mcmc.out$summary[(p+p2+2)]
  lvero=sum(dmixUGG(y,kappa1est,kappa2est,
                    gamaest,thetaest,tauest,q=q,log=T))
  EAIC=-2*lvero+2*(p+p2+3);EBIC=-2*lvero+(p+p2+3)*log(n)
  BARDEV = mean(colSums(t(-2*log(aux))))
  DEVBAR=-2*lvero
  DCI = DEVBAR + 2*(BARDEV-DEVBAR)
  WAIC=-2*(sum(lppd)-sum(pwaic))
  
  result=data.frame(EAIC,EBIC,DCI,WAIC,LPML)
  return(result)
}

# Bayesian Influence Anlysis
# Kullback-Leibler divergence                           
# if you want to highlight a point, use the argument ref to do so                          
kl=function(mod,ref=0.5){
  X=mod$X
  X2=mod$X2
  q=mod$q
  p=ncol(X)
  p2=ncol(X2)
  y=mod$y
  n=length(y)
  m=length(mod$mcmc.out$samples[,1])
  aux = matrix(0,m,n)
  if(mod$kappa.link=="logit"){
    for(j in 1:m){
      kappa1estj=ilogit(X%*%mod$mcmc.out$sample[j,1:p])
      kappa2estj=ilogit(X2%*%mod$mcmc.out$sample[j,(p+1):(p+p2)])
      gamaestj=mod$mcmc.out$sample[j,(p+p2+1)]
      thetaestj=mod$mcmc.out$sample[j,(p+p2+3)]
      tauestj=mod$mcmc.out$sample[j,(p+p2+2)]
      for(i in 1:n){
        aux[j,i] = dmixUGG(y[i],kappa1estj[i],kappa2estj[i],
                           gamaestj,thetaestj,tauestj,q)
      }
    }
  }
  if(mod$kappa.link=="probit"){
    for(j in 1:m){
      kappa1estj=iprobit(X%*%mod$mcmc.out$sample[j,1:p])
      kappa2estj=iprobit(X2%*%mod$mcmc.out$sample[j,(p+1):(p+p2)])
      gamaestj=mod$mcmc.out$sample[j,(p+p2+1)]
      thetaestj=mod$mcmc.out$sample[j,(p+p2+3)]
      tauestj=mod$mcmc.out$sample[j,(p+p2+2)]
      for(i in 1:n){
        aux[j,i] = dmixUGG(y[i],kappa1estj[i],kappa2estj[i],
                           gamaestj,thetaestj,tauestj,q)
      }
    }
  }
  if(mod$kappa.link=="cloglog"){
    for(j in 1:m){
      kappa1estj=icloglog(X%*%mod$mcmc.out$sample[j,1:p])
      kappa2estj=icloglog(X2%*%mod$mcmc.out$sample[j,(p+1):(p+p2)])
      gamaestj=mod$mcmc.out$sample[j,(p+p2+1)]
      thetaestj=mod$mcmc.out$sample[j,(p+p2+3)]
      tauestj=mod$mcmc.out$sample[j,(p+p2+2)]
      for(i in 1:n){
        aux[j,i] = dmixUGG(y[i],kappa1estj[i],kappa2estj[i],
                           gamaestj,thetaestj,tauestj,q)
      }
    }
  }
  
  cpo = 1/colMeans(1/aux)
  kpp=0
  for(i in 1:n){
    kpp[i]=colMeans(log(aux))[i]-log(cpo)[i]
  }
  gg<-ggplot() + geom_point(aes(x = 1:n, y = kpp)) +
    xlab("Subject index") + ylab("KL divergence") +
    geom_text(aes(x = 1:n, y = kpp, label=ifelse(kpp>ref,1:n,'')),
              hjust=1.5,vjust=1.5,size=6) +
    theme_bw() +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "grey"),
          text=element_text(size=25,family="serif"))
  print(gg)
}


# Estimation of parameters related to marginal quantile model
# The user can set the covariates to be considered in the fit,
# however, if is not set, the covariates will the same considered
# for component 1                          
marginal=function(mod,kappa3.formula=~1,data=NULL){
  X=mod$X
  X2=mod$X2
  q=mod$q
  p=ncol(X)
  p2=ncol(X2)
  
  if(kappa3.formula==~1){
    X3=X
  }else{
    if(is.null(data)){
      X3=model.matrix(kappa3.formula)
    }else{
      X3=model.matrix(kappa3.formula,data=data)
    }
  }
  
  m=length(mod$mcmc.out$samples[,1])
  gamaest=mod$mcmc.out$samples[,(p+p2+1)]
  thetaest=mod$mcmc.out$samples[,(p+p2+3)]
  tauest=mod$mcmc.out$samples[,(p+p2+2)]
  
  
  ff = function(par,kappa1e,kappa2e,gamaest, thetaest, tauest,q){
    dif=(pmixUGG(par, kappa1e,kappa2e,gamaest, thetaest, tauest,q) - q)^2
    return(dif)
  }
  
  aux2=matrix(,ncol=ncol(X3),nrow=m)
  
  if(mod$kappa.link=="logit"){
    for(j in 1:n){
      kappa1e[,j]<-ilogit(X[j,]%*%t(mod$mcmc.out$samples[,1:p]))
      kappa2e[,j]<-ilogit(X2[j,]%*%t(mod$mcmc.out$samples[,(p+1):(p+p2)]))
      for(i in 1:m){
        resu[i,j] <- optim(0.5, ff, method = "L-BFGS-B", kappa1e = kappa1e[i,j], 
                           kappa2e = kappa2e[i,j],gamaest = gamaest[i], thetaest = thetaest[i],
                           tauest = tauest[i], q = q,lower = c(0.0001), upper=c(0.9999))$par
      }
      
      print(j)
    }
    for(f in 1:m){
      mod=lm(logit(resu[f,])~X3[,-1])
      aux2[f,]=mod$coefficients
    }
    
  }
  if(mod$kappa.link=="probit"){
    for(j in 1:n){
      kappa1e[,j]<-iprobit(X[j,]%*%t(mod$mcmc.out$samples[,1:p]))
      kappa2e[,j]<-iprobit(X2[j,]%*%t(mod$mcmc.out$samples[,(p+1):(p+p2)]))
      for(i in 1:m){
        resu[i,j] <- optim(0.5, ff, method = "L-BFGS-B", kappa1e = kappa1e[i,j], 
                           kappa2e = kappa2e[i,j],gamaest = gamaest[i], thetaest = thetaest[i],
                           tauest = tauest[i], q = q,lower = c(0.0001), upper=c(0.9999))$par
      }
      
      print(j)
    }
    for(f in 1:m){
      mod=lm(probit(resu[f,])~X3[,-1])
      aux2[f,]=mod$coefficients
    }
    
  }
  if(mod$kappa.link=="cloglog"){
    for(j in 1:n){
      kappa1e[,j]<-icloglog(X[j,]%*%t(mod$mcmc.out$samples[,1:p]))
      kappa2e[,j]<-icloglog(X2[j,]%*%t(mod$mcmc.out$samples[,(p+1):(p+p2)]))
      for(i in 1:m){
        resu[i,j] <- optim(0.5, ff, method = "L-BFGS-B", kappa1e = kappa1e[i,j], 
                           kappa2e = kappa2e[i,j],gamaest = gamaest[i], thetaest = thetaest[i],
                           tauest = tauest[i], q = q,lower = c(0.0001), upper=c(0.9999))$par
      }
      
      print(j)
    }
    for(f in 1:m){
      mod=lm(cloglog(resu[f,])~X3[,-1])
      aux2[f,]=mod$coefficients
    }
    
  }
  
  colnames(aux2)<-colnames(X3)
  return(aux2)
}
                          
                        
