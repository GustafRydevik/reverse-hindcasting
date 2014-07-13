
###Setting up basic parameters

data.path<-file.path("./Data")
clean.data.path<-file.path("./Clean Data")
plates.path<-file.path(data.path,"Plates")
script.path<-"./Code"

####### Loading libraries #######
invisible(lapply(dir(file.path(script.path,"utilities"),full.names=T),source))
invisible(lapply(dir(file.path(script.path,"multitest library"),full.names=T),source))
My.device<-Gen.device("png",res=200)

autolib(lattice)
autolib(ggplot2)
autolib(reshape2)
autolib(plyr)
autolib(mvtnorm)

source(file.path(script.path,"HindcastingLib/TestCurves_Simonsen2009.R"))

source(file.path(script.path,"MH files/metropolisHastingsAdaptiveV2.R"))


obs.data<-rlnorm(100,meanlog=log(2),sdlog=log(1.5))

test.gen<-function(obs.data,sdlog.prior,meanlog.prior){
  
  function(pars){
    est.sdlog<-pars["est.sdlog"]
    est.meanlog<-pars["est.meanlog"]
    obs.ll<-dlnorm(obs.data,meanlog=est.meanlog,sdlog=est.sdlog,log=T)
    post.ll<-sum(obs.ll)+sdlog.prior(est.sdlog)+meanlog.prior(est.meanlog)
    return(post.ll)
  }
}


##Setting the priors for the curve and the lognormals.
unif.prior<-function(x,min.=0,max.=100){
  dunif(x,min=min.,max=max.,log=T)
}

gauss.prior<-function(x){
  dnorm(x,mean=0,sd=100)
}

gamma.prior<-function(sdlog,shape.=2,rate.=2){
  dgamma(sdlog,shape=shape.,rate=rate.,log=T)
}
rate.prior<-function(rate.par,inv.min=1,inv.max=20){
  dunif(1/rate.par,inv.min,inv.max,log=T)
}

##A jumping function for generating new proposal values based on current parameter values




jump.mvnorm<-function(pars,jump.cov=NULL){
  if(is.null(jump.cov)|length(jump.cov)==0){
    jump.cov<-diag(length(pars))*2.4^2/length(pars)/10
  }
  new.pars<-c(rmvnorm(1,
                    mean=pars,
                    sigma=jump.cov,method="chol"))
  names(new.pars)<-c("est.sdlog","est.meanlog")
  return(pars=new.pars)
}

test.init<-function(){
  init.pars<-rgamma(2,shape=2,rate=2)
  names(init.pars)<-c("est.sdlog","est.meanlog")
  return(init.pars)
}

###Initialising the parameter values
reverse.hindcast.init.actual<-function(){
  curve.pars<-c(0.15,0.3,0.05,15) ##Or change this to something more reasonable later...
  names(curve.pars)<-c("X.star.log","S.log","a.log","D.log")
  Ab.sdlog<-log(1.1)#rgamma(1,0.2,2)
  Na.sdlog<-log(1.1)#rgamma(1,0.2,2)
  #infection.times<-btv.test.data$infection.time#sample(seq_along(btv.2008.interpolated$y),61,prob=ceiling(btv.2008.interpolated$y*10)/10,replace=T)
  return(c(curve.pars=curve.pars,
           Ab.sdlog=Ab.sdlog,
           Na.sdlog=Na.sdlog)
         #Infection.times=infection.times
  )
}


reverse.hindcast.init<-function(){
  curve.pars<-c(log(runif(3)),log(runif(1,1,30))) ##Or change this to something more reasonable later...
  names(curve.pars)<-c("X.star.log","S.log","a.log","D.log")
  Ab.sdlog<-rgamma(1,0.2,2)
  Na.sdlog<-rgamma(1,0.2,2)
  delay.rate.log<-log(1/runif(1,2,8))
  delay.est.nonscaled<-rnorm(nrow(btv.test.data),0,5)
  #btv.test.data$infection.time#sample(seq_along(btv.2008.interpolated$y),61,prob=ceiling(btv.2008.interpolated$y*10)/10,replace=T)
  return(c(curve.pars=curve.pars,
           Ab.sdlog=Ab.sdlog,
           Na.sdlog=Na.sdlog,#,
           delay.rate.log=delay.rate.log,
           delay.est.nonscaled=delay.est.nonscaled
           #Infection.times=infection.times
  ))
}

prev.posterior.init<-function(Posterior){
  function(){
    par.means<-apply(Posterior[,,dim(Posterior)[3]-c(1000:0)],2,mean)
    curve.pars<-par.means[1:4]
    names(curve.pars)<-c("X.star.log","S.log","a.log","D.log")
    Ab.sdlog<-par.means[5]
    Na.sdlog<-par.means[6]
    delay.rate.log<-par.means[7]
    delay.est.nonscaled<-par.means[-c(1:7)]
    return(c(curve.pars=curve.pars,
             Ab.sdlog=Ab.sdlog,
             Na.sdlog=Na.sdlog,
             delay.rate.log=delay.rate.log,
             delay.est.nonscaled=delay.est.nonscaled#,
             #Infection.times=infection.times
    ))
  }}

Rprof()
temp<-metropolis(jump.fun=jump.mvnorm,
                 init.fun=test.init,
                 #prev.posterior.init(temp$Posterior),
                 ll.fun=test.gen(obs.data=obs.data,
                                                   sdlog.prior=gamma.prior,
                                                   meanlog.prior=gamma.prior),
                 npars=#nrow(btv.test.data)+ ##one infection time per obs
                   2,##two theta
                 nreps=50000,
                 jump.cov=#temp$jumping.kernel[[1]][[20]],
                   NULL,
                 Adapt=TRUE
)
plot(filter(ts(t(temp$Accepted)),rep(0.01,100)))
Rprof(NULL)