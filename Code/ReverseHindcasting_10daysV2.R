
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

##A version of hindcasting using MH, but assuming exact knowledge of infection times.
## This should mean that I'm only estimating the theta values and curve pars, which should thus converge to true values.

standard.dev<-log(1.1)
testdelay.mean<-3
sampling.time<-50
sample.size<-50

btv.2008.interpolated<-as.data.frame(predict(smooth.spline(btv.2008.padded$Day.since.19th.of.Sept,btv.2008.padded$Cases,df=nrow(btv.2008.df)/4),newdata=data.frame(Date.since.19th.of.Sept=1:176)))
btv.2008.interpolated[btv.2008.interpolated[,2]<0,2]<-0

#exact version
infection.times<-rep(1:50,times=round(btv.2008.interpolated[1:50,2]))
exposure.duration<-sampling.time-infection.times
igCurve.set<-igCurve(X.star=0.15,D=15,a=0.05,S=0.3)
bactLoad.set<-bactLoad(D=15)

btv.test.data<-data.frame(infection.time=infection.times,na=rlnorm(length(infection.times),
                                                                   meanlog=log(bactLoad.set(exposure.duration)),
                                                                   sdlog=standard.dev),
                          ab=rlnorm(length(infection.times),
                                    meanlog=log(igCurve.set(exposure.duration)),
                                    sdlog=standard.dev))


###This assigned an observation time after infection time, with an exponentially distributed wait.
## Should probably be discretised somehow
##Also, look up established methods for back projection.
btv.test.data$obs.time<-pmin(btv.test.data$infection.time+rexp(nrow(btv.test.data),rate=1/testdelay.mean),49)

#priors!!

#plogtheta.na
#plogmean.na






##Create a function that takes the current state of the MCMC, and calculates the loglikelihood


##Creating a function for calculating the log likelihood of paramaters, given test data, 
##epidemic curve and priors, that is log(P(params|AB,NA,P(Exposure.times),priors))

ll.reverse.hindcasting.gen<-function(ab.values,    #observed antibody test values
                                     na.values,        #observed na test values
                                     epi.curve,        #known distribution of infection times (i.e. observed epidemic curve
                                     curve.pars.prior,  ## function for calculating probability of parameters of the test curves
                                     sdlog.prior,       ##function for calculation  probability of standard deviation
                                     rate.prior,        ##function for calculating the prior prob of the rate parameter      ##is this prior needed?
                                     Infection.times.obs,
                                     Sampling.time=sampling.time
){
  
  
  ##Calculating the log likelihood P(AB|\theta,\mu)
  ll.ab<-function(Logmean.ab,sdlog.){
    dlnorm(ab.values,meanlog=Logmean.ab,sdlog=sdlog.,log=T)}
  
  ##Calculating the log likelihood P(NA|\theta,\mu)
  ll.na<-function(Logmean.na,sdlog.){
    dlnorm(na.values,meanlog=Logmean.na,sdlog.,log=T)}
  
  ##Calculating the log likelihood P(E.times|Epidemic.curve)
  ll.infectiontimes<-function(Infection.times){
    valid.infection<-Infection.times[Infection.times>0.5&Infection.times<(length(epi.curve)+0.5)]
    log(epi.curve[valid.infection]/sum(epi.curve))
  }
  
  ##Calculating a global fit, P({E.times}|Epidemic.curve)
  ll.infectiontimes.global<-function(Infection.times=infection.times.est){
    Infection.times[Infection.times<0]<-length(epi.curve)+1  ### estimated infection times before actual starting time is given outside the range value
    Infection.times[Infection.times>=length(epi.curve)+1]<-length(epi.curve)+1  ##Estimated infection times after time of sampling
    counts.per.day<-tabulate(floor(Infection.times)+1,nbins=length(epi.curve)+1) ## one extra bin for after
    ll<-dmultinom(counts.per.day,sum(counts.per.day),prob=c(epi.curve,0)/sum(epi.curve),log=T)
    return(ll)
  }
  
  ll.obs.delay<-function(delay.est,delay.rate){
    dexp(delay.est,rate=delay.rate,log=T)
  }
  
  ##This creates a function that accepts the current parameter values in the form of a vector,
  ##and then calculates the likelihood of data given the parameter values P(data|params)
  
  function(pars){
    
    ##Initializing the parameters
    curve.pars<-pars[grep("curve.pars",names(pars))]
    Ab.sdlog.scaled<-pars[grep("Ab.sdlog",names(pars))]
    Na.sdlog.scaled<-pars[grep("Na.sdlog",names(pars))]
    delay.rate.log<-pars[grep("delay.rate.log",names(pars))]
    delay.est.scaled<-pars[grep("delay.est.scaled",names(pars))]
    
    
    delay.est<-invlogit(delay.est.scaled)*Infection.times.obs
    Ab.sdlog<-exp(Ab.sdlog.scaled)
    Na.sdlog<-exp(Na.sdlog.scaled)
    delay.rate<-exp(delay.rate.log)
    ###force delay.est to always be positive!
    Infection.times.est<-Infection.times.obs-delay.est
    Exposure.duration<-Sampling.time-(Infection.times.est)
    
    
    
    igCurve.current<-igCurve(X.star=exp(curve.pars["curve.pars.X.star.log"])+1e-6,
                             S=exp(curve.pars["curve.pars.S.log"])+1e-6,
                             a=exp(curve.pars["curve.pars.a.log"])+1e-6,
                             D=exp(curve.pars["curve.pars.D.log"])+1e-6)
    bactLoad.current<-bactLoad(D=exp(curve.pars["curve.pars.D.log"])+1e-6)
    #Calculating the estimated mean, based on the Simonsen parametrization
    logmean.ab<-log(igCurve.current(Exposure.duration))
    logmean.na<-log(bactLoad.current(Exposure.duration))
    
    ##sum together the likelihood for all parts of the model at the current point
    post.ll<-sdlog.prior(Na.sdlog)+
      sdlog.prior((Ab.sdlog))+
      sum(curve.pars.prior(curve.pars))+
      #sum(delay.prior(Infection.times.obs-Infection.times.est))+
      rate.prior(delay.rate)+
      sum(ll.infectiontimes.global(Infection.times.est))+
      sum(ll.obs.delay(delay.est,delay.rate))+
      sum(ll.na(logmean.na,Na.sdlog))+
      sum(ll.ab(logmean.ab,Ab.sdlog)) 
    if(is.nan(post.ll)){browser()}
    return(post.ll)
  }
}


##Setting the priors for the curve and the lognormals.
curve.prior<-function(curve.pars){
  ll<-dlnorm(exp(curve.pars["curve.pars.X.star.log"]),meanlog=log(0.1),sdlog=log(1.4),log=T)+
    dlnorm(exp(curve.pars["curve.pars.D.log"]),meanlog=log(10),sdlog=log(1.4),log=T)+
    dlnorm(exp(curve.pars["curve.pars.a.log"]),meanlog=log(0.01),sdlog=log(1.4),log=T)+
    dlnorm(exp(curve.pars["curve.pars.S.log"]),meanlog=log(0.1),sdlog=log(1.4),log=T)
  return(ll)
}

gamma.prior<-function(sdlog,shape.=2,rate.=2){
  dgamma(sdlog,shape=shape.,rate=rate.,log=T)
}
rate.prior<-function(rate.par,inv.min=1,inv.max=20){
  dunif(1/rate.par,inv.min,inv.max,log=T)
}

##A jumping function for generating new proposal values based on current parameter values

reverse.hindcast.jump<-function(params,sd.=1){
  
  new.curve.pars<-rlnorm(length(params$curve.pars),
                         mean=log(params$curve.pars),sd=log(sd./10+1))
  names(new.curve.pars)<-names(params$curve.pars)
  new.ab.lognorm.theta<-rnorm(1,params$Ab.lognorm.theta,sd=sd.)
  new.na.lognorm.theta<-rnorm(1,params$Na.lognorm.theta,sd=sd.)
  new.na.lognorm.theta<-rnorm(1,params$Na.lognorm.theta,sd=sd.)
  
  
  return(params=c(curve.pars=new.curve.pars,
                  Ab.lognorm.theta=new.ab.lognorm.theta,
                  Na.lognorm.theta=new.na.lognorm.theta,
                  Infection.times=new.infection.times))
  
  
}



jump.mvnorm<-function(pars,jump.cov=NULL){
  if(is.null(jump.cov)|length(jump.cov)==0){
    jump.cov<-diag(length(pars))*2.4^2/length(pars)
  }
  new.pars<-rmvnorm(1,
                    mean=pars,
                    sigma=jump.cov,method="chol")
  new.pars<-c(curve.pars=new.pars[grepl("curve.pars",names(pars))],
              Ab.sdlog=new.pars[grepl("Ab.sdlog",names(pars))],
              Na.sdlog=new.pars[grepl("Na.sdlog",names(pars))],
              delay.rate.log=new.pars[grepl("delay.rate.log",names(pars))],
              delay.est.scaled=new.pars[grepl("delay.est.scaled",names(pars))]
  )  
  names(new.pars[grepl("curve.pars",names(pars))])<-paste("curve.pars.",c("X.star.log","S.log","a.log","D.log"),sep="")
  return(pars=new.pars)
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
  delay.est.scaled<-rnorm(nrow(btv.test.data),0,5)
  #btv.test.data$infection.time#sample(seq_along(btv.2008.interpolated$y),61,prob=ceiling(btv.2008.interpolated$y*10)/10,replace=T)
  return(c(curve.pars=curve.pars,
           Ab.sdlog=Ab.sdlog,
           Na.sdlog=Na.sdlog,#,
           delay.rate.log=delay.rate.log,
           delay.est.scaled=delay.est.scaled
           #Infection.times=infection.times
  ))
}


#this ideally needs to give a different starting value each time it is called.
prev.posterior.init<-function(Posterior,nchains=3){
  called<-0
  function(){
      called<<-called%%nchains+1
      
    par.means<-rowMeans(Posterior[called,,dim(Posterior)[3]-c(1000:0)],2)
    curve.pars<-par.means[1:4]
    names(curve.pars)<-c("X.star.log","S.log","a.log","D.log")
    Ab.sdlog<-par.means[5]
    Na.sdlog<-par.means[6]
    delay.rate.log<-par.means[7]
    delay.est.scaled<-par.means[-c(1:7)]
    return(c(curve.pars=curve.pars,
             Ab.sdlog=Ab.sdlog,
             Na.sdlog=Na.sdlog,
             delay.rate.log=delay.rate.log,
             delay.est.scaled=delay.est.scaled#,
             #Infection.times=infection.times
    ))
  }}

Rprof()
temp<-metropolis(jump.fun=jump.mvnorm,
                 init.fun=reverse.hindcast.init,
                 #prev.posterior.init(temp$Posterior),
                 ll.fun=ll.reverse.hindcasting.gen(ab.values=btv.test.data$ab,
                                                   na.values=btv.test.data$na,
                                                   epi.curve=btv.2008.interpolated$y,
                                                   curve.pars.prior=curve.prior,
                                                   sdlog.prior=gamma.prior,
                                                   rate.prior=rate.prior,
                                                   Infection.times.obs=btv.test.data$obs.time,
                                                   Sampling.time=sampling.time),
                 npars=#nrow(btv.test.data)+ ##one infection time per obs
                   2+##two theta
                   4+##four parameters for the simonsen curves
                   1+ ##rate parameter
                   nrow(btv.test.data), #delay estimates
                 nreps=2000,
                 jump.cov=#temp$jumping.kernel[[1]][[20]],
                   NULL,
                 Adapt=TRUE
)
Rprof(NULL)
temp2<-metropolis(jump.fun=jump.mvnorm,
                 init.fun=#reverse.hindcast.init,
                   prev.posterior.init(temp$Posterior),
                 ll.fun=ll.reverse.hindcasting.gen(ab.values=btv.test.data$ab,
                                                   na.values=btv.test.data$na,
                                                   epi.curve=btv.2008.interpolated$y,
                                                   curve.pars.prior=curve.prior,
                                                   sdlog.prior=gamma.prior,
                                                   rate.prior=rate.prior,
                                                   Infection.times.obs=btv.test.data$obs.time,
                                                   Sampling.time=sampling.time),
                 npars=#nrow(btv.test.data)+ ##one infection time per obs
                   2+##two theta
                   4+##four parameters for the simonsen curves
                   1+ ##rate parameter
                   nrow(btv.test.data), #delay estimates
                 nreps=500000,
                 jump.cov=temp$jumping.kernel[
                   which.max(
                     colMeans(tail(t(temp$Accepted),100))),
                   dim(temp$jumping.kernel)[2],,
                   ],
                 #NULL,
                 Adapt=TRUE,burn.in=500
)
save(temp,file="ObsKinetics_10Days.Rdata")
