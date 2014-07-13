
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
##Currently only estimating variance and two of the four curve parameters.
##
standard.dev<-log(1.1)
testdelay.mean<-14
sampling.time<-50
sample.size<-200
incubation.mean<-7
incubation.spread<-sqrt(2) ##Meaning 96% are within a factor two of the mean
btv.2008.interpolated<-as.data.frame(predict(smooth.spline(btv.2008.padded$Day.since.19th.of.Sept,btv.2008.padded$Cases,df=nrow(btv.2008.df)/4),newdata=data.frame(Date.since.19th.of.Sept=1:176)))
btv.2008.interpolated[btv.2008.interpolated[,2]<0,2]<-0

#observe.
obs.t1<-sample(1:sampling.time,sample.size,prob=round(btv.2008.interpolated[1:sampling.time,"y"]),replace=T)
obs.t2<-sampling.time
infection.times<-pmin(obs.t1-rlnorm(length(obs.t1),meanlog=log(incubation.mean),sdlog=log(incubation.spread)),sampling.time)
exposure.duration.t1<-obs.t1-infection.times
exposure.duration.t2<-obs.t2-infection.times
igCurve.set<-igCurve(X.star=0.15,D=15,a=0.05,S=0.3)
bactLoad.set<-bactLoad(D=15)

btv.test.data<-data.frame(infection.time=infection.times,na.t1=rlnorm(length(infection.times),
                                                                   meanlog=log(bactLoad.set(exposure.duration.t1)),
                                                                   sdlog=standard.dev),
                          ab.t1=rlnorm(length(infection.times),
                                    meanlog=log(igCurve.set(exposure.duration.t1)),
                                    sdlog=standard.dev),na.t2=rlnorm(length(infection.times),
                                                                        meanlog=log(bactLoad.set(exposure.duration.t2)),
                                                                        sdlog=standard.dev),
                          ab.t2=rlnorm(length(infection.times),
                                       meanlog=log(igCurve.set(exposure.duration.t2)),
                                       sdlog=standard.dev),
                          obs.t1=obs.t1,
                          obs.t2=obs.t2)


###This assigned an observation time after infection time, with an exponentially distributed wait.
## Should probably be discretised somehow
##Also, look up established methods for back projection.

#priors!!

#plogtheta.na
#plogmean.na






##Create a function that takes the current state of the MCMC, and calculates the loglikelihood


##Creating a function for calculating the log likelihood of paramaters, given test data, 
##epidemic curve and priors, that is log(P(params|AB,NA,P(Exposure.times),priors))



ll.reverse.hindcasting.gen<-function(ab.values.t1,      #observed antibody test values
                                     na.values.t1,     #observed na test values
                                     ab.values.t2,      #observed antibody test values
                                     na.values.t2,      #observed na test values
                                     epi.curve,         #known distribution of infection times (i.e. observed epidemic curve
                                     curve.pars.prior,  ## function for calculating probability of parameters of the test curves
                                     sdlog.prior,       ##function for calculation  probability of standard deviation
                                     rate.prior,        ##function for calculating the prior prob of the rate parameter      ##is this prior needed?
                                     Sampling.t1=obs.t1,
                                     Sampling.t2=obs.t2,
                                     infection.time.true=infection.times
){
  
  
  ##Calculating the log likelihood P(AB|\theta,\mu)
  ll.ab<-function(Logmean.ab,sdlog.,Ab.values=ab.values.t1){
    dlnorm(Ab.values,Logmean.ab,sdlog=sdlog.,log=T)}
  
  ##Calculating the log likelihood P(NA|\theta,\mu)
  ll.na<-function(Logmean.na,sdlog.,Na.values=na.values.t1){
    dlnorm(Na.values,Logmean.na,sdlog.,log=T)}
  

  ##Calculating the log likelihood P(E.times|Epidemic.curve)
  ll.infectiontimes<-function(Infection.times,Epi.curve=epi.curve){
    valid.infection<-Infection.times[Infection.times>0.5&Infection.times<(length(Epi.curve)+0.5)]
    log(Epi.curve[valid.infection]/sum(Epi.curve))
  }
  
  ##Calculating a global fit, P({E.times}|Epidemic.curve)
  ll.infectiontimes.global<-function(Infection.times=infection.times.est,Epi.curve=epi.curve){
    valid.infection<-Infection.times[Infection.times>0.5&Infection.times<(length(Epi.curve)+0.5)]
    counts.per.day<-data.frame(table(factor(round(Infection.times),levels=1:length(Epi.curve))))
    ll<-dmultinom(counts.per.day[,2],sum(counts.per.day[,2]),prob=Epi.curve/sum(Epi.curve),log=T)
    return(ll)
  }
  
  ll.obs.delay<-function(delay.est,delay.mean,delay.spread){
    dlnorm(delay.est,meanlog=log(delay.mean),sdlog=log(delay.spread),log=T)
  }
  
  ##This creates a function that accepts the current parameter values,
  ##and then calculates the likelihood of data given the parameter values P(data|params)
  

    
    function(pars){
      
      ##Initializing the parameters
      curve.pars<-exp(pars[grep("curve.pars",names(pars))])
      Ab.sdlog.scaled<-pars[grep("Ab.sdlog",names(pars))]
      Na.sdlog.scaled<-pars[grep("Na.sdlog",names(pars))]
      #delay.rate.log<-pars[grep("delay.rate.log",names(pars))]
      #delay.t1.scaled<-pars[grep("delay.t1.scaled",names(pars))]
      
      
      #delay.t1<-exp(delay.t1.scaled)
      Ab.sdlog<-exp(Ab.sdlog.scaled)
      Na.sdlog<-exp(Na.sdlog.scaled)
      delay.mean<-incubation.mean#exp(delay.rate.log)
      #delay.spread<-sqrt(2)
      ###force delay.est to always be positive!
      #Infection.times.est<-Sampling.t1-delay.t1
      delay.t2<-Sampling.t2-(infection.time.true)
      
      igCurve.current<-igCurve(X.star=(curve.pars["curve.pars.X.star.log"]),
                               S=(curve.pars["curve.pars.S.log"]),
                               a=(curve.pars["curve.pars.a.log"]),
                               D=15)#(curve.pars["curve.pars.D.log"]))
      bactLoad.current<-bactLoad(D=15)#(curve.pars["curve.pars.D.log"]))
    #Calculating the estimated mean, based on the Simonsen parametrization
    
    ##This is if we are estimating the infection times
#     logmean.ab.t1<-log(igCurve.current( delay.t1))
#     logmean.na.t1<-log(bactLoad.current( delay.t1))
#     logmean.ab.t2<-log(igCurve.current( delay.t2))
#     logmean.na.t2<-log(bactLoad.current( delay.t2))

#This is assuming exact knowledge about infection times.
    logmean.ab.t1<-log(igCurve.current( Sampling.t1-infection.time.true))
    logmean.na.t1<-log(bactLoad.current( Sampling.t1-infection.time.true))
    logmean.ab.t2<-log(igCurve.current( Sampling.t2-infection.time.true))
    logmean.na.t2<-log(bactLoad.current( Sampling.t2-infection.time.true))
    
    ##sum together the likelihood for all parts of the model at the current point
    post.ll<-sdlog.prior(Na.sdlog)+
      sdlog.prior((Ab.sdlog))+
      sum(curve.pars.prior(curve.pars))+
      #sum(delay.prior(Infection.times.obs-Infection.times.est))+
      #rate.prior(delay.rate)+
      #sum(ll.infectiontimes.global(Infection.times.est))+
      #sum(ll.obs.delay(delay.t1,delay.mean,delay.spread))+
      sum(ll.na(logmean.na.t1,Na.sdlog,na.values.t1))+
      sum(ll.ab(logmean.ab.t1,Ab.sdlog,ab.values.t1))+
      sum(ll.na(logmean.na.t2,Na.sdlog,na.values.t2))+
      sum(ll.ab(logmean.ab.t2,Ab.sdlog,ab.values.t2))  
    return(post.ll)
  }
}


##Setting the priors for the curve and the lognormals.
unif.prior<-function(curve.pars,min.=0,max.=100){
  dunif(exp(curve.pars),min=min.,max=max.,log=T)
}


curve.gamma.prior<-function(curve.pars,shape.=0.1,rate.=0.1)
  dgamma(curve.pars,shape=shape.,rate=rate.,log=T)


curve.lognorm.prior<-function(curve.pars,scale=c(0.1,0.1,0.1,10)){
  dlnorm(curve.pars[1],meanlog=log(scale[1]),sdlog=log(10)/2,log=T)+
  dlnorm(curve.pars[2],meanlog=log(scale[2]),sdlog=log(10)/2,log=T)+
  dlnorm(curve.pars[3],meanlog=log(scale[3]),sdlog=log(10)/2,log=T)#+
  #dlnorm(curve.pars[4],meanlog=log(scale[4]),sdlog=log(10)/2,log=T)
}

gamma.prior<-function(sdlog,shape.=2,rate.=4){
  dgamma(exp(sdlog)-1,shape=shape.,rate=rate.,log=T)
}
rate.prior<-function(rate.par,inv.min=1,inv.max=20){
  dunif(1/rate.par,inv.min,inv.max,log=T)
}

##A jumping function for generating new proposal values based on current parameter values




jump.mvnorm<-function(pars,jump.cov=NULL,scaling=1){
  if(is.null(jump.cov)|length(jump.cov)==0){
    jump.cov<-diag(length(pars))*2.4^2/length(pars)/scaling
  }
  new.pars<-rmvnorm(1,
                    mean=pars,
                    sigma=jump.cov,method="chol")
  names(new.pars)<-names(pars)
    return(pars=c(new.pars))
}



###Initialising the parameter values
reverse.hindcast.init.actual<-function(){
  curve.pars<-c(0.15,0.3,0.05,15) ##Or change this to something more reasonable later...
  names(curve.pars)<-c("X.star.log","S.log","a.log","D.log")
  curve.pars<-log(curve.pars)
  Ab.sdlog<-log(1.1)#rgamma(1,0.2,2)
  Na.sdlog<-log(1.1)#rgamma(1,0.2,2)
  delay.t1.scaled<-log(obs.t1-infection.times)
  #infection.times<-btv.test.data$infection.time#sample(seq_along(btv.2008.interpolated$y),61,prob=ceiling(btv.2008.interpolated$y*10)/10,replace=T)
  return(c(curve.pars=curve.pars,
          # Ab.sdlog=Ab.sdlog,
           #Na.sdlog=Na.sdlog)
          delay.t1.scaled=delay.t1.scaled
         #Infection.times=infection.times
  ))
}


reverse.hindcast.init<-function(){
  curve.pars<-c(log(runif(3)))#,log(runif(1,1,30))) ##Or change this to something more reasonable later...
  names(curve.pars)<-c("X.star.log","S.log","a.log")#,"D.log")
  Ab.sdlog<-rgamma(1,0.2,2)
  Na.sdlog<-rgamma(1,0.2,2)
  delay.rate.log<-log(1/runif(1,2,8))
  delay.t1.scaled<-rnorm(nrow(btv.test.data),0,5)
  #btv.test.data$infection.time#sample(seq_along(btv.2008.interpolated$y),61,prob=ceiling(btv.2008.interpolated$y*10)/10,replace=T)
  return(c(curve.pars=curve.pars,
           Ab.sdlog=Ab.sdlog,
           Na.sdlog=Na.sdlog#,
           #delay.rate.log=delay.rate.log,
           #delay.t1.scaled=delay.t1.scaled
           #Infection.times=infection.times
  ))
}


#this ideally needs to give a different starting value each time it is called.
prev.posterior.init<-function(Posterior,nchains=3){
  called<-0
  function(){
    called<<-called%%nchains+1
    
    par.means<-rowMeans(Posterior[called,,dim(Posterior)[3]-c(1000:0)],2)
    curve.pars<-par.means[1:3]
    names(curve.pars)<-c("X.star.log","S.log","a.log")#,"D.log")
    Ab.sdlog<-par.means[4]
    Na.sdlog<-par.means[5]
    #delay.rate.log<-par.means[5]
    #delay.t1.scaled<-par.means[-c(1:7)]
    return(c(curve.pars=curve.pars,
             Ab.sdlog=Ab.sdlog,
             Na.sdlog=Na.sdlog#,
             #delay.rate.log=delay.rate.log,
             #delay.t1.scaled=delay.t1.scaled#,
             #Infection.times=infection.times
    ))
  }}
##Setting the priors for the curve and the lognormals.





temp<-metropolis(jump.fun=jump.mvnorm,
                 init.fun=reverse.hindcast.init,#reverse.hindcast.init,
                   #prev.posterior.init(temp$Posterior),
                 ll.fun=ll.reverse.hindcasting.gen(ab.values.t1=btv.test.data$ab.t1,
                                                   na.values.t1=btv.test.data$na.t1,
                                                   ab.values.t2=btv.test.data$ab.t2,
                                                   na.values.t2=btv.test.data$na.t2,
                                                   epi.curve=btv.2008.interpolated$y,
                                                   curve.pars.prior=curve.lognorm.prior,
                                                   sdlog.prior=gamma.prior,
                                                   rate.prior=rate.prior,
                                                   Sampling.t1=btv.test.data$obs.t1,
                                                   Sampling.t2=btv.test.data$obs.t2),
                 npars=#nrow(btv.test.data)+ ##one infection time per obs
                   #2+##two theta
                   5##four parameters for the simonsen curves
                   #1+ ##rate parameter
                   ,#nrow(btv.test.data), #delay estimates
                 nreps=2000,
                 jump.cov=#temp$jumping.kernel[[1]][[20]],
                 NULL,
                 Adapt=TRUE
)

temp2<-metropolis(jump.fun=jump.mvnorm,
                  init.fun=#reverse.hindcast.init,
                    prev.posterior.init(temp$Posterior),
                  ll.fun=ll.reverse.hindcasting.gen(ab.values.t1=btv.test.data$ab.t1,
                                                    na.values.t1=btv.test.data$na.t1,
                                                    ab.values.t2=btv.test.data$ab.t2,
                                                    na.values.t2=btv.test.data$na.t2,
                                                    epi.curve=btv.2008.interpolated$y,
                                                    curve.pars.prior=curve.lognorm.prior,
                                                    sdlog.prior=gamma.prior,
                                                    rate.prior=rate.prior,
                                                    Sampling.t1=btv.test.data$obs.t1,
                                                    Sampling.t2=btv.test.data$obs.t2),
                  npars=#nrow(btv.test.data)+ ##one infection time per obs
                    #2+##two theta
                    5#+##four parameters for the simonsen curves
                    #1+ ##rate parameter
                    ,#nrow(btv.test.data), #delay estimates
                  nreps=10000,
                  jump.cov=
                     temp$jumping.kernel[
                     which.max(
                       colMeans(tail(t(temp$Accepted),100))),
                     dim(temp$jumping.kernel)[2],,
                     ],
#                  NULL,
                  Adapt=TRUE,burn.in=500
)

temp2<-metropolis(jump.fun=jump.mvnorm,
                 init.fun=#reverse.hindcast.init,
                   prev.posterior.init(temp2$Posterior),
                 ll.fun=ll.reverse.hindcasting.gen(ab.values.t1=btv.test.data$ab.t1,
                                                   na.values.t1=btv.test.data$na.t1,
                                                   ab.values.t2=btv.test.data$ab.t2,
                                                   na.values.t2=btv.test.data$na.t2,
                                                   epi.curve=btv.2008.interpolated$y,
                                                   curve.pars.prior=curve.lognorm.prior,
                                                   sdlog.prior=gamma.prior,
                                                   rate.prior=rate.prior,
                                                   Sampling.t1=btv.test.data$obs.t1,
                                                   Sampling.t2=btv.test.data$obs.t2),
                 npars=#nrow(btv.test.data)+ ##one infection time per obs
                  # 2+##two theta
                   4#+##four parameters for the simonsen curves
                   #1+ ##rate parameter
                   ,#nrow(btv.test.data), #delay estimates
                 nreps=20000,
                 jump.cov=
                   #temp2$jumping.kernel[
                  #   which.max(
                  #     colMeans(tail(t(temp2$Accepted),100))),
                   #  dim(temp2$jumping.kernel)[2],,
                    # ],   
                 NULL,
                 Adapt=TRUE
)

save(temp,file="ObsKinetics_3Days.Rdata")
