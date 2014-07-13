
###Setting up basic parameters

data.path<-file.path("./Original Data")
clean.data.path<-file.path("./Clean Data")
plates.path<-file.path(data.path,"Plates")
script.path<-"./ScriptFiles"

####### Loading libraries #######
lapply(dir(file.path(script.path,"utilities"),full.names=T),source)
My.device<-Gen.device("png",res=200)

autolib(lattice)
autolib(ggplot2)
autolib(reshape2)
autolib(plyr)
autolib(mvtnorm)

source("./TestCurves_Simonsen2009.R")

##A version of hindcasting using MH, but assuming exact knowledge of infection times.
## This should mean that I'm only estimating the theta values and curve pars, which should thus converge to true values.

btv.2008.interpolated<-as.data.frame(predict(smooth.spline(btv.2008.padded$Day.since.19th.of.Sept,btv.2008.padded$Cases,df=nrow(btv.2008.df)/4),newdata=data.frame(Date.since.19th.of.Sept=1:176)))


infection.times<-rep(1:50,times=round(btv.2008.interpolated[1:50,2]))
igCurve.set<-igCurve(X.star=0.15,D=15,a=0.05,S=0.3)
bactLoad.set<-bactLoad(D=15)

btv.test.data<-data.frame(infection.time=infection.times,na=rlnorm(length(infection.times),
                                                                   meanlog=log(bactLoad.set(infection.times)),
                                                                   sdlog=log(1.1)),
                          ab=rlnorm(length(infection.times),
                                    meanlog=log(igCurve.set(infection.times)),
                                    sdlog=log(1.1)))




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
                                 sdlog.prior,      ##function for calculation lognormal probability
                                Infection.times.true
                                 ){
  
  
  ##Calculating the log likelihood P(AB|\theta,\mu)
  ll.ab<-function(Logmean.ab,sdlog.,Ab.values=ab.values){
    dlnorm(Ab.values,Logmean.ab,sdlog=sdlog.,log=T)}
  
  ##Calculating the log likelihood P(NA|\theta,\mu)
  ll.na<-function(Logmean.na,sdlog.,Na.values=na.values){
    dlnorm(Na.values,Logmean.na,sdlog.,log=T)}
  
  ##Calculating the log likelihood P(E.times|Epidemic.curve)
  ll.infectiontimes<-function(Infection.times,Epi.curve=epi.curve){
    valid.infection<-Infection.times[Infection.times>0.5&Infection.times<(length(Epi.curve)+0.5)]
    log(Epi.curve[valid.infection]/sum(Epi.curve))
  }
  
  ##This creates a function that accepts the current parameter values,
  ##and then calculates the likelihood of data given the parameter values P(data|params)
  
 function(curve.pars,
          Ab.sdlog,
          Na.sdlog#,
          #Infection.times
          ){
  
   igCurve.current<-igCurve(X.star=curve.pars["X.star"],
                            S=curve.pars["S"],
                            a=curve.pars["a"],
                            D=curve.pars["D"])
   bactLoad.current<-bactLoad(D=curve.pars["D"])
  #Calculating the estimated mean, based on the Simonsen parametrization
  logmean.ab<-log(igCurve.current(Infection.times.true))
  logmean.na<-log(bactLoad.current(Infection.times.true))

  ##sum together the likelihood for all parts of the model at the current point
  post.ll<-sdlog.prior(Na.sdlog)+
    sdlog.prior((Ab.sdlog))+
    sum(curve.pars.prior(curve.pars))+
    #sum(curve.pars.prior(Na.curve.pars))+
    #sum(ll.infectiontimes(Infection.times))+
    sum(ll.na(logmean.na,Na.sdlog))+
    sum(ll.ab(logmean.ab,Ab.sdlog))  
  return(post.ll)
 }
}
  

##Setting the priors for the curve and the lognormals.
unif.prior<-function(curve.pars,min.=-100,max.=100){
  dunif(curve.pars,min=min.,max=max.,log=T)
}

gamma.prior<-function(sdlog,shape.=2,rate.=2){
  dgamma(sdlog,shape=shape.,rate=rate.,log=T)
}

##A jumping function for generating new proposal values based on current parameter values

reverse.hindcast.jump<-function(params,sd.=1){

                                new.curve.pars<-rlnorm(length(params$curve.pars),
                                                      mean=log(params$curve.pars),sd=log(sd./10+1))
                                names(new.curve.pars)<-names(params$curve.pars)
                                new.ab.lognorm.theta<-rnorm(1,params$Ab.lognorm.theta,sd=sd.)
                                new.na.lognorm.theta<-rnorm(1,params$Na.lognorm.theta,sd=sd.)
                                new.na.lognorm.theta<-rnorm(1,params$Na.lognorm.theta,sd=sd.)
         
                            
  return(params=list(curve.pars=new.curve.pars,
                     Ab.lognorm.theta=new.ab.lognorm.theta,
                     Na.lognorm.theta=new.na.lognorm.theta,
                     Infection.times=new.infection.times))
  
  
}



jump.mvnorm<-function(params,jump.cov=NULL){
  if(is.null(jump.cov)){
    jump.cov<-diag(x=unlist(params)/1000000,nrow=length(unlist(params)),ncol=length(unlist(params)))
  }
  new.params<-rmvnorm(1,
                      mean=unlist(params),
                      sigma=jump.cov,method="chol")
  new.params.list<-list(curve.pars=new.params[1:4],
          Ab.sdlog=new.params[5],
          Na.sdlog=new.params[6]#,
         # Infection.times=new.params[-c(1:6)]
                        )  
  names(new.params.list$curve.pars)<-c("X.star","S","a","D")
  return(params=new.params.list)
}



###Initialising the parameter values
reverse.hindcast.init.actual<-function(){
  curve.pars<-c(0.15,0.3,0.05,15) ##Or change this to something more reasonable later...
  names(curve.pars)<-c("X.star","S","a","D")
  Ab.sdlog<-log(1.1)#rgamma(1,0.2,2)
  Na.sdlog<-log(1.1)#rgamma(1,0.2,2)
  #infection.times<-btv.test.data$infection.time#sample(seq_along(btv.2008.interpolated$y),61,prob=ceiling(btv.2008.interpolated$y*10)/10,replace=T)
  return(list(curve.pars=curve.pars,
              Ab.sdlog=Ab.sdlog,
              Na.sdlog=Na.sdlog#,
              #Infection.times=infection.times
              ))
}


reverse.hindcast.init<-function(){
  curve.pars<-c(runif(3),runif(1,1,30)) ##Or change this to something more reasonable later...
  names(curve.pars)<-c("X.star","S","a","D")
  Ab.sdlog<-rgamma(1,0.2,2)
  Na.sdlog<-rgamma(1,0.2,2)
  #infection.times<-btv.test.data$infection.time#sample(seq_along(btv.2008.interpolated$y),61,prob=ceiling(btv.2008.interpolated$y*10)/10,replace=T)
  return(list(curve.pars=curve.pars,
              Ab.sdlog=Ab.sdlog,
              Na.sdlog=Na.sdlog#,
              #Infection.times=infection.times
  ))
}

prev.posterior.init<-function(Posterior){
  function(){
par.means<-apply(Posterior[,,dim(Posterior)[3]-c(1000:0)],2,mean)
curve.pars<-par.means[1:4]
names(curve.pars)<-c("X.star","S","a","D")
Ab.sdlog<-par.means[5]
Na.sdlog<-par.means[6]
return(list(curve.pars=curve.pars,
            Ab.sdlog=Ab.sdlog,
            Na.sdlog=Na.sdlog#,
            #Infection.times=infection.times
))
}}

temp<-metropolis(jump.fun=jump.mvnorm,
           init.fun=#reverse.hindcast.init,
           prev.posterior.init(temp$Posterior),
           ll.fun=ll.reverse.hindcasting.gen(ab.values=btv.test.data$ab,
                                             na.values=btv.test.data$na,
                                             epi.curve=btv.2008.interpolated$y,
                                             curve.pars.prior=unif.prior,
                                             sdlog.prior=gamma.prior,
                                             Infection.times.true=btv.test.data$infection.time),
           npars=#nrow(btv.test.data)+ ##one infection time per obs
             2+##two theta
             4,##four parameters for the simonsen curves
           nreps=20000,
           jump.cov=temp$jumping.kernel[[3]][[20]],
           Adapt=FALSE
           )



#old stuff
ll.norm<-function(muhat,sigmahat,data=sample.data){
  sum(dnorm(sample.data,muhat,sigmahat,log=T))
}
pmu<-function(mu){dnorm(mu,0,100,log=T)}
psigma<-function(sigma){dunif(sigma,0,10,log=T)}

post.norm<-function(mu,sigma,data=sample.data){
  ll.norm(mu,sigma,data)+pmu(mu)+psigma(sigma)
}

geninits<-function(){
  list(mu=runif(1,0,20),
       sigma=runif(1,0,10))
}

jump<-function(x,dist=.1,prior=NULL){
  x+rnorm(length(x),0,dist)
  prior()
}