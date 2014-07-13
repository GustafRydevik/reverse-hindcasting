
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

##A version of hindcasting using MH, but assuming exact knowledge of infection times.
## This should mean that I'm only estimating the theta values, which should thus converge to true values.

btv.2008.interpolated<-as.data.frame(predict(smooth.spline(btv.2008.padded$Day.since.19th.of.Sept,btv.2008.padded$Cases,df=nrow(btv.2008.df)/4),newdata=data.frame(Date.since.19th.of.Sept=1:176)))


infection.times<-rep(1:50,times=round(btv.2008.interpolated[1:50,2]))

btv.test.data<-data.frame(infection.time=infection.times,na=rlnorm(length(infection.times),
                                                                   meanlog=log(btv.spline.df[match(infection.times,btv.spline.df[,"DPI"]),"PCR"]),
                                                                   sdlog=log(1.1)),
                          ab=rlnorm(length(infection.times),
                                    meanlog=log(btv.spline.df[match(infection.times,btv.spline.df[,"DPI"]),"AB"]),
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
                                 theta.lognorm.prior,      ##function for calculation lognormal probability
                                Infection.times.true
                                 ){
  
  
  ##Calculating the log likelihood P(AB|\theta,\mu)
  ll.ab<-function(Logmean.ab,Theta.ab,Ab.values=ab.values){
    dlnorm(Ab.values,Logmean.ab,1/Theta.ab,log=T)}
  
  ##Calculating the log likelihood P(NA|\theta,\mu)
  ll.na<-function(Logmean.na,Theta.na,Na.values=na.values){
    dlnorm(Na.values,Logmean.na,1/Theta.na,log=T)}
  
  ##Calculating the log likelihood P(E.times|Epidemic.curve)
  ll.infectiontimes<-function(Infection.times,Epi.curve=epi.curve){
    Infection.times[Infection.times>length(Epi.curve)]<-length(Epi.curve)
    Infection.times[Infection.times<1]<-1
    log(Epi.curve[round(Infection.times)])
  }
  
  ##This creates a function that accepts the current parameter values,
  ##and then calculates the likelihood of data given the parameter values P(data|params)
  
 function(Ab.curve.pars,
          Na.curve.pars,
          Ab.lognorm.theta,
          Na.lognorm.theta,
          Infection.times
          ){
  
   
  #Calculating the estimated mean, based on a loess algorithm
   loess.ab<-loess(ab.values~Infection.times.true,weights=1/ab.values^2)  ##using actual infection times here
   loess.na<-loess(na.values~Infection.times.true,weights=1/ab.values^2)  ##using actual infection times here
   logmean.ab<-log(predict(loess.ab))
   logmean.na<-log(predict(loess.na))
   
#    lm.ab<-lm(ab.values~Infection.times+Infection.times^2+Infection.times^3+Infection.times^4,
#       weights=1/ab.values^2)
#    lm.na<-lm(na.values~Infection.times+Infection.times^2+Infection.times^3+Infection.times^4,
#              weights=1/na.values^2)
#    logmean.ab<-log(predict(lm.ab))
#    logmean.na<-log(predict(lm.na))
#    
#   logmean.ab<-log(Ab.curve.pars[1]+
#                     Ab.curve.pars[2]*Infection.times+
#                     Ab.curve.pars[3]*(Infection.times^2)+
#                     Ab.curve.pars[4]*(Infection.times^3)+
#                     Ab.curve.pars[5]*(Infection.times^4)
#   )
                  
#   logmean.na<-log(Na.curve.pars[1]+
#                     Na.curve.pars[2]*Infection.times+
#                     Na.curve.pars[3]*(Infection.times^2)+
#                     Na.curve.pars[4]*(Infection.times^3)+
#                     Na.curve.pars[5]*(Infection.times^4)
#   )
  
  
  ##sum together the likelihood for all parts of the model at the current point
  post.ll<-theta.lognorm.prior(Na.lognorm.theta)+
    theta.lognorm.prior((Ab.lognorm.theta))+
    #sum(curve.pars.prior(Ab.curve.pars))+
    #sum(curve.pars.prior(Na.curve.pars))+
    #sum(ll.infectiontimes(Infection.times))+
    sum(ll.na(logmean.na,Na.lognorm.theta))+
    sum(ll.ab(logmean.ab,Ab.lognorm.theta))
  
  
  return(post.ll)
 }
}
  

##Setting the priors for the curve and the lognormals.
unif.prior<-function(curve.pars,min.=-100,max.=100){
  dunif(curve.pars,min=min.,max=max.,log=T)
}

gamma.prior<-function(theta,shape.=2,rate.=2){
  dgamma(theta,shape=shape.,rate=rate.,log=T)
}

##A jumping function for generating new proposal values based on current parameter values
reverse.hindcast.jump<-function(params,jump.size=2){

                                new.ab.pars<-rnorm(length(params$Ab.curve.pars),
                                                      mean=params$Ab.curve.pars,jump.size)
                                new.na.pars<-rnorm(length(params$Na.curve.pars),
                                                   mean=params$Na.curve.pars,jump.size)
                                new.ab.lognorm.theta<-rnorm(1,params$Ab.lognorm.theta,jump.size)
                                new.na.lognorm.theta<-rnorm(1,params$Na.lognorm.theta,jump.size)
                                new.infection.times<-rnorm(length(infection.times),params$Infection.times,jump.size)
                                
                            
  return(params=list(Ab.curve.pars=new.ab.pars,
                     Na.curve.pars=new.na.pars,
                     Ab.lognorm.theta=new.ab.lognorm.theta,
                     Na.lognorm.theta=new.na.lognorm.theta,
                     Infection.times=new.infection.times))
  
  
}


jump.mvnorm<-function(params,jump.cov=diag(x=2,nrow=length(params),ncol=length(params))){
  new.params<-rmvnorm(n=length(unlist(params)),
                   mean=unlist(params),
                   sigma=jump.cov)
  list()
  new.ab.pars<-rnorm(length(params$Ab.curve.pars),
                     mean=params$Ab.curve.pars,jump.size)
  new.na.pars<-rnorm(length(params$Na.curve.pars),
                     mean=params$Na.curve.pars,jump.size)
  new.ab.lognorm.theta<-rnorm(1,params$Ab.lognorm.theta,jump.size)
  new.na.lognorm.theta<-rnorm(1,params$Na.lognorm.theta,jump.size)
  new.infection.times<-rnorm(length(infection.times),params$Infection.times,jump.size)
  
  
  return(params=list(Ab.curve.pars=new.ab.pars,
                     Na.curve.pars=new.na.pars,
                     Ab.lognorm.theta=new.ab.lognorm.theta,
                     Na.lognorm.theta=new.na.lognorm.theta,
                     Infection.times=new.infection.times))
  
  
}




###Initialising the parameter values
reverse.hindcast.init<-function(){
  ab.curve.pars<-runif(5,0,20)
  na.curve.pars<-runif(5,0,20)
  ab.lognorm.theta<-rgamma(1,2,2)
  na.lognorm.theta<-rgamma(1,2,2)
  infection.times<-sample(seq_along(btv.2008.interpolated$y),61,prob=ceiling(btv.2008.interpolated$y*10)/10,replace=T)
  return(list(Ab.curve.pars=ab.curve.pars,
              Na.curve.pars=na.curve.pars,
              Ab.lognorm.theta=ab.lognorm.theta,
              Na.lognorm.theta=na.lognorm.theta,
              Infection.times=infection.times))
}

###Testing the functions
ll.reverse.hindcasting<-ll.reverse.hindcasting.gen(ab.values=btv.test.data$ab,
                           na.values=btv.test.data$na,
                           epi.curve=btv.2008.interpolated$y,
                           curve.pars.prior=unif.prior,
                           theta.lognorm.prior=gamma.prior,
                           Infection.times.true=btv.test.data$infection.time)
                           
ll.reverse.hindcasting(Ab.curve.pars=c(1,2,3,4,5),
                       Na.curve.pars=c(1,2,3,4,5),
                       Ab.lognorm.theta=1/log(1.3),
                       Na.lognorm.theta=1/log(1.1),
                       Infection.times=btv.test.data$infection.time
                       )

reverse.hindcast.jump(list(Ab.curve.pars=c(1,2,3,4,5),
      Na.curve.pars=c(1,2,3,4,5),
      Ab.lognorm.theta=1/log(1.1),
      Na.lognorm.theta=1/log(1.1),
      Infection.times=sample(seq_along(btv.2008.interpolated$y),100,prob=ceiling(btv.2008.interpolated$y*10)/10,replace=T))
)


temp<-metropolis(jump.fun=reverse.hindcast.jump,
           init.fun=reverse.hindcast.init,
           ll.fun=ll.reverse.hindcasting.gen(ab.values=btv.test.data$ab,
                                             na.values=btv.test.data$na,
                                             epi.curve=btv.2008.interpolated$y,
                                             curve.pars.prior=unif.prior,
                                             theta.lognorm.prior=gamma.prior,
                                             Infection.times.true=btv.test.data$infection.time),
           npars=nrow(btv.test.data)+ ##one infection time per obs
             2+##two theta
             10,##two four-degree curves
           nreps=10000
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