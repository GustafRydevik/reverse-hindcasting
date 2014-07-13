
toy.ll.fun.gen<-function(observed.data){
  
  llnorm<-function(mean.est,sd.est){
    dnorm(observed.data,mean.est,sd.est,log=T)
  }
  
  gamma.prior<-function(x,par1,par2){
    dgamma(1/x,par1,par2,log=T)
  }
  lognorm.prior<-function(x,meanlog,sdlog){
    dlnorm(x,meanlog,sdlog,log=T)
  }
  
ll.fun<-function(pars){
  mean.est<-exp(pars["mean.est.scaled"])
  sd.est<-exp(pars["sd.est.scaled"])
  
  posterior.ll<-sum(llnorm(mean.est,sd.est))+
    gamma.prior(sd.est,0.001,0.001)
  +lognorm.prior(mean.est,meanlog=log(10),sdlog=log(sqrt(2)))
  return(posterior.ll)
}
  return(ll.fun)
}

toy.data<-rnorm(1000,10,0.1)


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






init.fun.toy<-function(){
  mean.est.scaled<-log(rnorm(1,10,2))
  sd.est.scaled<-log(rlnorm(1,log(0.1),log(sqrt(10))))
  return(c(mean.est.scaled=mean.est.scaled,sd.est.scaled=sd.est.scaled))
}


prev.posterior.init<-function(Posterior){
  called<-0
  nchains<-dim(Posterior)[1]
  function(){
    called<<-called%%nchains+1
    
    par.means<-rowMeans(Posterior[called,,dim(Posterior)[3]-c(1000:0)],2)
    names(par.means)<-colnames(Posterior)
    return(par.means)
    
  }}




temp<-metropolis(jump.fun=jump.mvnorm,
                 init.fun=init.fun.toy,
                 #prev.posterior.init(temp$Posterior),
                 ll.fun=toy.ll.fun.gen(toy.data),
                 npars=#nrow(btv.test.data)+ ##one infection time per obs
                   #2+##two theta
                   2,###four parameters for the simonsen curves
                 #1+ ##rate parameter
                 #nrow(btv.test.data), #delay estimates
                 nreps=10000,
                 jump.cov=#temp$jumping.kernel[[1]][[20]],
                   NULL,
                 Adapt=TRUE
)


temp2<-metropolis(jump.fun=jump.mvnorm,
                 init.fun=prev.posterior.init(temp$Posterior),
                 #prev.posterior.init(temp$Posterior),
                 ll.fun=toy.ll.fun.gen(toy.data),
                 npars=#nrow(btv.test.data)+ ##one infection time per obs
                   #2+##two theta
                   2,###four parameters for the simonsen curves
                 #1+ ##rate parameter
                 #nrow(btv.test.data), #delay estimates
                 nreps=10000,
                 jump.cov=with(temp,
                               jumping.kernel[which.max(rowMeans(Accepted)),
                                              dim(jumping.kernel)[2]
                                              ,,]
                               ),
                   #NULL,
                 Adapt=TRUE
)



temp2<-metropolis(jump.fun=jump.mvnorm,
                  init.fun=prev.posterior.init(temp2$Posterior),
                  #prev.posterior.init(temp$Posterior),
                  ll.fun=toy.ll.fun.gen(toy.data),
                  npars=#nrow(btv.test.data)+ ##one infection time per obs
                    #2+##two theta
                    2,###four parameters for the simonsen curves
                  #1+ ##rate parameter
                  #nrow(btv.test.data), #delay estimates
                  nreps=20000,
                  jump.cov=with(temp2,
                                jumping.kernel[which.max(rowMeans(Accepted)),
                                               dim(jumping.kernel)[2]
                                               ,,]
                  ),
                  #NULL,
                  Adapt=TRUE
)
