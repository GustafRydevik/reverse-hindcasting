
##based on http://mbjoseph.github.io/blog/2013/09/08/metropolis/
autolib(LaplacesDemon)
metropolis<-function(
  jump.fun,
  init.fun,
  ll.fun=function(x)dnorm(x,log=TRUE),
  nreps=10000,
  nchains=3,
  npars=2,
  jump.cov=NULL,
  Adapt=TRUE,
  epsilon=1e-6,
  burn.in=1000,
  cov.samplefreq=500){
  require(parallel)
  if(nchains>detectCores()){simpleError("More chains than available cores")}
  ################################################################
  ####This is a stand- alone funtion for calculating the MCMC ####
  ################################################################
  par.values<-list(jump.fun,
       init.fun,
       nreps,
       nchains,
       npars,
       jump.cov,
       Adapt,
       epsiolon,
       burn.in,
       cov.samplefreq)
  par.values.list<-lapply(seq_len(nchains),function(x){par.values})
  metropolis.fun<-function(par.values){
    attach(par.values)
    
    posterior<-array(dim=c(npars,nreps))
    accepted<-vector(length=c(nreps-1))
    likelihood.ratio<-vector(length=c(nreps-1))
    jump.state.array<-array(dim=c(nreps%/%cov.samplefreq,npars,npars))
    dimnames(jump.state.array)<-list(
                                     iteration=1:(nreps%/%cov.samplefreq),
                                     var1=names(init.fun()),
                                     var2=names(init.fun()))
    adaptive.tuning<-(2.4^2/npars)
    theta.post<-array(dim=c(npars,nreps)) 
    for(i in 1:nreps){
    if(!(i%%round(nreps/100))){cat("|")} 
    if(i==1){
      x.scale<-init.fun()
      rownames(theta.post)<-names(x.scale)
      theta.post[,1]<-x.scale
      x<-rep(1,length(x.scale))
      names(x)<-names(x.scale)
      prev.loglik<-ll.fun(x.scale)
    } 
    x.star<-jump.fun(x,current.jump.cov)
    theta.star<-x.scale*x.star
    star.loglik<-ll.fun(theta.star)
    
    log.lr<-star.loglik-prev.loglik
    lr<-exp(log.lr)
    #if(i==nreps){browser()}
    if(is.finite(lr)){
      accept<-(runif(1)<lr)}else{
        if(!is.na(lr)){
          accept<-(lr==Inf)}else{
            accept<-0}
      }
    accepted[i-1]<-accept
    likelihood.ratio[i-1]<-log.lr
    if(accept==1){
      x<-x.star
      theta.post[,i]<-theta.star
      prev.loglik<-star.loglik
    }else{
      theta.post[,i]<- theta.post[,i-1+(i==1)]}
    
    ###Adaptive MCMC here]
    if((i==burn.in)&Adapt){
      Xhat.now<-rowMeans(theta.post[,1:i,drop=FALSE])/x.scale
    }
    if((i>burn.in)&Adapt){  
      latest.ndx<-(i-burn.in):(i-1)
      
      if(is.null(current.jump.cov)|length(current.jump.cov)==0){
        new.jump.cov<-(cov(t(theta.post[,1:i,drop=FALSE]/x.scale))+epsilon*diag(npars))*adaptive.tuning
        current.jump.cov<-new.jump.cov
        if(i%%cov.samplefreq==0){
          jump.state.array[k,i%/%cov.samplefreq,,]<-current.jump.cov}
        next
      }
      
      X.now<-theta.post[,i,drop=F]/x.scale
      X.last<-theta.post[,i-1,drop=F]/x.scale
      Xhat.last<-Xhat.now
      Xhat.now<-(Xhat.last*(i-1)+X.now)/i ### weighted mean of old mean and another data point 
      
      cov.update<-(i-1)*Xhat.last%*%t(Xhat.last)-i*Xhat.now%*%t(Xhat.now)+X.now%*%t(X.now)
      #cov.update<-(X.now-Xhat.last)%*%t(X.now-Xhat.last)
      
      new.jump.cov<-(current.jump.cov*(i-2)/(i-1)+
                       adaptive.tuning/(i-1)*(cov.update+epsilon*diag(npars)))
      
      if(is.positive.semidefinite(new.jump.cov)){current.jump.cov<-new.jump.cov}
      if(i%%cov.samplefreq==0){jump.state.array[i%/%cov.samplefreq,,]<-current.jump.cov}
    }
      
    }
    return(list(Posterior=theta.post,
                Accepted=accepted,
                Likelihood.ratio=likelihood.ratio,
                jumping.kernel=jump.state.array))
  }
    
    #######################################
    ######## Running the calc. below  #####
    #######################################
      set.seed(123, "L'Ecuyer")  
      results<-mclapply(par.values.list,metropolis.fun)
     return(results)
     }



plot.mcmc<-function(mcmc.out)
{
  op=par(mfrow=c(2,2))
  plot(ts(mcmc.out),col=2)
  hist(mcmc.out,30,col=3)
  qqnorm(mcmc.out,col=4)
  abline(0,1,col=2)
  acf(mcmc.out,col=2,lag.max=100)
  par(op)
}
