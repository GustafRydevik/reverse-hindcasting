
##based on http://mbjoseph.github.io/blog/2013/09/08/metropolis/
autolib(LaplacesDemon)
metropolis<-function(
  jump.fun,
  init.fun,
  LikelihoodFun=function(x)dnorm(x,log=TRUE),
  nreps=10000,
  nchains=3,
  npars=NULL,
  jump.cov=NULL,
  Adapt=TRUE,
  epsilon=1e-6,
  burn.in=1000,
  cov.samplefreq=500,
  greedy=FALSE,
  ValidateFun=function(...){TRUE},
  store.ndx=NULL){
  
  if(is.null(store.ndx)){store.ndx<-1:nreps}
  if(is.null(npars)){npars<-length(init.fun())}
  posterior<-array(dim=c(nchains,npars,length(store.ndx)))
  accepted<-array(dim=c(nchains,length(store.ndx)-1))
  likelihood.ratio<-array(dim=c(nchains,length(store.ndx)-1))

  jump.state.array<-array(dim=c(nchains,nreps%/%cov.samplefreq,npars,npars))
  dimnames(jump.state.array)<-list(Chain=1:nchains,
                                   iteration=1:(nreps%/%cov.samplefreq),
                                   var1=names(init.fun()),
                                   var2=names(init.fun()))
  adaptive.tuning<-(2.4^2/npars)
  for(k in 1:nchains){
    cat(paste("\nchain ",k,"\n",sep=""))
    theta.post<-array(dim=c(npars,nreps)) 
    
    if(is.null(jump.cov)){jump.cov<-diag(npars)*2.4^2/npars}
    current.jump.cov<-jump.cov 
    for(i in 1:nreps){
      if(!(i%%round(nreps/100))){cat("|")} 
      if(i==1){
        x.scale<-init.fun()
        rownames(theta.post)<-names(x.scale)
        theta.post[,1]<-x.scale
        x<-rep(1,length(x.scale))
        x[x.scale==0]<-0
        names(x)<-names(x.scale)
        x.scale[x.scale==0]<-1
        prev.loglik<-LikelihoodFun(x*x.scale)
      } 
      
      invalid.proposal<-TRUE
      while(invalid.proposal){
        x.star<-jump.fun(x,current.jump.cov)
        theta.star<-x.scale*x.star
        if(any(is.na(theta.star))){browser()}
        invalid.proposal<-!ValidateFun(theta.star)
      }
      star.loglik<-LikelihoodFun(theta.star)
      
      log.lr<-star.loglik-prev.loglik
      
      lr<-exp(log.lr)
      #if(i==nreps){browser()}
      if(is.finite(lr)){
        accept<-(runif(1)<lr)}else{
          if(!is.na(lr)){
            accept<-(lr==Inf)}else{
              accept<-0}
        }
      if(i%in%store.ndx){accepted[k,i-1]<-accept
      likelihood.ratio[k,i-1]<-log.lr}
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
        if(greedy){Xhat.now<-(Xhat.last*min(99,(i-1))+X.now)/min(100,i)}
        cov.update<-(i-1)*Xhat.last%*%t(Xhat.last)-i*Xhat.now%*%t(Xhat.now)+X.now%*%t(X.now)
        #cov.update<-(X.now-Xhat.last)%*%t(X.now-Xhat.last)

        new.jump.cov<-(current.jump.cov*(i-2)/(i-1)+
                         adaptive.tuning/(i-1)*(cov.update+epsilon*diag(npars)))
        if(greedy){
          new.jump.cov<-(current.jump.cov*min(98,(i-2))/min(99,(i-1))+
                           adaptive.tuning/min(99,(i-1))*(cov.update+epsilon*diag(npars)))
        }
        if(is.positive.semidefinite(new.jump.cov)){current.jump.cov<-new.jump.cov}
        if(i%%cov.samplefreq==0){jump.state.array[k,i%/%cov.samplefreq,,]<-current.jump.cov}    
      }
      
    }
    posterior[k,,]<-theta.post[store.ndx,]
  }
  dimnames(posterior)[[2]]<-rownames(theta.post)
  return(list(Posterior=posterior,Accepted=accepted,Likelihood.ratio=likelihood.ratio,jumping.kernel=jump.state.array))
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
