
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
  Adapt=TRUE){
  posterior<-array(dim=c(nchains,npars,nreps))
  accepted<-array(dim=c(nchains,nreps-1))
  likelihood.ratio<-array(dim=c(nchains,nreps-1))
  jump.state.list<-vector(length=nchains,mode="list")
  for(k in 1:nchains){
    cat(paste("\nchain ",k,"\n",sep=""))
    x<-init.fun()
    prev.loglik<-do.call("ll.fun",x)
    
    theta.post<-array(dim=c(npars,nreps))
    rownames(theta.post)<-names(unlist(x))
    current.jump.cov<-jump.cov
    for(i in 1:npars){
      theta.post[i,1]<-unlist(x)[[i]]
    } 
    
    for(i in 2:nreps){
      if(!(i%%round(nreps/100))){cat("|")}
      theta.star<-jump.fun(x,current.jump.cov)
      
      star.loglik<-do.call("ll.fun",as.list(theta.star))
      
      log.lr<-star.loglik-prev.loglik
      lr<-exp(log.lr)
      #if(i==nreps){browser()}
      if(is.finite(lr)){
        accept<-(runif(1)<lr)}else{
          if(!is.na(lr)){
            accept<-(lr==Inf)}else{
              accept<-0}
        }
      accepted[k,i-1]<-accept
      likelihood.ratio[k,i-1]<-log.lr
      if(accept==1){
        x<-theta.star
        theta.post[,i]<-unlist(x)
        prev.loglik<-star.loglik
      }else{
        theta.post[,i]<- theta.post[,i-1]}
      if((i%%1000==0)&Adapt){   ###Changing the jumping kernel here - adaptive mcmc?
        recent.trials<-accepted[k,(i-999):(i-1)]
        
        if(is.null(current.jump.cov)|length(current.jump.cov)==0){
          new.jump.cov<-diag(x=unlist(x)/10000000,nrow=length(unlist(x)),ncol=length(unlist(x)))
          jump.state.list[[k]]<-c(jump.state.list[[k]],list(new.jump.cov))
          current.jump.cov<-new.jump.cov
          next
        }
        if(sum(recent.trials,na.rm=T)<10){
          new.jump.cov<-current.jump.cov*0.75}else{
            latest.ndx<-(i-1000+1):(i-1)
            latest.cov<-cov(t(theta.post[,latest.ndx]),use="complete.obs")
            if(!is.positive.semidefinite(latest.cov)){latest.cov<-as.positive.semidefinite(latest.cov)}
            new.jump.cov<-(current.jump.cov+latest.cov)/2
            
            if((mean(accepted[k,latest.ndx])<0.015)|(mean(accepted[k,latest.ndx])>0.5)){
              new.jump.cov<-new.jump.cov*qnorm(0.24/2)/qnorm(mean(mean(accepted[k,(i-999):(i-1)])/2))}
            if(!sum(!new.jump.cov==0)){browser()}}
        jump.state.list[[k]]<-c(jump.state.list[[k]],list(new.jump.cov))
        current.jump.cov<-new.jump.cov
      }
    }
    
    posterior[k,,]<-theta.post
  }
  dimnames(posterior)[[2]]<-rownames(theta.post)
  return(list(Posterior=posterior,Accepted=accepted,Likelihood.ratio=likelihood.ratio,jumping.kernel=jump.state.list))
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
