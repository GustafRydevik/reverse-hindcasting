
##based on http://mbjoseph.github.io/blog/2013/09/08/metropolis/
metropolis<-function(
  jump.fun,
  init.fun,
  nreps=10000,
  llfun=function(x)dnorm(x,log=TRUE),
  nchains=3,
  npars=2){

  posterior<-array(dim=c(nchains,npars,nreps))
  accepted<-array(dim=c(nchains,nreps-1))
  
  for(k in 1:nchains){
    x<-init.fun()
    prev.loglik<-do.call("llfun",x)
    
    theta.post<-array(dim=c(npars,nreps))
    rownames(theta.post)<-names(x)

    for(i in 1:npars){
      theta.post[i,1]<-x[[i]]
    } 
    
    for(i in 2:nreps){
      theta.star<-jump.fun(unlist(x))
      
      star.loglik<-do.call("llfun",as.list(theta.star))
      
      log.lr<-star.loglik-prev.loglik
      lr<-exp(log.lr)
      
      accept<-rbinom(1,1,prob=min(lr,1))
      if(!is.finite(star.loglik)){accept<-0}
      accepted[k,i-1]<-accept
      if(accept==1){
        theta.post[,i]<-x<-theta.star
        prev.loglik<-star.loglik
    }else{
      theta.post[,i]<- theta.post[,i-1]}
    }
  
    posterior[k,,]<-theta.post
  }
  return(list(Posterior=posterior,Accepted=accepted))
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


