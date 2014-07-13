test<-lapply(1:dim(temp2$Posterior)[1],function(x){
  x<-temp2$Posterior[x,,]
  x<-t(x)
  class(x)<-"mcmc"
  attr(x,"mcpar")<-c(1,dim(temp2$Posterior)[3],1)
  return(x)
})

prev.par<-par
par(mfrow=c(2,2))

matplot(t(exp(temp2$Posterior[,"curve.pars.X.star.log",])),type="l")
matplot(t(exp(temp2$Posterior[,"curve.pars.S.log",])),type="l")
matplot(t(exp(temp2$Posterior[,"curve.pars.D.log",])),type="l")
matplot(t(exp(temp2$Posterior[,"curve.pars.a.log",])),type="l")
par(prev.par)

plot(filter(ts(t(temp2$Accepted)),rep(0.01,100)))
matplot(t(temp2$jumping.kernel[,,1,1]),type="l")


library(rjags)
gelman.diag(mcmc.list(test)[,1:7])
curve.pars<-apply(temp2$Posterior[,1:4,],2,mean)
names(curve.pars)<-gsub("curve\\.pars\\.","",names(curve.pars))


plot(seq(1,sampling.time,by=0.1),igCurve.set(seq(1,sampling.time,by=0.1)),type="l",col="grey",lwd=2,ylim=c(0,5))


for(chain in 1:3){
t2<-do.call("cbind",list(temp2$Posterior[1,,],temp2$Posterior[2,,],temp2$Posterior[3,,]))
t2<-temp2$Posterior[chain,,]
dimnames(t2)[[1]][1:4]<-gsub("curve\\.pars\\.","",dimnames(t2)[[1]][1:4])

curve.matrix<-apply(t2[1:4,sample(1:dim(t2)[2],1000)],2,function(x){
  names(x)<-gsub("\\.log","",names(x))
  names(x)<-gsub("curve.pars.","",names(x))
  igCurve.tmp<-do.call("igCurve",as.list(exp(x)))
  curve.tmp<-igCurve.tmp(seq(1,sampling.time,by=0.1))
  return(curve.tmp)                    
})

curve.mean<-rowMeans(curve.matrix)
curve.05<-apply(curve.matrix,1,function(x)sort(x)[floor(0.05*length(x))])
curve.95<-apply(curve.matrix,1,function(x)sort(x)[floor(0.95*length(x))])
lines(seq(1,sampling.time,by=0.1),curve.mean,col="red",lwd=2)
lines(seq(1,sampling.time,by=0.1),curve.05,lty=2,col="red",lwd=2)
lines(seq(1,sampling.time,by=0.1),curve.95,lty=2,col="red",lwd=2)
}




est.pars<-apply(temp2$Posterior[,,1:10000],2,mean)
true.pars<-est.pars
est.pars[-c(1:7)]<-logit((btv.test.data$obs.time-btv.test.data$infection.time)/btv.test.data$obs.time)
true.pars[-(1:7)]<-logit((btv.test.data$obs.time-btv.test.data$infection.time)/btv.test.data$obs.time)
true.pars[1:7]<-c(log(c(X.star=0.15,S=0.3,a=0.05,D=15)),Ab.sdlog=log(1.1),Na.sdlog=log(1.1),log(1/10))
ll.fun(est.pars)
ll.fun(true.pars)

