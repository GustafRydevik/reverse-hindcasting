
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
autolib(matrixcalc)
source(file.path(script.path,"HindcastingLib/TestCurves_Simonsen2009.R"))

#source(file.path(script.path,"MH files/metropolisHastingsAdaptiveV3.R"))
#source(file.path(script.path,"Simonsen Model code/ObsKinetics_Simonsen_full_model.R"))


###Generating some data
simonsen.long<-dget(file=file.path(data.path,"simonsen_long.txt"))
simonsen.long<-simonsen.long[!is.nan(simonsen.long$IGG),]

simonsen.long$ID.seq<-as.numeric(factor(simonsen.long$ID))
simonsen.long<-simonsen.long[simonsen.long$Time>0,] #Dropping obs with Time<0
simonsen.long<-simonsen.long[!(rowSums(simonsen.long[,c("IGG","IGA","IGM")]<=0)),] #Dropping obs with any measurement < or equal to 0

n.igg<-3
n.pars<-4
Nobs<-nrow(simonsen.long)
N<-length(unique(simonsen.long$ID))
###  Running STAN here...

test.standata<-list(N=N,
                    Nobs=Nobs,
                         I=3,
                         SamplingTimes=simonsen.long$Time,
                         ID=simonsen.long$ID.seq,
                         TestData=as.matrix(simonsen.long[,c("IGG","IGA","IGM")])
                    )
simonsen.pseudo<-list(N=N,Nobs=Nobs,I=3,SamplingTimes=sample(simonsen.long$Time,length(simonsen.long$Time),replace=T),
                            ID=simonsen.long$ID.seq,
                            TestData=matrix(rlnorm(nrow(simonsen.long)*3, 
                                                   mean(as.matrix(simonsen.long[,c("IGG","IGA","IGM")])),
                                                      log(1.5)),ncol=3))

dput(simonsen.pseudo,file=file.path(data.path,"SimonsenPseudoData.txt"))
library(rstan)
library(parallel)
my.inits<-list(
  theta2IgLogmu=log(matrix(rep(c(Xstar=0.15,D=100,a=0.01,S=0.05),each=3),nrow=3)),
  theta1IgLogmu=log(aperm(array(rep(matrix(rep(c(Xstar=0.15,D=100,a=0.01,S=0.05),each=3),nrow=3),
                       N),dim=c(3,4,N)),c(3,1,2))),
  igCorr=aperm(array(rep(diag(4)+(1-diag(4))*0.0001,3),dim=c(4,4,3)),c(3,2,1)),
  igTau=array(rep(rep(1,4),3),dim=c(3,4)),
  omegaCorr=diag(3)+(1-diag(3))*0.0001,
  tau=rep(1,3)*log(1.5)
  )



set.seed(1000)
id.sampled<-sample(unique(simonsen.long$ID),50)
simonsen.small<-subset(simonsen.long,ID%in%id.sampled)
simonsen.small$ID.seq<-seq_along(id.sampled)[match(simonsen.small$ID,id.sampled)]
small.standata<-list(N=length(id.sampled),
                    Nobs=nrow(simonsen.small),
                    I=1,
                    SamplingTimes=simonsen.small$Time,
                    ID=simonsen.small$ID.seq,
                    TestData=c(simonsen.small[,c("IGM")])
)
niter<-500
warmup.iter=250
rng_seed<-1000:1003
simonsen.model <- stan(file.path(script.path,'stan implementation/simonsen2009_onetest.stan'), data=small.standata, chains = 0)
testing<-stan(fit=simonsen.model,
     data = small.standata,
     seed=rng_seed[1],
      warmup=25,
      iter = 50,
      chains = 2,refresh=-1,
      init = "random")
sflist <- 
  mclapply(1:4, mc.cores = 4, function(i){
                            stan(fit = simonsen.model,
                            data = small.standata,
                            seed=rng_seed[i],
                             warmup=warmup.iter,
                            iter = niter,
                            chains = 1,chain_id=i,refresh=-1,
                            init = "random")})



simonsen.posterior<-sflist2stanfit(sflist)
#save(simonsen.posterior,file=file.path(data.path,"fitted-stan-models/simonsen3Ig_50N_500iter.Rdata"))


pars.est<-monitor(extract(simonsen.posterior,pars="theta2IgLogmu", permuted = FALSE, inc_warmup = TRUE))
individual.pars<-monitor(extract(simonsen.posterior,pars="theta1IgLogmu",permuted=FALSE,inc_warmup=FALSE))

mean.pars<-c(pars.est[,"mean"])
names(mean.pars)<-c("X.star","D","a","S")
individual.matrix<-matrix(individual.pars[,"mean"],nrow=50,dimnames=list(1:50,c("X.star","D","a","S")))


individual.pred<-apply(individual.matrix,1,function(x){
  x<-exp(x)
  igCurve.est<-do.call("igCurve",as.list(x))
  path.est<-igCurve.est(1:500)
  return(path.est)}
)
mean.pred<-do.call("igCurve",as.list(exp(mean.pars)))(1:500)

plot(1:500,individual.pred[,1],type="l",ylim=c(0,5))
points(simonsen.small$Time,simonsen.small$IGG)

ggplot(data=melt(individual.pred,varnames=c("Time","ID")),aes(x=Time,y=value,group=ID))+
  geom_line(alpha=0.5,col="red")+
  geom_line(data=data.frame(value=mean.pred,Time=1:500,ID=1),size=2,col="darkred")+
  geom_point(data=melt(simonsen.small,id.vars=c("ID","ID.seq","Time"),measure.vars="IGM"),aes(x=Time,y=value))+
  geom_line(data=melt(simonsen.small,id.vars=c("ID","ID.seq","Time"),measure.vars="IGM"),aes(x=Time,y=value,group=ID),alpha=0.5)


dimnames(pars.est)[[1]]<-paste(rep(c("X.star","D","a","S"),each=3),c("igg","iga","igm"),sep=".")
igg<-3
plot(X1~time,type="l",data=data.frame(exp(tmp[[1]])[,igg,],iter=rep(1:4,each=800),time=rep(1:800,4))[1:800,])

lines(X1~time,type="l",data=data.frame(exp(tmp[[1]])[,igg,],iter=rep(1:4,each=800),time=rep(1:800,4))[1:800+800,],col=2)
lines(X1~time,type="l",data=data.frame(exp(tmp[[1]])[,igg,],iter=rep(1:4,each=800),time=rep(1:800,4))[1:800+1600,],col=3)
lines(X1~time,type="l",data=data.frame(exp(tmp[[1]])[,igg,],iter=rep(1:4,each=800),time=rep(1:800,4))[1:800+2400,],col=4)

tmp.last100.ndx<-c(sapply(c(8,16,24,32)*100,function(x)(x-100):(x)))
pars.individual.mean100<-exp(apply(tmp2[[1]][tmp.last100.ndx,,,],c(2:4),sample,1))


dimnames(pars.individual.mean100)[[3]]<-c("X.star","D","a","S")

individual.pred<-apply(pars.individual.mean100,c(1,2),function(x){
  igCurve.est<-do.call("igCurve",as.list(x))
  path.est<-igCurve.est(1:500)
  return(path.est)}
)


iggCurve.est<-do.call("igCurve",as.list(pars.mean100[3,]))
plot(IGG~Time,group=ID,data=simonsen.long,xlim=c(0,500))


xvals <- tapply(simonsen.long$Time,simonsen.long$ID,function(x) return(x))
yvals <- tapply(simonsen.long$IGG,simonsen.long$ID,function(x) return(x))

plot(x=1:max(unlist(xvals)),ylim=(c(0,max(unlist(yvals)))),type="n",xlim=c(0,500))
# thanks to @BenBolker for refining this next key line
mapply(lines,xvals,yvals ,col=c("red"),alpha=0.5,pch=1,lty=3,type="o",wed=0.5)
matplot(1:500,(individual.pred[,,1]),type="l",add=T,lwd=0.7)

points(1:500,iggCurve.est(1:500),type="l",col="red",lwd=2)



lines(1:500,iggCurve.est(1:500),type="l",col="black",lwd=5)

tmp<-data.frame(extract(simonsen.posterior,pars="lp__"),time=rep(1:500,4),chain=rep(1:4,each=500))
xyplot(lp__~time,group=chain,type="l",data=tmp,auto.key=T)

tmp<-extract(testing,pars="XStar",permuted=FALSE)
melt(tmp)
xyplot(value~iterations,group=chains,type="l",data=melt(tmp),auto.key=T)

lines(X1~time,type="l",data=tmp,col=2)
lines(X1~time,type="l",data=tmp,col=3)
lines(X1~time,type="l",data=data.frame(exp(tmp[[1]])[,igg,],iter=rep(1:4,each=500),time=rep(1:500,4))[1:500+1500,],col=4)

#                               list(list(theta1IgLogmu=log(aperm(array(rep(c(X.star=0.15,D=15,a=0.05,S=0.1, 
#                                                         X.star=0.3,D=15,a=0.1,S=0.2,
#                                                         X.star=0.45,D=15,a=0.15,S=0.3),N),dim=c(4,3,N)),c(3,2,1))),
#                                              theta2IgLogmu=log(matrix(c(X.star=0.15,D=15,a=0.05,S=0.1, 
#                                                                        X.star=0.3,D=15,a=0.1,S=0.2,
#                                                                        X.star=0.45,D=15,a=0.15,S=0.3),nrow=3,byrow=T)),
#                                              igCov=aperm(array(rep(diag(4),3),dim=c(4,4,3)),c(3,2,1)),
#                                              omega=diag(3))))
                          



temp<-metropolis(jump.fun=jump.mvnorm,
                  init.fun=
                    SimonsenInit,
                    #ResumeInit(temp$Posterior),
                  LikelihoodFun=SimonsenLLtest,
                  npars=NULL,
                  nreps=10000,
                  jump.cov=
                    #temp$jumping.kernel[1,20,,],
                  NULL,
                  ValidateFun=TestingCovariance,
                  Adapt=TRUE,greedy=T#,store.ndx=c(seq(1,50000-1000,length.out=10000),49001:50000)
)