
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


### loading true data
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

simonsen.truedata<-list(N=N,
                    Nobs=Nobs,
                         I=3,
                         SamplingTimes=simonsen.long$Time,
                         ID=simonsen.long$ID.seq,
                         TestData=as.matrix(simonsen.long[,c("IGG","IGA","IGM")])
                    )


##The parameters of the IGM, as estimated from the one-test model with 50 obs.
#parameter means
theta<-structure(c(0.16387042692861, 2.59639256030605, 0.0226677144229636, 
            1.17081744337851), .Names = c("X.star", "D", "a", "S"))
theta["D"]<-60  ##Increasing the time until pathogen disappears to 60 days.
 #individual parameter covariance matrix, log scale
igCov<-structure(c(0.37780789905664, 0.0737747515320145, -0.0316374343170642, 
            0.0646683272440163, 0.0737747515320145, 0.16532315340025, 0.0119353468806589, 
            -0.0134148279864887, -0.0316374343170642, 0.0119353468806589, 
            0.18364935997489, 0.0032842012043439, 0.0646683272440163, -0.0134148279864887, 
            0.0032842012043439, 0.12894598719025), .Dim = c(4L, 4L))
##measurement error, log scale
tau<-0.6

## This generates repeated measurement data based on the Simonsen model, and assuming equal measurement error in AB and NA
SimonsenDataGen<-function(N,times,times.seq=1:500,par.means,par.covariance,measurement.error,ntests=2){
  individual.pars<-exp(rmvnorm(N,log(par.means),par.covariance))
  periods.ndx<- c(1,round(rev(1/2^(1:(times-1))*max(times.seq))),max(times.seq))
  obs.times<-t(sapply(seq(N),function(x){
    obs.time<-sapply(1:times,function(x)sample(times.seq[periods.ndx[x]:periods.ndx[x+1]],1))
    return(obs.time)
  }))
  testmeans<-apply(cbind(individual.pars,obs.times),1,function(x){
    pars<-x[1:4]
    names(pars)<- c("X.star", "D", "a", "S")
    do.call("igCurve",as.list(pars))(x[-c(1:4)])
  })
  test.obs<-rlnorm(length(testmeans),log(c(t(testmeans))),measurement.error)
  na.mean<-sapply(1:N,function(x){
    ifelse(obs.times[x,]>individual.pars[x,"D"],1/10^6,(1-obs.times[x,]/individual.pars[x,"D"]))
  })
  na.obs<-na.mean
  na.obs[na.obs>0]<-rlnorm(sum(na.obs>0),log(na.obs[which(na.obs>0)]),measurement.error)
  test.results<-data.frame(ID=rep(seq(N),times),Time=c(obs.times),IG=test.obs,IG.true=c(t(testmeans)),na.obs=c(t(na.obs)))
  return(test.results)
}
set.seed(1000)
simonsen.abna<-SimonsenDataGen(50,6,par.means=theta,par.covariance = igCov,measurement.error=tau)
#ggplot(data=tmp,aes(x=Time,y=IG,group=ID))+geom_line()
N<-50
Nobs<-6*50
simonsen.pseudo<-list(N=N,Nobs=Nobs,SamplingTimes=simonsen.abna$Time,
                            ID=simonsen.abna$ID,
                            TestData=as.matrix(simonsen.abna[,c("IG","na.obs")])
                          )

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



niter<-10000
warmup.iter=round(niter/2)
rng_seed<-1000:1003
simonsen.model <- stan(file.path(script.path,'stan implementation/simonsen2009_na_ab.stan'), data=simonsen.pseudo, chains = 0)
testing<-stan(fit=simonsen.model,
     data = simonsen.pseudo,
     seed=rng_seed[1],
      warmup=25,
      iter = 50,
      chains = 2,refresh=-1,
      init = "random")
sflist <- 
  mclapply(1:4, mc.cores = 4, function(i){
                            stan(fit = simonsen.model,
                            data = simonsen.pseudo,
                            seed=rng_seed[i],
                             warmup=warmup.iter,
                            iter = niter,
                            chains = 1,chain_id=i,refresh=-1,
                            init = "random")})



simonsen.posterior<-sflist2stanfit(sflist)
save(simonsen.posterior,file=file.path(data.path,"fitted-stan-models/simonsen_igna_50N_2000iter.Rdata"))


tail.ndx<-tail(seq(dim(simonsen.posterior)[[1]]),300)
pars.est<-monitor(extract(simonsen.posterior,pars="theta2IgLogmu", permuted = FALSE, inc_warmup = TRUE)[tail.ndx,,])
individual.pars<-monitor(extract(simonsen.posterior,pars="theta1IgLogmu",permuted=FALSE,inc_warmup=FALSE)[tail.ndx,,])


mean.pars.chains<-apply(extract(simonsen.posterior,pars="theta2IgLogmu", permuted = FALSE, inc_warmup = TRUE)[tail.ndx,,],c(2,3),mean)
mean.pars<-c(pars.est[,"mean"])
names(mean.pars)<-c("X.star","D","a","S")
colnames(mean.pars.chains)<-c("X.star","D","a","S")
individual.matrix<-matrix(individual.pars[,"mean"],nrow=50,dimnames=list(1:50,c("X.star","D","a","S")))


individual.pred<-apply(individual.matrix,1,function(x){
  x<-exp(x)
  igCurve.est<-do.call("igCurve",as.list(x))
  path.est<-igCurve.est(1:500)
  return(path.est)}
)
chain.pred<-apply(mean.pars.chains,1,function(x){
  x<-exp(x)
  igCurve.est<-do.call("igCurve",as.list(x))
  path.est<-igCurve.est(1:500)
  return(path.est)}
)
mean.pred<-do.call("igCurve",as.list(exp(mean.pars)))(1:500)

##Plotting predicted curves and true curves over time
ggplot(data=melt(individual.pred,varnames=c("Time","ID")),aes(x=Time,y=value,group=ID))+
  geom_line(alpha=0.5,col="red")+
  geom_line(data=melt(chain.pred,varnames=c("Time","ID")),size=2,aes(col=ID))+
  geom_point(data=melt(simonsen.abna,id.vars=c("ID","Time"),measure.vars="IG"),aes(x=Time,y=value))+
  geom_line(data=melt(simonsen.abna,id.vars=c("ID","Time"),measure.vars="IG"),aes(x=Time,y=value,group=ID),alpha=0.5)



#plotting predicted vs true values.
ggplot(data=merge(simonsen.abna,melt(individual.pred,varnames=c("Time","ID")),
             by.x=c("Time","ID"),by.y=c("Time","ID")),
       aes(x=log(value),y=log(IG)))+
  geom_point()+
  geom_line(data=data.frame(value=range(individual.pred),IG=range(individual.pred)))

ggplot(melt(extract(simonsen.posterior,pars="lp__",permute=FALSE)),
       aes(x=iterations,y=value,group=chains,color=chains))+geom_line()