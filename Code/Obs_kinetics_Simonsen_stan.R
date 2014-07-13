
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

source(file.path(script.path,"MH files/metropolisHastingsAdaptiveV3.R"))
#source(file.path(script.path,"Simonsen Model code/ObsKinetics_Simonsen_full_model.R"))


###Generating some data
ig.values.simonsen<-dget(file=file.path(data.path,"simonesn_ig_array.txt"))
obs.array.simonsen<-dget(file=file.path(data.path,"simonsen_obstime_array.txt"))
n.igg<-3
n.pars<-4
N<-nrow(obs.array.simonsen)


###  Running STAN here...

test.standata<-list(N=N,
                         T=4,
                         I=3,
                         SamplingTimes=obs.array.simonsen,
                         TestData=ig.values.simonsen)#make this generic!
niter<-100
warmup.iter=50
library(rstan)
simonsen.results.stan<-stan(file = file.path(script.path,'stan implementation/simonsen2009.stan'),
                            data = test.standata,
                            verbose = TRUE,
                             warmup=warmup.iter,
                            iter = niter,
                            chains = 4,
                            init = "random")


rng_seed <- as.numeric(Sys.time())
simonsen.model <- stan(file.path(script.path,'stan implementation/simonsen2009.stan'), data=test.standata, chains = 0)
sflist <- 
  mclapply(1:4, mc.cores = 4, 
           function(i) stan(fit=simonsen.model, data=test.standata, seed = rng_seed[i], 
                            chains = 1, chain_id = i, refresh = -1,iter=niter,warmup=warmup.iter))


simonsen.posterior<-sflist2stanfit(sflist)
tmp<-extract(simonsen.posterior,pars="theta2IgLogmu")
igg<-3
plot(X1~time,type="l",data=data.frame(exp(tmp[[1]])[,igg,],iter=rep(1:4,each=500),time=rep(1:500,4))[1:500,])
lines(X1~time,type="l",data=data.frame(exp(tmp[[1]])[,igg,],iter=rep(1:4,each=500),time=rep(1:500,4))[1:500+500,],col=2)
lines(X1~time,type="l",data=data.frame(exp(tmp[[1]])[,igg,],iter=rep(1:4,each=500),time=rep(1:500,4))[1:500+1000,],col=3)
lines(X1~time,type="l",data=data.frame(exp(tmp[[1]])[,igg,],iter=rep(1:4,each=500),time=rep(1:500,4))[1:500+1500,],col=4)


tmp<-data.frame(extract(simonsen.posterior,pars="lp__"),time=rep(1:500,4),chain=rep(1:4,each=500))
xyplot(lp__~time,group=chain,type="l",data=tmp,auto.key=T)

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