
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
source(file.path(script.path,"Simonsen Model code/ObsKinetics_Simonsen_full_model.R"))


convert.from.triangle<-function(triangular.elements,dim){    ###Transforming vector elements of lower triangular matrices, to an array of matrices
  n.tri<-(dim^2+dim)/2
  result.array<-aperm(sapply(1:(length(triangular.elements)/n.tri),function(i){   
    orig.matrix<-diag(dim)
    orig.matrix[lower.tri(orig.matrix,diag=T)]<-triangular.elements[n.tri*(i-1)+seq_len(n.tri)]
    orig.matrix<-orig.matrix+
      t(orig.matrix-diag(diag(orig.matrix)))
    return(orig.matrix)
  },simplify="array"),c(3,1,2))
  if(dim(result.array)[1]==1){result.array<-result.array[1,,]}
  return(result.array)
}


###Generating some data
N<-200
igg.curve<-igCurve(X.star=0.15,D=15,a=0.05,S=0.1)  ## Change these to something more realistic later...
iga.curve<-igCurve(X.star=0.3,D=15,a=0.1,S=0.2)
igm.curve<-igCurve(X.star=0.45,D=15,a=0.15,S=0.3)
obs.array<-array(data=sort(sample(1:50,N*4,replace=T)),dim=c(N,4))
ig.values.sim<-array(data=rlnorm(N*4*3,log(c(igg.curve(c(obs.array)),
                            iga.curve(c(obs.array)),
                            igm.curve(c(obs.array)))),log(1.3)),dim=c(N,4,3))


n.igg<-3
n.pars<-4



##Setting the priors
Theta2.igg.prior<-function(...){return(0)}  #Priors for the mean of curve parameters
Theta2.igm.prior<-function(...){return(0)}
Theta2.iga.prior<-function(...){return(0)}
omega.prior<-function(...){return(0)}



##A jumping function for generating new proposal values based on current parameter values

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

###initialization functions
SimonsenInit<-function(){
  theta2.ig.mu<-array(data=log(c(c(X.star=0.15,D=15,a=0.05,S=0.3),
                                 c(X.star=0.3,D=15,a=0.1,S=0.3),
                                 c(X.star=0.45,D=15,a=0.15,S=0.3))),
                      dim=c(n.igg,n.pars),dimnames=list(ig=c("igg","iga","igm"),
                                                        pars=c("X.star","D","a","S")))
  
  theta2.ig.sigma<-aperm(sapply(1:n.igg,function(x){diag(4)},simplify="array"))
  dimnames(theta2.ig.sigma)<-list(ig=c("igg","iga","igma"),
                                  pars1=c("X.star","D","a","S"),pars2=c("X.star","D","a","S"))
  omega<-diag(3)
  dimnames(omega)<-list(ig1=c("igg","iga","igm"),
                        ig2=c("igg","iga","igm"))
  
  theta1.ig.mu<-aperm(sapply(1:n.igg,function(x){rmnorm(N,mean=theta2.ig.mu[x,],theta2.ig.sigma[x,,])},simplify="array"),c(1,3,2))
  dimnames(theta1.ig.mu)<-list(id=(1:N),
                               ig=c("igg","iga","igm"),
                               pars=c("X.star","D","a","S"))
  
  parameter.vector<-c(c(theta1.ig.mu),
                    c(theta2.ig.mu),
                    c(apply(theta2.ig.sigma,1,function(x)x[lower.tri(x,diag=T)])),
                    omega[lower.tri(omega,diag=T)])
  
  
  tri.ndx.npars<-which(lower.tri(diag(n.pars),diag=T),arr.ind=T)
  tri.ndx.igg<-which(lower.tri(diag(n.igg),diag=T),arr.ind=T)
  names(parameter.vector)<-c(  ##This gives names to each element of the looong vector, so we can look it up later. 
    do.call("paste",  
            args=c(list("theta1.ig.mu"),
                   expand.grid(dimnames(theta1.ig.mu),stringsAsFactors=FALSE)[c(2,3,1)],
                   list(sep="."))),
    do.call("paste",
            args=c(list("theta2.ig.mu"),
                   expand.grid(dimnames(theta2.ig.mu),stringsAsFactors=FALSE),
                   list(sep="."))),
    do.call("paste",
            args=c(list("theta2.ig.sigma"),
                   list(rep(dimnames(theta2.ig.sigma)[[1]],each=nrow(tri.ndx.npars))),
                   list(dimnames(theta2.ig.sigma)[[2]][tri.ndx.npars[,1]]),
                   list(dimnames(theta2.ig.sigma)[[3]][tri.ndx.npars[,2]]),
                   list(sep="."))),
    do.call("paste",
            args=c(list("omega"),
                   list(dimnames(omega)[[1]][tri.ndx.igg[,1]]),
                   list(dimnames(omega)[[2]][tri.ndx.igg[,2]]),
                   list(sep=".")))
  )
  
  #btv.test.data$infection.time#sample(seq_along(btv.2008.interpolated$y),61,prob=ceiling(btv.2008.interpolated$y*10)/10,replace=T)
  return(parameter.vector)
}

ResumeInit<-function(Posterior,nchains=3){
  called<-0
  function(){
    called<<-called%%nchains+1
    
    par.means<-rowMeans(Posterior[called,,dim(Posterior)[3]-c(1000:0)],2)
    return(par.means)
  }}



##Likelihood function
SimonsenLLtest<-Simonsen.likelihood.gen(ig.values.array=ig.values.sim,     #observed antibody test values 
                                          #Each of length N do dim=c(n,4,3)          
                                          obs.array=obs.array,            #times of observation
                                          # each of length N 
                                          Theta2.igg.prior=Theta2.igg.prior,  #Priors for the mean of curve parameters
                                          Theta2.igm.prior=Theta2.igm.prior,
                                          Theta2.iga.prior=Theta2.iga.prior,
                                          omega.prior=omega.prior,       ##Prior for the correlation of observed test values at a given time 
                                          N,n.times=4,n.igg=3,n.pars=4)


##Function for validating pos. definite 

TestingCovariance<-function(pars){
  theta2.ig.sigma.triangle<-pars[grep("theta2.ig.sigma",names(pars))] ##3*(4+3+2+1) pars.. How do we force a proper corr matrix?
  theta2.ig.sigma<-convert.from.triangle(theta2.ig.sigma.triangle,n.pars)
  omega.triangle<-pars[grep("omega",names(pars))] ###3*2*1 pars: cov matrix for the test values
  omega<-convert.from.triangle(omega.triangle,n.igg)
  all.semidefinite<-is.positive.semi.definite(omega)&
    all(apply(theta2.ig.sigma,1,is.positive.semi.definite))
  return(all.semidefinite)
}

###  Running the MCMC here...


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

temp<-metropolis(jump.fun=jump.mvnorm,
           init.fun=
             #SimonsenInit,
             ResumeInit(temp$Posterior),
           LikelihoodFun=SimonsenLLtest,
           npars=NULL,
           nreps=100000,
           jump.cov=
             temp$jumping.kernel[1,20,,],
             #NULL,
           ValidateFun=TestingCovariance,
           Adapt=TRUE,greedy=F,,store.ndx=c(seq(1,100000-1000,length.out=10000),99001:100000)
)

temp<-metropolis(jump.fun=jump.mvnorm,
                  init.fun=
                    #SimonsenInit,
                    ResumeInit(temp$Posterior),
                  LikelihoodFun=SimonsenLLtest,
                  npars=NULL,
                  nreps=150000,
                  jump.cov=
                    temp$jumping.kernel[1,20,,],
                  #NULL,
                  ValidateFun=TestingCovariance,
                  Adapt=TRUE,greedy=F,store.ndx=c(seq(1,150000-1000,length.out=10000),149001:150000)
)

temp<-metropolis(jump.fun=jump.mvnorm,
                  init.fun=
                    #SimonsenInit,
                    ResumeInit(temp3$Posterior),
                  LikelihoodFun=SimonsenLLtest,
                  npars=NULL,
                  nreps=150000,
                  jump.cov=
                    temp$jumping.kernel[1,20,,],
                  #NULL,
                  ValidateFun=TestingCovariance,
                  Adapt=TRUE,greedy=F,store.ndx==c(seq(1,150000-1000,length.out=10000),149001:150000)
)