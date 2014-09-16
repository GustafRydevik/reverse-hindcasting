
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



load(file=file.path(data.path,"fitted-stan-models/simonsen_igna_50N_2000iter.Rdata"))
igna.posterior<-simonsen.posterior

##Stanfit.object should have all variables from the simonsen stan files
##model.data should have columns Time, ID, and IG
simonsen.plot<-function(stanfit.object,model.data=NULL){
  tail.ndx<-tail(seq(dim(stanfit.object)[[1]]),300)
  pars.est<-monitor(extract(stanfit.object,pars="theta2IgLogmu", permuted = FALSE, inc_warmup = TRUE)[tail.ndx,,])
  individual.pars<-monitor(extract(stanfit.object,pars="theta1IgLogmu",permuted=FALSE,inc_warmup=FALSE)[tail.ndx,,])
  
  
  mean.pars.chains<-apply(extract(stanfit.object,pars="theta2IgLogmu", permuted = FALSE, inc_warmup = TRUE)[tail.ndx,,],c(2,3),mean)
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
    geom_point(data=melt(model.data,id.vars=c("ID","Time"),measure.vars="IG"),aes(x=Time,y=value))+
    geom_line(data=melt(model.data,id.vars=c("ID","Time"),measure.vars="IG"),aes(x=Time,y=value,group=ID),alpha=0.5)
  
  
  #plotting predicted vs true values.
  ggplot(data=merge(model.data,melt(individual.pred,varnames=c("Time","ID")),
                    by.x=c("Time","ID"),by.y=c("Time","ID")),
         aes(x=log(value),y=log(IG)))+
    geom_point()+
    geom_line(data=data.frame(value=range(individual.pred),IG=range(individual.pred)))
  
  ggplot(melt(extract(stanfit.object,pars="lp__",permute=FALSE)),
         aes(x=iterations,y=value,group=chains,color=chains))+geom_line()
}