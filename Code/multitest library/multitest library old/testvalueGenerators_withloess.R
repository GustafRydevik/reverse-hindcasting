#######Deterministic models for NA/antibody interaction
###using differential equations

###Parameters for the lottka-voltera diff. equations below

Pars<-c(alpha.na=30,beta.na=0.5,alpha.ab=5,beta.ab=0.1)
State<-c(nat=1,abt=0)
Time=seq(0,10,by=0.001)


###Lotke-Voltera equation based on Giles' idea.
autolib(deSolve)
ab.na.diffmod<-function(Time,State,Pars){
  with(as.list(c(State,Pars)),{
      dna=nat*(alpha.na-beta.na*abt)
        dab=alpha.ab*nat-beta.ab*abt
        return(list(c(dna,dab)))
      })}
inf.levels.lv <- as.data.frame(ode(func = ab.na.diffmod, y = State, parms = Pars, times = Time))
inf.levels.lv$abt<-2^(inf.levels.lv$abt/12)

###Deterministic nucleic acid function, pulling values from Giles' diff model
naDetFun<-(function(Inf.levels=inf.levels.lv,timescale=NULL){
  function(infection.time){
   sapply(infection.time,function(x){
     if(x>max(Inf.levels$time)){return(Inf.levels[which.max(Inf.levels$time>=x),"nat"])}
     if(x<=max(Inf.levels$time)){Inf.levels[min(which(Inf.levels$time>=x)),"nat"]}}
        )}}
   )()



###Deterministic antibody function, pulling values from Giles' diff model
abDetFun<-(function(Inf.levels=inf.levels.lv,timescale=NULL){
  function(infection.time){
   abt.response<-sapply(infection.time,function(x){
     if(x>max(Inf.levels$time)){return(Inf.levels[which.max(Inf.levels$time>=x),"abt"])}
     if(x<=max(Inf.levels$time)){Inf.levels[min(which(Inf.levels$time>=x)),"abt"]}}
          )
   abt.response[infection.time==0]<-rlnorm(sum(infection.time==0),mean=log(mean(Inf.levels$abt[100])))
   return(abt.response)
	}}
   )()


###Old AB decay function just based on a lognormal distribution of  responses with added half life
AntibodyDecay<-function(response.mean,half.life=1){
function(infection.time,number.infections){
initial.response<-rlnorm(length(infection.time),mean=log(response.mean),sd=log(response.mean)/10)*(infection.time!=0)
decayed.response<- 2^(-(infection.time/half.life))*initial.response+
                       rlnorm(length(infection.time),
                              mean=log(response.mean/2^8),
                              sd=log(response.mean)/10)
return(decayed.response)
}}

###Old NA decay function just based on a gamma distribution of  responses with added half life
NaDecay<-function(infection.duration){
function(infection.time){
realised.infection.duration<-rgamma(length(infection.duration),shape=1,
                                    scale=infection.duration)
test.result<-(infection.time<=realised.infection.duration)*(infection.time!=0)
return(test.result)
}}
 

###Realistic scenarios ###

##Function for optimize
# Pseudo-code:
#   
#   1) set some start parameters 
#   2) Have a generic L-V model function
#   3) Have Some data to fit
#   4) Have a function that generates LV, looks at the time points
#      closest to data time points, and calculate Squared minimum distance
#   5) wrap all in optimize


###Start parameters, generic

#### Generic L-V function
####Wrapper to use in optim
lv.mse<-function(start.pars,
         start.state,
         test1.data,
         test2.data,
         t1.time.data,
         t2.time.data,
         est.timepoints
         ){

  ##L-V function
test.diffmod<-function(Time,State,Pars){
  with(as.list(c(State,Pars)),{
      dt1=t1t*(alpha.t1-beta.t1*t2t)
        dt2=alpha.t1*t1t-beta.t2*t2t
        return(list(c(dt1,dt2)))
      })}

  ##Solving the ode
test.levels <- as.data.frame(ode(func = test.diffmod, 
                                 y = start.state, 
                                 parms =start.pars, times = est.timepoints))
#Calculating MSE
                                 
t1.comparison.ndx<-sapply(t1.time.data,function(x)which.min(abs(x-est.timepoints)))
t2.comparison.ndx<-sapply(t2.time.data,function(x)which.min(abs(x-est.timepoints)))                                
test.mse<-mean(c((test1.data-test.levels[t1.comparison.ndx,"t1t"])^2,
                (test2.data-test.levels[t2.comparison.ndx,"t2t"])^2),na.rm=T)        
return(test.mse)}
                                                 


## Anaplasma
##Stuen 2001 data; AB to Ehrlichia equi
anaplasma.ab<-data.frame(deer.id=rep(1:3,each=7),
           time=rep(c(0,7,14,21,28,56,84),3),
           deer.ab.values=c(Deer.one=c(NA,NA,160,640,1280,160,40),
           deer.two=c(NA,NA,NA,160,160,160,40),
           deer.three=c(NA,NA,NA,NA,NA,NA,NA))
                         )


reindeer.granulocytes<-data.frame(
reindeer.id=c("adult.one", "adult.one", "adult.one", "adult.one", "adult.one", 
"adult.one", "adult.one", "adult.one", "adult.one", "adult.one", 
"adult.one", "adult.two", "adult.two", "adult.two", "adult.two", 
"adult.two", "adult.two", "adult.two", "adult.two", "adult.two", 
"adult.two", "adult.two", "calf", "calf", "calf", "calf", "calf", 
"calf", "calf", "calf", "calf", "calf", "calf", "mean", "mean", 
"mean", "mean", "mean", "mean", "mean", "mean", "mean", "mean", 
"mean"),
days=c(0,5,6,7,8,9,10,14,16,18,21),
percent.granulocytes=
c(adult.one=c(0,10,11,50,46,37,38,16,15,0,0),
adult.two=c(0,10,50,84,21,18,18,12,12,10,0),
calf=c(0,10,40,32,24,25,33,10,12,15,11),
mean=c(0,10,33,55,30,26,29,12,13,10,10)))


### using optim for anaplasma
####The proper way to do this would probably be to fit a nonparametric curve to the data of interest (ie anaplasma),
  ###calculate dt's and use a linear regression to fit parameters of the LV/situation.
  
#   
# start.pars<-c(alpha.t1=4,beta.t1=0.01,alpha.t2=80,beta.t2=30)
# State<-c(t1t=1,t2t=0)
# Time=seq(0,100,by=0.1)
# 
# test<-optim(par=start.pars,fn=lv.mse,
#          start.state=State,
#          test2.data=anaplasma.ab$deer.ab.values,
#          test1.data=reindeer.granulocytes$percent.granulocytes,
#          t2.time.data=anaplasma.ab$time,
#          t1.time.data=reindeer.granulocytes$days,
#          est.timepoints=seq(0,100,by=0.1))
# 


## Bluetongue

btv.df<-read.table(file.path(data.path,"Test response data","BTV.csv"),header=T,sep=",")
btv.df$AB[is.na(btv.df$AB)]<-0
dpi.seq<-seq(head(btv.df$DPI,1),tail(btv.df$DPI,1),by=0.5)
btv.loess.df<-data.frame(DPI=dpi.seq,
                         AB=predict(loess(btv.df$AB~btv.df$DPI,surface="direct"),dpi.seq)+0.29,  ###adding constant to avoid negative interpolation results.
                         PCR=1/predict(loess(btv.df$PCR~btv.df$DPI,surface="direct"),dpi.seq))
btv.loess.df$PCR[btv.loess.df$PCR<0]<-0

naDetFunBTV<-(function(Inf.levels=btv.loess.df,timescale=NULL){
  function(infection.time){
    sapply(infection.time,function(x){
      if(x>max(Inf.levels$DPI)){return(Inf.levels[which.max(Inf.levels$DPI>=x),"PCR"])}
      if(x<=max(Inf.levels$DPI)){Inf.levels[min(which(Inf.levels$DPI>=x)),"PCR"]}}
           )}}
           )()




abDetFunBTV<-(function(Inf.levels=btv.loess.df,timescale=NULL){
  function(infection.time){
    abt.response<-sapply(infection.time,function(x){
      if(x>max(Inf.levels$DPI)){return(Inf.levels[which.max(Inf.levels$DPI>=x),"AB"])}
      if(x<=max(Inf.levels$DPI)){Inf.levels[min(which(Inf.levels$DPI>=x)),"AB"]}}
                         )
    abt.response[infection.time==0]<-rlnorm(sum(infection.time==0),mean=log(mean(Inf.levels$AB[56])))
    return(abt.response)
  }}
           )()

inf.levels.btv <- btv.loess.df
names(inf.levels.btv)<-c("time","abt","nat")


##The subset of data after the peak of btv PCR values
btv.loess.df.tail<-btv.loess.df[btv.loess.df$DPI>btv.loess.df[which.max(btv.loess.df$PCR),"DPI"],]

btv.loess.df.50pc.dec<-btv.loess.df
btv.loess.df.50pc.dec[btv.loess.df$DPI>btv.loess.df[which.max(btv.loess.df$PCR),"DPI"],"PCR"]<-
  btv.loess.df.tail$PCR-(max(btv.loess.df$PCR)-btv.loess.df.tail$PCR)

naDetFunBTV.50pc<-(function(Inf.levels=btv.loess.df.50pc.dec,timescale=NULL){
  function(infection.time){
    sapply(infection.time,function(x){
      if(x>max(Inf.levels$DPI)){return(Inf.levels[which.max(Inf.levels$DPI>=x),"PCR"])}
      if(x<=max(Inf.levels$DPI)){Inf.levels[min(which(Inf.levels$DPI>=x)),"PCR"]}}
    )}}
)()





btv.loess.df.75pc.dec<-btv.loess.df
btv.loess.df.75pc.dec[btv.loess.df$DPI>btv.loess.df[which.max(btv.loess.df$PCR),"DPI"],"PCR"]<-
  btv.loess.df.tail$PCR-(max(btv.loess.df$PCR)-btv.loess.df.tail$PCR)*2

naDetFunBTV.75pc<-(function(Inf.levels=btv.loess.df.75pc.dec,timescale=NULL){
  function(infection.time){
    sapply(infection.time,function(x){
      if(x>max(Inf.levels$DPI)){return(Inf.levels[which.max(Inf.levels$DPI>=x),"PCR"])}
      if(x<=max(Inf.levels$DPI)){Inf.levels[min(which(Inf.levels$DPI>=x)),"PCR"]}}
    )}}
)()


btv.loess.df.50pc.higher<-btv.loess.df
btv.loess.df.50pc.higher[btv.loess.df$DPI>btv.loess.df[which.max(btv.loess.df$PCR),"DPI"],"PCR"]<-
  btv.loess.df.tail$PCR+(max(btv.loess.df$PCR)-btv.loess.df.tail$PCR)/2

naDetFunBTV.50pc.higher<-(function(Inf.levels=btv.loess.df.50pc.higher,timescale=NULL){
  function(infection.time){
    sapply(infection.time,function(x){
      if(x>max(Inf.levels$DPI)){return(Inf.levels[which.max(Inf.levels$DPI>=x),"PCR"])}
      if(x<=max(Inf.levels$DPI)){Inf.levels[min(which(Inf.levels$DPI>=x)),"PCR"]}}
    )}}
)()

btv.loess.df.75pc.higher<-btv.loess.df
btv.loess.df.75pc.higher[btv.loess.df$DPI>btv.loess.df[which.max(btv.loess.df$PCR),"DPI"],"PCR"]<-
  btv.loess.df.tail$PCR+(max(btv.loess.df$PCR)-btv.loess.df.tail$PCR)*0.75

naDetFunBTV.75pc.higher<-(function(Inf.levels=btv.loess.df.75pc.higher,timescale=NULL){
  function(infection.time){
    sapply(infection.time,function(x){
      if(x>max(Inf.levels$DPI)){return(Inf.levels[which.max(Inf.levels$DPI>=x),"PCR"])}
      if(x<=max(Inf.levels$DPI)){Inf.levels[min(which(Inf.levels$DPI>=x)),"PCR"]}}
    )}}
)()


###AB mod
dpi.seq.fine<-seq(head(btv.df$DPI,1),tail(btv.df$DPI,1),by=0.01)
btv.loess.df.fine<-data.frame(DPI=dpi.seq.fine,
                         AB=predict(loess(btv.df$AB~btv.df$DPI,surface="direct"),dpi.seq.fine)+0.29,  ###adding constant to avoid negative interpolation results.
                         PCR=1/predict(loess(btv.df$PCR~btv.df$DPI,surface="direct"),dpi.seq.fine))


comp<-0.25
temp<-round(seq(1,which.max(btv.loess.df.fine$AB),
                length.out=floor(nrow(btv.loess.df)*which.max(btv.loess.df$AB)/nrow(btv.loess.df)*comp)))
temp2<-round(seq(which.max(btv.loess.df.fine$AB)+1,nrow(btv.loess.df.fine),length.out=
  floor(nrow(btv.loess.df)*(1-which.max(btv.loess.df$AB)/nrow(btv.loess.df)*comp)+1)))
btv.loess.df.comp.25<-data.frame(AB=btv.loess.df.fine[c(temp,temp2),"AB"],DPI=btv.loess.df$DPI)

abDetFunBTV.25<-(function(Inf.levels=btv.loess.df.comp.25,timescale=NULL){
  function(infection.time){
    abt.response<-sapply(infection.time,function(x){
      if(x>max(Inf.levels$DPI)){return(Inf.levels[which.max(Inf.levels$DPI>=x),"AB"])}
      if(x<=max(Inf.levels$DPI)){Inf.levels[min(which(Inf.levels$DPI>=x)),"AB"]}}
    )
    abt.response[infection.time==0]<-rlnorm(sum(infection.time==0),mean=log(mean(Inf.levels$AB[56])))
    return(abt.response)
  }}
)()


comp<-0.5
temp<-round(seq(1,which.max(btv.loess.df.fine$AB),
          length.out=round(nrow(btv.loess.df)*which.max(btv.loess.df$AB)/nrow(btv.loess.df)*comp)))
temp2<-round(seq(which.max(btv.loess.df.fine$AB)+1,nrow(btv.loess.df.fine),length.out=
    floor(nrow(btv.loess.df)*(1-which.max(btv.loess.df$AB)/nrow(btv.loess.df)*comp))))
btv.loess.df.comp.50<-data.frame(AB=btv.loess.df.fine[c(temp,temp2),"AB"],DPI=btv.loess.df$DPI)

abDetFunBTV.50<-(function(Inf.levels=btv.loess.df.comp.50,timescale=NULL){
  function(infection.time){
    abt.response<-sapply(infection.time,function(x){
      if(x>max(Inf.levels$DPI)){return(Inf.levels[which.max(Inf.levels$DPI>=x),"AB"])}
      if(x<=max(Inf.levels$DPI)){Inf.levels[min(which(Inf.levels$DPI>=x)),"AB"]}}
    )
    abt.response[infection.time==0]<-rlnorm(sum(infection.time==0),mean=log(mean(Inf.levels$AB[56])))
    return(abt.response)
  }}
)()



comp<-0.75
temp<-round(seq(1,which.max(btv.loess.df.fine$AB),
                length.out=floor(nrow(btv.loess.df)*which.max(btv.loess.df$AB)/nrow(btv.loess.df)*comp)))
temp2<-round(seq(which.max(btv.loess.df.fine$AB)+1,nrow(btv.loess.df.fine),length.out=
  floor(nrow(btv.loess.df)*(1-which.max(btv.loess.df$AB)/nrow(btv.loess.df)*comp)+1)))
btv.loess.df.comp.75<-data.frame(AB=btv.loess.df.fine[c(temp,temp2),"AB"],DPI=btv.loess.df$DPI)

abDetFunBTV.75<-(function(Inf.levels=btv.loess.df.comp.75,timescale=NULL){
  function(infection.time){
    abt.response<-sapply(infection.time,function(x){
      if(x>max(Inf.levels$DPI)){return(Inf.levels[which.max(Inf.levels$DPI>=x),"AB"])}
      if(x<=max(Inf.levels$DPI)){Inf.levels[min(which(Inf.levels$DPI>=x)),"AB"]}}
    )
    abt.response[infection.time==0]<-rlnorm(sum(infection.time==0),mean=log(mean(Inf.levels$AB[56])))
    return(abt.response)
  }}
)()


comp<-1.25
temp<-round(seq(1,which.max(btv.loess.df.fine$AB),
                length.out=floor(nrow(btv.loess.df)*which.max(btv.loess.df$AB)/nrow(btv.loess.df)*comp)))
temp2<-round(seq(which.max(btv.loess.df.fine$AB)+1,nrow(btv.loess.df.fine),length.out=
  floor(nrow(btv.loess.df)*(1-which.max(btv.loess.df$AB)/nrow(btv.loess.df)*comp)+1)))
btv.loess.df.comp.125<-data.frame(AB=btv.loess.df.fine[c(temp,temp2),"AB"],DPI=btv.loess.df$DPI)

abDetFunBTV.125<-(function(Inf.levels=btv.loess.df.comp.125,timescale=NULL){
  function(infection.time){
    abt.response<-sapply(infection.time,function(x){
      if(x>max(Inf.levels$DPI)){return(Inf.levels[which.max(Inf.levels$DPI>=x),"AB"])}
      if(x<=max(Inf.levels$DPI)){Inf.levels[min(which(Inf.levels$DPI>=x)),"AB"]}}
    )
    abt.response[infection.time==0]<-rlnorm(sum(infection.time==0),mean=log(mean(Inf.levels$AB[56])))
    return(abt.response)
  }}
)()


###EHEC

library(reshape)
ehec.ab.df<-read.table(file.path(data.path,"Test response data","EHEC_ab.csv"),header=T,sep=",")
ehec.ab.df<-ehec.ab.df[1:4,1:10]
ehec.ab.df<-melt(ehec.ab.df,measure.vars=paste("C",1:9,sep=""),variable_name="cow")

ehec.shedding.df<-read.table(file.path(data.path,"Test response data","EHEC_shedding.csv"),header=T,sep=",")


dpi.seq<-seq(head(ehec.ab.df$DPI,1),tail(ehec.ab.df$DPI,1),by=0.5)
ehec.loess.df<-data.frame(DPI=dpi.seq,
                         AB=predict(loess(ehec.ab.df$value~ehec.ab.df$DPI,surface="direct"),dpi.seq)+0.1,  ###adding constant to avoid negative interpolation results.
                         Shedding=predict(loess(ehec.shedding.df$Shedding~ehec.shedding.df$DPI,surface="direct"),dpi.seq))

naDetFun.ehec<-(function(Inf.levels=ehec.loess.df,timescale=NULL){
  function(infection.time){
    sapply(infection.time,function(x){
      if(x>max(Inf.levels$DPI)){return(Inf.levels[which.max(Inf.levels$DPI>=x),"Shedding"])}
      if(x<=max(Inf.levels$DPI)){Inf.levels[min(which(Inf.levels$DPI>=x)),"Shedding"]}}
    )}}
)()




abDetFun.ehec<-(function(Inf.levels=ehec.loess.df,timescale=NULL){
  function(infection.time){
    abt.response<-sapply(infection.time,function(x){
      if(x>max(Inf.levels$DPI)){return(Inf.levels[which.max(Inf.levels$DPI>=x),"AB"])}
      if(x<=max(Inf.levels$DPI)){Inf.levels[min(which(Inf.levels$DPI>=x)),"AB"]}}
    )
    abt.response[infection.time==0]<-rlnorm(sum(infection.time==0),mean=log(mean(Inf.levels$AB[56])))
    return(abt.response)
  }}
)()

inf.levels.ehec <- ehec.loess.df
names(inf.levels.ehec)<-c("time","abt","nat")



###Measles

library(reshape)
measles.ab.df<-read.table(file.path(data.path,"Test response data","Binnendijk 2003 - Measles throat swab.csv"),header=T,sep=",") ### this is also PCR!!!
measles.pcr.df<-read.table(file.path(data.path,"Test response data","beard 2007 - measles pcr.csv"),header=T,sep=",")

measles.ab.df$Y.log<-log(measles.ab.df$Y)
dpi.seq<-seq(-10,30,by=0.1)
measles.ab.loess<-data.frame(DPI=dpi.seq,igM=
  exp(predict(smooth.spline(y=subset(measles.ab.df,X>-5&X<10)$Y.log,x=subset(measles.ab.df,X>-5&X<10)$X,df=7),x=dpi.seq)$y
))

plot(Y~X,data=measles.ab.df,log="y")
lines(igM~DPI,data=measles.ab.loess,xlim=c(-5,12),col="red")


measles.pcr.loess<-data.frame(DPI=dpi.seq,ct=
  predict(smooth.spline(measles.pcr.df$X,measles.pcr.df$Y,df=7),x=dpi.seq)$y)

plot(Y~X,data=measles.pcr.df)
lines(ct~DPI,data=measles.pcr.loess,col="red")

plot(measles.pcr.loess$ct[-(1:50)],measles.ab.loess$igM[-(1:50)],type="l",log="y")




influenza.ab.df<-read.table(file.path(data.path,"Test response data","hancioglu 2007 - influenza ab.csv"),header=T,sep=",")
influenza.pcr.df<-read.table(file.path(data.path,"Test response data","baccam 2006 influenza VL.csv"),header=T,sep=",")

influenza.ab.df$A.log<-log(influenza.ab.df$A)
dpi.seq<-seq(0,15,by=0.1)
influenza.ab.loess<-data.frame(DPI=dpi.seq,ab=
  exp(predict(smooth.spline(y=infuenza.ab.df$A.log,x=influenza.ab.df$Days,df=25),x=dpi.seq)$y))

plot(A~Days,data=influenza.ab.df,xlim=c(0,30),log="y")
lines(ab~DPI,data=influenza.ab.loess,col="red")


influenza.pcr.loess<-data.frame(DPI=dpi.seq,vl=
  predict(smooth.spline(y=infuenza.pcr.df$viral.titre,x=infuenza.pcr.df$days,df=25),x=dpi.seq)$y)

influenza.pcr.loess$vl[influenza.pcr.loess$vl<0]<-0

plot(viral.titre~days,data=infuenza.pcr.df,xlim=c(0,30))
lines(vl~DPI,data=influenza.pcr.loess,col="red")

plot(influenza.pcr.loess$vl,influenza.ab.loess$ab,type="l",log="y")