#######Deterministic models for NA/antibody interaction
###using differential equations

###Parameters for the lottka-voltera diff. equations below
# 
# Pars<-c(alpha.na=30,beta.na=0.5,alpha.ab=5,beta.ab=0.1)
# State<-c(nat=1,abt=0)
# Time=seq(0,10,by=0.001)
# 
# 
# ###Lotke-Voltera equation based on Giles' idea.
# autolib(deSolve)
# ab.na.diffmod<-function(Time,State,Pars){
#   with(as.list(c(State,Pars)),{
#       dna=nat*(alpha.na-beta.na*abt)
#         dab=alpha.ab*nat-beta.ab*abt
#         return(list(c(dna,dab)))
#       })}
# inf.levels.lv <- as.data.frame(ode(func = ab.na.diffmod, y = State, parms = Pars, times = Time))
# inf.levels.lv$abt<-2^(inf.levels.lv$abt/12)
# 
# ###Deterministic nucleic acid function, pulling values from Giles' diff model
# naDetFun<-(function(Inf.levels=inf.levels.lv,timescale=NULL){
#   function(infection.time){
#    sapply(infection.time,function(x){
#      if(x>max(Inf.levels$time)){return(Inf.levels[which.max(Inf.levels$time>=x),"nat"])}
#      if(x<=max(Inf.levels$time)){Inf.levels[min(which(Inf.levels$time>=x)),"nat"]}}
#         )}}
#    )()
# 
# 
# 
# ###Deterministic antibody function, pulling values from Giles' diff model
# abDetFun<-(function(Inf.levels=inf.levels.lv,timescale=NULL){
#   function(infection.time){
#    abt.response<-sapply(infection.time,function(x){
#      if(x>max(Inf.levels$time)){return(Inf.levels[which.max(Inf.levels$time>=x),"abt"])}
#      if(x<=max(Inf.levels$time)){Inf.levels[min(which(Inf.levels$time>=x)),"abt"]}}
#           )
#    abt.response[infection.time==0]<-rlnorm(sum(infection.time==0),mean=log(mean(Inf.levels$abt[100])))
#    return(abt.response)
# 	}}
#    )()


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
 

## Bluetongue
btv.df<-read.table(file.path(data.path,"Test response data","BTV.csv"),header=T,sep=",")
btv.df$AB[is.na(btv.df$AB)]<-0
dpi.seq<-seq(1,120,by=0.5)
#dropping obs 1-6 in ab below because they're all zero and skew the fit.
btv.spline.df<-data.frame(DPI=dpi.seq,
                          AB=predict(smooth.spline(c(0,btv.df$DPI),c(0.5,btv.df$AB),df=7),x=dpi.seq)$y,  ###adding constant to avoid negative interpolation results.
                          PCR=scale(predict(smooth.spline(btv.df$DPI,2^(-btv.df$PCR),df=7),x=dpi.seq)$y,center=F))
btv.spline.df$PCR[btv.spline.df$PCR<0]<-min(btv.spline.df$PCR[btv.spline.df$PCR>0])
plot(btv.df$DPI,btv.df$AB,xlim=c(0,120))
lines(btv.spline.df$DPI,btv.spline.df$AB)
btv.spline.df$AB

naDetFunBTV<-(function(Inf.levels=btv.spline.df,timescale=NULL){
  function(infection.time){
    sapply(infection.time,function(x){
      if(x>max(Inf.levels$DPI)){return(Inf.levels[which.max(Inf.levels$DPI>=x),"PCR"])}
      if(x<=max(Inf.levels$DPI)){Inf.levels[min(which(Inf.levels$DPI>=x)),"PCR"]}}
    )}}
)() 




abDetFunBTV<-(function(Inf.levels=btv.spline.df,timescale=NULL){
  function(infection.time){
    abt.response<-sapply(infection.time,function(x){
      if(x>max(Inf.levels$DPI)){return(Inf.levels[which.max(Inf.levels$DPI>=x),"AB"])}
      if(x<=max(Inf.levels$DPI)){Inf.levels[min(which(Inf.levels$DPI>=x)),"AB"]}}
    )
    #abt.response[infection.time==0]<-rlnorm(sum(infection.time==0),mean=log(mean(Inf.levels$AB[56])))
    return(abt.response)
  }}
)()

inf.levels.btv <- btv.spline.df
names(inf.levels.btv)<-c("time","abt","nat")


##The subset of data after the peak of btv PCR values
btv.spline.df.tail<-btv.spline.df[btv.spline.df$DPI>btv.spline.df[which.max(btv.spline.df$PCR),"DPI"],]

btv.spline.df.50pc.dec<-btv.spline.df
btv.spline.df.50pc.dec[btv.spline.df$DPI>btv.spline.df[which.max(btv.spline.df$PCR),"DPI"],"PCR"]<-
  btv.spline.df.tail$PCR-(max(btv.spline.df$PCR)-btv.spline.df.tail$PCR)

naDetFunBTV.50pc<-(function(Inf.levels=btv.spline.df.50pc.dec,timescale=NULL){
  function(infection.time){
    sapply(infection.time,function(x){
      if(x>max(Inf.levels$DPI)){return(Inf.levels[which.max(Inf.levels$DPI>=x),"PCR"])}
      if(x<=max(Inf.levels$DPI)){Inf.levels[min(which(Inf.levels$DPI>=x)),"PCR"]}}
    )}}
)()





btv.spline.df.75pc.dec<-btv.spline.df
btv.spline.df.75pc.dec[btv.spline.df$DPI>btv.spline.df[which.max(btv.spline.df$PCR),"DPI"],"PCR"]<-
  btv.spline.df.tail$PCR-(max(btv.spline.df$PCR)-btv.spline.df.tail$PCR)*2

naDetFunBTV.75pc<-(function(Inf.levels=btv.spline.df.75pc.dec,timescale=NULL){
  function(infection.time){
    sapply(infection.time,function(x){
      if(x>max(Inf.levels$DPI)){return(Inf.levels[which.max(Inf.levels$DPI>=x),"PCR"])}
      if(x<=max(Inf.levels$DPI)){Inf.levels[min(which(Inf.levels$DPI>=x)),"PCR"]}}
    )}}
)()


btv.spline.df.50pc.higher<-btv.spline.df
btv.spline.df.50pc.higher[btv.spline.df$DPI>btv.spline.df[which.max(btv.spline.df$PCR),"DPI"],"PCR"]<-
  btv.spline.df.tail$PCR+(max(btv.spline.df$PCR)-btv.spline.df.tail$PCR)/2

naDetFunBTV.50pc.higher<-(function(Inf.levels=btv.spline.df.50pc.higher,timescale=NULL){
  function(infection.time){
    sapply(infection.time,function(x){
      if(x>max(Inf.levels$DPI)){return(Inf.levels[which.max(Inf.levels$DPI>=x),"PCR"])}
      if(x<=max(Inf.levels$DPI)){Inf.levels[min(which(Inf.levels$DPI>=x)),"PCR"]}}
    )}}
)()

btv.spline.df.75pc.higher<-btv.spline.df
btv.spline.df.75pc.higher[btv.spline.df$DPI>btv.spline.df[which.max(btv.spline.df$PCR),"DPI"],"PCR"]<-
  btv.spline.df.tail$PCR+(max(btv.spline.df$PCR)-btv.spline.df.tail$PCR)*0.75

naDetFunBTV.75pc.higher<-(function(Inf.levels=btv.spline.df.75pc.higher,timescale=NULL){
  function(infection.time){
    sapply(infection.time,function(x){
      if(x>max(Inf.levels$DPI)){return(Inf.levels[which.max(Inf.levels$DPI>=x),"PCR"])}
      if(x<=max(Inf.levels$DPI)){Inf.levels[min(which(Inf.levels$DPI>=x)),"PCR"]}}
    )}}
)()


###AB mod
dpi.seq.fine<-seq(head(btv.df$DPI,1),tail(btv.df$DPI,1),by=0.01)
#dropping obs 1-6 in ab below because they're all zero and skew the fit.
btv.spline.df.fine<-data.frame(DPI=dpi.seq.fine,
                               AB=predict(smooth.spline(btv.df$DPI[-(1:6)],btv.df$AB[-(1:6)],df=7),x=dpi.seq.fine)$y, 
                               PCR=1/predict(smooth.spline(btv.df$DPI,btv.df$PCR,df=7),x=dpi.seq.fine)$y)
btv.spline.df.fine$AB[btv.spline.df.fine$AB<=0]<-min(subset(btv.spline.df.fine,AB>0)$AB)

comp<-0.25
temp<-round(seq(1,which.max(btv.spline.df.fine$AB),
                length.out=floor(nrow(btv.spline.df)*which.max(btv.spline.df$AB)/nrow(btv.spline.df)*comp)))
temp2<-round(seq(which.max(btv.spline.df.fine$AB)+1,nrow(btv.spline.df.fine),length.out=
  floor(nrow(btv.spline.df)*(1-which.max(btv.spline.df$AB)/nrow(btv.spline.df)*comp)+1)))
btv.spline.df.comp.25<-data.frame(AB=btv.spline.df.fine[c(temp,temp2),"AB"],DPI=btv.spline.df$DPI)

###Trim neg values
btv.spline.df.comp.25$AB[btv.spline.df.comp.25$AB<=0]<-min(subset(btv.spline.df.comp.25,AB>0)$AB)

abDetFunBTV.25<-(function(Inf.levels=btv.spline.df.comp.25,timescale=NULL){
  function(infection.time){
    abt.response<-sapply(infection.time,function(x){
      if(x>max(Inf.levels$DPI)){return(Inf.levels[which.max(Inf.levels$DPI>=x),"AB"])}
      if(x<=max(Inf.levels$DPI)){Inf.levels[min(which(Inf.levels$DPI>=x)),"AB"]}}
    )
    #abt.response[infection.time==0]<-rlnorm(sum(infection.time==0),mean=log(mean(Inf.levels$AB[56])))
    return(abt.response)
  }}
)()


comp<-0.5
temp<-round(seq(1,which.max(btv.spline.df.fine$AB),
                length.out=round(nrow(btv.spline.df)*which.max(btv.spline.df$AB)/nrow(btv.spline.df)*comp)))
temp2<-round(seq(which.max(btv.spline.df.fine$AB)+1,nrow(btv.spline.df.fine),length.out=
  floor(nrow(btv.spline.df)*(1-which.max(btv.spline.df$AB)/nrow(btv.spline.df)*comp))))
btv.spline.df.comp.50<-data.frame(AB=btv.spline.df.fine[c(temp,temp2),"AB"],DPI=btv.spline.df$DPI)

btv.spline.df.comp.50$AB[btv.spline.df.comp.50$AB<=0]<-min(subset(btv.spline.df.comp.50,AB>0)$AB)

abDetFunBTV.50<-(function(Inf.levels=btv.spline.df.comp.50,timescale=NULL){
  function(infection.time){
    abt.response<-sapply(infection.time,function(x){
      if(x>max(Inf.levels$DPI)){return(Inf.levels[which.max(Inf.levels$DPI>=x),"AB"])}
      if(x<=max(Inf.levels$DPI)){Inf.levels[min(which(Inf.levels$DPI>=x)),"AB"]}}
    )
    #abt.response[infection.time==0]<-rlnorm(sum(infection.time==0),mean=log(mean(Inf.levels$AB[56])))
    return(abt.response)
  }}
)()



comp<-0.75
temp<-round(seq(1,which.max(btv.spline.df.fine$AB),
                length.out=floor(nrow(btv.spline.df)*which.max(btv.spline.df$AB)/nrow(btv.spline.df)*comp)))
temp2<-round(seq(which.max(btv.spline.df.fine$AB)+1,nrow(btv.spline.df.fine),length.out=
  floor(nrow(btv.spline.df)*(1-which.max(btv.spline.df$AB)/nrow(btv.spline.df)*comp)+1)))
btv.spline.df.comp.75<-data.frame(AB=btv.spline.df.fine[c(temp,temp2),"AB"],DPI=btv.spline.df$DPI)

btv.spline.df.comp.75$AB[btv.spline.df.comp.75$AB<=0]<-min(subset(btv.spline.df.comp.75,AB>0)$AB)

abDetFunBTV.75<-(function(Inf.levels=btv.spline.df.comp.75,timescale=NULL){
  function(infection.time){
    abt.response<-sapply(infection.time,function(x){
      if(x>max(Inf.levels$DPI)){return(Inf.levels[which.max(Inf.levels$DPI>=x),"AB"])}
      if(x<=max(Inf.levels$DPI)){Inf.levels[min(which(Inf.levels$DPI>=x)),"AB"]}}
    )
    #abt.response[infection.time==0]<-rlnorm(sum(infection.time==0),mean=log(mean(Inf.levels$AB[56])))
    return(abt.response)
  }}
)()


comp<-1.25
temp<-round(seq(1,which.max(btv.spline.df.fine$AB),
                length.out=floor(nrow(btv.spline.df)*which.max(btv.spline.df$AB)/nrow(btv.spline.df)*comp)))
temp2<-round(seq(which.max(btv.spline.df.fine$AB)+1,nrow(btv.spline.df.fine),length.out=
  floor(nrow(btv.spline.df)*(1-which.max(btv.spline.df$AB)/nrow(btv.spline.df)*comp)+1)))
btv.spline.df.comp.125<-data.frame(AB=btv.spline.df.fine[c(temp,temp2),"AB"],DPI=btv.spline.df$DPI)

btv.spline.df.comp.125$AB[btv.spline.df.comp.125$AB<=0]<-min(subset(btv.spline.df.comp.125,AB>0)$AB)

abDetFunBTV.125<-(function(Inf.levels=btv.spline.df.comp.125,timescale=NULL){
  function(infection.time){
    abt.response<-sapply(infection.time,function(x){
      if(x>max(Inf.levels$DPI)){return(Inf.levels[which.max(Inf.levels$DPI>=x),"AB"])}
      if(x<=max(Inf.levels$DPI)){Inf.levels[min(which(Inf.levels$DPI>=x)),"AB"]}}
    )
    #abt.response[infection.time==0]<-rlnorm(sum(infection.time==0),mean=log(mean(Inf.levels$AB[56])))
    return(abt.response)
  }}
)()

###Measles

library(reshape)
measles.ab.df<-read.table(file.path(data.path,"Test response data","Binnendijk 2003 - Measles throat swab.csv"),header=T,sep=",") ### this is also PCR!!!
measles.pcr.df<-read.table(file.path(data.path,"Test response data","beard 2007 - measles pcr.csv"),header=T,sep=",")
measles.pcr.df$Y<- 1/(-measles.pcr.df$Y)  ### what's the best unit of measurement here?
measles.ab.df$Y.log<-log(measles.ab.df$Y)
dpi.seq<-seq(-10,100,by=0.1)
measles.ab.spline<-data.frame(DPI=dpi.seq,AB=
  exp(predict(smooth.spline(y=subset(measles.ab.df,X>-5&X<11)$Y.log,x=subset(measles.ab.df,X>-5&X<11)$X,df=7),x=dpi.seq)$y
))



measles.pcr.spline<-data.frame(DPI=dpi.seq,PCR=
  predict(smooth.spline(measles.pcr.df$X,measles.pcr.df$Y,df=7),x=dpi.seq)$y)


measles.spline.df<-merge(measles.ab.spline,measles.pcr.spline)


abDetFunMeasles<-(function(Inf.levels=measles.ab.spline,timescale=NULL){
  function(infection.time){
    abt.response<-sapply(infection.time,function(x){
      if(x>max(Inf.levels$DPI)){return(Inf.levels[which.max(Inf.levels$DPI>=x),"AB"])}
      if(x<=max(Inf.levels$DPI)){Inf.levels[min(which(Inf.levels$DPI>=x)),"AB"]}}
    )
    #abt.response[infection.time==0]<-rlnorm(sum(infection.time==0),mean=log(mean(Inf.levels$AB[56])))
    return(abt.response)
  }}
)()

naDetFunMeasles<-(function(Inf.levels=measles.pcr.spline,timescale=NULL){
  function(infection.time){
    sapply(infection.time,function(x){
      if(x>max(Inf.levels$DPI)){return(Inf.levels[which.max(Inf.levels$DPI>=x),"PCR"])}
      if(x<=max(Inf.levels$DPI)){Inf.levels[min(which(Inf.levels$DPI>=x)),"PCR"]}}
    )}}
)()



######
##### Influensa
#####

influenza.ab.df<-read.table(file.path(data.path,"Test response data","hancioglu 2007 - influenza ab.csv"),header=T,sep=",")
influenza.pcr.df<-read.table(file.path(data.path,"Test response data","baccam 2006 influenza VL.csv"),header=T,sep=",")

influenza.ab.df$A.log<-log(influenza.ab.df$A)
dpi.seq<-seq(0,15,by=0.1)
influenza.ab.spline<-data.frame(DPI=dpi.seq,ab=
  exp(predict(smooth.spline(y=influenza.ab.df$A.log,x=influenza.ab.df$Days,df=25),x=dpi.seq)$y))

#plot(A~Days,data=influenza.ab.df,xlim=c(0,30),log="y")
#lines(ab~DPI,data=influenza.ab.spline,col="red")


influenza.pcr.spline<-data.frame(DPI=dpi.seq,vl=
  predict(smooth.spline(y=influenza.pcr.df$viral.titre,x=influenza.pcr.df$days,df=25),x=dpi.seq)$y)

influenza.pcr.spline$vl[influenza.pcr.spline$vl<0]<-0

#plot(viral.titre~days,data=influenza.pcr.df,xlim=c(0,30))
#lines(vl~DPI,data=influenza.pcr.spline,col="red")

#plot(influenza.pcr.spline$vl,influenza.ab.spline$ab,type="l",log="y")

##########
##########



####Pertusssis
#file.rename("Teunis etal 2002 Pertussis.csv","Teunis etal 2002 Pertussis_bad.csv")
#temp<-scan(file.path(data.path,"Test response data","Teunis etal 2002 Pertussis_bad.csv"),what="char")
#temp[-1]<-sapply(strsplit(temp[-1],"\\."),function(x)paste(x[1],".",x[2],",",x[3],".",x[4],sep=""))
#write(temp,file=file.path(data.path,"Test response data","Teunis etal 2002 Pertussis.csv"))
pertussis.ab.df<-read.table(file.path(data.path,"Test response data","Teunis etal 2002 Pertussis.csv"),sep=",",header=T)
pertussis.ab.df<-as.data.frame(10^(pertussis.ab.df))
names(pertussis.ab.df)<-c("days","EU_per_ml")
#pertussis.ab.df<-read.table(file.path(data.path,"Test response data","Hallander_etal_AB_2009.csv"),header=T,sep=",")
pertussis.pcr.df<-read.table(file.path(data.path,"Test response data","Bidet etal 2008 pertussis antibiotics pcr.csv"),header=T,sep=",")
pertussis.pcr.df$patient<-cumsum(c(1,diff(pertussis.pcr.df$X)<1))
pertussis.pcr.df<-subset(pertussis.pcr.df,X<=17) ###last X value is an outlier
autolib(lme4)
pertussis.pcr.df.round<-round(pertussis.pcr.df)
pertussis.pcr.df.round$centered<-unlist((with(pertussis.pcr.df.round,by(1/2^Y,list(X),function(x)x-mean(x))))) ##creating a centered PCR result
pertussis.pcr.sd<-exp(mean(
  c(by(log(1/pertussis.pcr.df.round$Y),pertussis.pcr.df.round$X,FUN=sd))
  ,na.rm=T)
                      ) ##~sd==14%
pertussis.pcr.predicted<-fixef(lmer(Y~X+I(X^2)+(1|patient),data=pertussis.pcr.df))%*%
  t(matrix(c(rep(1,31),seq(0,15,length.out=31),seq(0,15,length.out=31)^2),ncol=3))


#pertussis.ab.df$Y.log<-log(pertussis.ab.df$EU_per_ml)
dpi.seq<-seq(0,200,by=0.5)
pertussis.ab.spline<-data.frame(DPI=dpi.seq,ab=
  (predict(smooth.spline(y=pertussis.ab.df$EU_per_ml,x=pertussis.ab.df$days,df=20),x=dpi.seq)$y))

plot(EU_per_ml~days,data=pertussis.ab.df,xlim=c(0,150))
lines(ab~DPI,data=pertussis.ab.spline,col="red")


pertussis.pcr.spline<-data.frame(DPI=dpi.seq,bl=
  predict(smooth.spline(y=1/pertussis.pcr.df$Y,x=pertussis.pcr.df$X,df=8),x=dpi.seq)$y)
pertussis.pcr.spline$bl[pertussis.pcr.spline$bl<.02]<-0.02*head(pertussis.pcr.spline[pertussis.pcr.spline$bl<.02,"DPI"],1)/(pertussis.pcr.spline[pertussis.pcr.spline$bl<.02,"DPI"])
plot(1/Y~X,data=pertussis.pcr.df,xlim=c(0,150),ylim=c(0,0.1))
lines(bl~DPI,data=pertussis.pcr.spline,col="red")




pertussis.spline.df<-merge(pertussis.ab.spline,pertussis.pcr.spline)


abDetFunPertussis<-(function(Inf.levels=pertussis.ab.spline,timescale=NULL){
  function(infection.time){
    abt.response<-sapply(infection.time,function(x){
      if(x>max(Inf.levels$DPI)){return(Inf.levels[which.max(Inf.levels$DPI>=x),"ab"])}
      if(x<=max(Inf.levels$DPI)){Inf.levels[min(which(Inf.levels$DPI>=x)),"ab"]}}
    )
    #abt.response[infection.time==0]<-rlnorm(sum(infection.time==0),mean=Inf.levels$ab[56])
    return(abt.response)
  }}
)()

naDetFunPertussis<-(function(Inf.levels=pertussis.pcr.spline,timescale=NULL){
  function(infection.time){
    sapply(infection.time,function(x){
      if(x>max(Inf.levels$DPI)){return(Inf.levels[which.max(Inf.levels$DPI>=x),"bl"])}
      if(x<=max(Inf.levels$DPI)){Inf.levels[min(which(Inf.levels$DPI>=x)),"bl"]}}
    )}}
)()


inf.levels.pertussis <- pertussis.spline.df
names(inf.levels.pertussis)<-c("time","abt","nat")
