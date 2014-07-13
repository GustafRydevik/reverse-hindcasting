#####
##Various functions for generating a set of 
# timed events according to different distributions.
#
###When do we need to define these as functions returning functions?

#### a uniform event generator, simulating a fully stable, endemic infection.
stableTimeFun<-(function(){function(n.infections,start.time,end.time,...){
                         runif(n.infections,start.time,end.time)}})()   

##Intervention breakpoint generator
intervention.at7<-(function(){function(n.infections,start.time,end.time,...){
  runif(n.infections,start.time,min(end.time,7))}})()

introduction.at7<-(function(){function(n.infections,start.time,end.time,...){
  runif(n.infections,7,end.time)}})()

##epidemic introduction generator; with a bit of low-level noise
epidemic.introduction<-(function(){function(n.infections,start.time,end.time,epidemic.start=6,slope=3){
  time.length<-end.time-epidemic.start
  scale.factor<-time.length
  dist.unif<-runif(n.infections,0,1)
  dist.exp<-log((exp(slope)-1)*(dist.unif-1/(1-exp(slope))))/slope
  dist.exp<-dist.exp*scale.factor+epidemic.start
  #low.level.noise<-rpois(1,n.infections/100)
  #dist.exp[sample(seq_along(dist.exp), low.level.noise)]<-runif(low.level.noise,start.time,end.time)
  return(dist.exp)
  }})()




#### a increasing event generator, simulating an increasing, endemic infection.
increasingTimeFun<-(function(){function(n.infections,start.time,end.time){
  autolib(triangle)
      rtriangle(n.infections,a=start.time,b=end.time,c=end.time)}})()   

#### a decreasing event generator, simulating an decreasing, endemic infection.
decreasingTimeFun<-(function(){function(n.infections,start.time,end.time){
  autolib(triangle)
      rtriangle(n.infections,a=start.time,b=end.time,c=start.time)}})()   


#### an event generator with all events in the beginning
earlyTimeFun<-(function(){function(n.infections,start.time,end.time){
  rep(start.time,n.infections)}})()   

#### an event generator with all events when 10% of the time remains before the end
lateTimeFun<-(function(){function(n.infections,start.time,end.time){
  start.time+(end.time-start.time)*0.9}})()   

expTimeFun<-(function(){function(n.infections,start.time=0,end.time,slope=4){
  time.length<-end.time-start.time
  scale.factor<-time.length
  dist.unif<-runif(n.infections,0,1)
  dist.exp<-log((exp(slope)-1)*(dist.unif-1/(1-exp(slope))))/slope
  dist.exp<-dist.exp*scale.factor+start.time
  return(dist.exp)
  }})()



invexpTimeFun<-(function(){function(n.infections,start.time=0,end.time,slope=4){
  time.length<-end.time-start.time
  scale.factor<-time.length
  dist.unif<-runif(n.infections,0,1)
  dist.exp<-log((exp(slope)-1)*(dist.unif-1/(1-exp(slope))))/slope
  dist.exp<-(1-dist.exp)*scale.factor
  dist.exp<-dist.exp+start.time
  return(dist.exp)
	}})()

###Truncated exponential distribution
trunc.dexp<-function(n.infections,start.time=0,end.time=1,epidemic.start=0,slope){
  obs.time<-end.time-start.time
  above.any<-T
  results<-rep(epidemic.start+1,n.infections)
  above.ndx<-1:n.infections
  while(above.any){
    results[above.ndx]<-rexp(length(above.ndx),slope)
    above.ndx<-which(results>epidemic.start)
    above.any<-any(results>epidemic.start)}
  return(results)
}

######


pertussis.wi.df<-read.table(file.path(data.path,"Test response data","Wisconsin_pertussis_labcases.csv"),
                                 header=T,sep=",")
pertussis.wi.df<-data.frame(pertussis.wi.df[1:which(diff(pertussis.wi.df[,1])<0),],
                            conf.cases=pertussis.wi.df[-(1:which(diff(pertussis.wi.df[,1])<0)),2])
pertussis.wi.df<-round(pertussis.wi.df)

##with spline, interpolates the observed incidence curve, roughly guesstimating the underlying force of infection

pertussis.interpolated<-as.data.frame(predict(smooth.spline(1:(max(pertussis.wi.df$two.week.block)*14),rep(pertussis.wi.df$No.cases,each=14)/14,df=nrow(pertussis.wi.df)*2)))
pertussis.interpolated$y[pertussis.interpolated$y<0]<-0
pertussis.part.cases.ndx<-c(by(1:(nrow(pertussis.interpolated)-1),
                               factor(cumsum(diff(round( ###finding the groups where round error adds up to one
                                 cumsum(pertussis.interpolated$y-floor(pertussis.interpolated$y)### summing up rounding error
                                 ))))),function(x)sample(x,1))) ### assigning an extra case randomly when rounding error>1
pertussis.interpolated.round<-round(pertussis.interpolated)
pertussis.interpolated.round$y[pertussis.part.cases.ndx]<-pertussis.interpolated.round$y[pertussis.part.cases.ndx]+1


pertussis.wi.Fun<-with(list(pertussis.interpolated.round,pertussis.interpolated),(function(){function(n.infections=55,start.time=1,Actual=F,
                                                                  end.time=120,spline=FALSE,FIX=T,
                                                                  only.confirmed=FALSE,...){
  round.sample<-pertussis.interpolated.round[sample(rep(start.time:end.time,times=pertussis.interpolated.round$y[start.time:end.time]),size=n.infections,replace=T),"x"]
  round.sample<-c(by(round.sample,round.sample,FUN=length))
  round.sample<-data.frame(x=as.numeric(names(round.sample)),y=round.sample)[order(as.numeric(names(round.sample))),]
  
  mean.sample<-sample(start.time:end.time,prob=pertussis.interpolated$y[start.time:end.time],size=n.infections,replace=T)
  mean.sample.aggregate<-c(by(mean.sample,mean.sample,FUN=length))
  mean.sample.aggregate<-data.frame(x=as.numeric(names(mean.sample.aggregate)),y=mean.sample.aggregate)[order(as.numeric(names(mean.sample.aggregate))),]
  infection.times<-mean.sample
    if(Actual) infection.times<-with(subset(pertussis.interpolated.round,x%in%c(start.time:end.time)),rep(x,times=y))
  infection.times<-end.time-infection.times
  
  return(infection.times)
}})()
)



###BTV epidemic simulator

btv.2008.df<-read.table(file.path(data.path,"Test response data","btv_epireport_typepos.csv"),
                            header=T,sep=",")
btv.2008.df<-round(btv.2008.df)

##with spline, interpolates the observed incidence curve, roughly guesstimating the underlying force of infection
btv.2008.padded<-data.frame(Day.since.19th.of.Sept=seq(1:max(btv.2008.df$Day.since.19th.of.Sept)),Cases=0)
btv.2008.padded[btv.2008.padded[,1]%in%c(btv.2008.df$Day.since.19th.of.Sept),2]<-btv.2008.df[,"Cases"]
btv.2008.interpolated<-as.data.frame(predict(smooth.spline(btv.2008.padded$Day.since.19th.of.Sept,btv.2008.padded$Cases,df=nrow(btv.2008.df)/4),newdata=data.frame(Date.since.19th.of.Sept=1:176)))
btv.2008.interpolated$y[btv.2008.interpolated$y<0]<-0
plot(btv.2008.interpolated)


btv.2008.Fun<-with(list(btv.2008.interpolated),(function(){function(n.infections=55,start.time=1,
                                                                                                      end.time=50,spline=FALSE,Actual=F,
                                                                                                      only.confirmed=FALSE,Time.since.Sampling=T,...){
  mean.sample<-sample(start.time:end.time,prob=btv.2008.interpolated$y[start.time:end.time],size=n.infections,replace=T)
  mean.sample.aggregate<-c(by(mean.sample,mean.sample,FUN=length))
  mean.sample.aggregate<-data.frame(x=as.numeric(names(mean.sample.aggregate)),y=mean.sample.aggregate)[order(as.numeric(names(mean.sample.aggregate))),]
  infection.times<-mean.sample
  if(Actual) infection.times<-with(subset(btv.2008.padded,Day.since.19th.of.Sept%in%c(start.time:end.time)),rep(Day.since.19th.of.Sept,times=Cases))
  if(Time.since.Sampling){infection.times<-end.time-infection.times}
  return(infection.times)
}})()
)