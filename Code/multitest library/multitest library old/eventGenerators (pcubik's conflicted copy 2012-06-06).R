#####
##Various functions for generating a set of 
# timed events according to different distributions.
#


#### a uniform event generator, simulating a fully stable, endemic infection.
stableTimeFun<-(function(){function(n.infections,start.time,end.time){
                         runif(n.infections,start.time,end.time)}})()   

#### a increasing event generator, simulating an increasing, endemic infection.
increasingTimeFun<-(function(){function(n.infections,start.time,end.time){
  autolib(triangle)
                         rtriangle(n.infections,a=start.time,b=end.time,c=end.time)}})()   

#### a decreasing event generator, simulating an increasing, endemic infection.
decreasingTimeFun<-(function(){function(n.infections,start.time,end.time){
  autolib(triangle)
                         rtriangle(n.infections,a=start.time,b=end.time,c=start.time)}})()   


#### an event generator with all events in the beginning
earlyTimeFun<-(function(){function(n.infections,start.time,end.time){
  start.time}})()   

#### an event generator with all events when 10% of the time remains before the end
lateTimeFun<-(function(){function(n.infections,start.time,end.time){
  start.time+(end.time-start.time)*0.9}})()   

expTimeFun<-(function(){function(n.infections,start.time=0,end.time,slope=2){
  time.length<-end.time-start.time
  scale.factor<-time.length
  dist.unif<-runif(n.infections,0,1)
  dist.exp<-log((exp(slope)-1)*(dist.unif-1/(1-exp(slope))))/slope
  dist.exp<-dist.exp*scale.factor+start.time
  return(dist.exp)
  }})()



invexpTimeFun<-(function(){function(n.infections,start.time=0,end.time,slope=2){
  time.length<-end.time-start.time
  scale.factor<-time.length
  dist.unif<-runif(n.infections,0,1)
  dist.exp<-log((exp(slope)-1)*(dist.unif-1/(1-exp(slope))))/slope
  dist.exp<-(1-dist.exp)*scale.factor
  dist.exp<-dist.exp+start.time
  return(dist.exp)
	}})()