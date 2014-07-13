###Function by Gustaf Rydevik to generate simulated  population antibody
###and nucleic acid data,
### for some disease, including both infected and non-infected subjects.
### April-May 2011

LabdataGenerator<-function(n=100,start.time=0,end.time=NULL, incidence=NULL,SerumDecayFun,NaDecayFun,timeFun,errorFun,reinfect=FALSE,SD=log(1.1),FIX=FALSE,pop.actual=FALSE,...){
  ##Generate basic population
  obs.time<-end.time-start.time
  pop<-vector(length=n,mode="numeric")
  total.infected<-rpois(1,n*obs.time*incidence)
  if(FIX){total.infected<-n*obs.time*incidence}
  if(pop.actual){total.infected<-n}
  ## Distribute a number of infections in the population, proportional to
  ##incidence and time
  if(!FIX)infected.ndx<-sort(sample(1:n,total.infected,replace=T))
  if(FIX)infected.ndx<-1:total.infected 
  if(!reinfect)infected.ndx<-unique(infected.ndx)
  number.infections<-aggregate(infected.ndx,
                               list(n.infected=infected.ndx),length)[,2]
  infected.ndx<-sort(unique(infected.ndx))
  
  ##Decide when those infections occurred (currently uniformily over
  ##the time period)
  if(reinfect){infection.time<-sapply(number.infections,
                                       function(x)max(timeFun(n.infections=x,start.time=start.time,end.time=end.time,...)))}
  if(!reinfect){infection.time<-timeFun(n.infections=length(infected.ndx),start.time=start.time,end.time=end.time,...)}
  

  ##Generate output data frame
  antibody.data<-data.frame(index=infected.ndx,number.infections,
                            infection.time=infection.time,
                            antibody.levels=NA,
                            na.result=0)
  antibody.data<-rbind(antibody.data,data.frame(
                               index=Ndx<-which(!((1:n)%in%infected.ndx)),
                               number.infections=rep(0,length(Ndx)),
                               infection.time=rep(0,length(Ndx)),
                               antibody.levels=rep(0,length(Ndx)),
                               na.result=rep(0,length(Ndx))))
  

  ##Generate an antibody response, based on when the infection occcured. Counts time as time since infection.
  antibody.data$antibody.levels<-SerumDecayFun(infection.time=
                                               antibody.data$infection.time)#,
                                               #the line below can be added if we care about
                                               #previous infections
                                               #number.infections=antibody.data$number.infections)


  ##Generate presence or absence of nucleic acids
  antibody.data$na.result<-NaDecayFun(infection.time=
                                      antibody.data$infection.time)

##Add observational and individual variation.
  antibody.data[c("antibody.levels","na.result")]<-errorFun(antibody.data[c("antibody.levels","na.result")],standard.deviation=SD)

  ##Return end result
  return(antibody.data)
}
