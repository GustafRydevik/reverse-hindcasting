##Author: Gustaf Rydevik
##Purpose: Proof of concept for reverse hindcasting

###Setting up basic parameters

data.path<-file.path("./Original Data")
clean.data.path<-file.path("./Clean Data")
plates.path<-file.path(data.path,"Plates")
script.path<-"./ScriptFiles"

####### Loading libraries #######
lapply(dir(file.path(script.path,"utilities"),full.names=T),source)
My.device<-Gen.device("png",res=200)

autolib(lattice)
autolib(ggplot2)
autolib(reshape2)
autolib(plyr)


#Code


### Parametrize an indicator curve


###Parameters for the lottka-voltera diff. equations below

Pars<-c(alpha.na=0.3,beta.na=0.17,alpha.ab=0.09,beta.ab=0.06)
State<-c(nat=2.65,abt=0.04)


##Attempt to fit parameters to an actual BTV curve
btv.df<-read.table(file.path(data.path,"BTV.csv"),header=T,sep=",")
btv.df$AB[is.na(btv.df$AB)]<-0
Time<-dpi.seq<-seq(1,120,by=0.5)
#dropping obs 1-6 in ab below because they're all zero and skew the fit.
btv.spline.df<-data.frame(DPI=dpi.seq,
                          AB=predict(smooth.spline(c(0,btv.df$DPI),c(0.5,btv.df$AB),df=7),x=dpi.seq)$y,  ###adding constant to avoid negative interpolation results.
                          PCR=1/predict(smooth.spline(btv.df$DPI,btv.df$PCR,df=7),x=dpi.seq)$y)

btv.spline.df$PCR<-btv.spline.df$PCR*100
###Lotke-Voltera equation based on Giles' idea.
autolib(deSolve)
ab.na.diffmod<-function(Time,State,Pars){
  with(as.list(c(State,Pars)),{
    dna=nat*(alpha.na-beta.na*abt)
    dab=alpha.ab*nat-beta.ab*abt
    return(list(c(dna,dab)))
  })}

t0<-proc.time()
optim.pars<-optim(Pars,
  fn=function(Pars){ode.result<-ode(func = ab.na.diffmod, y = State, parms = Pars, times = Time)
                    distance<-(ode.result[,"nat"]-btv.spline.df$PCR)^2+(ode.result[,"abt"]-btv.spline.df$AB)^2
                    distance<-sqrt(mean(distance))
                    return(distance)})

t1<-proc.time()

ode.result<-ode(func = ab.na.diffmod, y = State, parms = optim.pars$par, times = Time)
plot(PCR~AB,data=btv.spline.df)
lines(PCR~AB,data=btv.spline.df)
lines(nat~abt,data=ode.result)
plot(nat~abt,data=ode.result,type="l",ylim=c(0,7),xlim=c(0,5))
points(PCR~AB,data=btv.spline.df)


###Create data set, using a known epidemic curve 


###BTV epidemic simulator

btv.2008.df<-read.table(file.path(data.path,"btv_epireport_typepos.csv"),
                        header=T,sep=",")
btv.2008.df<-round(btv.2008.df)

##with spline, interpolates the observed incidence curve, roughly guesstimating the underlying force of infection
btv.2008.padded<-data.frame(Day.since.19th.of.Sept=seq(1:max(btv.2008.df$Day.since.19th.of.Sept)),Cases=0)
btv.2008.padded[btv.2008.padded[,1]%in%c(btv.2008.df$Day.since.19th.of.Sept),2]<-btv.2008.df[,"Cases"]
btv.2008.interpolated<-as.data.frame(predict(smooth.spline(btv.2008.padded$Day.since.19th.of.Sept,btv.2008.padded$Cases,df=nrow(btv.2008.df)/4),newdata=data.frame(Date.since.19th.of.Sept=1:176)))
btv.2008.interpolated$y[btv.2008.interpolated$y<0]<-0
plot(btv.2008.interpolated)

infection.times<-rep(1:50,times=round(btv.2008.interpolated[1:50,2]))
btv.test.data<-data.frame(infection.time=infection.times,na=rlnorm(length(infection.times),
                  meanlog=log(ode.result[match(infection.times,ode.result[,"time"]),"nat"]),
                  sdlog=log(1.6)),
                      ab=rlnorm(length(infection.times),
                                  meanlog=log(ode.result[match(infection.times,ode.result[,"time"]),"abt"]),
                                  sdlog=log(1.6)))
                      

btv.test.data<-data.frame(infection.time=infection.times,na=rlnorm(length(infection.times),
                                                                   meanlog=log(btv.spline.df[match(infection.times,btv.spline.df[,"DPI"]),"PCR"]),
                                                                   sdlog=log(1.6)),
                          ab=rlnorm(length(infection.times),
                                    meanlog=log(btv.spline.df[match(infection.times,btv.spline.df[,"DPI"]),"AB"]),
                                    sdlog=log(1.6)))


##use BUGS to estimate parameters for indicator curve

autolib(rjags)
reverse.hindcasting.test<-jags.model(file.path(script.path,"BUGS code", "ReverseHindcasting.bugs"),
                         data=list(N=nrow(btv.test.data),
                                   ab.values=(btv.test.data$ab),
                                   na.values=(btv.test.data$na),
                                   btv.probs=rev(btv.2008.interpolated[1:50,2])/sum(btv.2008.interpolated[1:50,2])
                         )
)
update(reverse.hindcasting.test,5000)
test<-jags.samples(reverse.hindcasting.test,variable.names=c("a.ab",
                                                             "b.ab",
                                                             #"c.ab",
                                                             #"d.ab",
                                                             "a.na",
                                                             "b.na",
                                                             #"c.na",
                                                             #"d.na",
                                                             "scale",
                                                             "mean.ab",
                                                             "mean.na",
                                                             #"intercept.na","intercept.ab",
                                                             "est.time","theta.ab","theta.na"),n.iter=5000)
#test$est.time[,,]
plot(test$a.ab[,,1],type="l")


plot(rowMeans(test$est.time),mean(test$scale)*dbeta(rowMeans(test$est.time)/mean(test$scale),mean(test$a.na),mean(test$b.na)),xlim=c(0,60),ylim=c(0,10),ylab="AB",xlab="days")
lines(btv.spline.df[,"DPI"],btv.spline.df[,"PCR"],col="red")
points(btv.test.data$infection.time,btv.test.data$na,col="blue")

plot(rowMeans(test$est.time),mean(test$scale)*dbeta(rowMeans(test$est.time)/mean(test$scale),mean(test$a.ab),mean(test$b.ab)),xlim=c(0,60),ylim=c(0,6),ylab="NA",xlab="days")
lines(btv.spline.df[,"DPI"],btv.spline.df[,"AB"],col="red")
points(btv.test.data$infection.time,btv.test.data$ab,col="blue")


plot(rowMeans(test$est.time[,,1]),rowMeans(test$mean.na[,,1]),ylim=c(0,5),xlim=c(0,60))
lines(ode.result[,"time"],ode.result[,"abt"],col="red") 
points(btv.test.data$infection.time,btv.test.data$ab,col="red")
time<-seq(0,42,by=1)
with(test,plot(rowMeans(est.time),mean(intercept.ab)+mean(a.ab)*rowMeans(est.time)+mean(b.ab)*rowMeans(est.time)^2+mean(c.ab)*rowMeans(est.time)^3,type="l",ylab="AB",ylim=c(-2,6),xlim=c(0,60))) 
lines(btv.spline.df[,"DPI"],btv.spline.df[,"AB"],col="red")
with(test,plot(rowMeans(est.time),mean(intercept.na)+mean(a.na)*rowMeans(est.time)+mean(b.na)*rowMeans(est.time)^2+mean(c.na)*rowMeans(est.time)^3,type="l",ylab="NA",ylim=c(0,4),xlim=c(0,60)))
lines(btv.spline.df[,"DPI"],btv.spline.df[,"PCR"],col="red")
##Estimate it using linear regression to start with. Do that with the curve above as well! 
###Take the likely distribution of time since infection into account. Bugs code should be simple
#AB~rlnorm(mean.AB,theta)
#mean.ab=AB.FUN.est(estimated.infection.time)
#NA~rlnorm(NA.FUN.est(estimated.infection.time),theta)
#mean.na=NA.FUN.est(estimated.infection.time)
#estimated.infection.time~multinomial(btv.epidemic)
### AND priors...
