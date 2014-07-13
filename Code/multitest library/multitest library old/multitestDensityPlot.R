###Plotting function for visualizing marginal densities of multiple test data.
### 
multitestDensityPlot<-function(LabDataGenerator.data,show.real=TRUE,time.cutoff=0,na.cutoff=8,ab.neg=TRUE,ab.pos=TRUE,ab.all=TRUE,legend=TRUE,...){
d0<-density(log(LabDataGenerator.data$antibody.levels[LabDataGenerator.data$infection.time<=time.cutoff]))
d1<-density(log(LabDataGenerator.data$antibody.levels[LabDataGenerator.data$infection.time>time.cutoff]))
d0.obs<-density(log(LabDataGenerator.data$antibody.levels[LabDataGenerator.data$na.result<na.cutoff]))
d1.obs<-density(log(LabDataGenerator.data$antibody.levels[LabDataGenerator.data$na.result>na.cutoff]))
d<-density(log(LabDataGenerator.data$antibody.levels))

plot(d$x,d$y,type="l",lty=1,col="white",...)

axis(1,col="lightgrey",tck=1)
axis(2,col="lightgrey",tck=1)

if(show.real){lines(d1$x,d1$y,col="red")
              lines(d0$x,d0$y,col="black")
             }

if(ab.pos){lines(d1.obs$x,d1.obs$y,col="red",lty=2,...)}
if(ab.all){lines(d$x,d$y,col="green",...)}
if(ab.neg){lines(d0.obs$x,d0.obs$y,type="l",lty=2,...)}
if(legend){legend("topleft",legend=c("True ab of noninfected",
                         "Obs ab of na negative",
                       "True ab of infected",
                         "Obs ab of na positive",
                         "Observed overall ab distribution"),
	col=c("black","black","red","red","green"),lty=c(1,2,1,2,1))}
}
