
kinetics.data<-data.frame(time=seq(1,sampling.time,by=0.1),
           curve.true=igCurve.set(seq(1,sampling.time,by=0.1)),
           curve.mean=curve.mean,
           curve.05=curve.05,
           curve.95=curve.95
           )
plot(seq(1,sampling.time,by=0.1),igCurve.set(seq(1,sampling.time,by=0.1)),type="l",col="grey",lwd=2,ylim=c(0,2.5))
lines(seq(1,sampling.time,by=0.1),curve.mean,col="red",lwd=2)
lines(seq(1,sampling.time,by=0.1),curve.05,lty=2,col="red",lwd=2)
lines(seq(1,sampling.time,by=0.1),curve.95,lty=2,col="red",lwd=2)



ggplot(kinetics.data, aes(x=time,y=curve.true))+
  geom_line(col="black",size=2)+
  geom_smooth(aes(y=curve.mean,ymin=curve.05,ymax=curve.95),col="chartreuse4",size=1.5,lty=2)+
  theme_minimal()+labs(y="Test level",x="Days")+theme(text=element_text(size=40))



yval.epistart<-pert.long.est.exp.plotdata$density[min(which(pert.long.est.exp.plotdata$day>pert.long.obs-quantile(pert.long.allchains[,"EpiStart"],0.5)))]*scaling
p.est.vs.true.pert<-ggplot(pert.long.est.exp.plotdata,aes(x=day,y=density))+
  geom_bar(data=data.frame(day=pertussis.interpolated.round$x[seq(1,nrow(pertussis.interpolated.round),by=2)],density=hist(rep(pertussis.interpolated.round$x,times=pertussis.interpolated.round$y),plot=F,breaks=pert.breaks)$density),
           aes(y=density*260),stat="identity",colour=sruc.col,fill=sruc.col,alpha=0.5)+
  geom_smooth(size=3,alpha=1,aes(y=density*scaling,ymin=density.0.05*scaling,ymax=density.95*scaling),col="darkred",data=pert.long.est.exp.plotdata,,stat="identity")+
  #theme_gray(28)+
  theme(legend.position = c(0.85, 0.85))+
  ylab("Density of infected cases")+density.scales+
  scale_x_continuous(name="weeks after start of epidemic",breaks=seq(7,pertussis.endtime,by=7),labels=ifelse(rep(c(1,0),length.out=floor(pertussis.endtime/7)),(1):floor(pertussis.endtime/7),""),limits=c(0,pertussis.endtime))+
  geom_errorbarh(aes(y=yval.epistart,xmax=pert.long.obs-quantile(pert.long.allchains[,"EpiStart"],0.05), xmin=pert.long.obs-quantile(pert.long.allchains[,"EpiStart"],0.95)), 
                 width=0,colour="black",alpha=0.5, lty=1,size=2,height=0.09) + 
  geom_point(aes(y=yval.epistart,
                 x=pert.long.obs-quantile(pert.long.allchains[,"EpiStart"],0.5)),size=4,colour="black",alpha=0.5)+
  theme_minimal()+theme(text=element_text(size=40))+
  labs(title="Wisconsin Epidemic of Whooping Cough",size=25)+
  scale_y_continuous(limits=c(0,4),name="Cases")

print(p.est.vs.true.pert)