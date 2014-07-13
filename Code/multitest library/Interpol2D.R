Interpol.2D.path<-function(xpos,ypos,tpos,npoints){
xpos<-xpos[order(tpos)]
ypos<-ypos[order(tpos)]
points.dist<-sqrt(diff(xpos)^2+diff(ypos)^2)
path.cumlength<-cumsum(points.dist)
equaldist.cumlength<-seq(0,tail(path.cumlength,1),
length.out=npoints-1)
stddist<-equaldist.cumlength[2]
Interpol.points<-list(x=
   vector(length=npoints,mode="numeric"),
   y=vector(length=npoints,mode="numeric"))
Interpol.points$x[1]=xpos[1]
Interpol.points$y[1]=ypos[1]
for(i in 2:(npoints-5)){
start.point=c(Interpol.points$x[i-1], Interpol.points$y[i-1])
ndx<-which(path.cumlength>equaldist.cumlength[i-1]&
path.cumlength<=equaldist.cumlength[i+5])
x.lm<-xpos[ndx+1]-start.point[1]
y.lm<-ypos[ndx+1]-start.point[2]
##Slope<-coef(lm(y-1~x))#deal with wrong dir?
##alt?:
cos.values <- x.lm/sqrt(x.lm^2+y.lm^2)
Angle=mean(((y.lm<0)*2*pi)+acos(cos.values)*sign(y.lm))
Ydist=sin(Angle)*stddist
Xdist=cos(Angle)*stddist
##Xdist=sqrt(stddist^2/(1+Slope^2))
##Ydist=Xdist*Slope
Next.point<-start.point+c(Xdist,Ydist)
##if(any(is.nan(Next.point))){browser()}
Interpol.points$x[i]=Next.point[1]
Interpol.points$y[i]=Next.point[2]                  
}
Interpol.points<-as.data.frame(Interpol.points)
Interpol.points<-Interpol.points[!is.na(Interpol.points$x),]
Interpol.points<-Interpol.points[!is.na(Interpol.points$y),]
Interpol.points<-Interpol.points[!duplicated(Interpol.points),]
return(Interpol.points)
}
