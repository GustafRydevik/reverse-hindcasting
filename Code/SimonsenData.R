require(sas7bdat)
simonsen.df<-read.sas7bdat(file.path(data.path,"mixed_elisa.sas7bdat"))
View(simonsen.df)
simonsen.df$DEBUT<-as.Date("1960-01-01")+simonsen.df$DEBUT
simonsen.df$DATO1<-as.Date("1960-01-01")+simonsen.df$DATO1
simonsen.df$DATO2<-as.Date("1960-01-01")+simonsen.df$DATO2
simonsen.df$DATO3<-as.Date("1960-01-01")+simonsen.df$DATO3
simonsen.df$DATO4<-as.Date("1960-01-01")+simonsen.df$DATO4

obs.array.simonsen<-as.matrix(simonsen.df[c("TID1","TID2","TID3","TID4")])

##"tmp can be used to Validate that I am storing the data correctly.
tmp<-array(rbind(matrix(paste(paste("N",1:10,sep=""),rep(paste(".IG",1:3,sep=""),each=10),".T1",sep=""),nrow=10,ncol=3),
                 matrix(paste(paste("N",1:10,sep=""),rep(paste(".IG",1:3,sep=""),each=10),".T2",sep=""),nrow=10,ncol=3),
                 matrix(paste(paste("N",1:10,sep=""),rep(paste(".IG",1:3,sep=""),each=10),".T3",sep=""),nrow=10,ncol=3),
                 matrix(paste(paste("N",1:10,sep=""),rep(paste(".IG",1:3,sep=""),each=10),".T4",sep=""),nrow=10,ncol=3)),
           dim=c(10,4,3))
simonsen.long<-data.frame(rbind(as.matrix(simonsen.df[c("IGG1","IGA1","IGM1")]),
                                as.matrix(simonsen.df[c("IGG2","IGA2","IGM2")]),
                                as.matrix(simonsen.df[c("IGG3","IGA3","IGM3")]),
                                as.matrix(simonsen.df[c("IGG4","IGA4","IGM4")])),
      Time=c(simonsen.df$TID1,simonsen.df$TID2,simonsen.df$TID3,simonsen.df$TID4),
      ID=simonsen.df$NR_)
names(simonsen.long)[1:3]<-c("IGG","IGA","IGM")


dput(simonsen.long,file=file.path(data.path,"simonsen_long.txt"))


