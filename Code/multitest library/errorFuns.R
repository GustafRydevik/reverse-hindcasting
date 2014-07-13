errorFun.lognorm<-function(xy.df,standard.deviation=log(1.1)){
  obs.error<-data.frame(rlnorm(nrow(xy.df),meanlog=0,sdlog=standard.deviation),
                        rlnorm(nrow(xy.df),meanlog=0,sdlog=standard.deviation))
return.df<-data.frame(x=(obs.error[,1])*(xy.df[,1]),
                      y=(obs.error[,2])*(xy.df[,2]))
return(return.df)
}


norm.error.heteroscedastic<-function(x){sapply(
                              x,function(x)rnorm(1,mean=x,sd=x/10))}


errorFun.exact<-function(xy.df,standard.deviation=log(1.1)){
  return(xy.df)
}
