###Function by Gustaf Rydevik to generate simulated  population antibody
###and nucleic acid data,
### for some disease, including both infected and non-infected subjects.
### April-May 2011

LabdataGeneratorSimple<-function(SerumDecayFun,NaDecayFun,
                           timeFun,timeFun.args=NULL, 
                           errorFun,errorFun.args=NULL,...){

  #timeFun generates a set of infection times for a number of individuals

infection.times<-do.call("timeFun",c(timeFun.args,...))


##Generating expected mean test results
antibody.expected<-SerumDecayFun(infection.times,...)
na.expected<-NaDecayFun(infection.times,...)
  ##Generate output data frame
  antibody.data<-data.frame(index=seq_along(infection.times),
                            infection.time=infection.times,
                            antibody.mean=antibody.expected,
                            na.mean=na.expected)


##Add observational and individual variation.
antibody.data<-data.frame(antibody.data,
            do.call("errorFun",c(xy.df=list(antibody.data[c("antibody.mean","na.mean")]),errorFun.args))
                          )
             
names(antibody.data)[tail(seq_along(names(antibody.data)),2)]<-c("antibody.levels","na.results")
  ##Return end result
  return(antibody.data)
}
