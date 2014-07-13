
#Simulates data based on the expected path an infected
#animal should take through the test-space
DensityEstimator <- function(expected.path,errorFun,reps=1){
        expected.path.long<-do.call("data.frame",lapply(expected.path,rep,reps))
  			generated.data<-errorFun(expected.path.long)
       	return(generated.data)
}




