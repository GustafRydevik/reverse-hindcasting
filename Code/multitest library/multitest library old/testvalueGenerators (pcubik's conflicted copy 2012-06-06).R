#######Deterministic models for NA/antibody interaction
###using differential equations

###Parameters for the lottka-voltera diff. equations below
  
Pars<-c(alpha.na=30,beta.na=0.5,alpha.ab=5,beta.ab=0.1)
State<-c(nat=1,abt=0)
Time=seq(0,10,by=0.001)


###Lotke-Voltera equation based on Giles' idea.
autolib(deSolve)
ab.na.diffmod<-function(Time,State,Pars){
  with(as.list(c(State,Pars)),{
      dna=nat*(alpha.na-beta.na*abt)
        dab=alpha.ab*nat-beta.ab*abt
        return(list(c(dna,dab)))
      })}
inf.levels <- as.data.frame(ode(func = ab.na.diffmod, y = State, parms = Pars, times = Time))
inf.levels$abt<-2^(inf.levels$abt/12)

###Deterministic nucleic acid function, pulling values from Giles' diff model
naDetFun<-(function(Inf.levels=inf.levels,timescale=NULL){
  function(infection.time){
   sapply(infection.time,function(x){
          Inf.levels[min(which(Inf.levels$time>=x)),"nat"]}
        )}}
   )()



###Deterministic antibody function, pulling values from Giles' diff model
abDetFun<-(function(Inf.levels=inf.levels,timescale=NULL){
  function(infection.time){
   abt.response<-sapply(infection.time,function(x){
          Inf.levels[min(which(Inf.levels$time>=x)),"abt"]}
          )
   abt.response[infection.time==0]<-rlnorm(sum(infection.time==0),mean=log(mean(Inf.levels$abt[100])))
   return(abt.response)
	}}
   )()


###Old AB decay function just based on a lognormal distribution of  responses with added half life
AntibodyDecay<-function(response.mean,half.life=1){
function(infection.time,number.infections){
initial.response<-rlnorm(length(infection.time),mean=log(response.mean),sd=log(response.mean)/10)*(infection.time!=0)
decayed.response<- 2^(-(infection.time/half.life))*initial.response+
                       rlnorm(length(infection.time),
                              mean=log(response.mean/2^8),
                              sd=log(response.mean)/10)
return(decayed.response)
}}

###Old NA decay function just based on a gamma distribution of  responses with added half life
NaDecay<-function(infection.duration){
function(infection.time){
realised.infection.duration<-rgamma(length(infection.duration),shape=1,
                                    scale=infection.duration)
test.result<-(infection.time<=realised.infection.duration)*(infection.time!=0)
return(test.result)
}}
 
