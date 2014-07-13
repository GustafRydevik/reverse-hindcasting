
data {
  int <lower=0> N; // number of individuals
  int <lower=0> T; // number of sampling times per individual
  int <lower=0> I; // Number of different Igg antibodies
  real SamplingTimes[N,T]; // Times of sampling for each individual
  real TestData[N,T,I];      // Array of measured test data
  int Timelength;
  real PriorStart;
}


parameters {
  real<lower=0,upper=1> mu;
  real<lower=0> psi;
  real<lower=0,upper=PriorStart> eta;
  vector<lower=0,upper=eta>[N] it; //it =infection time
  real<lower=log(1.02),upper=log(4)> logsd1;
  real<lower=log(1.02),upper=log(4)> logsd2;
}

transformed parameters {
  real<lower=0> alpha;
  real<lower=0> beta;
  vector<lower=0,upper=1>[N] itrel;
  vector<lower=0>[N] TestMean1;
  vector<lower=0>[N] TestMean2;
  
  alpha<-mu*psi+0.1;
  beta<-psi-mu*psi+0.1;
  itrel<-it/eta;
  
  
  
  for(i in 1:N){
    int ndxTime;  //temporary index tracker 
    real itt;    ///temporary infection time it
    
    itt<-it[i];
    ndxTime <- 1; 
    while ((ndxTime + 1< floor(itt)) && (ndxTime< Timelength)){
      ndxTime <- ndxTime + 1;
    }
    
    
    TestMean1[i]<- t1Look[ndxTime];
    TestMean2[i]<- t2Look[ndxTime];
  }
  
  
}
model {
  
  ### This models a situation where infections give complete immunity
  #assumes an exponential distribution for time since infection. 
  
  itrel~beta(alpha,beta);
  
  // might want to add something like  later EpiStart<-max(InfTime);
  
  
  
  Test1Data ~ lognormal(log(TestMean1),logsd1); 
  Test2Data ~ lognormal(log(TestMean2),logsd2);
  
  
  logsd1~uniform(log(1.02),log(4));
  logsd2~uniform(log(1.02),log(4));
  psi~exponential(0.15);
  mu~uniform(0.01,0.99);
  eta~uniform(1,PriorStart);
}



# Simonsen.likelihood.gen<-function(ig.values.array,     #observed antibody test values 
#                                   #Each of length N do dim=c(n,4,3)          
#                                   obs.array,            #times of observation
#                                   # each of length N 
#                                   Theta2.igg.prior,  #Priors for the mean of curve parameters
#                                   Theta2.igm.prior,
#                                   Theta2.iga.prior,
#                                   omega.prior,       ##Prior for the correlation of observed test values at a given time 
#                                   N,n.times=4,n.igg=3,n.pars=4){
#   
#   
#   ##Calculating the log likelihood for individual pars given global mean and variance
#   ll.theta.individual<-function(Theta.individual,Theta.mu,Theta.sigma){
#     dmnorm(x=log(Theta.individual),mean=log(Theta.mu),varcov=Theta.sigma,log=T)}
#   
#   ##Calculating the log likelihood for ab values given expected mean at a certain time and variance matrix
#   ll.theta.values<-function(ig.values,ig.values.mean,omega){
#     dmnorm(x=log(ig.values),mean=log(ig.values.mean),varcov=omega,log=T)}
#   
#   
#   function(pars){  
#     
#     ##Initializing the parameters
#     theta1.ig.mu<-array(exp(pars[grep("theta1.ig.mu",names(pars))]),dim=c(N,3,4)) ##N*3*4 parameters
#     
#     #^^ the above needs to be coerced into an array of cov matrices!
#     theta2.ig.mu<-array(exp(pars[grep("theta2.ig.mu",names(pars))]),dim=c(3,4)) ##3*4 mean curve parameters    
#     theta2.ig.sigma.triangle<-pars[grep("theta2.ig.sigma",names(pars))] ##3*(4+3+2+1) pars.. How do we force a proper corr matrix?
#     
#     theta2.ig.sigma<-convert.from.triangle(theta2.ig.sigma.triangle,n.pars)
#     
#     omega.triangle<-pars[grep("omega",names(pars))] ###3*2*1 pars: cov matrix for the test values
#     omega<-convert.from.triangle(omega.triangle,n.igg)
#     
#     ig.mean<-array(dim=c(N,n.times,n.igg),dimnames=c("ID","Igtype","time"))
#     for(i in 1:N){
#       for(ig in 1:n.igg){
#         igCurve.current<-igCurve(X.star=(theta1.ig.mu[i,ig,1]),
#                                  S=(theta1.ig.mu[i,ig,2]),
#                                  a=(theta1.ig.mu[i,ig,3]),
#                                  D=(theta1.ig.mu[i,ig,4]))  
#         ig.mean[i,,ig]<-igCurve.current(obs.array[i,])
#       }}
#     
#     
#     
#     ##sum together the likelihood for all parts of the model at the current point
#     post.ll<-omega.prior(omega)+
#       Theta2.igg.prior(Theta2.ig.mu[1,],Theta2.ig.sigma[1,,])+
#       Theta2.igm.prior(Theta2.ig.mu[2,],Theta2.ig.sigma[2,,])+
#       Theta2.igg.prior(Theta2.ig.mu[3,],Theta2.ig.sigma[3,,])+
#       sum(sapply(1:n.igg,function(x){
#         ll.theta.individual(theta1.ig.mu[,x,],
#                             theta2.ig.mu[x,],
#                             theta2.ig.sigma[x,,])}))+
#       sum(ll.theta.values(apply(ig.values.array,3,rbind),
#                           apply(ig.mean,3,rbind),
#                           omega))
#     
#     return(post.ll)
#   }
# }
# 
# 
# 
# ###Testing
# igg.curve<-igCurve(X.star=0.15,D=15,a=0.05,S=0.3)
# iga.curve<-igCurve(X.star=0.3,D=15,a=0.1,S=0.3)
# igm.curve<-igCurve(X.star=0.45,D=15,a=0.15,S=0.3)
# obs.array<-array(data=sort(sample(1:50,40,replace=T)),dim=c(10,4))
# ig.values.sim<-array(data=c(igg.curve(c(obs.array)),
#                             iga.curve(c(obs.array)),
#                             igm.curve(c(obs.array))),dim=c(10,4,3))
# 
# Theta2.igg.prior<-function(...){return(0)}  #Priors for the mean of curve parameters
# Theta2.igm.prior<-function(...){return(0)}
# Theta2.iga.prior<-function(...){return(0)}
# omega.prior<-function(...){return(0)}
# 
# n.igg<-3
# n.pars<-4
# N<-10
# 
# theta2.ig.mu<-array(data=log(c(c(X.star=0.15,D=15,a=0.05,S=0.3),
#                                c(X.star=0.3,D=15,a=0.1,S=0.3),
#                                c(X.star=0.45,D=15,a=0.15,S=0.3))),
#                     dim=c(n.igg,n.pars),dimnames=list(ig=c("igg","iga","igm"),
#                                                       pars=c("X.star","D","a","S")))
# 
# theta2.ig.sigma<-aperm(sapply(1:n.igg,function(x){diag(4)},simplify="array"))
# dimnames(theta2.ig.sigma)<-list(ig=c("igg","iga","igma"),
#                                 pars1=c("X.star","D","a","S"),pars2=c("X.star","D","a","S"))
# omega<-diag(3)
# dimnames(omega)<-list(ig1=c("igg","iga","igm"),
#                       ig2=c("igg","iga","igm"))
# 
# theta1.ig.mu<-aperm(sapply(1:n.igg,function(x){rmnorm(10,mean=theta2.ig.mu[x,],theta2.ig.sigma[x,,])},simplify="array"),c(1,3,2))
# dimnames(theta1.ig.mu)<-list(id=(1:N),
#                              ig=c("igg","iga","igm"),
#                              pars=c("X.star","D","a","S"))
# 
# testpar.vector<-c(c(theta1.ig.mu),
#                   c(theta2.ig.mu),
#                   c(apply(theta2.ig.sigma,1,function(x)x[lower.tri(x,diag=T)])),
#                   omega[lower.tri(omega,diag=T)])
# 
# 
# tri.ndx.npars<-which(lower.tri(diag(n.pars),diag=T),arr.ind=T)
# tri.ndx.igg<-which(lower.tri(diag(n.igg),diag=T),arr.ind=T)
# names(testpar.vector)<-c(  ##This gives names to each element of the looong vector, so we can look it up later. 
#   do.call("paste",  
#           args=c(list("theta1.ig.mu"),
#                  expand.grid(dimnames(theta1.ig.mu),stringsAsFactors=FALSE)[c(2,3,1)],
#                  list(sep="."))),
#   do.call("paste",
#           args=c(list("theta2.ig.mu"),
#                  expand.grid(dimnames(theta2.ig.mu),stringsAsFactors=FALSE),
#                  list(sep="."))),
#   do.call("paste",
#           args=c(list("theta2.ig.sigma"),
#                  list(rep(dimnames(theta2.ig.sigma)[[1]],each=nrow(tri.ndx.npars))),
#                  list(dimnames(theta2.ig.sigma)[[2]][tri.ndx.npars[,1]]),
#                  list(dimnames(theta2.ig.sigma)[[3]][tri.ndx.npars[,2]]),
#                  list(sep="."))),
#   do.call("paste",
#           args=c(list("omega"),
#                  list(dimnames(omega)[[1]][tri.ndx.igg[,1]]),
#                  list(dimnames(omega)[[2]][tri.ndx.igg[,2]]),
#                  list(sep=".")))
# )
# 
# simonsen.ll.test<-Simonsen.likelihood.gen(ig.values.array=ig.values.sim,     #observed antibody test values 
#                                           #Each of length N do dim=c(n,4,3)          
#                                           obs.array=obs.array,            #times of observation
#                                           # each of length N 
#                                           Theta2.igg.prior=Theta2.igg.prior,  #Priors for the mean of curve parameters
#                                           Theta2.igm.prior=Theta2.igm.prior,
#                                           Theta2.iga.prior=Theta2.iga.prior,
#                                           omega.prior=omega.prior,       ##Prior for the correlation of observed test values at a given time 
#                                           N,n.times=4,n.igg=3,n.pars=4)
# 
# simonsen.ll.test(testpar.vector)


