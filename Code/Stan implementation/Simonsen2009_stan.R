
data {
  int <lower=0> N; // number of observations
  real Test1Data[N];      // estimated treatment effects
  real Test2Data[N];      // estimated treatment effects
  int Timelength;
  real t1Look[Timelength];      // looking up test curve for test 1
  real t2Look[Timelength];      // looking up test curve for test 2
  real timeLook[Timelength];    // time indexing vector
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


