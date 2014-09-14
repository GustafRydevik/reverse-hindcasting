
data {
  int <lower=0> Nobs; // number of observations
  int <lower=0> N; // number of individuals
  int <lower=0> ID[Nobs]; // ID numbers
  real <lower=0> SamplingTimes[Nobs]; // Times of sampling for each individual
  matrix<lower=0>[Nobs,2] TestData;      // matrix measured test data
}

transformed data{
matrix[Nobs,2] logTestData;
for(n in 1:Nobs){
  logTestData[n,1]<-log(TestData[n,1]);
  logTestData[n,2]<-log(TestData[n,2]);
}}



parameters {

  matrix[N,4] theta1IgLogmu;  // Array of parameters for each individual, on logscale
  real theta2IgLogmu[4]; // Array of paramater means for each IGG, in logscale
  corr_matrix[4] igCorr; //Correlation matrix  between parameters for each IGG
  vector<lower=0> [4] igPars; //variance for parameters  -- See p. 38 in the STAN manual.
  real<lower=0> tauIg; //measurement error for IG
  real<lower=0> tauNa; //measurement error for NA
}

transformed parameters {
  matrix<lower=0>[Nobs,2] estimatedIGG; // An array to store estimated Test mean values in
  real logIGG[Nobs,2];                 //matrix to store their log values in
  real <lower=0> XStar;
  real <lower=0> S;
  real <lower=0> a;
  real <lower=0> D;
  real <lower=0>  meanTmp;
  real <lower=0>  meanNa;
  {  /// Calculating the estimated IGG values acc. to the used model 
    for(n in 1:Nobs){
    int index;
    index <- ID[n]; 
      XStar<-exp(theta1IgLogmu[ID[n],1]);
      D<-exp(theta1IgLogmu[ID[n],2]);
      a<-exp(theta1IgLogmu[ID[n],3]);
      S<-exp(theta1IgLogmu[ID[n],4]);
       if (SamplingTimes[n]<D) {
          meanTmp<-(XStar+(S+a*S*(D-SamplingTimes[n])-S*(1+a*D)*exp(-a*SamplingTimes[n]))/(D*square(a)));
          meanNa<-(1-SamplingTimes[n]/D);
          
      } else {
        meanTmp<-XStar+((S*exp(a*D)-S*(1+a*D))*exp(-a*SamplingTimes[n]))/(D*square(a));
        meanNa<-1/10^6;

      }
    estimatedIGG[n,1]<-meanTmp;
    estimatedIGG[n,2]<-meanNa;
    logIGG[n,1]<-log(estimatedIGG[n,1]);
    logIGG[n,2]<-log(estimatedIGG[n,2]);

    }}

}
  

model {
       matrix[4,4] igCov;  

       igCov<-diag_matrix(igPars)*igCorr*diag_matrix(igPars);
       igPars~cauchy(0,0.1);  ///prior for the within-parameter scale
       igCorr~lkj_corr(1.5);   ///prior for the within-parameter correlation
       
       tauIg ~ cauchy(0,0.0001); ///prior for the between-measurement scale
       tauNa ~ cauchy(0,0.0001); ///prior for the between-measurement scale

          

theta2IgLogmu[1]~cauchy(0,0.1); //normal(log(1),log(sqrt(10)));
theta2IgLogmu[2]~cauchy(0,0.1); //normal(log(10),log(sqrt(10)));
theta2IgLogmu[3]~cauchy(0,0.1); //normal(log(0.1),log(sqrt(10)));
theta2IgLogmu[4]~cauchy(0,0.1); //normal(log(1),log(sqrt(10)));

    for(n in 1:N){
      to_vector(row(theta1IgLogmu,n))~multi_normal(to_vector(theta2IgLogmu),igCov);
    }


  
      for(n in 1:Nobs) {
  logTestData[n,1]~normal(logIGG[n,1],tauIg);
  logTestData[n,2]~normal(logIGG[n,2],tauNa);
      }
}
  
  
