
data {
  int <lower=0> Nobs; // number of observations
  int <lower=0> N; // number of individuals
  int <lower=0> ID[Nobs]; // ID numbers
  int <lower=0> I; // Number of different Igg antibodies
  real <lower=0> SamplingTimes[Nobs]; // Times of sampling for each individual
  matrix<lower=0>[Nobs,I] TestData;      // matrix measured test data
}

transformed data{
matrix[Nobs,I] logTestData;
for(i in 1:I){
for(n in 1:Nobs){
  logTestData[n,i]<-log(TestData[n,i]);
}}

}

parameters {

  matrix[I,4] theta1IgLogmu[N];  // Array of parameters for each individual and igg, on logscale
  matrix[I,4] theta2IgLogmu; // Array of paramater means for each IGG, in logscale
  corr_matrix[4] igCorr[I]; //Covariance matrix  between parameters for each IGG
  vector<lower=0>[4] igTau[I]; //scale for the correlation  -- See p. 38 in the STAN manual.
  corr_matrix[I] omegaCorr; // Correlation matrix between the  IGG values
  vector<lower=0>[I] tau; //Scale for the correlation
}

transformed parameters {
  matrix<lower=0>[Nobs,I] estimatedIGG; // An array to store estimated IGG mean values in
  matrix[Nobs,I] logIGG;                 //matrix to store their log values in
  real <lower=0> XStar;
  real <lower=0> S;
  real <lower=0> a;
  real <lower=0> D;
  real <lower=0>  meanTmp;
  matrix[4,4] igCov[I];  
  matrix[I,I] omega;


  for( i in 1:I){  /// Calculating the estimated IGG values acc. to the used model 
    for(n in 1:Nobs){
    int index;
    index <- ID[n]; 
      XStar<-exp(theta1IgLogmu[ID[n]][I,1]);
      D<-exp(theta1IgLogmu[ID[n]][I,2]);
      a<-exp(theta1IgLogmu[ID[n]][I,3]);
      S<-exp(theta1IgLogmu[ID[n]][I,4]);
       if (SamplingTimes[n]<D) {
          meanTmp<-(XStar+(S+a*S*(D-SamplingTimes[n])-S*(1+a*D)*exp(-a*SamplingTimes[n]))/(D*square(a)));
          
      } else {
        meanTmp<-XStar+((S*exp(a*D)-S*(1+a*D))*exp(-a*SamplingTimes[n]))/(D*square(a));

      }
    estimatedIGG[n,i]<-meanTmp;
    logIGG[n,i]<-log(estimatedIGG[n,i]);
    }}


  for(i in 1:I){
       igCov[i]<-diag_matrix(igTau[i])*igCorr[i]*diag_matrix(igTau[i]);
    }

       omega <- diag_matrix(tau) * omegaCorr * diag_matrix(tau);


}
  

model {




       omegaCorr ~ lkj_corr(1);  /////prior for the within-measurement correlation
      
  for(i in 1:I){
  tau[i] ~ cauchy(0,0.05); ///prior for the within-measurement scale
  igTau[i]~cauchy(0,0.5);  ///prior for the within-parameter scale
  igCorr[i]~lkj_corr(1.5);   ///prior for the within-parameter correlation

theta2IgLogmu[i,1]~cauchy(0,0.1); //normal(log(1),log(sqrt(10)));
theta2IgLogmu[i,2]~cauchy(0,0.1); //normal(log(10),log(sqrt(10)));
theta2IgLogmu[i,3]~cauchy(0,0.1); //normal(log(0.1),log(sqrt(10)));
theta2IgLogmu[i,4]~cauchy(0,0.1);//normal(log(1),log(sqrt(10)));

    for(n in 1:N){
      to_vector(row(theta1IgLogmu[n],i))~multi_normal(to_vector(row(theta2IgLogmu,i)),igCov[i]);
    }

}
  
      for(n in 1:Nobs) {
  to_vector(row(logTestData,n))~multi_normal(to_vector(logIGG[n]),omega);
      }
}
  
  
