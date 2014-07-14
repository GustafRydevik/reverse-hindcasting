
data {
  int <lower=0> Nobs; // number of observations
  int <lower=0> N; // number of individuals
  int <lower=0> ID[Nobs]; // ID numbers
  int <lower=0> I; // Number of different Igg antibodies
  real <lower=0> SamplingTimes[Nobs]; // Times of sampling for each individual
  matrix[Nobs,I] TestData;      // matrix measured test data
}

transformed data{
matrix[Nobs,I] logTestData;
vector[4] parsVector;
vector[I] omegaVector;
matrix[I,I] idOmega;
matrix[4,4] idPars;
parsVector<-rep_vector(1,4);
omegaVector<-rep_vector(1,I);
idOmega<-diag_matrix(omegaVector);
idPars<-diag_matrix(parsVector);

for(i in 1:I){
for(n in 1:Nobs){
  logTestData[n,i]<-log(TestData[n,i]);
}}

}

parameters {

  matrix[I,4] theta1IgLogmu[N];  // Array of parameters for each individual and igg, on logscale
  matrix[I,4] theta2IgLogmu; // Array of paramater means for each IGG, in logscale
  cov_matrix[4] igCov[I] ; //Covariance matrix  between parameters for each IGG
  cov_matrix[3] omega; // Covariance matrix between the three IGG values
}

transformed parameters {
  matrix<lower=0>[Nobs,I] estimatedIGG; // An array to store estimated IGG mean values in
  matrix[Nobs,I] logIGG;                 //matrix to store their log values in
  real <lower=0> XStar;
  real <lower=0> S;
  real <lower=0> a;
  real <lower=0> D;
  real <lower=0>  meanTmp;
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


}
  

model {
  for(i in 1:I){
    for(n in 1:N){
      to_vector(row(theta1IgLogmu[n],i))~multi_normal(to_vector(row(theta2IgLogmu,i)),igCov[i]);
    }


   igCov[i]~wishart(4,idPars);
}
  
      for(n in 1:Nobs) {
  to_vector(row(logTestData,n))~multi_normal(to_vector(logIGG[n]),omega);
      }
   omega~wishart(3,idOmega);
    
//Needs to add priors and starting values!

}
  
  
