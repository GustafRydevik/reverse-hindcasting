
data {
  int <lower=0> N; // number of individuals
  int <lower=0> T; // number of sampling times per individual
  int <lower=0> I; // Number of different Igg antibodies
  real <lower=0> SamplingTimes[N,T]; // Times of sampling for each individual
  matrix[T,I] TestData[N];      // Array of measured test data
}

transformed data{
  matrix[T,I] logTestData[N];
vector[4] parsVector;
vector[I] omegaVector;
matrix[I,I] idOmega;
matrix[4,4] idPars;

parsVector<-rep_vector(1,4);
omegaVector<-rep_vector(1,I);
idOmega<-diag_matrix(omegaVector);
idPars<-diag_matrix(parsVector);

for(i in 1:I){
for(t in 1:T){
for(n in 1:N){
  logTestData[n,t,i]<-log(TestData[n,t,i]);
}}}

}

parameters {

  matrix[I,4] theta1IgLogmu[N];  // Array of parameters for each individual and igg, on logscale
  matrix[I,4] theta2IgLogmu; // Array of paramater means for each IGG, in logscale
  cov_matrix[4] igCov[I] ; //Covariance matrix  between parameters for each IGG
  cov_matrix[3] omega; // Covariance matrix between the three IGG values

}

transformed parameters {
  real <lower=0> estimatedIGG[N,T,I]; // An array to store estimated IGG mean values in
  matrix[T,I] logIGG[N];                 //array to store their log values in
  real <lower=0> XStar;
  real <lower=0> S;
  real <lower=0> a;
  real <lower=0> D;
  real <lower=0>  meanTmp;
  for( i in 1:I){  /// Calculating the estimated IGG values acc. to the used model 
    for(t in 1:T){
    for(n in 1:N){

      XStar<-exp(theta1IgLogmu[n,i,1]);
      D<-exp(theta1IgLogmu[n,i,2]);
      a<-exp(theta1IgLogmu[n,i,3]);
      S<-exp(theta1IgLogmu[n,i,4]);
       if (SamplingTimes[n,t]<D) {
          meanTmp<-(XStar+(S+a*S*(D-SamplingTimes[n,t])-S*(1+a*D)*exp(-a*SamplingTimes[n,t]))/(D*square(a)));
          
      } else {
        meanTmp<-XStar+((S*exp(a*D)-S*(1+a*D))*exp(-a*SamplingTimes[n,t]))/(D*square(a));

      }
    estimatedIGG[n,t,i]<-meanTmp;
    logIGG[n,t,i]<-log(estimatedIGG[n,t,i]);

    }}}


}
  

model {

  for(i in 1:I){
    for(n in 1:N){
      to_vector(theta1IgLogmu[n][i])~multi_normal(to_vector(theta2IgLogmu[i]),igCov[i]);
    }
   igCov[i]~wishart(4,idPars);
}
  
    for(t in 1:T){
      for(n in 1:N) {
  to_vector(logTestData[n][t])~multi_normal(to_vector(logIGG[n][t]),omega);
      }}
   omega~wishart(3,idOmega);
    
//Needs to add priors and starting values!

}
  
  
