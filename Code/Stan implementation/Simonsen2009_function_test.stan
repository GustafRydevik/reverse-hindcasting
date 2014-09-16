functions {
real smoothCeiling(real t,real D){

real tTrunc;
real smooth;
real tscaled;
real pi;
smooth<-1;
pi<-3.14159;
tscaled<-(t-D+2.0)/smooth;
//shifting d by two above to avoid overrun

tTrunc<-(0.5-atan(tscaled)/pi)*t+(D)*(atan(tscaled)/pi+0.5)-1/pi;

return tTrunc;
}



real igKinetic(real t, real XStar, real S,real a, real D) {
real igValue;
real tTrunc;
tTrunc<-smoothCeiling(t,D);

          igValue<-(XStar+(S+a*S*(D-tTrunc)-S*(1+a*D)*exp(-a*t))/(D*square(a)));
return igValue;
}
}

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

   
  vector[4] theta1IgLogmu_t[I,N]; // Array of parameters for each individual and igg, on logscale

  matrix[I,4] theta2IgLogmu; // Array of paramater means for each IGG, in logscale
  corr_matrix[4] igCorr[I]; //Covariance matrix  between parameters for each IGG
  real<lower=0> igTau[I,4]; //scale for the correlation  -- See p. 38 in the STAN manual.
  corr_matrix[I] omegaCorr; // Correlation matrix between the  IGG values
  real<lower=0> tau[I]; //Scale for the correlation
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
   XStar<-exp(theta1IgLogmu_t[I,ID[n],1]);
      D<-exp(theta1IgLogmu_t[I,ID[n],2]);
      a<-exp(theta1IgLogmu_t[I,ID[n],3]);
      S<-exp(theta1IgLogmu_t[I,ID[n],4]);

     meanTmp<-igKinetic(SamplingTimes[n], XStar,S,a,D);
    estimatedIGG[n,i]<-meanTmp;
    logIGG[n,i]<-log(estimatedIGG[n,i]);
    }}


  for(i in 1:I){
       igCov[i]<- quad_form_diag(igCorr[i],to_vector(igTau[i]));
    }

       omega <-  quad_form_diag(omegaCorr,to_vector(tau));


}
  

model {

  omegaCorr ~ lkj_corr(1);  /////prior for the within-measurement correlation      
  tau ~ cauchy(0,1); ///prior for the within-measurement scale

  for(i in 1:I){
  igTau[I] ~ cauchy(0,1);  ///prior for the within-parameter scale
  igCorr[i]~lkj_corr(1.5);   ///prior for the within-parameter correlation
}

  for(i in 1:I){
to_vector(theta2IgLogmu[i])~cauchy(0,1); //normal(log(1),log(sqrt(10)));
}

  for(i in 1:I){
    theta1IgLogmu_t[i]~multi_normal(theta2IgLogmu[i],igCov[i]);
}
  
      for(n in 1:Nobs) {
  logTestData[n]~multi_normal(to_vector(logIGG[n]),omega);
      }
}
  
  
