
##As used by Simonsen etal in 
#### "Estimation of incidences of infectious diseases
####### based on antibody measurements"
###### Statist. Med. 2009; 28:1882–1895

igCurve<-function(X.star,S,a,D){
  function(t){
    if(a<1e-8){a<-1e-8}  ### This is a fudge for numerical errror problems...
    Ig.t<-ifelse(t<D,
    X.star+(S+a*S*(D-t)-S*(1+a*D)*exp(-a*t))/(D*a^2),
      X.star+(S*exp(a*(D-t))-S*(1+a*D)*exp(-a*t))/(D*a^2)
    )
    return(Ig.t)
  }}


bactLoad<-function(D){
  function(t){
    ifelse(t<D,
           1-t/D,
           0.00001)
  }
}


bactLoad.exp<-function(D){
function(t){
  exp(-t*log(D)/D)
  }
}


###Likelihood functions...
###1-level version

###Theta=c(X.star,S,a,D)
##P(theta|obs)=P(obs|theta)P(theta)
###P(Obs|theta)=lognorm(igCurve,bactload|\sigma,theta)


