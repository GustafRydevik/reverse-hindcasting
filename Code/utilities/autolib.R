autolib<-function(Package,mirror="http://stat.ethz.ch/CRAN",...){
if(!suppressWarnings(require(as.character(substitute(Package)),
            character.only=TRUE,quietly=TRUE,...))){
  Package<-as.character(substitute(Package))
  install.packages(Package,repos=mirror,...)
	library(Package,character.only=TRUE)	
	}		
}

