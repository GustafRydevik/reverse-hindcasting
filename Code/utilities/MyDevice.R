###Function by Gustaf Rydevik, 2009-12-03 gustaf.rydevik@gmail.com
## Created to facilitate easy changes in the file format of generated graphs.
## Gen.device() generates a device function that is a copy of an existing function, but
## with (possibly) new defaults.
## Wanted.device can be the name of any device you choose: png(),jpeg(),postcript(),etc.
## if fileEnding is missing, the function uses the Wanted.device name as file ending.
## I then use My.device in the rest of the script file, meaning that I only have to change file 
##format in one location  (in the argument of Gen.device()) to do so for all susequent graphs.

Gen.device<-function(Wanted.device="png",fileEnding=NULL,...){
    dots<-list(...)
    ending<-Wanted.device
    Wanted.device<-get(Wanted.device)
    if(!is.null(fileEnding)) ending<-fileEnding
        generated.device<-function(File,...){
        dots2<-list(...)
        File<-paste(File,ending,sep=".")
        dots[which(names(dots)%in%names(dots2))]<-NULL
        if(ending!="pdf"){do.call(Wanted.device,c(filename=File,dots,dots2))}
        if(ending=="pdf"){do.call(Wanted.device,c(file=File,dots,dots2))}
    }
    return(generated.device)
}


##example
# My.device(File="test")
# plot(rnorm(1999))
# dev.off()

