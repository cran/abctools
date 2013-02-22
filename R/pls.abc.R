pls.abc <-function (obs,param,sumstats, obspar=NULL, abcmethod=abc,transfile = "Routput_test", bc=FALSE,err.only=TRUE,errfn=rsse,...) {

# For the function to work, the script "transformer" should be in the PATH.

if(!is.matrix(obs)|is.data.frame(obs)){
        obs<-matrix(obs,nrow=1)
}
if(!is.matrix(param)|is.data.frame(param)){
        param<-as.matrix(param)
}
if(!is.matrix(sumstats)|is.data.frame(sumstats)){
        sumstats<-as.matrix(sumstats)
}
if(!is.null(obspar)|is.data.frame(obspar)){
        if(!is.matrix(obspar)){
                obspar<-matrix(obspar,byrow=T,ncol=ncol(param))
        }
        if(nrow(obs)!=nrow(obspar)){
                stop("Please supply observed statistics and observed parameter matrices with the same number of rows!\n")
        }
}

# dump out observed stats:

write.table(obs,file="obsfile.txt",quote=FALSE,col.names=TRUE,row.names=FALSE)

# dump out simulation stats:

write.table(sumstats,file="sumfile.txt",quote=FALSE,col.names=TRUE,row.names=FALSE)

cmd1 <- paste("transformer", transfile, "sumfile.txt", "newstats.txt")
cmd2 <- paste("transformer", transfile, "obsfile.txt", "newobsstats.txt")
if (bc) {
        cmd1 <- paste(cmd1, "boxcox")
        cmd2 <- paste(cmd2, "boxcox")
}
cat("transforming simulations...")
system(cmd1)
cat("transforming observations...")
system(cmd2)

Gs <- data.matrix(read.table("newstats.txt", header = T))
Gobs <- data.matrix(read.table("newobsstats.txt", header = T))

if (!is.matrix(Gobs)) {
        Gobs <- matrix(Gobs, nrow = 1)
}

if(!is.null(obspar)){ 
        if(!is.matrix(obspar)){ 
                obspar<-matrix(obspar,byrow=T,ncol=ncol(param))
        }
        if(nrow(Gobs)!=nrow(obspar)){ 
                stop("Please supply observed statistics and observed parameter matrices with the same number of rows\n")
        }
}

nsims<-nrow(Gs)
ndatasets <- nrow(Gobs)
npar<-ncol(param)

cat("Doing ABC with:",npar,"parameters, ",ncol(Gs)," PLS components...\n")

vals<-err<-NULL

for (i in 1:ndatasets) {
        cat("dataset...", i, "\n")
	valsi<-abcmethod(Gobs[i, ],param,Gs,...)
        if(is.null(valsi$adj.values)){
                valsi<-valsi$unadj.values
        }
        else{
                valsi<-valsi$adj.values
        }
	if(is.null(obspar)|!(err.only)){
		vals<-abind(vals,valsi,along=3)
	}
	if(!is.null(obspar)){
		err[i]<-errfn(valsi,obspar[i,])
	}
}

l<-list()

if(is.null(obspar)){
        l$post.sample<-vals
}
else{
        l$err<-matrix(err)
	if(!err.only){
	        l$post.sample<-vals		
	}
}

return(l)

}
