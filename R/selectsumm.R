selectsumm <-function (obs,param, sumstats, obspar=NULL, ssmethod, verbose=TRUE,final.dens=FALSE, ...) {

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

argl <- as.list(match.call())

smargind <-match(names(argl),"ssmethod")
smargind <-which(!is.na(smargind))
sm<-as.character(argl[[smargind]])

ssargind <-match(names(argl),"sumsubs")
ssargind <-which(!is.na(ssargind))

if (length(ssargind)==0){
        sumsubs<-1:ncol(sumstats)
}
else{
	sumsubs<-eval(argl[[ssargind]])
}

sumstats<-sumstats[,sumsubs]
obs<-obs[,sumsubs]

if(!is.matrix(obs)){
        obs<-matrix(obs,nrow=1)
}
if(!is.matrix(sumstats)){
        sumstats<-as.matrix(sumstats)
}

ndatasets <- nrow(obs)
nstats<-ncol(obs)
   
err <- critvals <- best<-vals<-NULL

record<-matrix(0,ndatasets,length(sumsubs))

if(sm=="pls.abc"){
	tmp<-pls.abc(obs,param,sumstats, obspar=NULL,...)
	record<-NULL
	if(!is.null(tmp$err)){    
        	err<-tmp$err
	}
	if(!is.null(tmp$post.sample)){
        	vals<-tmp$post.sample
		final.dens<-TRUE
	}
}
else{
	for (i in 1:ndatasets) {
		if(verbose){
        		cat("dataset...", i, "\n")
		}
		resi<-ssmethod(obs[i,],param, sumstats, obspar[i,], verbose = verbose,final.dens=final.dens, ...)
		if(!is.null(resi$best)){
			record[i,resi$best]<-1
		}
        	if(!is.null(resi$err)){    
    		    err <- rbind(err, resi$err)
        	}
        	if(!is.null(resi$critvals)){    
        		critvals <- rbind(critvals, resi$critvals)
        	}
        	if (final.dens) {
			if(ncol(resi$post.sample)==1){
                		vals <- cbind(vals, resi$post.sample)
			}
			else{
                		vals <- abind(vals, resi$post.sample,along=3)
			}
        	}
	}
}

l<-list()

if (length(ssargind)!=0){
	l$sumsubs<-sumsubs
}

if(!is.null(critvals)){    
        l$critvals<-critvals
}
if(!is.null(err)){    
        l$err<-err
}
if(final.dens){
        l$post.sample<-vals
}

if(!is.null(resi$order)){
	l$order<-resi$order
}
if(!is.null(resi$sainfo)){
	l$sainfo<-resi$sainfo
}

if(sum(record)!=0){
	l$best<-record
}

return(l)

}
