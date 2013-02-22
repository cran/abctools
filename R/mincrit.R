mincrit <-function (obs,param, sumstats, obspar=NULL, abcmethod=abc,crit=nn.ent,sumsubs=1:ncol(sumstats), limit = length(sumsubs), do.only = NULL, verbose = TRUE, do.crit = TRUE, do.err=FALSE,final.dens=FALSE,errfn=rsse, ...) 
{

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


sumstats<-sumstats[,sumsubs]
obs<-obs[,sumsubs]

data <- matrix(obs, nrow = 1)
npar<-ncol(param)
nr<-nrow(sumstats)
nstats<-length(sumsubs)

q<- (!is.null(obspar))&(do.err)

if(!q){
    	do.err<-FALSE
}

cm <- combmat(nstats, limit)
nc <- nrow(cm)
if (is.null(do.only)) {
	do.only <- 1:nc
}

vals<-err<-NULL

critvals <- matrix(0, 1, length(do.only))
best<-1

for (i in 1:length(do.only)) {
        I <- which(cm[do.only[i], ] == 1)
        if (verbose) {
        	cat("doing statistics: ", sumsubs[I], "   (", i, "/", length(do.only), ")\n")
        }

	valsi<-abcmethod(data[I], param, sumstats[,I], ...)
	if(is.null(valsi$adj.values)){
                valsi<-valsi$unadj.values
        }
        else{
                valsi<-valsi$adj.values
        }
        if(do.err){
                err[i]<-errfn(valsi,obspar,apply(param,2,var))
        }	
        if (do.crit) {
        	cat("doing criterion calculation...\n")
#        	cat("dimension of theta is:", dim(valsi), "\n")
        	critvals[i] <- crit(valsi)
  		stick<-(critvals[best]<=critvals[i])
		best<-ifelse(stick,best,i)
		if(final.dens){
			if((i==1)|!stick){
				vals<-valsi
			}
		}
        }
}
err<-err[best]

best <- which(cm[do.only[best], ] == 1)

l<-list()

if(do.crit){
	l$critvals<-critvals
	best<-do.only[best]
	l$best<-best
}
if(do.err){
	l$err<-err
}
if(final.dens){
	l$post.sample<-vals
}
l$order<-do.only
l$sumsubs<-sumsubs

return(l)
}

