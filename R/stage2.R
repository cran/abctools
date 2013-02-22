stage2 <-function (obs,param, sumstats, obspar=NULL, chosen, dsets = 100, sumsubs = 1:ncol(sumstats), limit = length(sumsubs), do.only = NULL, do.err=FALSE,final.dens=FALSE,verbose=TRUE,...) {

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

aargind <- match(names(argl),"abcmethod")

tolargind<-which(!is.na(aargind))

if (length(aargind)==0){
        stop("'abcmethod' is missing")
}

abcmethod<-eval(argl[[aargind]])

sumstats<-sumstats[,sumsubs]
obs<-obs[,sumsubs]


ndatasets <- nrow(obs)

err <- best<-vals<-NULL

nr<-nrow(param)
npar<-ncol(param)

eps2 <- dsets/nr
if(verbose){
 	cat("chosen subset is:", chosen, "\n")
}

cm <- combmat(length(sumsubs), limit)
if (is.null(do.only)) {
        do.only <- 1:nrow(cm)
}
if (nrow(cm) < chosen) {
        stop("value of chosen is too big for value of limit!")
}
if (length(chosen) == 1) {
        chosen <- which(cm[chosen, ] == 1)
}
l.true <- abc(obs[chosen], param, sumstats[, chosen], tol = eps2, ...)
closest <- which(l.true$region)
obss <- sumstats[closest,]
obst <- param[closest,]
    
tmp<-selectsumm(obss,param, sumstats, obst, sumsubs = sumsubs,limit = limit,do.only = do.only, do.crit = FALSE, do.err=TRUE,final.dens=FALSE, ...) 

ave<- rowMeans(tmp$err)
best<-which.min(ave)  
best <- do.only[best]

if (final.dens) {
        cat("getting final posterior sample...\n")
        Ifin <- which(cm[best, ] == 1)
        valsI <- abcmethod(obs[Ifin], param, sumstats[, Ifin], ...)
        if(is.null(valsI$adj.values)){
                post.sample<-valsI$unadj.values
        }
        else{
                post.sample<-valsI$adj.values
        }
}

l<-list()

l$best<-best
l$closest<-which(closest==TRUE)
l$order<-do.only
if(do.err){
	l$err<-ave
}

if (final.dens) {
        l$post.sample <- post.sample
}
l$sumsubs<-sumsubs


return(l)

}

