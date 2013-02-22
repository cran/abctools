AS.select <-function(obs,param,sumstats,obspar=NULL,abcmethod=abc,grid=10,inturn=FALSE,limit=ncol(sumstats),allow.none=TRUE,do.err=FALSE,final.dens=FALSE,errfn=rsse,...){

argl <- as.list(match.call())

tolargind <- match(names(argl),"tol")

tolargind<-which(!is.na(tolargind))

if (length(tolargind)==0){
	stop("'tol' is missing")
}

eps<-eval(argl[[tolargind]])

supp<-range(param)

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


nstats<-ncol(sumstats)

limit<-min(nstats,limit)

order<-sample(1:nstats,limit,F)
cat("order is:",order,"\n")

I<-numeric(0)

data<-obs

for(j in 1:length(order)){
	cat("I is:",I,", testing:", c(I,order[j]),"\n")
	
	if((length(I)==0)){
		if(allow.none){
			index1<-sample(param,ceiling(length(param)*eps),T)
			index2<-abcmethod(data[c(I,order[j])],param,sumstats[,c(I,order[j])],...)
			if(is.null(index2$adj.values)){
				index2<-index2$unadj.values
			}
			else{
				index2<-index2$adj.values
			}
			add<-AS.test(grid,index1,index2,supp)			
		}
		else{
			add<-TRUE
		}
	}
	else{
		index1<-abcmethod(data[I],param,sumstats[,I],...)
		if(is.null(index1$adj.values)){ 
                	index1<-index1$unadj.values
                }
                else{
                        index1<-index1$adj.values
                }
		index2<-abcmethod(data[c(I,order[j])],param,sumstats[,c(I,order[j])],...)
		if(is.null(index2$adj.values)){ 
                	index2<-index2$unadj.values
                }
                else{
                	index2<-index2$adj.values
                }
		add<-AS.test(grid,index1,index2,supp)			
	}


	if(add){
							#now do cyclic dropping procedure.........

		I<-c(I,order[j])				#(present) new "in" and "out" subsets
								#note that I should stay in the order that the stats were added
		bad<-NULL					#vector to hold removed statistics
		if(length(I)>1){
			if(inturn){					#########inturn means drop after testing (sequential dropping)
				for(i in 1:(length(I)-1)){		#cycle through previously added stats
					subset2<-I
					subset1<-setdiff(subset2,I[i])

					if(length(subset1)==0){		#this could happen for inturn=T
						if(allow.none){
							i1<-sample(param,ceiling(length(param)*eps),T)	
							i2<-abcmethod(data[subset2],param,sumstats[,subset2],...)
							if(is.null(i2$adj.values)){ 
								i2<-i2$unadj.values
							}
							else{
								i2<-i2$adj.values
							}
							add2<-AS.test(grid,i1,i2,supp)
						}
						else{
							add2<-TRUE
						}
					}
					else{
						i1<-abcmethod(data[subset1],param,sumstats[,subset1],...)
						if(is.null(i1$adj.values)){
                                                        i1<-i1$unadj.values
                                                }
                                                else{
                                                	i1<-i1$adj.values
                                                }
						i2<-abcmethod(data[subset2],param,sumstats[,subset2],...)
						if(is.null(i2$adj.values)){
                                                	i2<-i2$unadj.values
                                                }
                                                else{
                                                	i2<-i2$adj.values
                                                }
						add2<-AS.test(grid,i1,i2,supp)
					}
					if(add2){
						I<-subset2			#add to subset of rejected summaries
					}
					else{
						I<-subset1
					}
				}
		#	I<-setdiff(I,bad)				#get final "in" and "out".
			}
			else{						#########drop all bad stats afterwards
				for(i in 1:(length(I)-1)){              #cycle through previously added stats
                                	subset2<-I
                                	subset1<-setdiff(subset2,I[i])
					if(length(subset1)==0){
						if(allow.none){
							i1<-sample(param,ceiling(length(param)*eps),T)	
							i2<-abcmethod(data[subset2],param,sumstats[,subset2],...)
							if(is.null(i2$adj.values)){
                                                                i2<-i2$unadj.values
                                                        }
                                                        else{
                                                                i2<-i2$adj.values
                                                        }
							add2<-AS.test(grid,i1,i2,supp)
						}
						else{
							add2<-TRUE
						}
					}
					else{
						i1<-abcmethod(data[subset1],param,sumstats[,subset1],...)
						if(is.null(i1$adj.values)){
                                                        i1<-i1$unadj.values
                                                }
                                                else{
                                                	i1<-i1$adj.values
                                                }
						i2<-abcmethod(data[subset2],param,sumstats[,subset2],...)
						if(is.null(i2$adj.values)){
                                                	i2<-i2$unadj.values
                                                }
                                                else{
                                                        i2<-i2$adj.values
                                                }
	                                	add2<-AS.test(grid,i1,i2,supp)
			        	}
                                	if(!add2){
                                        	bad<-c(bad,I[i])
                                	}
                        	}
                	I<-setdiff(I,bad)                               #get final "in" and "out".
			}
		}
	} #end if add
} #end for

l<-list()

best<-I

#print(I)

if((length(I)==0)){
	vals<-sample(param,ceiling(length(param)*eps),T)
}
else{
	vals<-abcmethod(obs[I], param, sumstats[,I], ...)
	if(is.null(vals$adj.values)){
		vals<-vals$unadj.values
	}
	else{
		vals<-vals$adj.values
	}

	if(do.err){
		err<-errfn(vals,obspar,apply(param,2,var))
	}
}

if(do.err){
        l$err<-err
}
if(final.dens){
        l$post.sample<-matrix(vals,ncol=1)
}
l$best<-best

return(l)

}
