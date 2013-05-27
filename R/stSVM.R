#' Cross validation for smoothed t-statistic to select significant top ranked differential expressed genes
#'
#'
#' @param x: a p x n matrix of expression measurements with p samples and n genes.
#' @param y: a factor of length p comprising the class labels.
#' @param DEBUG: show debugging information in screen or not.
#' @param Cs: soft-margin tuning parameter of the SVM. Defaults to \code{10^c(-3:3)}.
#' @param Gsub: an adjacency matrix that represents the underlying biological network.
#' @param pt.pvalue: cut off p value of permutation test
#' @param pt.step:  permutation test steps
#' @param op: optimal on top op% ranked genes
#' @param a: constant value of random walk kernel
#' @param p: random walk step(s) of random walk kernel
#' @param allF: Using all features (TRUE) or only these genes mapped to prior information (FALSE). 
#' @seed  seed: seed for random sampling.
#' @return a list with the results of the cross-validation.  \item{feat}{the selected features}  \item{fit}{the fitted SVM model} \item{auc}{the prediction auc of the fitted SVM model} \item{lables}{the test labels}  
#' @references Y.Cun and H. Frohlcih, (2013), Data Integration via Network Smoothed T-Statistics, submitted    
#' @export
#' @note The optimal number of features is found by using the span estimate. See Chapelle, O., Vapnik, V., Bousquet, O., and Mukherjee, S. (2002). Choosing multiple parameters for support vector machines. \emph{Machine Learning}, 46(1), 131-159.
#'
#' @author  Yupeng Cun \email{yupeng.cun@gmail.com}
#' @examples
#' 
#' data(expr)
#' x <- expr$genes
#' y <- expr$y
#'
#' library(netClass)
#' 
#' r.stsvm <- stsvm.cv(x=x,x.mi=NULL,y=y, folds=5,Gsub=ad.matrix,  repeats=3, parallel = TRUE, cores=4,DEBUG=TRUE,pt.pvalue=0.05,op=0.9,pt.step=1000,a=1,p=2,allF=TRUE,Cs=10^(-4:4))
#' mean(r.stsvm$auc)


stsvm.cv <- function(x=x, x.mi=NULL,y=y, folds=5,Gsub=matrix(1,100,100),repeats=3, parallel = FALSE, cores=2,DEBUG=TRUE, pt.pvalue=0.05,op=0.85,pt.step=1000,a=1,p=2,allF=TRUE,seed=1234,Cs=10^c(-3:3))
{
	multicore <- ("package:multicore" %in% search())
  
	if(multicore == TRUE && parallel == TRUE)  
	{
    
		if(is.null(cores))       
			cores <- multicore:::detectCores()    
		options(cores = cores - 1)     
		
		cat("There's ", cores," cores detected, and ",getOption("cores")," will used\n")    
		parallel <- TRUE  
	}  
	else  
	{
		if(parallel == TRUE)       
		cat('You nedd Package \'multicore\'\n')
		cat('Computation will performed sequential crossvalidation.\n', sep='')
		parallel <- FALSE
	}

 	n     <- length(y)
	folds <- trunc(folds)

	if (folds < 2) stop("folds =< 2.\n")
	if (folds > n) stop("folds > the number of observations.\n")
	
  	dk.tf = pt.pvalue
  	pred.type="class"
  	aa=pt.step
  	mapping=NULL
	ex.sum=x
  	
	if(allF)
	{
		
		int.all = colnames(x)
		cat("All ",length(int.all),"features.\n")
		ad.all <- matrix(0, ncol=length(int.all), nrow=length(int.all))
		colnames(ad.all) <- colnames(x)
		rownames(ad.all) <- colnames(x)
		inta = intersect(int.all, colnames(Gsub))
		ad.all[inta,inta]= Gsub[inta,inta]		
		Gsub <- ad.all
	}	
  
	
	int = intersect(colnames(Gsub),colnames(x))
	x=x[,int]
	Gsub=Gsub[int,int]
	
	#Calculating random walk kernel matrix of network.
	dk <- calc.diffusionKernelp(L=Gsub, is.adjacency=TRUE,p=p,a=a)
	
	#if(fs.method=="TRG+SVM")dk <- calc.diffusionKernelp(L=Gsub[int,int], is.adjacency=TRUE,p=p,a=a)
	
	cuts  <- cv.repeats <- list()  
				
	set.seed(1234)	
	for(r in 1:repeats)
	{		
		###random sampling 
		stratified=FALSE
		if(stratified){
			perm0 = sample(which(y == levels(y)[1]))
			perm1 = sample(which(y == levels(y)[2]))
			perm = c(perm0, perm1)
		}
		else
			perm = sample(1:n)		
		#perm <- sample(1:n) #Sampling a random integer between 1:n
		repeat.models <- NULL 
    
		for(k in 1:folds) #randomly divide the training set in to 10 folds
		{
			tst <- perm[seq(k, n, by=folds)]  #      
			trn <- setdiff(1:n, tst)            
			cuts[[k]] <- list(trn=trn, tst=tst)    
		}    
		
		#pb <- txtProgressBar(min = 0, max = folds, style = 3)   #txtProgressBar id funtion in packages "utils", Text progress bar in the R console
    
		if(DEBUG) cat('Starting classification of repeat:',r,'\n')
    
		if(parallel)	repeat.models <- mclapply(1:folds, classify.stsvm, cuts=cuts, ex.sum=ex.sum, x=x, y=y,cv.repeat=r, DEBUG=DEBUG, Gsub=Gsub,op=op,aa=aa,dk=dk,dk.tf=dk.tf,p=p,a=a,seed=seed,Cs=Cs)
		else		repeat.models <-   lapply(1:folds, classify.stsvm, cuts=cuts, ex.sum=ex.sum,x=x, y=y,cv.repeat=r, DEBUG=DEBUG, Gsub=Gsub,op=op,aa=aa,dk=dk, dk.tf=dk.tf,p=p,a=a,seed=seed,Cs=Cs)
	
		#close(pb)
    
		if(length(repeat.models) != folds)
		{      
			geterrmessage()      
			stop("One or more processes did not return. May be due to lack of memory.\n")
		}
	
		if(DEBUG)	cat('All models of repeat:',r,'have been trained.\n')
    
		cv.repeats[[r]] <- repeat.models
  
  	} 
	 
  	auc <- sapply(cv.repeats, function(cv.repeat) sapply(cv.repeat, function(model) model$auc))  
  	colnames(auc) <- paste("Repeat",1:repeats,sep="")  
  	rownames(auc) <- paste("Fold",1:folds,sep="")    	
  	
  	fits <- lapply(cv.repeats,function(cv.repeat) lapply(cv.repeat, function(model) model$model))  
  	names(fits) <- paste("Repeat",1:repeats,sep="")  
  	fits <- lapply(fits, function(x)  {names(x) = paste("Fold", 1:folds, sep = ""); x })  
  	
  	feat <- lapply(cv.repeats,function(cv.repeat) lapply(cv.repeat, function(model) model$feat))  
  	names(feat) <- paste("Repeat",1:repeats,sep="")  
  	feat <- lapply(feat, function(x)  {names(x) = paste("Fold", 1:folds, sep = ""); x })
  	
  	res <- list(feat=feat, auc=auc,fits=fits, labels=y)  
  	class(res) <- 'netClassResult' 
  	#res<-cv.repeats
  	return(res)  
	#return ( cv.repeats)
}

# Training and predicting using stSVM classification methods
#-----------------------------  
# x: a p x n matrix of expression measurements with p samples and n genes.
# y: a factor of length p comprising the class labels.
# DEBUG: show debugging information in screen or not.
# Cs: soft-margin tuning parameter of the SVM. Defaults to \code{10^c(-3:3)}.
# Gsub: an adjacency matrix that represents the underlying biological network.
# pt.pvalue: cut off p value of permutation test
# aa:  permutation test steps
# op: optimal on top op% ranked genes
#-----------------------------  
# Return a list with the results of the training model and predcit AUC   


classify.stsvm <- function(fold, cuts,ex.sum, x, p,a, y, cv.repeat, DEBUG=DEBUG,Gsub=Gsub, op=op,aa=aa,dk=dk,dk.tf=dk.tf,seed=seed,Cs=Cs)
{
	gc() 		
	
	if(DEBUG) cat("stSVM ==>> starting Fold:",fold,"\n")	
	
  
	## get training and test indices
	trn <- cuts[[fold]]$trn
	tst <- cuts[[fold]]$tst		
	#print(tst)
	label 	<- sign(as.numeric(y[tst]) - 1.5) # because y is a factor			
	
	## train and test the model						
		
		
		xmmiTrn 	<- x[trn,]#,x.mi[trn,])#sig.mi])
		xmmiTst 	<- x[tst,]#,x.mi[tst,]#sig.mi])
		
		sca=TRUE		
		rank.method="dk"
		
		# ranked by smoothed t-statistics  				
		ranks = getGraphRank(x=x[trn,], y=y[trn], Gsub=dk,sca=sca)	
		
						
		topGenes = names(ranks[which(ranks > quantile(ranks,op))])				
		#ranksI=sort(ranks,decreasing=T)
		#topGenes = names(ranksI[1:op])
		#cat("==================\nrange(ranks[topGenes]):",range(ranks[topGenes]),"\n")
		#topRanks =  getGraphRank1(x=x[trn,topGenes], y=y[trn], Gsub=dk,sca=sca)
		topRanks =ranks[topGenes]
		
		#topRan
		#topRanks = rank(topRanks)
		#topRanks=ranks
		if(DEBUG) cat("range of topRanks:",range(topRanks),"\n=======================\n")
		
		ytrn = y[trn]
		xtrn =xmmiTrn[,topGenes]			

		R = aa
		set.seed(seed)
		p_ranks=matrix(0,ncol=length(topRanks))[1,]
		names(p_ranks)=topGenes
		#cat("length(topRanks) ",length(topRanks),"\n")
		
		if(DEBUG) cat("\nStart permutation test for significant topRanks. \n")	
		
		for(i in 1 :R )
		{							
			ys = sample(ytrn)															
			rankt <-   getGraphRank(x=xtrn, y=ys, Gsub=dk, sca=sca)			
			for(j in topGenes) if(rankt[j] >= topRanks[j])p_ranks[j] = p_ranks[j]+1					
		}	
		cat("head(p_ranks)   ", head(p_ranks),"\n")
		p_ranks = p_ranks/R
		#cat("head(p_ranks)/R ", head(p_ranks),"\n")
		#cat("range(rankt) ", range(rankt),"\n")
		#cat("head(rankt) ", head(rankt),"\n")
		p_ranks = p.adjust(p_ranks, method="BH")
		names(p_ranks)=names(topRanks)
		
		sigRank <- names(p_ranks[which(p_ranks < dk.tf)])		
		if(length(sigRank)==0)sigRank=names(sort(topRanks))[1:10]							
		
		feat  = sigRank					
		
		train 	 <- svm.fit(x=x[trn,feat], y=y[trn], Cs=Cs, scale="scale", DEBUG=FALSE) 
		test    <- svm.predict(fit=train, newdata=x[tst,], type="response")
		
		auc	<- calc.auc(test, label)				
			
		#auc3=spec
		if(DEBUG) cat("fold:",fold,"\tlength(test):",length(test),"\tlength(label)",length(label),"\n")			
	
		if(DEBUG) cat("\n\n=> the best AUC	is ", auc, "\t Best features length:  ", length(feat),"\n\n")			
		if(DEBUG) {cat('Finished fold:',fold,'\n\n')}
		#else {setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)}
		gc()	
		res=list(fold=fold, model=train, auc=auc, feat= feat)
	
	
	return(res)
}


## 
# Computing the Random Walk Kernel matrix of network
#-----------------------------  
# L: an adjacency matrix that represents the underlying biological network.
# a: constant value of random walk kernel
# p: #(p) random walk step(s) of random walk kernel
#-----------------------------  

calc.diffusionKernelp = function(L, is.adjacency=TRUE, p=3,a=2)
{
  if(missing(L)) stop('You have to provide either a laplace or a adjacency matrix to calculate the diffusion kernel.')
  print("thresh")
  
  if(is.adjacency)
  {    
    dnames <- dimnames(L)
    L = graph.adjacency(L, mode='undirected', diag=FALSE)
    L = graph.laplacian(L, normalized=TRUE)
    dimnames(L) = dnames
  }
  
  n=ncol(L)
  I=diag(n)      
      
  if( p==1) R = a*I-L
    
  else
  {
    R=a*I -L
    for(ii in seq(from=2,to=p, by=1))
    {
      R = R %*% (a*I -L)
      ii = ii+1
    }
  }
  
  #KernelDiag <- sqrt(diag(R) + 1e-10)
  #R.norm 	<-  R/(KernelDiag %*% t(KernelDiag))
  #R			<-  R.norm  
  
  colnames(R) =colnames(L)
  rownames(R) =rownames(L)
  R    
}

##### random walk kernel matrix smoothing t-statistic
#------------------------------------
# x: a p x n matrix of expression measurements with p samples and n genes.
# y: a factor of length p comprising the class labels. 
# Gsub: Random Walk Kernel matrix of network
# -----------------------------------------
getGraphRank = function(x=x, y=y, Gsub=Gsub,sca=TRUE)
{    
    int= intersect(colnames(x),colnames(Gsub))  
    x=x[,int]
    Gsub=Gsub[int,int]   
    #print(length(int)) 
    if(sca==TRUE) x = scale(x)	
    #calculate the t-scorce of each probe.
    xtt = rep(0, ncol(x))
    yy=sign(as.numeric(y)-1.5)
    yy[yy==-1]=0

    for(i in 1: ncol(x))  xtt[i]= t.test(x[,i], yy,paired=T)$statistic    	   
    
    names(xtt)= colnames(x)	
    exprs = abs(xtt)
    exprs = exprs/ max(exprs)       	
    ranks   <- t(exprs[int]) %*% Gsub[int,int]
    r=ranks[1,]
    r= r/max(r)
    #print(dim(ranks))
   # print(dim(Gsub))
    names(r) = colnames(x)	 	   
   # print("finised GrapRank")
    return(r)
}				


