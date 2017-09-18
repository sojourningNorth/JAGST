


#' JAGST self-contained test
#'
#' JAGST self-contained test
#'
#' @param array expression array. Genes along rows.
#' @param design design matrix (intercept included by default)
#' @param inds indices of transcripts to be tested
#'
#' @export
#' @examples
#' dat <- matrix(rnorm(100*100),nrow=100,ncol=100)
#' design <- rep(0:1,each=50)
#' inds <- sample(100,10)
#' JAGSTself(arr,des,ind)
JAGSTself <- function(array,design,inds){
    require(limma)
        
    des <- cbind(1,design)
	ob <- limma::lmFit(array,design=des)
	obmod <- limma::eBayes(ob)
	statsCor <- obmod$t[,2]
	# could use unmoderated for self-contained test
	
	cov_mat <- cor(t(array[inds,])) # QR instead

	if(kappa(cov_mat)>1500) {
		increa <- svd(cov_mat)$d[1]/1500 
		cov_mat <- cov_mat + diag(increa,dim(cov_mat)[1])
		}	

	test_stat_cor <- t(statsCor[inds])%*%solve(cov_mat)%*%(statsCor[inds])
	
	p_self_cor <- pchisq(test_stat_cor[1],df=length(inds),lower.tail=F)
	return(p_self_cor)	
}





#' JAGST competitive test
#'
#' JAGST competitive test
#'
#' @param array expression array. Genes along rows.
#' @param design design matrix (intercept included by default)
#' @param inds indices of transcripts to be tested
#' @param its number of iterations for null generation.  This number times num_nest is the total number of samples for the null.
#' @param num_nest number of samples within each null iteration.  Cheap computationally, but decreasing returns efficiency gain.
#' @export
#' @examples
#' dat <- matrix(rnorm(100*100),nrow=100,ncol=100)
#' design <- rep(0:1,each=50)
#' inds <- sample(100,10)
#' JAGSTcomp(arr,des,ind,its,num_nest)
JAGSTcomp <- function(array,design,inds,its=200,num_nest=2000){
    require(limma)
    require(selectiveInference)

    len_inds <- length(inds)
	dim_x <- dim(array)[1]
    
    des <- cbind(1,design)
	ob <- limma::lmFit(array,design=des)
	obmod <- limma::eBayes(ob)
	statsCor <- obmod$t[,2]

	cov_mat <- cor(t(array[inds,])) # QR instead

	if(kappa(cov_mat)>1500) {
		increa <- svd(cov_mat)$d[1]/1500 
		cov_mat <- cov_mat + diag(increa,dim(cov_mat)[1])
		}

	test_stat_cor <- t(statsCor[inds])%*%solve(cov_mat)%*%(statsCor[inds])

	mat_trans <- t(array)
	stat_nulls <- list()
	for(ee in 1:its){
	null_sam_inds <- sample(1:dim_x,len_inds,replace=F)
	
	larfit = lar(mat_trans[,null_sam_inds],design) # is there a way to stop the algorithm early in general since aic will likely be achieved early in the calculation
	  out.aic <- larInf(larfit,type="aic")
	  ncps_cor <- round(out.aic$sign * out.aic$vmat %*% out.aic$y/(out.aic$sigma * sqrt(rowSums(out.aic$vmat^2))),3) # test stats
		
	  stat_nulls[[ee]] <- rchisq(2000,df=length(inds),ncp=sum(ncps_cor^2)) ## in these simulations would generally be assuming same size of correlated and uncorrelated test set so should not necessarily need to change df parameter but doing so here to prevent future issues
		}
	
	comp_test_p <- calc_emp_p(unlist(stat_nulls),test_stat_cor[1])
	return(comp_test_p)	
	}








#' ROAST directed test
#'
#' ROAST directed test
#'
#' @param array expression array. Genes along rows.
#' @param design design matrix (intercept included by default)
#' @param inds indices of transcripts to be tested
#' @export
#' @examples
#' dat <- matrix(rnorm(100*100),nrow=100,ncol=100)
#' design <- rep(0:1,each=50)
#' inds <- sample(100,10)
#' roastDir(arr,des,ind)
roastDir <- function(array,design,inds){
	require(limma)
    len_inds <- length(inds)
	dim_x <- dim(array)[1]
    
    des <- cbind(1,design)
	ob <- limma::lmFit(array,design=des)
	obmod <- limma::eBayes(ob)
	statsCor <- obmod$t[,2]
	
	test_stat <- sum(statsCor[inds])
	cov_mat <- cor(t(array[inds,])) # QR instead
	
	emp_dir <- 2*pnorm(-abs(test_stat),mean=0,sd=sqrt(sum(cov_mat)))
		
	return(emp_dir)
	}





#' ROAST undirected/mixed test
#'
#' ROAST undirected/mixed test
#'
#' @param array expression array. Genes along rows.
#' @param design design matrix (intercept included by default)
#' @param inds indices of transcripts to be tested
#' @param samp_size number of samples off of which empirical p-value is calculated (not permutation based, but sampling from parametric distribution)
#' @export
#' @examples
#' dat <- matrix(rnorm(100*100),nrow=100,ncol=100)
#' design <- rep(0:1,each=50)
#' inds <- sample(100,10)
#' roastMix(arr,des,ind,samp_size)
roastMix <- function(array,design,inds,samp_size=2000){
	require(limma)
    len_inds <- length(inds)
	dim_x <- dim(array)[1]
    
    des <- cbind(1,design)
	ob <- limma::lmFit(array,design=des)
	obmod <- limma::eBayes(ob)
	statsCor <- obmod$t[,2]

	test_stat <- sum((statsCor[inds])^2)
	cov_mat <- cor(t(array[inds,])) # QR instead

	eg <- eigen(cov_mat)
	mea <- rep(0,len_inds)
	un_norm_vecs <- eg$vectors
	deltas2 <- ((t(un_norm_vecs)%*%mea)/sqrt(eg$values))^2
	
	tmp <- list()
	for (i in seq_along(eg$value)){
		tmp[[i]] <- (eg$values[i])*rchisq(samp_size,df=1,ncp=deltas2[i])
		}
	vec <- rowSums(matrix(unlist(tmp),nrow=samp_size))
	emp <- calc_emp_p(vec,test_stat)
	return(emp)
	
	}


#' calculate empirical p-value
#'
#' Not for export
#' @param dist distribution
#' @param num test statistic
calc_emp_p <- function(dist,num){
	L <- length(dist)	
	plac <- which(order(c(dist,num))==(L+1))
	pval <- 1-plac/(L+2)
	return(pval)	
	}
	
	
#' calculate map reals to (0,1)
#'
#' Not for export
#' @param x argument to be evaluated by expit	
expit <- function(x){
	exp(x)/(1+exp(x))
	}

	
	
	
	