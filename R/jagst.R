


#' JAGST self-contained test
#'
#' JAGST self-contained test
#' @param array expression array
#' design design matrix (intercept included by default)
#' inds indices of transcripts to be tested
#' @keywords 
#' @export
#' @examples
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

	test_stat_cor <- t(statsCor[inds])%*%solve(cov_mat[inds,inds])%*%(statsCor[inds])
	
	p_self_cor <- pchisq(test_stat_cor[1],df=length(inds),lower.tail=F)
	return(p_self_cor)	
}





#' JAGST competitive test
#'
#' JAGST competitive test
#' @param array expression array
#' design design matrix (intercept included by default)
#' inds indices of transcripts to be tested
#' its number of iterations for null generation.  This number times num_nest is the total number of samples for the null.
#' num_nest number of samples within each null iteration.  Cheap computationally, but decreasing returns efficiency gain.
#' @keywords 
#' @export
#' @examples
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

	test_stat_cor <- t(statsCor[inds])%*%solve(cov_mat[inds,inds])%*%(statsCor[inds])

	mat_trans <- t(array)
	stat_nulls <- list()
	for(ee in 1:its){
	null_sam_inds <- sample(1:dim_x,len_inds,replace=F)
	
	larfit = lar(mat_trans[,null_sam_inds],design) # is there a way to stop the algorithm early in general since aic will likely be achieved early in the calculation
	  out.aic <- larInf(larfit,type="aic")
	  ncps_cor <- round(out.aic$sign * out.aic$vmat %*% out.aic$y/(out.aic$sigma * sqrt(rowSums(out.aic$vmat^2))),3) # test stats
		
	  stat_nulls[[ee]] <- rchisq(2000,df=length(inds),ncp=sum(ncps_cor^2)) ## in these simulations would generally be assuming same size of correlated and uncorrelated test set so should not necessarily need to change df parameter but doing so here to prevent future issues
		}
	
	comp_test_p <- calc_emp_p(unlist(stat_nulls,test_stat_cor[1])
	return(comp_test_p)	
	}




#' calculate empirical p-value
#'
#' Not for export
#' @param 
#' @keywords 
calc_emp_p <- function(dist,num){
	L <- length(dist)	
	plac <- which(order(c(dist,num))==(L+1))
	pval <- 1-plac/(L+2)
	return(pval)	
	}
	
	
	
	