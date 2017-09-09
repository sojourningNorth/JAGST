


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
	# unmoderated for self-contained test

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
#' @keywords 
#' @export
#' @examples
#' JAGSTcomp(arr,des,ind)
 
JAGSTcomp <- function(array,design,inds){
    require(limma)
    require(selectiveInference)
        
    des <- cbind(1,design)
	ob <- limma::lmFit(array,design=des)
	obmod <- limma::eBayes(ob)
	statsCor <- obmod$t[,2]
	# unmoderated for self-contained test

	test_stat_cor <- t(statsCor[inds])%*%solve(cov_mat[inds,inds])%*%(statsCor[inds])
	te_null_cor <- JAGST::easier_null(len_inds,dim_x,its,num_nest,cancer_cor)

	
	
	
	
	p_self_cor <- pchisq(test_stat_cor[1],df=length(inds),lower.tail=F)
	return(p_self_cor)
}

















#' Calculate null distribution for JAGST competitive test
#'
#' null distribution competitive test
#' @param array expression array
#' design design matrix (intercept included by default)
#' inds indices of transcripts to be tested
#' @keywords 
easier_null <- function(len_inds,dim_x,its,num_nest,cancer_cor){
	stat_nulls <- list()
	for(ee in 1:its){
	null_sam_inds <- sample(1:dim_x,len_inds,replace=F)
	
	larfit = lar(mat_trans[,null_sam_inds],cancer_cor) # is there a way to stop the algorithm early in general since aic will likely be achieved early in the calculation
	  out.aic <- larInf(larfit,type="aic")
	  ncps_cor <- round(out.aic$sign * out.aic$vmat %*% out.aic$y/(out.aic$sigma * sqrt(rowSums(out.aic$vmat^2))),3) # test stats
		
	  stat_nulls[[ee]] <- rchisq(num_nest,df=length(inds1),ncp=sum(ncps_cor^2)) ## in these simulations would generally be assuming same size of correlated and uncorrelated test set so should not necessarily need to change df parameter but doing so here to prevent future issues
		}
		return(unlist(stat_nulls))	
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