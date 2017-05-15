

#' AAC score
#'
#' @param mat a numeric matrix. AAC score is calculated by columns 
#' @param cor_method pass to [stats::cor()]
#' @param min_cor minimal absolute correlation
#' 
#' @details 
#' AAC score for a given item is the area above the curve of the curmulative 
#' distribution of the absolute correlation to all other items.
#'
#' @return A vector of AAC scores
#' @export
#' @import stats
#'
#' @examples
#' set.seed(12345)
#' require(matrixStats)
#' nr1 = 100
#' mat1 = matrix(rnorm(100*nr1), nrow = nr1)
#' 
#' nr2 = 10
#' require(mvtnorm)
#' sigma = matrix(0.8, nrow = nr2, ncol = nr2); diag(sigma) = 1
#' mat2 = t(rmvnorm(100, mean = rep(0, nr2), sigma = sigma))
#' 
#' nr3 = 50
#' sigma = matrix(0.5, nrow = nr3, ncol = nr3); diag(sigma) = 1
#' mat3 = t(rmvnorm(100, mean = rep(0, nr3), sigma = sigma))
#' 
#' mat = t(rbind(mat1, mat2, mat3))
#' AAC_cor(t(mat))
AAC_cor = function(mat, cor_method = "pearson", min_cor = 0.2) {
	n = ncol(mat)
	v = numeric(n)
	for(i in 1:n) {
		suppressWarnings(cor_v <- abs(cor(mat[, i, drop = FALSE], mat[, -i, drop = FALSE], method = cor_method)))
		if(sum(is.na(cor_v))/length(cor_v) >= 0.75) {
			v[i] = 1
		} else {
			f = ecdf(cor_v)
			cor_v = seq(min_cor, 1, length = 100)
			n2 = length(cor_v)
			v[i] = sum((cor_v[2:n2] - cor_v[1:(n2-1)])*f(cor_v[-n2]))
		}
		# if(interactive()) {
		# 	cat(strrep("\b", 100))
		# 	cat("AAC of the ecdf of correlations for ", i, "/", n, sep = "")
		# }
	}
	# if(interactive()) cat("\n")
	v = (1 - min_cor) - v
	return(v)
}

entropy = function(p) {
	n = length(p)
	p = p[p > 0]
	p2 = rep(1/n, n)
	-sum(p*log2(p))/abs(sum(p2*log2(p2)))
}

PAC = function(consensus_mat, x1 = 0.1, x2 = 0.9) {
	Fn = ecdf(consensus_mat[lower.tri(consensus_mat)])
	Fn(x2) - Fn(x1)
}


#' Cophenetic correlation coefficient in all methods
#'
#' @param res_list a `run_all` object
#' @param top_method a vector of top methods
#' @param partition_method a vector of partition methods
#'
#' @export
#' @import ggplot2
mean_cophcor = function(res_list, top_method = res_list$top_method, partition_method = res_list$partition_method) {
	
	df = data.frame(top_method = character(0), partition_method = character(0), k = numeric(0), mean_cophcor = numeric(0))
	for(i in seq_along(top_method)) {
	    for(j in seq_along(partition_method)) {    
	    	res = get_single_run(res_list, top_method = top_method[i], partition_method = partition_method[j])

	        for(ik in seq_along(res$k)) {
	        	k = res$k[ik]
	        	df = rbind(df, data.frame(top_method = top_method[i], partition_method = partition_method[j], k = k, 
	        		cophcor = cophcor(res$object_list[[ik]]$consensus)))
	        }
	    }
	}

	gp = ggplot(df, aes(x = k, y = cophcor)) + geom_line() + geom_point() + facet_grid(top_method ~ partition_method)
	print(gp)
}
